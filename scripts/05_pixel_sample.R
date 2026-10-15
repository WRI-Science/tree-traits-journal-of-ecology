#!/usr/bin/env Rscript
# 05 — Stratified random pixel sample by FRG
#
# Question: Which pixels should Fig 2 / Fig 6 / Table S2 use under BA weighting?
# Why: Main text samples pixels, not unique tm_ids. Hypergeometric reservoir
#      sampling while streaming TreeMap + FRG keeps equal quota per FRG in
#      constant memory (do not bind every eligible pixel).
# Inputs:  study TreeMap raster, FRG raster, LF16 lookup, derived/tm_id_cwmts.csv
# Outputs: OUT_ROOT/shared/pixel_sample.csv
#          OUT_ROOT/shared/pixel_sample_manifest.csv
# Run:     Rscript scripts/05_pixel_sample.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(terra)

# ============================================================
# ⚙️ Settings ----
# ============================================================
N_PER_FRG <- as.integer(Sys.getenv("TREE_TRAITS_PIXEL_SUBSAMPLE_PER_FRG", "10000"))
SAMPLE_SEED <- as.integer(Sys.getenv("TREE_TRAITS_PIXEL_SUBSAMPLE_SEED", "123"))
ROW_CHUNK <- as.integer(Sys.getenv("TREE_TRAITS_BA_ROW_CHUNK", "1024"))
NICE_VALUE <- as.integer(Sys.getenv("TREE_TRAITS_BA_NICE", "10"))
FRG_LEVELS <- c("I", "III", "IV", "V")

# ============================================================
# 📁 File paths ----
# ============================================================
OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
CWM_CSV <- file.path(OUT_ROOT, "derived", "tm_id_cwmts.csv")
TM_RAST <- Sys.getenv("TREE_TRAITS_STUDY_TREEMAP_TIF", "")
if (!nzchar(TM_RAST)) {
  stop(
    "Set TREE_TRAITS_STUDY_TREEMAP_TIF to the path of: Study-region TreeMap raster GeoTIFF (derived from TreeMap 2016)"
  )
}
FRG_RAST <- Sys.getenv("TREE_TRAITS_FRG_TIF", "")
if (!nzchar(FRG_RAST)) {
  stop(
    "Set TREE_TRAITS_FRG_TIF to the path of: LANDFIRE FRG raster masked to study region (https://www.landfire.gov/fire-regime/frg)"
  )
}
FRG_LOOKUP_CSV <- Sys.getenv("TREE_TRAITS_FRG_LOOKUP_CSV", "")
if (!nzchar(FRG_LOOKUP_CSV)) {
  stop(
    "Set TREE_TRAITS_FRG_LOOKUP_CSV to the path of: LANDFIRE FRG attribute lookup CSV (LF 2016)"
  )
}
SAMPLE_CSV <- file.path(OUT_ROOT, "shared", "pixel_sample.csv")
MANIFEST_CSV <- file.path(OUT_ROOT, "shared", "pixel_sample_manifest.csv")
VRT_PATH <- file.path(OUT_ROOT, "shared", "tmp_frg_on_tm.vrt")

# ============================================================
# 🔧 Helpers ----
# ============================================================
log_pt <- function(...) {
  message(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/Los_Angeles"),
    " ",
    paste0(...)
  )
  flush(stderr())
}

# ============================================================
# 1. 🧹 Session ----
# ============================================================
if (requireNamespace("tools", quietly = TRUE)) {
  try(tools::psnice(value = NICE_VALUE), silent = TRUE)
}
force <- identical(Sys.getenv("TREE_TRAITS_FORCE_RESAMPLE", ""), "1")
if (file.exists(SAMPLE_CSV) && !force) {
  log_pt("Existing pixel sample: ", SAMPLE_CSV)
  quit(save = "no", status = 0L)
}
need <- c(CWM_CSV, TM_RAST, FRG_RAST, FRG_LOOKUP_CSV)
missing <- need[!file.exists(need)]
if (length(missing)) {
  stop("Missing inputs:\n", paste(missing, collapse = "\n"))
}
dir.create(dirname(SAMPLE_CSV), recursive = TRUE, showWarnings = FALSE)

valid_ids <- unique(fread(CWM_CSV, select = "tm_id")$tm_id)
max_id <- max(valid_ids)
id_ok <- rep(FALSE, max_id)
id_ok[valid_ids] <- TRUE
log_pt("Valid tm_ids for sample: ", length(valid_ids))

lookup <- fread(FRG_LOOKUP_CSV)
lookup_clean <- lookup[
  !is.na(FRG_NEW) & FRG_NEW != "-9999" & FRG_NEW != "-1111" &
    FRG_NEW != -9999 & FRG_NEW != -1111
]
frg_codes <- unique(lookup_clean[, .(FRG_NEW)])
setorder(frg_codes, FRG_NEW)
frg_codes[, FRG_code := .I]
frg_codes[, collapsed_frg := fcase(
  FRG_NEW %chin% c("I-A", "I-B", "I-C"), "I",
  FRG_NEW %chin% c("III-A", "III-B"), "III",
  FRG_NEW %chin% c("IV-A", "IV-B"), "IV",
  FRG_NEW %chin% c("V-A", "V-B"), "V",
  default = NA_character_
)]

# ============================================================
# 2. 🗺️ Align FRG onto TreeMap grid ----
# ============================================================
tm <- rast(TM_RAST)
if (file.exists(VRT_PATH)) {
  unlink(VRT_PATH)
}
te <- as.vector(ext(tm))
warp_status <- system2(
  "gdalwarp",
  c(
    "-r", "near", "-of", "VRT",
    "-te", te[1], te[3], te[2], te[4],
    "-tr", "30", "30",
    "-dstnodata", "65535",
    FRG_RAST, VRT_PATH
  ),
  stdout = TRUE, stderr = TRUE
)
if (!is.null(attr(warp_status, "status")) && !identical(attr(warp_status, "status"), 0L)) {
  stop("gdalwarp failed:\n", paste(warp_status, collapse = "\n"))
}
frg_on_tm <- rast(VRT_PATH)

# ============================================================
# 3. 🎲 Reservoir-sample eligible pixels by pixel FRG ----
# ============================================================
# Pixel FRG (not tm_id majority). Hypergeometric batch reservoir keeps a
# uniform sample of size N_PER_FRG without storing every forest pixel.
n_rows <- nrow(tm)
row_starts <- seq(1L, n_rows, by = ROW_CHUNK)
n_chunks <- length(row_starts)
readStart(tm)
readStart(frg_on_tm)
on.exit({
  try(readStop(tm), silent = TRUE)
  try(readStop(frg_on_tm), silent = TRUE)
}, add = TRUE)

code_lu <- rep(NA_character_, max(frg_codes$FRG_code))
code_lu[frg_codes$FRG_code] <- frg_codes$collapsed_frg

reservoirs <- setNames(replicate(length(FRG_LEVELS), integer(0), simplify = FALSE), FRG_LEVELS)
counts <- setNames(numeric(length(FRG_LEVELS)), FRG_LEVELS)

update_reservoir <- function(frg, tm_ids) {
  m <- length(tm_ids)
  if (m == 0L) {
    return(invisible(NULL))
  }
  n_old <- counts[[frg]]
  counts[[frg]] <<- n_old + m
  if (n_old >= N_PER_FRG) {
    k_new <- rhyper(1L, m = m, n = n_old, k = N_PER_FRG)
    if (k_new > 0L) {
      take <- tm_ids[sample.int(m, k_new)]
      dest <- sample.int(N_PER_FRG, k_new)
      reservoirs[[frg]][dest] <<- take
    }
    return(invisible(NULL))
  }
  if (n_old + m <= N_PER_FRG) {
    reservoirs[[frg]] <<- c(reservoirs[[frg]], tm_ids)
    return(invisible(NULL))
  }
  need <- N_PER_FRG - n_old
  reservoirs[[frg]] <<- c(reservoirs[[frg]], tm_ids[seq_len(need)])
  rest <- tm_ids[(need + 1L):m]
  m_rest <- length(rest)
  k_new <- rhyper(1L, m = m_rest, n = N_PER_FRG, k = N_PER_FRG)
  if (k_new > 0L) {
    take <- rest[sample.int(m_rest, k_new)]
    dest <- sample.int(N_PER_FRG, k_new)
    reservoirs[[frg]][dest] <<- take
  }
}

log_pt("Streaming pixels for stratified sample, seed=", SAMPLE_SEED)
set.seed(SAMPLE_SEED)
t0 <- proc.time()[["elapsed"]]
for (i in seq_along(row_starts)) {
  r0 <- row_starts[[i]]
  nr <- min(ROW_CHUNK, n_rows - r0 + 1L)
  v_tm <- as.integer(readValues(tm, row = r0, nrows = nr))
  v_frg <- as.integer(readValues(frg_on_tm, row = r0, nrows = nr))
  ok <- !is.na(v_tm) & v_tm >= 1L & v_tm <= max_id & id_ok[pmax(v_tm, 1L)] &
    !is.na(v_frg) & v_frg >= 1L & v_frg <= length(code_lu)
  if (any(ok)) {
    frg_lab <- code_lu[v_frg[ok]]
    tm_ok <- v_tm[ok]
    keep <- !is.na(frg_lab) & frg_lab %chin% FRG_LEVELS
    if (any(keep)) {
      frg_keep <- frg_lab[keep]
      tm_keep <- tm_ok[keep]
      for (frg in FRG_LEVELS) {
        update_reservoir(frg, tm_keep[frg_keep == frg])
      }
    }
  }
  if (i == 1L || i %% 5L == 0L || i == n_chunks) {
    elapsed <- proc.time()[["elapsed"]] - t0
    eta <- elapsed / i * (n_chunks - i)
    log_pt(
      "  chunk ", i, "/", n_chunks,
      " elapsed=", round(elapsed), "s eta=", round(eta), "s",
      " rss_seen=", paste(paste0(names(counts), "=", counts), collapse = ",")
    )
  }
}

n_sampled <- vapply(reservoirs, length, integer(1L))
rows <- rbindlist(lapply(FRG_LEVELS, function(frg) {
  data.table(collapsed_frg = frg, tm_id = reservoirs[[frg]])
}))
manifest <- data.table(
  collapsed_frg = FRG_LEVELS,
  n_avail = unname(counts[FRG_LEVELS]),
  n_sampled = unname(n_sampled[FRG_LEVELS]),
  target_n_per_frg = N_PER_FRG,
  seed = SAMPLE_SEED
)
fwrite(rows, SAMPLE_CSV)
fwrite(manifest, MANIFEST_CSV)
log_pt("Wrote ", SAMPLE_CSV, " n=", nrow(rows))
print(manifest)
