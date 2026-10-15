#!/usr/bin/env Rscript
# 03 — Rebuild community-weighted mean traits (basal area only)
#
# Question: What are BA-weighted CWMTs and resilience indices for each TreeMap tm_id?
# Why: TreeMap is imputed FIA; basal area (TPA × 0.005454 × DIA²) is the mass-ratio
#      weight used for all community figures. CWMTs are a property of tm_id (the donor
#      plot), so we rebuild there — not the ~234 GB pixel table. Majority FRG among
#      study-region pixels labels each tm_id for Fig 2 / Fig 6 / supplements.
# Inputs:  TreeMap2016_tree_table.csv
#          traits_95_west_conus.csv
#          study-region TreeMap raster + FRG raster + LF16_FRG_200.csv
# Outputs: OUT_ROOT/derived/tm_id_cwmts.csv
#          OUT_ROOT/derived/tm_id_frg_majority.csv
#          OUT_ROOT/derived/rebuild_manifest.csv
# Run:     Rscript scripts/03_rebuild_ba_cwmts.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(terra)

# ============================================================
# ⚙️ Settings ----
# ============================================================
# FIA nested-plot expansion lives on each tree record; do not overwrite it.
STATUS_LIVE <- 1L
# FIA basal area (ft²) = 0.005454 × DIA² (DIA in inches). Times TPA_UNADJ →
# ft²/acre. Any positive constant times DIA² × TPA gives the same CWM.
BA_FT2_PER_IN2 <- 0.005454
ROW_CHUNK <- as.integer(Sys.getenv("TREE_TRAITS_BA_ROW_CHUNK", "1024"))
NICE_VALUE <- as.integer(Sys.getenv("TREE_TRAITS_BA_NICE", "10"))
TRAIT_COLS <- c(
  "avg_bark_percent", "height_m", "pruning",
  "longevity_yrs", "seed_mass_mg", "serotiny", "resprout"
)
SCALE_MINMAX_COLS <- c("avg_bark_percent", "resprout", "height_m", "serotiny", "pruning")
FRG_LEVELS <- c("I", "III", "IV", "V")

# ============================================================
# 📁 File paths ----
# ============================================================
# Repo root (publish package): parent of scripts/ when run via Rscript --file=
REPO_ROOT <- Sys.getenv("TREE_TRAITS_REPO_ROOT", "")
if (!nzchar(REPO_ROOT)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  REPO_ROOT <- if (length(file_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), winslash = "/")
  } else {
    normalizePath(getwd(), winslash = "/")
  }
}

OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
TREE_CSV <- Sys.getenv("TREE_TRAITS_TREEMAP_TREE_CSV", "")
if (!nzchar(TREE_CSV)) {
  stop(
    "Set TREE_TRAITS_TREEMAP_TREE_CSV to the path of: TreeMap 2016 tree table CSV (Riley et al. 2021; doi:10.2737/RDS-2021-0074)"
  )
}
SPECIES_CSV <- Sys.getenv("TREE_TRAITS_SPECIES_CSV", file.path(REPO_ROOT, "data", "traits_95_west_conus.csv"))
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
OUTPUT_DIR <- file.path(OUT_ROOT, "derived")

# ============================================================
# 🔧 Helpers ----
# ============================================================
log_pt <- function(...) {
  message(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/Los_Angeles"),
    " ",
    paste0(...)
  )
}

minmax01 <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (!is.finite(rng) || rng == 0) {
    return(rep(0, length(x)))
  }
  (x - min(x, na.rm = TRUE)) / rng
}

weighted_cwm <- function(trees, wcol, trait_cols) {
  w <- trees[[wcol]]
  keep <- is.finite(w) & w > 0
  dt <- trees[keep]
  if (nrow(dt) == 0L) {
    return(dt[0])
  }
  dt[, w := get(wcol)]
  out <- dt[, c(
    list(n_trees_used = .N, weight_sum = sum(w)),
    lapply(.SD, function(x) sum(w * x) / sum(w))
  ), by = tm_id, .SDcols = trait_cols]
  out
}

scale_and_resilience <- function(cwm) {
  dt <- copy(cwm)
  dt[, (SCALE_MINMAX_COLS) := lapply(.SD, minmax01), .SDcols = SCALE_MINMAX_COLS]
  # Match published community table: log longevity, invert (long-lived → low recovery).
  dt[, longevity_yrs := {
    x <- log(longevity_yrs)
    1 - minmax01(x)
  }]
  # Log seed mass, min-max, not inverted.
  dt[, seed_mass_mg := {
    x <- log(seed_mass_mg)
    minmax01(x)
  }]
  dt[, resistance := (avg_bark_percent + height_m + pruning) / 3]
  dt[, recovery := (serotiny + resprout + longevity_yrs + seed_mass_mg) / 4]
  dt[, resilience_equivalence := 1 - ((1 - resistance) * (1 - recovery))]
  dt[, resilience_complementarity := resistance * recovery]
  dt[, resilience_compensation_50 := resistance + 0.5 * (1 - resistance) * recovery]
  dt
}

# ============================================================
# 1. 🧹 Session ----
# ============================================================
if (requireNamespace("tools", quietly = TRUE)) {
  try(tools::psnice(value = NICE_VALUE), silent = TRUE)
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

need <- c(TREE_CSV, SPECIES_CSV, TM_RAST, FRG_RAST, FRG_LOOKUP_CSV)
missing <- need[!file.exists(need)]
if (length(missing)) {
  stop("Missing inputs:\n", paste(missing, collapse = "\n"))
}

# ============================================================
# 2. 🌲 Species traits ----
# ============================================================
traits <- fread(SPECIES_CSV)
need_trait <- c("FIA_CODE", TRAIT_COLS)
miss_tr <- setdiff(need_trait, names(traits))
if (length(miss_tr)) {
  stop("Species CSV missing columns: ", paste(miss_tr, collapse = ", "))
}
traits[, FIA_CODE := as.integer(FIA_CODE)]
traits <- traits[, c("FIA_CODE", TRAIT_COLS), with = FALSE]
if (anyDuplicated(traits$FIA_CODE)) {
  stop("Duplicate FIA_CODE in species traits.")
}
log_pt("Species traits n = ", nrow(traits))

# ============================================================
# 3. 📋 Live TreeMap trees + basal-area weight ----
# ============================================================
log_pt("Reading tree table: ", TREE_CSV)
trees <- fread(
  TREE_CSV,
  select = c("tm_id", "STATUSCD", "TPA_UNADJ", "SPCD", "DIA")
)
trees <- trees[STATUSCD == STATUS_LIVE]
n_live <- nrow(trees)
trees <- trees[SPCD %in% traits$FIA_CODE]
n_matched <- nrow(trees)
trees <- merge(trees, traits, by.x = "SPCD", by.y = "FIA_CODE", all.x = FALSE)
# BA weight is NA when DIA is missing/non-positive so that plot is dropped from BA CWM, not silently zeroed.
trees[, w_ba := fifelse(
  is.finite(DIA) & DIA > 0 & is.finite(TPA_UNADJ) & TPA_UNADJ > 0,
  TPA_UNADJ * BA_FT2_PER_IN2 * (DIA^2),
  NA_real_
)]
log_pt(
  "Live trees ", n_live, "; matched to 95-species list ", n_matched,
  " (", round(100 * n_matched / n_live, 1), "%)"
)

cwm_ba <- weighted_cwm(trees, "w_ba", TRAIT_COLS)
log_pt("CWM tm_ids (BA): ", nrow(cwm_ba))

# ============================================================
# 4. 🗺️ Majority FRG per tm_id on the study-region raster ----
# ============================================================
# Raster values are sequential FRG_code from distinct FRG_NEW (same as the
# published extract→merge), not raw LANDFIRE BPS codes.
lookup <- fread(FRG_LOOKUP_CSV)
lookup_clean <- lookup[
  !is.na(FRG_NEW) &
    FRG_NEW != "-9999" & FRG_NEW != "-1111" &
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

log_pt("Opening rasters (file-backed). Aligning FRG onto TreeMap grid via VRT.")
tm <- rast(TM_RAST)
frg <- rast(FRG_RAST)
NAflag(frg) <- 65535

vrt_path <- file.path(OUTPUT_DIR, "tmp_frg_on_tm.vrt")
if (file.exists(vrt_path)) {
  unlink(vrt_path)
}
te <- as.vector(ext(tm))
# terra ext vector is xmin, xmax, ymin, ymax; gdalwarp -te is xmin ymin xmax ymax.
warp_status <- system2(
  "gdalwarp",
  c(
    "-r", "near",
    "-of", "VRT",
    "-te", te[1], te[3], te[2], te[4],
    "-tr", "30", "30",
    "-dstnodata", "65535",
    FRG_RAST,
    vrt_path
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!identical(attr(warp_status, "status"), NULL) && !identical(attr(warp_status, "status"), 0L)) {
  stop("gdalwarp VRT failed:\n", paste(warp_status, collapse = "\n"))
}
if (!file.exists(vrt_path)) {
  stop("gdalwarp VRT failed: ", vrt_path)
}
frg_on_tm <- rast(vrt_path)
if (nrow(frg_on_tm) != nrow(tm) || ncol(frg_on_tm) != ncol(tm)) {
  stop(
    "Aligned FRG grid mismatch: ", ncol(frg_on_tm), "x", nrow(frg_on_tm),
    " vs TreeMap ", ncol(tm), "x", nrow(tm)
  )
}

n_rows <- nrow(tm)
row_starts <- seq(1L, n_rows, by = ROW_CHUNK)
n_chunks <- length(row_starts)
log_pt("FRG majority: ", n_chunks, " row chunks of ", ROW_CHUNK, " (", n_rows, " rows)")

readStart(tm)
readStart(frg_on_tm)
on.exit({
  try(readStop(tm), silent = TRUE)
  try(readStop(frg_on_tm), silent = TRUE)
}, add = TRUE)

counts <- data.table(tm_id = integer(), frg = integer(), N = integer())
t0 <- proc.time()[["elapsed"]]
for (i in seq_along(row_starts)) {
  r0 <- row_starts[[i]]
  nr <- min(ROW_CHUNK, n_rows - r0 + 1L)
  v_tm <- as.integer(readValues(tm, row = r0, nrows = nr))
  v_frg <- as.integer(readValues(frg_on_tm, row = r0, nrows = nr))
  ok_pix <- !is.na(v_tm) & v_tm > 0L &
    !is.na(v_frg) & v_frg > 0L & v_frg < 65535L
  if (!any(ok_pix)) {
    if (i == 1L || i %% 20L == 0L || i == n_chunks) {
      log_pt("  chunk ", i, "/", n_chunks, " (empty)")
    }
    next
  }
  stripe <- data.table(tm_id = v_tm[ok_pix], frg = v_frg[ok_pix])[, .N, by = .(tm_id, frg)]
  counts <- rbindlist(list(counts, stripe), use.names = TRUE)
  counts <- counts[, .(N = sum(N)), by = .(tm_id, frg)]
  if (i == 1L || i %% 20L == 0L || i == n_chunks) {
    elapsed <- proc.time()[["elapsed"]] - t0
    eta <- elapsed / i * (n_chunks - i)
    log_pt(
      "  chunk ", i, "/", n_chunks,
      " unique tm_id so far=", uniqueN(counts$tm_id),
      " elapsed=", round(elapsed), "s eta=", round(eta), "s"
    )
  }
}

frg_maj <- counts[counts[, .I[which.max(N)], by = tm_id]$V1]
tot <- counts[, .(n_pixels = sum(N)), by = tm_id]
frg_maj <- merge(frg_maj, tot, by = "tm_id")
frg_maj[, purity := N / n_pixels]
frg_maj <- merge(frg_maj, frg_codes, by.x = "frg", by.y = "FRG_code", all.x = TRUE)
setnames(frg_maj, c("frg", "N"), c("FRG_code", "n_pixels_majority"))
fwrite(frg_maj, file.path(OUTPUT_DIR, "tm_id_frg_majority.csv"))
log_pt(
  "Study-region tm_ids with FRG: ", nrow(frg_maj),
  "; with collapsed I/III/IV/V: ",
  frg_maj[collapsed_frg %chin% FRG_LEVELS, .N]
)

# ============================================================
# 5. 🔗 Join, scale, write BA-only CWMTs ----
# ============================================================
keep_ids <- intersect(cwm_ba$tm_id, frg_maj$tm_id)
log_pt("Study tm_ids with BA CWM + FRG: ", length(keep_ids))
if (length(keep_ids) < 10L) {
  stop("Too few tm_ids after join (n = ", length(keep_ids), ").")
}

raw <- cwm_ba[tm_id %in% keep_ids]
scaled <- scale_and_resilience(raw[, c("tm_id", TRAIT_COLS), with = FALSE])
out <- merge(
  raw[, c("tm_id", "n_trees_used", "weight_sum", TRAIT_COLS), with = FALSE],
  scaled[, c(
    "tm_id", "resistance", "recovery",
    "resilience_equivalence", "resilience_complementarity",
    "resilience_compensation_50"
  ), with = FALSE],
  by = "tm_id"
)
out <- merge(
  out,
  frg_maj[, .(tm_id, collapsed_frg, FRG_NEW, n_pixels, purity)],
  by = "tm_id",
  all.x = TRUE
)

fwrite(out, file.path(OUTPUT_DIR, "tm_id_cwmts.csv"))

manifest <- data.table(
  abundance_weight = "ba",
  n_live_trees = n_live,
  n_live_matched_species = n_matched,
  n_cwm_ba = nrow(cwm_ba),
  n_study_tm_id_with_frg = nrow(frg_maj),
  n_tm_id_written = nrow(out),
  n_frg_maintext = out[collapsed_frg %chin% FRG_LEVELS, .N]
)
fwrite(manifest, file.path(OUTPUT_DIR, "rebuild_manifest.csv"))
log_pt("Wrote ", file.path(OUTPUT_DIR, "tm_id_cwmts.csv"), " rows=", nrow(out))
log_pt("Done 03.")
