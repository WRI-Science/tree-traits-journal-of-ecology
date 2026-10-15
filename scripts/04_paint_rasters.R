#!/usr/bin/env Rscript
# 04 — Paint TreeMap pixels with BA-weighted tm_id CWMTs
#
# Question: What do wall-to-wall CWMT / resistance rasters look like under basal area?
# Why: CWMTs live on tm_id. Painting the study-region TreeMap raster avoids rebuilding
#      the ~234 GB pixel CSV. One pass writes all map layers.
# Inputs:  OUT_ROOT/derived/tm_id_cwmts.csv
#          study-region TreeMap raster
# Outputs: OUT_ROOT/rasters/*.tif
# Run:     Rscript scripts/04_paint_rasters.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(terra)

# ============================================================
# ⚙️ Settings ----
# ============================================================
ROW_CHUNK <- as.integer(Sys.getenv("TREE_TRAITS_BA_ROW_CHUNK", "1024"))
NICE_VALUE <- as.integer(Sys.getenv("TREE_TRAITS_BA_NICE", "10"))
PAINT_COLS <- c(
  "avg_bark_percent", "height_m", "pruning",
  "longevity_yrs", "seed_mass_mg", "serotiny", "resprout",
  "resistance", "recovery",
  "resilience_equivalence", "resilience_complementarity",
  "resilience_compensation_50"
)

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
RASTER_DIR <- file.path(OUT_ROOT, "rasters")

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
if (!file.exists(CWM_CSV)) {
  stop("Missing tm_id CWMTs. Run 03_rebuild_ba_cwmts.R first: ", CWM_CSV)
}
if (!file.exists(TM_RAST)) {
  stop("Missing TreeMap raster: ", TM_RAST)
}

dir.create(RASTER_DIR, recursive = TRUE, showWarnings = FALSE)
fnames_check <- file.path(RASTER_DIR, paste0(PAINT_COLS, ".tif"))
force_paint <- identical(Sys.getenv("TREE_TRAITS_FORCE_PAINT", ""), "1")
if (!force_paint && all(file.exists(fnames_check)) && all(file.info(fnames_check)$size > 1e6)) {
  log_pt("Existing painted rasters: ", RASTER_DIR)
  quit(save = "no", status = 0L)
}

# ============================================================
# 2. 📋 Lookup tables (integer index, not named character) ----
# ============================================================
cwm <- fread(CWM_CSV, select = c("tm_id", PAINT_COLS))
if (nrow(cwm) == 0L) {
  stop("No rows in ", CWM_CSV)
}
max_id <- max(cwm$tm_id)
lookups <- lapply(PAINT_COLS, function(col) {
  lu <- rep(NA_real_, max_id)
  lu[cwm$tm_id] <- cwm[[col]]
  lu
})
names(lookups) <- PAINT_COLS
log_pt("BA paint: lookup tm_ids = ", nrow(cwm), " max tm_id = ", max_id)

# ============================================================
# 3. 🖌️ Paint all layers in one read pass ----
# ============================================================
tm <- rast(TM_RAST)

outs <- lapply(PAINT_COLS, function(col) {
  rast(tm)
})
names(outs) <- PAINT_COLS
fnames <- file.path(RASTER_DIR, paste0(PAINT_COLS, ".tif"))
names(fnames) <- PAINT_COLS

readStart(tm)
on.exit(try(readStop(tm), silent = TRUE), add = TRUE)
# In terra 1.9, writeStart modifies the SpatRaster by side-effect and
# returns a block plan. Empty template rasters often get n=1 (the whole
# grid). Ignore that plan and write ROW_CHUNK rows at a time.
n_rows <- nrow(tm)
row_starts <- seq(1L, n_rows, by = ROW_CHUNK)
n_chunks <- length(row_starts)
for (col in PAINT_COLS) {
  writeStart(
    outs[[col]],
    filename = fnames[[col]],
    overwrite = TRUE,
    n = 32,
    datatype = "FLT4S",
    gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
  )
}
on.exit({
  for (col in PAINT_COLS) {
    try(writeStop(outs[[col]]), silent = TRUE)
  }
}, add = TRUE)

log_pt("Painting ", length(PAINT_COLS), " layers, ", n_chunks, " chunks of ", ROW_CHUNK, " rows")
t0 <- proc.time()[["elapsed"]]
for (i in seq_along(row_starts)) {
  r0 <- row_starts[[i]]
  nr <- min(ROW_CHUNK, n_rows - r0 + 1L)
  v_tm <- as.integer(readValues(tm, row = r0, nrows = nr))
  ok <- !is.na(v_tm) & v_tm >= 1L & v_tm <= max_id
  idx <- integer(length(v_tm))
  idx[ok] <- v_tm[ok]
  for (col in PAINT_COLS) {
    v_out <- rep(NA_real_, length(v_tm))
    if (any(ok)) {
      v_out[ok] <- lookups[[col]][idx[ok]]
    }
    writeValues(outs[[col]], v_out, start = r0, nrows = nr)
  }
  if (i == 1L || i %% 5L == 0L || i == n_chunks) {
    elapsed <- proc.time()[["elapsed"]] - t0
    eta <- elapsed / i * (n_chunks - i)
    log_pt("  chunk ", i, "/", n_chunks, " elapsed=", round(elapsed), "s eta=", round(eta), "s")
  }
}

for (col in PAINT_COLS) {
  writeStop(outs[[col]])
}
try(readStop(tm), silent = TRUE)
# Drop on.exit handlers so a second writeStop cannot error after success.
on.exit(NULL)
log_pt("Wrote rasters under ", RASTER_DIR)
