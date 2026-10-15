#!/usr/bin/env Rscript
# 00 — Run the basal-area manuscript pipeline end to end
#
# Question: Reproduce all main-text and supplement figures/tables under BA weighting.
# Why: Reviewer-facing single pipeline; no count/TPA abundance-weight toggle.
# Inputs:  data/traits_95_west_conus.csv + TreeMap/FRG paths via env vars (see README)
# Outputs: TREE_TRAITS_OUT_ROOT (default: <repo>/output)
# Run:     Rscript scripts/00_run_all.R
#          (from repository root)

# ============================================================
# 📁 File paths ----
# ============================================================
repo_root <- Sys.getenv("TREE_TRAITS_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  repo_root <- if (length(file_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), winslash = "/")
  } else {
    normalizePath(getwd(), winslash = "/")
  }
}

out_root <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(repo_root, "output"))

scripts <- c(
  "01_individual_trait_space.R",
  "02_supplement_individual_traits_by_group.R",
  "03_rebuild_ba_cwmts.R",
  "04_paint_rasters.R",
  "05_pixel_sample.R",
  "06_community_pca.R",
  "07_cmtv_by_frg.R",
  "08_cwmt_trait_maps.R",
  "09_resistance_recovery_bivariate_maps.R",
  "10_resilience_formulation_maps.R",
  "11_frg_boxplots_equivalence.R",
  "12_supplement_tables.R",
  "13_supplement_frg_sensitivity.R"
)

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

run_script <- function(scr) {
  fp <- file.path(repo_root, "scripts", scr)
  if (!file.exists(fp)) {
    stop("Missing script: ", fp)
  }
  log_pt("========== ", scr, " ==========")
  status <- system2("Rscript", fp, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Failed: ", scr, " (exit ", status, ")")
  }
}

# ============================================================
# 1. ▶️ Run 01–13 ----
# ============================================================
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TREE_TRAITS_OUT_ROOT = out_root)
Sys.setenv(TREE_TRAITS_REPO_ROOT = repo_root)
log_pt("repo root: ", repo_root)
log_pt("OUT_ROOT: ", out_root)
for (scr in scripts) {
  run_script(scr)
}
log_pt("BA manuscript pipeline complete.")
