#!/usr/bin/env Rscript
# 12 — Supplement tables S1 (plot PERMANOVA) and S2 (pixel ANOVA)
#
# Question: Does multivariate CWMT composition differ by plot-level FRG
#           (S2/S3), and do resistance / recovery / resilience differ by
#           pixel FRG in the main-text sample (Fig 6)?
# Why: CWMTs are BA-weighted at tm_id. Plot FRG is raster majority/purity
#      (same filter as 13_supplement_frg_sensitivity.R). PERMANOVA scales
#      the seven traits so units do not dominate Euclidean distance.
# Inputs:  OUT_ROOT/derived/tm_id_cwmts.csv
#          OUT_ROOT/derived/tm_id_frg_majority.csv
#          OUT_ROOT/shared/pixel_sample.csv
# Outputs: OUT_ROOT/tables/table_S1_permanova_plot_frg.csv
#          OUT_ROOT/tables/table_S2_anova_pixel_sample.csv
# Run:     Rscript scripts/12_supplement_tables.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(vegan)

# ============================================================
# ⚙️ Settings ----
# ============================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
N_PERM <- as.integer(Sys.getenv("TREE_TRAITS_N_PERM", "199"))
PERM_SEED <- 123L
FRG_LEVELS <- c("I", "III", "IV", "V")
TRAIT_COLS <- c(
  "avg_bark_percent", "height_m", "pruning",
  "longevity_yrs", "seed_mass_mg", "serotiny", "resprout"
)
AOV_COLS <- c(
  "resistance", "recovery",
  "resilience_equivalence", "resilience_complementarity",
  "resilience_compensation_50"
)
# Same cutoffs as supplement Fig S2 / S3.
SCENARIOS <- list(
  list(scenario = "Majority 66", slug = "majority_66", min_share = 0.66),
  list(scenario = "Strict purity", slug = "strict_purity", min_share = 1.0)
)

# ============================================================
# 📁 File paths ----
# ============================================================
OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
OUTPUT_DIR <- file.path(OUT_ROOT, "tables")
SAMPLE_CSV <- file.path(OUT_ROOT, "shared", "pixel_sample.csv")
CWM_CSV <- file.path(OUT_ROOT, "derived", "tm_id_cwmts.csv")
FRG_MAJ_CSV <- file.path(OUT_ROOT, "derived", "tm_id_frg_majority.csv")

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

aov_summary_dt <- function(fit, response_label) {
  tab <- as.data.frame(summary(fit)[[1L]], stringsAsFactors = FALSE)
  tab$response <- response_label
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  as.data.table(tab)
}

# Euclidean PERMANOVA on z-scored traits (equal trait weight; same seven
# CWMTs as Fig 2 / S2). Returns the Model row plus n.
permanova_scaled <- function(dt, scenario, slug) {
  X <- scale(as.matrix(dt[, TRAIT_COLS, with = FALSE]))
  grp <- data.frame(collapsed_frg = factor(dt$collapsed_frg, levels = FRG_LEVELS))
  set.seed(PERM_SEED)
  ad <- adonis2(
    X ~ collapsed_frg,
    data = grp,
    permutations = N_PERM,
    method = "euclidean"
  )
  ad_dt <- as.data.table(as.data.frame(ad), keep.rownames = "term")
  ad_dt[, scenario := scenario]
  ad_dt[, slug := slug]
  ad_dt[, n := nrow(dt)]
  model <- ad_dt[term == "collapsed_frg" | term == "Model"]
  if (nrow(model) != 1L) {
    stop("Unexpected adonis2 terms: ", paste(ad_dt$term, collapse = ", "))
  }
  list(
    summary = data.table(
      Scenario = scenario,
      n = nrow(dt),
      R2 = model$R2,
      F = model$F,
      p = model$`Pr(>F)`
    ),
    adonis = ad_dt
  )
}

# ============================================================
# 1. 🧪 Table S1 — plot-level PERMANOVA (Fig S2 / S3) ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cwm <- fread(CWM_CSV)
drop <- intersect(c("collapsed_frg", "FRG_NEW", "n_pixels", "purity"), names(cwm))
if (length(drop)) {
  cwm[, (drop) := NULL]
}
maj <- fread(FRG_MAJ_CSV)
plots <- merge(cwm, maj, by = "tm_id")
plots <- plots[collapsed_frg %chin% FRG_LEVELS]
plots <- plots[complete.cases(plots[, TRAIT_COLS, with = FALSE])]
log_pt("plot-level candidate tm_ids = ", nrow(plots))

s1_rows <- list()
s1_adonis <- list()
for (sc in SCENARIOS) {
  sub <- plots[purity >= sc$min_share]
  log_pt("S1 ", sc$scenario, ": n = ", nrow(sub), " (min_share = ", sc$min_share, ")")
  if (nrow(sub) < 10L) {
    stop("Too few plots for ", sc$scenario)
  }
  res <- permanova_scaled(sub, sc$scenario, sc$slug)
  s1_rows[[sc$slug]] <- res$summary
  s1_adonis[[sc$slug]] <- res$adonis
}

table_s1 <- rbindlist(s1_rows)
fwrite(table_s1, file.path(OUTPUT_DIR, "table_S1_permanova_plot_frg.csv"))
fwrite(
  rbindlist(s1_adonis, fill = TRUE),
  file.path(OUTPUT_DIR, "table_S1_permanova_plot_frg_adonis.csv")
)
log_pt("Table S1")
print(table_s1)

# ============================================================
# 2. 📐 Table S2 — pixel-sample ANOVA (Fig 6) ----
# ============================================================
ids <- fread(SAMPLE_CSV)
dt <- merge(ids, cwm, by = "tm_id")
dt <- dt[complete.cases(dt[, AOV_COLS, with = FALSE])]
log_pt("Table S2 pixel n = ", nrow(dt))
fwrite(
  dt[, .(n_rows = .N), by = collapsed_frg],
  file.path(OUTPUT_DIR, "table_S2_n_rows_by_frg.csv")
)

df_aov <- as.data.frame(dt[, c("collapsed_frg", AOV_COLS), with = FALSE])
fits <- list(
  resistance = aov(resistance ~ collapsed_frg, data = df_aov),
  recovery = aov(recovery ~ collapsed_frg, data = df_aov),
  resilience_equivalence = aov(resilience_equivalence ~ collapsed_frg, data = df_aov),
  resilience_complementarity = aov(resilience_complementarity ~ collapsed_frg, data = df_aov),
  resilience_compensation_50 = aov(resilience_compensation_50 ~ collapsed_frg, data = df_aov)
)
anova_long <- rbindlist(
  mapply(function(fit, nm) aov_summary_dt(fit, nm), fits, names(fits), SIMPLIFY = FALSE),
  use.names = TRUE, fill = TRUE
)
fwrite(anova_long, file.path(OUTPUT_DIR, "table_S2_anova_pixel_sample_full.csv"))
table_s2 <- anova_long[term == "collapsed_frg"]
fwrite(table_s2, file.path(OUTPUT_DIR, "table_S2_anova_pixel_sample.csv"))
log_pt("Table S2 (FRG term)")
print(table_s2)
log_pt("Wrote tables to ", OUTPUT_DIR)
