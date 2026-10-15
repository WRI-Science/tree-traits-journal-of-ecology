#!/usr/bin/env Rscript
# 07 — CWMT pseudo-boxplots by FRG (Fig 2 companion; basal area)
#
# Question: How do native-unit CWMTs vary across FRG under BA weighting?
# Why: Same pixel sample as community PCA.
# Inputs:  OUT_ROOT/shared/pixel_sample.csv
#          OUT_ROOT/derived/tm_id_cwmts.csv
# Outputs: OUT_ROOT/figure_02/cmtv_boxplots_by_frg.png
# Run:     Rscript scripts/07_cmtv_by_frg.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(ggplot2)

# ============================================================
# ⚙️ Settings ----
# ============================================================
TRAIT_COLS <- c(
  "avg_bark_percent", "height_m", "pruning",
  "longevity_yrs", "seed_mass_mg", "serotiny", "resprout"
)
FRG_LEVELS <- c("I", "III", "IV", "V")
FRG_COLORS <- c("I" = "red", "III" = "orange", "IV" = "yellow", "V" = "green")
TRAIT_LABELS <- c(
  avg_bark_percent = "b.) Bark thickness (%)",
  height_m = "c.) Height (m)",
  pruning = "d.) Self-pruning",
  longevity_yrs = "e.) Longevity (yrs)",
  seed_mass_mg = "f.) Seed mass (mg)",
  serotiny = "g.) Serotiny",
  resprout = "h.) Resprout"
)

# ============================================================
# 📁 File paths ----
# ============================================================
OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
OUTPUT_DIR <- file.path(OUT_ROOT, "figure_02")
SAMPLE_CSV <- file.path(OUT_ROOT, "shared", "pixel_sample.csv")
CWM_CSV <- file.path(OUT_ROOT, "derived", "tm_id_cwmts.csv")

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

load_sample <- function() {
  ids <- fread(SAMPLE_CSV)
  cwm <- fread(CWM_CSV)
  drop <- intersect(c("collapsed_frg", "FRG_NEW", "n_pixels", "purity"), names(cwm))
  if (length(drop)) {
    cwm[, (drop) := NULL]
  }
  merge(ids, cwm, by = "tm_id")
}

# ============================================================
# 1. 📈 Boxplots ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dt <- load_sample()
dt <- dt[complete.cases(dt[, TRAIT_COLS, with = FALSE])]
log_pt("CMTV n = ", nrow(dt))

dt_long <- melt(
  dt[, c("collapsed_frg", TRAIT_COLS), with = FALSE],
  id.vars = "collapsed_frg", variable.name = "trait", value.name = "value"
)
dt_long[, collapsed_frg := factor(as.character(collapsed_frg), levels = FRG_LEVELS)]
dt_long[, trait := factor(
  as.character(trait), levels = TRAIT_COLS, labels = unname(TRAIT_LABELS[TRAIT_COLS])
)]
box_by_frg <- dt_long[, {
  b <- boxplot.stats(value, coef = 1.5)
  list(
    ymin = b$stats[1L], lower = b$stats[2L], middle = b$stats[3L],
    upper = b$stats[4L], ymax = b$stats[5L]
  )
}, by = .(collapsed_frg, trait)]
means_by_frg <- dt_long[, .(mean_value = mean(value, na.rm = TRUE)), by = .(collapsed_frg, trait)]
p <- ggplot(box_by_frg, aes(collapsed_frg, fill = collapsed_frg)) +
  geom_boxplot(
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", width = 0.6, color = "black"
  ) +
  geom_point(
    data = means_by_frg, aes(x = collapsed_frg, y = mean_value),
    inherit.aes = FALSE, shape = 16, size = 1.5, color = "black"
  ) +
  facet_wrap(~trait, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = FRG_COLORS) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, face = "bold", size = 11),
    legend.position = "none"
  ) +
  labs(x = "Fire Regime Group", y = NULL)
ggsave(
  file.path(OUTPUT_DIR, "cmtv_boxplots_by_frg.png"), p,
  width = 4.5, height = 7.5, dpi = 300, units = "in"
)
log_pt("Wrote CMTV boxplots")
