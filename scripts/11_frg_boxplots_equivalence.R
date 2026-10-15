#!/usr/bin/env Rscript
# 11 — Resistance, recovery, equivalence by FRG (Fig 6; basal area)
#
# Question: Do FRG boxplots of resistance / recovery / equivalence hold under BA weighting?
# Inputs:  OUT_ROOT/shared/pixel_sample.csv
#          OUT_ROOT/derived/tm_id_cwmts.csv
# Outputs: OUT_ROOT/figure_06/boxplot_*.png (1.5 in × 1 in)
#          OUT_ROOT/figure_06/figure_06_composite.png (stacked a–c)
# Run:     Rscript scripts/11_frg_boxplots_equivalence.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(ggplot2)
library(cowplot)

# ============================================================
# ⚙️ Settings ----
# ============================================================
FRG_LEVELS <- c("I", "III", "IV", "V")
FRG_COLORS <- c("I" = "red", "III" = "orange", "IV" = "yellow", "V" = "green")
# Stacked Fig 6: one panel each, 1.5 in wide × 1 in tall.
PANEL_W <- 1.5
PANEL_H <- 1
METRICS <- c("resistance", "recovery", "resilience_equivalence")
METRIC_LABELS <- c(
  resistance = "a.) Resistance",
  recovery = "b.) Recovery",
  resilience_equivalence = "c.) Resilience"
)

# ============================================================
# 📁 File paths ----
# ============================================================
OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
OUTPUT_DIR <- file.path(OUT_ROOT, "figure_06")
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

# ============================================================
# 1. 📈 Boxplots ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
ids <- fread(SAMPLE_CSV)
cwm <- fread(CWM_CSV)
drop <- intersect(c("collapsed_frg", "FRG_NEW", "n_pixels", "purity"), names(cwm))
if (length(drop)) {
  cwm[, (drop) := NULL]
}
dt <- merge(ids, cwm, by = "tm_id")
log_pt("Fig 6 n = ", nrow(dt))

panels <- vector("list", length(METRICS))
for (i in seq_along(METRICS)) {
  m <- METRICS[[i]]
  is_last <- i == length(METRICS)
  dt_m <- dt[, .(collapsed_frg, value = get(m))]
  dt_m[, collapsed_frg := factor(as.character(collapsed_frg), levels = FRG_LEVELS)]
  box <- dt_m[, {
    b <- boxplot.stats(value, coef = 1.5)
    list(
      ymin = b$stats[1L], lower = b$stats[2L], middle = b$stats[3L],
      upper = b$stats[4L], ymax = b$stats[5L]
    )
  }, by = collapsed_frg]
  means <- dt_m[, .(mean_value = mean(value, na.rm = TRUE)), by = collapsed_frg]
  p <- ggplot(box, aes(collapsed_frg, fill = collapsed_frg)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity", width = 0.6, color = "black"
    ) +
    geom_point(
      data = means, aes(x = collapsed_frg, y = mean_value),
      inherit.aes = FALSE, shape = 16, size = 1.2, color = "black"
    ) +
    scale_fill_manual(values = FRG_COLORS) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(b = 1)),
      axis.title.x = if (is_last) element_text() else element_blank(),
      plot.margin = margin(2, 4, 1, 2)
    ) +
    labs(
      x = if (is_last) "FRG" else NULL,
      y = NULL,
      title = METRIC_LABELS[[m]]
    )
  ggsave(
    file.path(OUTPUT_DIR, paste0("boxplot_", m, ".png")),
    p, width = PANEL_W, height = PANEL_H, dpi = 300, units = "in", bg = "white"
  )
  panels[[i]] <- p
}

log_pt("Building Fig 6 composite")
composite <- plot_grid(plotlist = panels, ncol = 1, align = "v")
ggsave(
  file.path(OUTPUT_DIR, "figure_06_composite.png"),
  composite,
  width = PANEL_W, height = PANEL_H * length(METRICS),
  dpi = 300, units = "in", bg = "white"
)
log_pt("Wrote Fig 6 to ", OUTPUT_DIR)
