#!/usr/bin/env Rscript
# 06 — Community PCA (Fig 2; basal-area CWMTs)
#
# Question: How do community-weighted mean traits separate across FRG under BA weighting?
# Why: Same stratified pixel identities as 05_pixel_sample.R; traits from tm_id CWMTs.
# Inputs:  OUT_ROOT/shared/pixel_sample.csv
#          OUT_ROOT/derived/tm_id_cwmts.csv
# Outputs: OUT_ROOT/figure_02/community_pca_*.png
# Run:     Rscript scripts/06_community_pca.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(ggplot2)

# ============================================================
# ⚙️ Settings ----
# ============================================================
PCA_SEED <- 77L
PCA_WIDTH <- as.numeric(Sys.getenv("TREE_TRAITS_PCA_WIDTH", "3.85"))
PCA_HEIGHT <- as.numeric(Sys.getenv("TREE_TRAITS_PCA_HEIGHT", "2.75"))
PCA_DPI <- as.integer(Sys.getenv("TREE_TRAITS_PCA_DPI", "300"))
TRAIT_COLS <- c(
  "avg_bark_percent", "height_m", "pruning",
  "longevity_yrs", "seed_mass_mg", "serotiny", "resprout"
)
FRG_LEVELS <- c("I", "III", "IV", "V")
FRG_COLORS <- c("I" = "red", "III" = "orange", "IV" = "yellow", "V" = "green")
TRAIT_LABELS <- c(
  avg_bark_percent = "% Bark", height_m = "Height", pruning = "Pruning",
  longevity_yrs = "Longevity", seed_mass_mg = "Seed Mass",
  serotiny = "Serotiny", resprout = "Resprout"
)
STRUCTURAL_TRAITS <- c("avg_bark_percent", "pruning", "height_m")

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
  if (!file.exists(SAMPLE_CSV)) {
    stop("Run 05_pixel_sample.R first: ", SAMPLE_CSV)
  }
  ids <- fread(SAMPLE_CSV)
  cwm <- fread(CWM_CSV)
  drop <- intersect(c("collapsed_frg", "FRG_NEW", "n_pixels", "purity"), names(cwm))
  if (length(drop)) {
    cwm[, (drop) := NULL]
  }
  merge(ids, cwm, by = "tm_id")
}

build_community_pca_plot <- function(scores, centroids, loadings, pc1_var, pc2_var, show_labels) {
  p <- ggplot() +
    stat_ellipse(
      data = scores, aes(PC1, PC2, group = collapsed_frg),
      color = "black", linewidth = 0.5, level = 0.5
    ) +
    stat_ellipse(
      data = scores, aes(PC1, PC2, color = collapsed_frg),
      linewidth = 0.4, level = 0.5
    ) +
    geom_point(
      data = centroids, aes(PC1, PC2, fill = collapsed_frg),
      shape = 21, size = 3, stroke = 0.4, color = "black"
    ) +
    geom_segment(
      data = loadings,
      aes(x = 0, y = 0, xend = PC1, yend = PC2, color = arrow_color),
      arrow = arrow(length = unit(0.12, "cm")),
      linewidth = 0.35, inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(FRG_COLORS, "#376387" = "#376387", "#b36600" = "#b36600"),
      labels = c("I" = "FRG I", "III" = "FRG III", "IV" = "FRG IV", "V" = "FRG V")
    ) +
    scale_fill_manual(values = FRG_COLORS) +
    theme_classic(base_size = 10) +
    labs(
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)"),
      color = "FRG / Trait Type", fill = "FRG"
    )
  if (show_labels && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = loadings,
      aes(x = PC1, y = PC2, label = trait_label, color = arrow_color),
      inherit.aes = FALSE, size = 2.8, show.legend = FALSE
    )
  }
  p
}

# ============================================================
# 1. 📥 Sample ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
traits_sub <- load_sample()
traits_sub <- traits_sub[complete.cases(traits_sub[, TRAIT_COLS, with = FALSE])]
log_pt("PCA n = ", nrow(traits_sub))

# ============================================================
# 2. 📊 PCA ----
# ============================================================
set.seed(PCA_SEED)
pca <- prcomp(traits_sub[, TRAIT_COLS, with = FALSE], center = TRUE, scale. = TRUE)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)
scores <- as.data.table(pca$x)
scores[, collapsed_frg := traits_sub$collapsed_frg]
centroids <- scores[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = collapsed_frg]
loadings <- as.data.table(pca$rotation, keep.rownames = "trait")
loadings[, trait_label := unname(TRAIT_LABELS[trait])]
loadings[, arrow_color := ifelse(trait %in% STRUCTURAL_TRAITS, "#376387", "#b36600")]

ggsave(
  file.path(OUTPUT_DIR, "community_pca_clean.png"),
  build_community_pca_plot(scores, centroids, loadings, pc1_var, pc2_var, FALSE),
  width = PCA_WIDTH, height = PCA_HEIGHT, dpi = PCA_DPI, units = "in"
)
ggsave(
  file.path(OUTPUT_DIR, "community_pca_labeled.png"),
  build_community_pca_plot(scores, centroids, loadings, pc1_var, pc2_var, TRUE),
  width = PCA_WIDTH, height = PCA_HEIGHT, dpi = PCA_DPI, units = "in"
)
fwrite(loadings[, .(trait, trait_label, PC1, PC2)], file.path(OUTPUT_DIR, "pca_loadings.csv"))
log_pt("Wrote PCA to ", OUTPUT_DIR)
