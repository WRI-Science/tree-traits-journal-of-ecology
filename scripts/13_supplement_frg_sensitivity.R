#!/usr/bin/env Rscript
# 13 — Supplement S2/S3: FRG assignment sensitivity (majority_66 / strict_purity)
#
# Question: Are PCA / CMTV patterns robust to majority_66 vs strict_purity
#           when communities are BA-weighted?
# Why: Pixel FRG still varies for a tm_id. Use raster majority/purity from
#      03_rebuild_ba_cwmts.R (no 234 GB CSV). This is FRG-assignment sensitivity
#      for the paper — not abundance-weight (count/TPA) sensitivity.
# Inputs:  OUT_ROOT/derived/tm_id_frg_majority.csv
#          OUT_ROOT/derived/tm_id_cwmts.csv
# Outputs: OUT_ROOT/supplement_S2/{majority_66,strict_purity}/
# Run:     Rscript scripts/13_supplement_frg_sensitivity.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(data.table)
library(ggplot2)

# ============================================================
# ⚙️ Settings ----
# ============================================================
PCA_WIDTH <- as.numeric(Sys.getenv("TREE_TRAITS_PCA_WIDTH", "3.85"))
PCA_HEIGHT <- as.numeric(Sys.getenv("TREE_TRAITS_PCA_HEIGHT", "2.75"))
PCA_DPI <- as.integer(Sys.getenv("TREE_TRAITS_PCA_DPI", "300"))
CMTV_WIDTH <- as.numeric(Sys.getenv("TREE_TRAITS_CMTV_WIDTH", "4.5"))
CMTV_HEIGHT <- as.numeric(Sys.getenv("TREE_TRAITS_CMTV_HEIGHT", "7.5"))
CMTV_DPI <- as.integer(Sys.getenv("TREE_TRAITS_CMTV_DPI", "300"))
PCA_SEED <- 77L
SCENARIOS <- list(
  list(name = "majority_66", min_share = 0.66),
  list(name = "strict_purity", min_share = 1.0)
)
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
OUTPUT_ROOT <- file.path(OUT_ROOT, "supplement_S2")
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

run_pca_panels <- function(traits_sub, out_dir, scenario_name) {
  set.seed(PCA_SEED)
  pca <- prcomp(traits_sub[, TRAIT_COLS, with = FALSE], center = TRUE, scale. = TRUE)
  scores <- as.data.table(pca$x)
  scores[, collapsed_frg := traits_sub$collapsed_frg]
  centroids <- scores[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = collapsed_frg]
  loadings <- as.data.table(pca$rotation, keep.rownames = "trait")
  loadings[, arrow_color := ifelse(
    trait %in% c("avg_bark_percent", "pruning", "height_m"), "#376387", "#b36600"
  )]
  pc1 <- round(summary(pca)$importance[2, 1] * 100, 1)
  pc2 <- round(summary(pca)$importance[2, 2] * 100, 1)
  build_plot <- function(show_labels) {
    p <- ggplot() +
      stat_ellipse(data = scores, aes(PC1, PC2, group = collapsed_frg),
                   color = "black", linewidth = 0.5, level = 0.5) +
      stat_ellipse(data = scores, aes(PC1, PC2, color = collapsed_frg),
                   linewidth = 0.4, level = 0.5) +
      geom_point(data = centroids, aes(PC1, PC2, fill = collapsed_frg),
                 shape = 21, size = 3, stroke = 0.4, color = "black") +
      geom_segment(
        data = loadings,
        aes(x = 0, y = 0, xend = PC1, yend = PC2, color = arrow_color),
        arrow = arrow(length = unit(0.12, "cm")),
        linewidth = 0.35, inherit.aes = FALSE
      ) +
      scale_color_manual(values = c(FRG_COLORS, "#376387" = "#376387", "#b36600" = "#b36600")) +
      scale_fill_manual(values = FRG_COLORS) +
      theme_classic(base_size = 10) +
      labs(x = paste0("PC1 (", pc1, "%)"), y = paste0("PC2 (", pc2, "%)"),
           title = scenario_name)
    if (show_labels && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = loadings, aes(x = PC1, y = PC2, label = trait, color = arrow_color),
        inherit.aes = FALSE, size = 2.8, show.legend = FALSE
      )
    }
    p
  }
  ggsave(file.path(out_dir, "community_pca_clean.png"), build_plot(FALSE),
         width = PCA_WIDTH, height = PCA_HEIGHT, dpi = PCA_DPI, units = "in")
  ggsave(file.path(out_dir, "community_pca_labeled.png"), build_plot(TRUE),
         width = PCA_WIDTH, height = PCA_HEIGHT, dpi = PCA_DPI, units = "in")
  fwrite(as.data.table(pca$rotation, keep.rownames = "trait"),
         file.path(out_dir, "pca_loadings.csv"))
}

run_cmtv_panel <- function(traits_sub, out_dir) {
  dt_long <- melt(
    traits_sub[, c("collapsed_frg", TRAIT_COLS), with = FALSE],
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
    file.path(out_dir, "cmtv_boxplots_by_frg.png"), p,
    width = CMTV_WIDTH, height = CMTV_HEIGHT, dpi = CMTV_DPI, units = "in"
  )
}

# ============================================================
# 1. 📥 tm_id CWMTs + raster purity ----
# ============================================================
cwm <- fread(CWM_CSV)
drop <- intersect(c("collapsed_frg", "FRG_NEW", "n_pixels", "purity"), names(cwm))
if (length(drop)) {
  cwm[, (drop) := NULL]
}
maj <- fread(FRG_MAJ_CSV)
dt <- merge(cwm, maj, by = "tm_id")
dt <- dt[collapsed_frg %chin% FRG_LEVELS]
dt <- dt[complete.cases(dt[, TRAIT_COLS, with = FALSE])]
log_pt("S2 candidate tm_ids = ", nrow(dt))

# ============================================================
# 2. 🔬 Scenarios ----
# ============================================================
for (sc in SCENARIOS) {
  out_dir <- file.path(OUTPUT_ROOT, sc$name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  traits_sub <- dt[purity >= sc$min_share]
  log_pt("S2 ", sc$name, ": n = ", nrow(traits_sub), " (min_share = ", sc$min_share, ")")
  if (nrow(traits_sub) < 10L) {
    stop("Too few plots for ", sc$name)
  }
  run_pca_panels(traits_sub, out_dir, sc$name)
  run_cmtv_panel(traits_sub, out_dir)
}
log_pt("Done. ", OUTPUT_ROOT)
