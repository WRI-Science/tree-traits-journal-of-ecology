#!/usr/bin/env Rscript
# 01 — Individual trait space and emergent functional groups (Fig 1)
#
# Question: How do species cluster in trait space, and what are mean trait values per group?
# Why: Species-level ordination defines emergent functional groups (PAM on Gower distance).
#      k is chosen from silhouette (max mean width over k = 2..K_MAX) and gap (firstSEmax);
#      we use the silhouette optimum for the final partition. Gap tables are written for
#      the process figure. Override with TREE_TRAITS_N_EMERGENT_GROUPS if needed.
# Inputs:  data/traits_95_west_conus.csv (or TREE_TRAITS_SPECIES_CSV)
# Outputs: OUT_ROOT/figure_01/species_functional_groups.csv
#          OUT_ROOT/figure_01/individual_trait_pca_clean.png
#          OUT_ROOT/figure_01/individual_trait_pca_labeled.png  (reference only)
#          OUT_ROOT/figure_01/group_mean_traits_bars.png
#          OUT_ROOT/figure_01/genus_by_cluster.csv
#          OUT_ROOT/figure_01/genus_cluster_permutation.csv
# Run:     Rscript scripts/01_individual_trait_space.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(cluster)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# ⚙️ Settings ----
# ============================================================
K_MIN <- 2L
K_MAX_SILHOUETTE <- as.integer(Sys.getenv("TREE_TRAITS_K_MAX_SILHOUETTE", "25"))
K_MAX_GAP <- as.integer(Sys.getenv("TREE_TRAITS_K_MAX_GAP", "10"))
GAP_B <- as.integer(Sys.getenv("TREE_TRAITS_GAP_B", "500"))
# Seed for clusGap bootstrap reference distributions — fixed so gap tables are reproducible
# across runs; 123 matches other pipeline RNG defaults.
GAP_SEED <- as.integer(Sys.getenv("TREE_TRAITS_GAP_SEED", "123"))
K_OVERRIDE <- Sys.getenv("TREE_TRAITS_N_EMERGENT_GROUPS", "")
PCA_WIDTH <- as.numeric(Sys.getenv("TREE_TRAITS_PCA_WIDTH", "3.85"))
PCA_HEIGHT <- as.numeric(Sys.getenv("TREE_TRAITS_PCA_HEIGHT", "2.75"))
PCA_DPI <- as.integer(Sys.getenv("TREE_TRAITS_PCA_DPI", "300"))
# Genus × cluster association: permute cluster labels with genus fixed (reviewer ask).
GENUS_PERM_N <- as.integer(Sys.getenv("TREE_TRAITS_GENUS_PERM_N", "9999"))
GENUS_PERM_SEED <- as.integer(Sys.getenv("TREE_TRAITS_GENUS_PERM_SEED", "123"))
GROUP_COLOR_PALETTE <- c("#3a5311", "#728c69", "#028a04", "#74b72e", "#b2d3c2")

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
SPECIES_CSV <- Sys.getenv("TREE_TRAITS_SPECIES_CSV", file.path(REPO_ROOT, "data", "traits_95_west_conus.csv"))
OUTPUT_DIR <- file.path(OUT_ROOT, "figure_01")

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

select_emergent_k <- function(gower_dist, k_min, k_max_sil, k_max_gap, gap_b, gap_seed) {
  n <- attr(gower_dist, "Size")
  k_hi <- min(k_max_sil, n - 1L)
  if (k_hi < k_min) {
    stop("Too few species (n = ", n, ") to cluster with k >= ", k_min, ".")
  }

  sil_tbl <- data.frame(k = k_min:k_hi, avg_silhouette = NA_real_)
  for (i in seq_len(nrow(sil_tbl))) {
    fit <- cluster::pam(gower_dist, diss = TRUE, k = sil_tbl$k[i])
    sil_tbl$avg_silhouette[i] <- fit$silinfo$avg.width
  }
  k_sil <- sil_tbl$k[which.max(sil_tbl$avg_silhouette)]

  k_gap_hi <- min(k_max_gap, n - 1L)
  # Gap statistic draws bootstrap nulls — seed so diagnostics match across machines.
  set.seed(gap_seed)
  gap_stat <- cluster::clusGap(
    as.matrix(gower_dist),
    FUN = function(x, k) list(cluster = cluster::pam(x, k, diss = TRUE)$clustering),
    K.max = k_gap_hi,
    B = gap_b,
    verbose = FALSE
  )
  k_gap <- cluster::maxSE(
    gap_stat$Tab[, "gap"],
    gap_stat$Tab[, "SE.sim"],
    method = "firstSEmax"
  )
  gap_tbl <- data.frame(
    k = seq_len(k_gap_hi),
    gap = gap_stat$Tab[, "gap"],
    SE = gap_stat$Tab[, "SE.sim"]
  )

  list(
    k_silhouette = k_sil,
    k_gap = k_gap,
    silhouette_table = sil_tbl,
    gap_table = gap_tbl
  )
}

# ============================================================
# 1. 📥 Load species traits ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SPECIES_CSV)) {
  stop("Species traits CSV not found: ", SPECIES_CSV)
}

traits_raw <- read.csv(SPECIES_CSV, stringsAsFactors = FALSE)
species_names <- traits_raw$SCIENTIFIC_NAME
traits <- traits_raw[, -(1:2)]

# ============================================================
# 2. 🔄 Transform traits for Gower / PCA ----
# ============================================================
# log10 compresses right-skewed longevity and seed mass so they don't dominate Gower.
traits$longevity_yrs <- log10(traits$longevity_yrs)
traits$seed_mass_mg <- log10(traits$seed_mass_mg)

# Min–max rescale to [0, 1] so mixed binary/continuous traits share a common range for Gower.
traits_rescaled <- traits %>%
  mutate(across(everything(), ~ (. - min(.)) / (max(.) - min(.))))
# Longevity flipped so long lifespans align with low recovery in downstream maps.
traits_rescaled$longevity_yrs <- 1 - traits_rescaled$longevity_yrs

# ============================================================
# 3. 🔍 Choose k and assign PAM groups ----
# ============================================================
gower_dist <- cluster::daisy(traits_rescaled, metric = "gower")
k_diag <- select_emergent_k(
  gower_dist, K_MIN, K_MAX_SILHOUETTE, K_MAX_GAP, GAP_B, GAP_SEED
)
N_EMERGENT_GROUPS <- if (nzchar(K_OVERRIDE)) {
  as.integer(K_OVERRIDE)
} else {
  k_diag$k_silhouette
}
if (is.na(N_EMERGENT_GROUPS) || N_EMERGENT_GROUPS < K_MIN) {
  stop("Invalid k = ", N_EMERGENT_GROUPS, ". Set TREE_TRAITS_N_EMERGENT_GROUPS >= ", K_MIN, ".")
}
if (N_EMERGENT_GROUPS > length(GROUP_COLOR_PALETTE)) {
  stop(
    "k = ", N_EMERGENT_GROUPS, " exceeds palette length (",
    length(GROUP_COLOR_PALETTE), "). Extend GROUP_COLOR_PALETTE."
  )
}

GROUP_LABELS <- LETTERS[seq_len(N_EMERGENT_GROUPS)]
GROUP_COLORS <- GROUP_COLOR_PALETTE[seq_len(N_EMERGENT_GROUPS)]

write.csv(k_diag$silhouette_table, file.path(OUTPUT_DIR, "silhouette_by_k.csv"), row.names = FALSE)
write.csv(k_diag$gap_table, file.path(OUTPUT_DIR, "gap_statistic_by_k.csv"), row.names = FALSE)
log_pt(
  "Cluster diagnostics: k_silhouette = ", k_diag$k_silhouette,
  ", k_gap_firstSEmax = ", k_diag$k_gap,
  "; using k = ", N_EMERGENT_GROUPS,
  if (nzchar(K_OVERRIDE)) " (TREE_TRAITS_N_EMERGENT_GROUPS override)" else " (silhouette optimum)"
)

pam_fit <- cluster::pam(gower_dist, diss = TRUE, k = N_EMERGENT_GROUPS)
group_ids <- as.integer(pam_fit$clustering)

groups_out <- data.frame(
  scientific_name = species_names,
  emergent_group = group_ids,
  emergent_group_label = GROUP_LABELS[group_ids],
  stringsAsFactors = FALSE
)
write.csv(
  groups_out,
  file.path(OUTPUT_DIR, "species_functional_groups.csv"),
  row.names = FALSE
)
log_pt(
  "Wrote species assignments (k = ", N_EMERGENT_GROUPS, "): ",
  file.path(OUTPUT_DIR, "species_functional_groups.csv")
)

# ============================================================
# 4. 📊 Species PCA panels ----
# ============================================================
# Traits already [0,1]; center only — scale.=FALSE avoids re-weighting after min–max.
pca_result <- prcomp(traits_rescaled, center = TRUE, scale. = FALSE)
pca_scores <- pca_result$x
explained_variance <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
pc1_var <- round(explained_variance[1], 1)
pc2_var <- round(explained_variance[2], 1)

pca_df <- data.frame(
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2],
  group = factor(GROUP_LABELS[group_ids], levels = GROUP_LABELS)
)

loadings_df <- data.frame(
  trait = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[, 1],
  PC2 = pca_result$rotation[, 2],
  stringsAsFactors = FALSE
)
trait_label_map <- c(
  avg_bark_percent = "% Bark",
  longevity_yrs = "Longevity",
  resprout = "Resprout",
  seed_mass_mg = "Seed Mass",
  height_m = "Height",
  serotiny = "Serotiny",
  pruning = "Pruning"
)
loadings_df$trait_label <- unname(trait_label_map[loadings_df$trait])
loadings_df$arrow_color <- ifelse(
  loadings_df$trait_label %in% c("% Bark", "Pruning", "Height"),
  "#376387", "#b36600"
)

build_pca_plot <- function(show_labels) {
  p <- ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = group), size = 1.5) +
    stat_ellipse(aes(color = group), type = "t", linewidth = 0.5, level = 0.5) +
    geom_segment(
      data = loadings_df,
      aes(x = 0, y = 0, xend = PC1, yend = PC2, color = arrow_color),
      arrow = arrow(length = unit(0.12, "cm")),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    scale_color_manual(
      values = c(
        setNames(GROUP_COLORS, GROUP_LABELS),
        "#376387" = "#376387",
        "#b36600" = "#b36600"
      )
    ) +
    theme_classic(base_size = 10) +
    labs(
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    )
  if (show_labels) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("Install ggrepel for labeled PCA, or use clean output only.")
    }
    p <- p + ggrepel::geom_text_repel(
      data = loadings_df,
      aes(x = PC1, y = PC2, label = trait_label, color = arrow_color),
      inherit.aes = FALSE,
      size = 2.8,
      show.legend = FALSE
    ) +
      labs(title = "Individual trait PCA (reference labels)")
  }
  p
}

ggsave(
  file.path(OUTPUT_DIR, "individual_trait_pca_clean.png"),
  build_pca_plot(FALSE),
  width = PCA_WIDTH, height = PCA_HEIGHT, dpi = PCA_DPI, units = "in"
)
ggsave(
  file.path(OUTPUT_DIR, "individual_trait_pca_labeled.png"),
  build_pca_plot(TRUE),
  width = PCA_WIDTH, height = PCA_HEIGHT, dpi = PCA_DPI, units = "in"
)
log_pt("Wrote PCA panels to ", OUTPUT_DIR)

# ============================================================
# 5. 📈 Group mean trait bars (raw units) ----
# ============================================================
traits_bars <- traits_raw
traits_bars$Group <- factor(GROUP_LABELS[group_ids], levels = GROUP_LABELS)
trait_order <- c(
  "avg_bark_percent", "height_m", "pruning",
  "longevity_yrs", "seed_mass_mg", "serotiny", "resprout"
)
trait_labels <- c(
  avg_bark_percent = "b.) Bark thickness (%)",
  height_m = "c.) Height (m)",
  pruning = "d.) Self-pruning",
  longevity_yrs = "e.) Longevity (yrs)",
  seed_mass_mg = "f.) Seed mass (mg)",
  serotiny = "g.) Serotiny",
  resprout = "h.) Resprout"
)

summary_group <- traits_bars %>%
  select(Group, all_of(trait_order)) %>%
  pivot_longer(-Group, names_to = "Trait", values_to = "Value") %>%
  group_by(Group, Trait) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
    .groups = "drop"
  )
summary_group$Trait <- factor(
  summary_group$Trait,
  levels = trait_order,
  labels = trait_labels[trait_order]
)

p_means <- ggplot(summary_group, aes(Group, mean_value, fill = Group)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(~Trait, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = GROUP_COLORS) +
  theme_classic(base_size = 16) +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, face = "bold", size = 14),
    axis.text = element_text(size = 13)
  )

ggsave(
  file.path(OUTPUT_DIR, "group_mean_traits_bars.png"),
  p_means, width = 5, height = 10, dpi = 300
)
log_pt("Wrote group mean bars: ", file.path(OUTPUT_DIR, "group_mean_traits_bars.png"))

# ============================================================
# 6. 🧪 Genus × cluster permutation (taxonomy, not a phylogeny) ----
# ============================================================
# Reviewer ask: quantify the qualitative genus-level remark. Genus from the
# binomial; permute cluster labels with genus fixed. Not Blomberg's K / a tree.
genus_vec <- vapply(
  strsplit(species_names, "\\s+"),
  function(x) x[[1L]],
  character(1L)
)
cluster_vec <- GROUP_LABELS[group_ids]
genus_n <- as.integer(table(genus_vec)[genus_vec])

genus_by_cluster <- as.data.frame(table(
  genus = genus_vec,
  cluster = factor(cluster_vec, levels = GROUP_LABELS)
), stringsAsFactors = FALSE)
names(genus_by_cluster)[3] <- "n"
genus_by_cluster <- genus_by_cluster[genus_by_cluster$n > 0L, ]
genus_by_cluster$n_in_genus <- as.integer(table(genus_vec)[genus_by_cluster$genus])
write.csv(
  genus_by_cluster[order(genus_by_cluster$n_in_genus, genus_by_cluster$genus, decreasing = TRUE), ],
  file.path(OUTPUT_DIR, "genus_by_cluster.csv"),
  row.names = FALSE
)

# χ² on the genus × cluster contingency (association).
chi2_stat <- function(genus, cluster) {
  suppressWarnings(as.numeric(chisq.test(table(genus, cluster))$statistic))
}

# Among species pairs that share a genus, fraction in the same cluster.
congeneric_concordance <- function(genus, cluster) {
  n <- length(genus)
  n_pairs <- 0L
  n_same <- 0L
  for (i in seq_len(n - 1L)) {
    same_g <- which(genus[(i + 1L):n] == genus[i])
    if (!length(same_g)) {
      next
    }
    n_pairs <- n_pairs + length(same_g)
    n_same <- n_same + sum(cluster[i + same_g] == cluster[i])
  }
  c(concordance = n_same / n_pairs, n_pairs = n_pairs, n_same = n_same)
}

# Multi-species genera that fall entirely in one cluster.
frac_pure_multi <- function(genus, cluster, genus_n) {
  multi <- unique(genus[genus_n >= 2L])
  if (!length(multi)) {
    return(c(frac_pure = NA_real_, n_multi = 0, n_pure = 0))
  }
  n_pure <- sum(vapply(multi, function(g) {
    length(unique(cluster[genus == g])) == 1L
  }, logical(1L)))
  c(frac_pure = n_pure / length(multi), n_multi = length(multi), n_pure = n_pure)
}

obs_chi2 <- chi2_stat(genus_vec, cluster_vec)
obs_conc <- congeneric_concordance(genus_vec, cluster_vec)
obs_pure <- frac_pure_multi(genus_vec, cluster_vec, genus_n)

set.seed(GENUS_PERM_SEED)
null_chi2 <- numeric(GENUS_PERM_N)
null_conc <- numeric(GENUS_PERM_N)
null_pure <- numeric(GENUS_PERM_N)
for (b in seq_len(GENUS_PERM_N)) {
  cl_perm <- sample(cluster_vec)
  null_chi2[b] <- chi2_stat(genus_vec, cl_perm)
  null_conc[b] <- congeneric_concordance(genus_vec, cl_perm)[["concordance"]]
  null_pure[b] <- frac_pure_multi(genus_vec, cl_perm, genus_n)[["frac_pure"]]
}

# One-sided: greater association / concordance / purity than null.
p_chi2 <- (1 + sum(null_chi2 >= obs_chi2, na.rm = TRUE)) / (1 + GENUS_PERM_N)
p_conc <- (1 + sum(null_conc >= obs_conc[["concordance"]], na.rm = TRUE)) / (1 + GENUS_PERM_N)
p_pure <- (1 + sum(null_pure >= obs_pure[["frac_pure"]], na.rm = TRUE)) / (1 + GENUS_PERM_N)

perm_out <- data.frame(
  test = c(
    "chi2_genus_x_cluster",
    "congeneric_pair_concordance",
    "frac_pure_multi_species_genera"
  ),
  statistic = c(
    obs_chi2,
    unname(obs_conc[["concordance"]]),
    unname(obs_pure[["frac_pure"]])
  ),
  null_mean = c(
    mean(null_chi2, na.rm = TRUE),
    mean(null_conc, na.rm = TRUE),
    mean(null_pure, na.rm = TRUE)
  ),
  p_permute = c(p_chi2, p_conc, p_pure),
  n_perm = GENUS_PERM_N,
  seed = GENUS_PERM_SEED,
  n_species = length(species_names),
  n_genera = length(unique(genus_vec)),
  n_multi_species_genera = unname(obs_pure[["n_multi"]]),
  n_pure_multi_species_genera = unname(obs_pure[["n_pure"]]),
  n_congeneric_pairs = unname(obs_conc[["n_pairs"]]),
  n_congeneric_pairs_same_cluster = unname(obs_conc[["n_same"]]),
  stringsAsFactors = FALSE
)
write.csv(
  perm_out,
  file.path(OUTPUT_DIR, "genus_cluster_permutation.csv"),
  row.names = FALSE
)
log_pt(
  "Genus × cluster permutation (n=", GENUS_PERM_N, "): ",
  "chi2 p=", signif(p_chi2, 3),
  "; concordance=", round(obs_conc[["concordance"]], 3),
  " (null mean ", round(mean(null_conc), 3), ", p=", signif(p_conc, 3), ")",
  "; pure multi-spp genera=", obs_pure[["n_pure"]], "/", obs_pure[["n_multi"]],
  " (p=", signif(p_pure, 3), ")"
)
log_pt("Wrote ", file.path(OUTPUT_DIR, "genus_cluster_permutation.csv"))
