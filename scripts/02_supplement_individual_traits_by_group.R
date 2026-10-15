#!/usr/bin/env Rscript
# 02 — Supplement S1: Individual trait values by emergent functional group
#
# Question: How do raw trait values distribute across species within each emergent group?
# Why: Complements Fig 1 with per-species bars (continuous) and a grid (binary traits).
# Inputs:  traits_95_west_conus.csv
#          OUT_ROOT/figure_01/species_functional_groups.csv (from 01)
# Outputs: OUT_ROOT/supplement_S1/continuous_traits_by_species.png
#          OUT_ROOT/supplement_S1/binary_traits_grid.png
# Run:     Rscript scripts/02_supplement_individual_traits_by_group.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# ⚙️ Settings ----
# ============================================================
GROUP_COLOR_PALETTE <- c("#3a5311", "#728c69", "#028a04", "#74b72e", "#b2d3c2")
BINARY_TRAITS <- c("serotiny", "pruning", "resprout")
CONTINUOUS_TRAITS <- c("avg_bark_percent", "longevity_yrs", "seed_mass_mg", "height_m")

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
GROUPS_CSV <- file.path(OUT_ROOT, "figure_01", "species_functional_groups.csv")
OUTPUT_DIR <- file.path(OUT_ROOT, "supplement_S1")

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

# ============================================================
# 1. 📥 Load traits + group assignments ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SPECIES_CSV)) {
  stop("Species CSV not found: ", SPECIES_CSV)
}
if (!file.exists(GROUPS_CSV)) {
  stop("Run 01_individual_trait_space.R first. Missing: ", GROUPS_CSV)
}

traits <- read.csv(SPECIES_CSV, stringsAsFactors = FALSE)
groups <- read.csv(GROUPS_CSV, stringsAsFactors = FALSE)
traits <- merge(traits, groups[, c("scientific_name", "emergent_group", "emergent_group_label")],
                by.x = "SCIENTIFIC_NAME", by.y = "scientific_name", all.x = TRUE)
if (any(is.na(traits$emergent_group))) {
  stop("Some species lack emergent group assignments.")
}

n_groups <- max(groups$emergent_group, na.rm = TRUE)
group_levels <- groups %>%
  distinct(emergent_group, emergent_group_label) %>%
  arrange(emergent_group) %>%
  pull(emergent_group_label)
if (length(group_levels) != n_groups) {
  group_levels <- LETTERS[seq_len(n_groups)]
}
GROUP_COLORS <- GROUP_COLOR_PALETTE[seq_len(n_groups)]
names(GROUP_COLORS) <- group_levels

species_order <- traits %>%
  arrange(emergent_group, SCIENTIFIC_NAME) %>%
  pull(SCIENTIFIC_NAME) %>%
  unique()

traits$emergent_group <- factor(
  traits$emergent_group_label,
  levels = group_levels
)
traits$SCIENTIFIC_NAME <- factor(traits$SCIENTIFIC_NAME, levels = species_order)

# ============================================================
# 2. 📊 Continuous + binary trait panels ----
# ============================================================
traits_long <- traits %>%
  select(SCIENTIFIC_NAME, emergent_group, all_of(c(CONTINUOUS_TRAITS, BINARY_TRAITS))) %>%
  pivot_longer(-c(SCIENTIFIC_NAME, emergent_group), names_to = "Trait", values_to = "Value") %>%
  mutate(
    Trait = dplyr::recode(
      Trait,
      avg_bark_percent = "a.) Bark (%)",
      longevity_yrs = "b.) Longevity (yrs)",
      seed_mass_mg = "c.) Seed mass (mg)",
      height_m = "d.) Height (m)",
      serotiny = "e.) Serotiny",
      pruning = "f.) Self-pruning",
      resprout = "g.) Resprout"
    )
  )

BINARY_TRAIT_LABELS <- c("e.) Serotiny", "f.) Self-pruning", "g.) Resprout")

p_cont <- traits_long %>%
  filter(!Trait %in% BINARY_TRAIT_LABELS) %>%
  ggplot(aes(SCIENTIFIC_NAME, Value, fill = emergent_group)) +
  geom_col(width = 0.7) +
  facet_wrap(~Trait, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = GROUP_COLORS) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL)

p_bin <- traits_long %>%
  filter(Trait %in% BINARY_TRAIT_LABELS) %>%
  mutate(
    Trait = factor(Trait, levels = BINARY_TRAIT_LABELS)
  ) %>%
  ggplot(aes(SCIENTIFIC_NAME, Trait)) +
  geom_tile(fill = "white", color = "gray80") +
  geom_text(
    aes(label = ifelse(Value == 1, "\u2726", NA_character_)),
    size = 3.5,
    na.rm = TRUE
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic", size = 7),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(file.path(OUTPUT_DIR, "continuous_traits_by_species.png"), p_cont,
       width = 8, height = 10, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "binary_traits_grid.png"), p_bin,
       width = 12, height = 3, dpi = 300)
log_pt("Wrote S1 panels to ", OUTPUT_DIR)
