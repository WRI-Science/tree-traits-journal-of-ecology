#!/usr/bin/env Rscript
# 08 — Mapped CWMTs from painted BA rasters (Fig 3 composite)
#
# Question: Spatial pattern of each CWMT under basal-area weighting.
# Why: Mask to the 11 western states (not the bbox). Resistance traits on the
#      top row, recovery on the bottom. Longevity and seed mass use a log fill.
# Inputs:  OUT_ROOT/rasters/<trait>.tif
# Outputs: OUT_ROOT/figure_03/map_<trait>.png
#          OUT_ROOT/figure_03/figure_03_composite.png
# Run:     Rscript scripts/08_cwmt_trait_maps.R

# ============================================================
# 📦 Libraries ----
# ============================================================
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(tidyterra)
library(viridis)
library(scales)
library(patchwork)
library(cowplot)
library(grid)

# ============================================================
# ⚙️ Settings ----
# ============================================================
WESTERN_STATES <- c(
  "California", "Oregon", "Washington", "Nevada", "Idaho",
  "Montana", "Wyoming", "Colorado", "Utah", "Arizona", "New Mexico"
)
# Resistance (top row, 3 + spacer) then recovery (bottom row, 4).
# Index traits are 0–1 by definition. Seed mass is log10 with 1–1000 mg
# so a few very large CWMs (p99 ~ 2800 mg) do not wash out the map.
TRAIT_SPECS <- list(
  list(
    col = "avg_bark_percent", letter = "a",
    title = "Bark thickness", fill = "% bark",
    log10 = FALSE, limits = NULL
  ),
  list(
    col = "height_m", letter = "b",
    title = "Height", fill = "m",
    log10 = FALSE, limits = NULL
  ),
  list(
    col = "pruning", letter = "c",
    title = "Self-pruning", fill = "index",
    log10 = FALSE, limits = c(0, 1)
  ),
  list(
    col = "longevity_yrs", letter = "d",
    title = "Longevity (log)", fill = "years",
    log10 = TRUE, limits = c(100, 2000),
    breaks = c(100, 300, 1000)
  ),
  list(
    col = "seed_mass_mg", letter = "e",
    title = "Seed mass (log)", fill = "mg",
    log10 = TRUE, limits = c(1, 1000),
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")
  ),
  list(
    col = "serotiny", letter = "f",
    title = "Serotiny", fill = "index",
    log10 = FALSE, limits = c(0, 1)
  ),
  list(
    col = "resprout", letter = "g",
    title = "Resprout", fill = "index",
    log10 = FALSE, limits = c(0, 1)
  )
)
# Display-grid size after the 30 m crop. Cached under figure_03/display/.
DISPLAY_N <- 5e5

# ============================================================
# 📁 File paths ----
# ============================================================
OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
RASTER_DIR <- file.path(OUT_ROOT, "rasters")
OUTPUT_DIR <- file.path(OUT_ROOT, "figure_03")
DISPLAY_DIR <- file.path(OUTPUT_DIR, "display")

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

row_label <- function(txt) {
  ggplot() +
    theme_void() +
    annotate(
      "text", x = 0.5, y = 0.5, label = txt, angle = 90,
      fontface = "bold", size = 4.2
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
}

# ============================================================
# 1. 🗺️ Western-state mask (polygon, not bbox) ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DISPLAY_DIR, recursive = TRUE, showWarnings = FALSE)
boundaries_sf <- ne_states(country = "United States of America", returnclass = "sf")
boundaries_sf <- subset(boundaries_sf, name %in% WESTERN_STATES)

# ============================================================
# 2. 🖼️ One panel per trait ----
# ============================================================
panels <- vector("list", length(TRAIT_SPECS))
b <- NULL
outside <- NULL
bb <- NULL
for (i in seq_along(TRAIT_SPECS)) {
  spec <- TRAIT_SPECS[[i]]
  src_tif <- file.path(RASTER_DIR, paste0(spec$col, ".tif"))
  if (!file.exists(src_tif)) {
    stop("Paint rasters first (04_paint_rasters.R): ", src_tif)
  }
  log_pt("Plotting ", spec$col)
  display_tif <- file.path(DISPLAY_DIR, paste0(spec$col, ".tif"))
  r <- rast(src_tif)
  if (is.null(b)) {
    b <- st_make_valid(st_transform(boundaries_sf, crs(r)))
    b_union <- st_make_valid(st_union(b))
    bb <- st_bbox(b)
    universe <- st_buffer(st_as_sfc(bb, crs = st_crs(b)), dist = 5e4)
    outside <- st_sf(geometry = st_make_valid(st_difference(universe, b_union)))
  }
  if (file.exists(display_tif)) {
    log_pt("Display cache ", spec$col)
    r_plot <- rast(display_tif)
  } else {
    plot_ext <- ext(as.vector(bb))
    r_plot <- crop(r, plot_ext)
    # Downsample, then mask to the union of the 11 states. A bbox crop still
    # includes the Dakotas/Nebraska sliver east of MT/WY (CO's eastern
    # border is farther east, so the rectangle is not the state shapes).
    r_plot <- spatSample(r_plot, size = DISPLAY_N, method = "regular", as.raster = TRUE)
    r_plot <- crop(r_plot, vect(st_union(b)), mask = TRUE, touches = FALSE)
    writeRaster(
      r_plot, display_tif, overwrite = TRUE,
      wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"))
    )
    r_plot <- rast(display_tif)
  }
  fill_args <- list(
    option = "plasma",
    na.value = "white",
    oob = squish,
    name = spec$fill
  )
  if (isTRUE(spec$log10)) {
    fill_args$trans <- "log10"
  }
  if (!is.null(spec$limits)) {
    fill_args$limits <- spec$limits
  }
  if (identical(spec$col, "seed_mass_mg")) {
    v <- values(r_plot, mat = FALSE)
    v <- v[is.finite(v) & v > 0]
    fill_args$limits <- as.numeric(quantile(v, c(0.02, 0.98)))
  }
  if (!is.null(spec$breaks)) {
    fill_args$breaks <- spec$breaks
  }
  if (!is.null(spec$labels)) {
    fill_args$labels <- spec$labels
  }
  fill_scale <- do.call(scale_fill_viridis_c, fill_args)
  p <- ggplot() +
    geom_spatraster(data = r_plot, maxcell = Inf) +
    geom_sf(data = outside, fill = "white", color = NA) +
    geom_sf(data = b, fill = NA, color = "black", linewidth = 0.2) +
    fill_scale +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE, datum = NA
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 10, hjust = 0, face = "bold"),
      plot.margin = margin(2, 2, 0, 2),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt")
    ) +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        barwidth = unit(1.6, "in"),
        barheight = unit(0.12, "in"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(title = paste0(spec$letter, ".) ", spec$title))
  ggsave(
    file.path(OUTPUT_DIR, paste0("map_", spec$col, ".png")),
    p, width = 4.2, height = 4.4, dpi = 300, units = "in", bg = "white"
  )
  panels[[i]] <- p
  rm(r, r_plot, p)
  gc()
}

# ============================================================
# 3. 🧩 Composite: resistance (3) over recovery (4) ----
# ============================================================
log_pt("Building Fig 3 composite")
# 4-column grid so every map is the same size. Empty cell at top-right
# (resistance has 3 traits; recovery has 4).
top_maps <- plot_grid(
  panels[[1]], panels[[2]], panels[[3]], NULL,
  nrow = 1, rel_widths = c(1, 1, 1, 1), align = "hv"
)
bot_maps <- plot_grid(
  panels[[4]], panels[[5]], panels[[6]], panels[[7]],
  nrow = 1, rel_widths = c(1, 1, 1, 1), align = "hv"
)
maps <- plot_grid(top_maps, bot_maps, ncol = 1, rel_heights = c(1, 1))
labs <- plot_grid(
  row_label("Resistance"), row_label("Recovery"),
  ncol = 1, rel_heights = c(1, 1)
)
composite <- plot_grid(labs, maps, nrow = 1, rel_widths = c(0.045, 1))
ggsave(
  file.path(OUTPUT_DIR, "figure_03_composite.png"),
  composite,
  width = 13, height = 8.2, dpi = 300, units = "in", bg = "white"
)
log_pt("Done. ", OUTPUT_DIR)
