#!/usr/bin/env Rscript
# 10 — Resilience formulation maps (Fig 5; basal area)
#
# Question: Spatial pattern of equivalence, complementarity, compensation (β=0.5).
# Why: Mask to the 11 western-state polygons (not the bbox). A white overlay
#      covers any cell that still straddles the border line.
# Inputs:  OUT_ROOT/rasters/resilience_*.tif
# Outputs: OUT_ROOT/figure_05/
# Run:     Rscript scripts/10_resilience_formulation_maps.R

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

# ============================================================
# ⚙️ Settings ----
# ============================================================
DISPLAY_N <- 5e5
SCENARIOS <- list(
  list(
    file = "resilience_equivalence",
    slug = "equivalence",
    label = "a.) Resilience equivalence"
  ),
  list(
    file = "resilience_compensation_50",
    slug = "compensation_beta_0_50",
    label = "b.) Resilience compensation (β = 0.50)"
  ),
  list(
    file = "resilience_complementarity",
    slug = "complementarity",
    label = "c.) Resilience complementarity"
  )
)
WESTERN_STATES <- c(
  "California", "Oregon", "Washington", "Nevada", "Idaho",
  "Montana", "Wyoming", "Colorado", "Utah", "Arizona", "New Mexico"
)

# ============================================================
# 📁 File paths ----
# ============================================================
OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
RASTER_DIR <- file.path(OUT_ROOT, "rasters")
OUTPUT_DIR <- file.path(OUT_ROOT, "figure_05")
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

# Downsample, then mask to the state union. Cached so layout tweaks are cheap.
downsample_mask <- function(src, display_tif, plot_ext, vb) {
  if (file.exists(display_tif) && !identical(Sys.getenv("TREE_TRAITS_FORCE_FIG5", ""), "1")) {
    log_pt("Display cache ", basename(display_tif))
    return(rast(display_tif))
  }
  r <- rast(src)
  r_plot <- crop(r, plot_ext)
  r_plot <- spatSample(r_plot, size = DISPLAY_N, method = "regular", as.raster = TRUE)
  # Rasterize the state union onto the display grid, then mask. Vector mask
  # alone still left Canada / Mexico / ocean cells in the PNG.
  state_grid <- rasterize(vb, r_plot, touches = FALSE)
  r_plot <- mask(r_plot, state_grid)
  writeRaster(
    r_plot, display_tif, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
  rast(display_tif)
}

map_theme <- function() {
  theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ============================================================
# 1. 🗺️ Western-state polygons (not the bbox) ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DISPLAY_DIR, recursive = TRUE, showWarnings = FALSE)

force <- identical(Sys.getenv("TREE_TRAITS_FORCE_FIG5", ""), "1")
out_pngs <- vapply(SCENARIOS, function(sc) {
  file.path(OUTPUT_DIR, paste0("map_resilience_", sc$slug, ".png"))
}, character(1L))
if (!force && all(file.exists(out_pngs)) && all(file.info(out_pngs)$size > 1e5)) {
  log_pt("Existing Fig 5 maps: ", OUTPUT_DIR, " (set TREE_TRAITS_FORCE_FIG5=1 to rebuild)")
  quit(save = "no", status = 0L)
}

r0_path <- file.path(RASTER_DIR, paste0(SCENARIOS[[1]]$file, ".tif"))
if (!file.exists(r0_path)) {
  stop("Paint rasters first (04_paint_rasters.R): ", RASTER_DIR)
}
r0 <- rast(r0_path)
b <- ne_states(country = "United States of America", returnclass = "sf")
b <- subset(b, name %in% WESTERN_STATES)
b <- st_make_valid(st_transform(b, crs(r0)))
b_union <- st_make_valid(st_union(b))
vb <- vect(b_union)
bb <- st_bbox(b)
plot_ext <- ext(as.vector(bb))
# Buffer the frame so the white patch covers panel corners (Canada / Mexico /
# Dakotas slivers), not just the exact bbox.
universe <- st_as_sfc(bb, crs = st_crs(b))
universe <- st_buffer(universe, dist = 1e5)
outside <- st_sf(geometry = st_make_valid(st_difference(universe, b_union)))

# ============================================================
# 2. 🖼️ Three panels ----
# ============================================================
for (sc in SCENARIOS) {
  src_tif <- file.path(RASTER_DIR, paste0(sc$file, ".tif"))
  if (!file.exists(src_tif)) {
    stop("Paint rasters first: ", src_tif)
  }
  log_pt(sc$label)
  display_tif <- file.path(DISPLAY_DIR, paste0(sc$file, ".tif"))
  if (force && file.exists(display_tif)) {
    unlink(display_tif)
  }
  r_plot <- downsample_mask(src_tif, display_tif, plot_ext, vb)
  p <- ggplot() +
    geom_spatraster(data = r_plot, maxcell = Inf) +
    geom_sf(data = outside, fill = "white", color = NA) +
    geom_sf(data = b, fill = NA, color = "black", linewidth = 0.35) +
    scale_fill_viridis_c(
      option = "plasma", limits = c(0, 1), oob = squish,
      na.value = "white", name = "Resilience"
    ) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE, datum = NA
    ) +
    map_theme() +
    labs(title = sc$label)
  ggsave(
    file.path(OUTPUT_DIR, paste0("map_resilience_", sc$slug, ".png")),
    p, width = 5.5, height = 5.2, dpi = 300, units = "in", bg = "white"
  )
}
log_pt("Done. ", OUTPUT_DIR)
