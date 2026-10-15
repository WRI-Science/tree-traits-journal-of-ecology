#!/usr/bin/env Rscript
# 09 — Resistance, recovery, and bivariate maps (Fig 4; basal area)
#
# Question: Spatial resistance / recovery under BA weighting.
# Why: Tertile bivariate is built from the painted BA layers (not the published
#      count-weighted bivariate). Pixels are masked to the 11 western-state
#      polygons (not the bounding box); a white overlay covers any cell that
#      still straddles the line.
# Inputs:  OUT_ROOT/rasters/resistance.tif, recovery.tif
# Outputs: OUT_ROOT/figure_04/
# Run:     Rscript scripts/09_resistance_recovery_bivariate_maps.R

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
library(pals)
library(cowplot)

# ============================================================
# ⚙️ Settings ----
# ============================================================
DISPLAY_N <- 5e5
WESTERN_STATES <- c(
  "California", "Oregon", "Washington", "Nevada", "Idaho",
  "Montana", "Wyoming", "Colorado", "Utah", "Arizona", "New Mexico"
)

# ============================================================
# 📁 File paths ----
# ============================================================
OUT_ROOT <- Sys.getenv("TREE_TRAITS_OUT_ROOT", file.path(getwd(), "output"))
RASTER_DIR <- file.path(OUT_ROOT, "rasters")
OUTPUT_DIR <- file.path(OUT_ROOT, "figure_04")
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

# Downsample, mask to the state union, and write a small GeoTIFF so layout
# tweaks do not re-read the 30 m grid. The white ggplot overlay is the
# visual guarantee: Albers cells that straddle a border still show otherwise.
downsample_mask <- function(src, display_tif, plot_ext, vb) {
  if (file.exists(display_tif)) {
    log_pt("Display cache ", basename(display_tif))
    return(rast(display_tif))
  }
  r <- rast(src)
  r_plot <- crop(r, plot_ext)
  r_plot <- spatSample(r_plot, size = DISPLAY_N, method = "regular", as.raster = TRUE)
  r_plot <- crop(r_plot, vb, mask = TRUE, touches = FALSE)
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
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      plot.margin = margin(4, 2, 2, 2)
    )
}

# ============================================================
# 1. 🗺️ Western-state polygons (not the bbox) ----
# ============================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DISPLAY_DIR, recursive = TRUE, showWarnings = FALSE)
res_tif <- file.path(RASTER_DIR, "resistance.tif")
rec_tif <- file.path(RASTER_DIR, "recovery.tif")
if (!file.exists(res_tif) || !file.exists(rec_tif)) {
  stop("Paint rasters first (04_paint_rasters.R): ", RASTER_DIR)
}

r0 <- rast(res_tif)
b <- ne_states(country = "United States of America", returnclass = "sf")
b <- subset(b, name %in% WESTERN_STATES)
b <- st_make_valid(st_transform(b, crs(r0)))
b_union <- st_make_valid(st_union(b))
vb <- vect(b_union)
bb <- st_bbox(b)
plot_ext <- ext(as.vector(bb))
# Buffer the frame so the white patch covers panel corners, not just the
# exact bbox (coord_sf can still show a sliver of Dakotas / Mexico).
universe <- st_as_sfc(bb, crs = st_crs(b))
universe <- st_buffer(universe, dist = 5e4)
outside <- st_sf(geometry = st_make_valid(st_difference(universe, b_union)))

# ============================================================
# 2. 🖼️ Three panels ----
# ============================================================
log_pt("Resistance")
r_res <- downsample_mask(
  res_tif, file.path(DISPLAY_DIR, "resistance.tif"), plot_ext, vb
)
p_res <- ggplot() +
  geom_spatraster(data = r_res, maxcell = Inf) +
  geom_sf(data = outside, fill = "white", color = NA) +
  geom_sf(data = b, fill = NA, color = "black", linewidth = 0.35) +
  scale_fill_viridis_c(
    option = "plasma", limits = c(0, 1), oob = squish,
    na.value = "white", guide = "none"
  ) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE, datum = NA) +
  map_theme() +
  labs(title = "a.) Resistance")
ggsave(
  file.path(OUTPUT_DIR, "map_resistance.png"),
  p_res, width = 4.2, height = 5.2, dpi = 300, units = "in", bg = "white"
)

log_pt("Recovery")
r_rec <- downsample_mask(
  rec_tif, file.path(DISPLAY_DIR, "recovery.tif"), plot_ext, vb
)
p_rec <- ggplot() +
  geom_spatraster(data = r_rec, maxcell = Inf) +
  geom_sf(data = outside, fill = "white", color = NA) +
  geom_sf(data = b, fill = NA, color = "black", linewidth = 0.35) +
  scale_fill_viridis_c(
    option = "plasma", limits = c(0, 1), oob = squish,
    na.value = "white", guide = "none"
  ) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE, datum = NA) +
  map_theme() +
  labs(title = "b.) Recovery")
ggsave(
  file.path(OUTPUT_DIR, "map_recovery.png"),
  p_rec, width = 4.2, height = 5.2, dpi = 300, units = "in", bg = "white"
)

# Tertiles on the display grid (avoid ifel on the 30 m rasters).
log_pt("Bivariate")
bivar_display <- file.path(DISPLAY_DIR, "resistance_recovery_bivariate.tif")
bivar_raster <- file.path(RASTER_DIR, "resistance_recovery_bivariate_display.tif")
if (file.exists(bivar_display)) {
  r_bivar <- rast(bivar_display)
} else if (file.exists(bivar_raster)) {
  r_bivar <- crop(rast(bivar_raster), vb, mask = TRUE, touches = FALSE)
  writeRaster(
    r_bivar, bivar_display, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
} else {
  # Fixed cutpoints 0.33 / 0.66 on [0,1] indices — same encoding as published maps.
  log_pt("Building bivariate tertiles from display resistance/recovery")
  tertile <- function(r) ifel(r <= 0.33, 1L, ifel(r <= 0.66, 2L, 3L))
  r_bivar <- (tertile(r_res) - 1L) * 3L + tertile(r_rec)
  writeRaster(
    r_bivar, bivar_display, overwrite = TRUE,
    wopt = list(datatype = "INT1U", gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
  writeRaster(
    r_bivar, bivar_raster, overwrite = TRUE,
    wopt = list(datatype = "INT1U", gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
}
bivar_colors <- setNames(rev(brewer.seqseq2(9)), as.character(1:9))
p_bivar <- ggplot() +
  geom_spatraster(data = r_bivar, aes(fill = factor(after_stat(value))), maxcell = Inf) +
  geom_sf(data = outside, fill = "white", color = NA) +
  geom_sf(data = b, fill = NA, color = "black", linewidth = 0.35) +
  scale_fill_manual(values = bivar_colors, na.value = "white", guide = "none") +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE, datum = NA) +
  map_theme() +
  labs(title = "c.) Resistance x Recovery")
ggsave(
  file.path(OUTPUT_DIR, "map_resistance_recovery_bivariate.png"),
  p_bivar, width = 4.2, height = 5.2, dpi = 300, units = "in", bg = "white"
)

# ============================================================
# 3. 🧩 Composite (full panels; do not trim the east edge) ----
# ============================================================
log_pt("Building Fig 4 composite")
composite <- plot_grid(p_res, p_rec, p_bivar, nrow = 1, align = "hv")
ggsave(
  file.path(OUTPUT_DIR, "figure_04_composite.png"),
  composite, width = 12.6, height = 5.2, dpi = 300, units = "in", bg = "white"
)
log_pt("Done. ", OUTPUT_DIR)
