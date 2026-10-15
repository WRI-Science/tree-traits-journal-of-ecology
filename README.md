# tree-traits-journal-of-ecology

Analysis code for a *Journal of Ecology* manuscript on western CONUS tree traits and fire regime groups (FRG). Community-weighted mean traits (CWMTs) use **basal area** abundance weighting only (`TPA_UNADJ × 0.005454 × DIA²`).

This repository contains **scripts** and the **species trait table**. Large geospatial inputs (TreeMap, LANDFIRE) are not bundled; configure local paths before running the full pipeline.

**License:** MIT

---

## Contents

```text
├── data/
│   └── traits_95_west_conus.csv   # archived with this release (see Data availability)
├── scripts/                       # numbered pipeline 00–13
├── DEPENDENCIES.md
├── packages.txt
└── README.md
```

---

## Data availability

**Species traits (this repository).**  
`data/traits_95_west_conus.csv` is archived with the analysis code at Zenodo: https://doi.org/10.5281/zenodo.23401676 (release `v1.0.1`).

**TreeMap 2016.**  
Riley, K.L.; Grenfell, I.C.; Finney, M.A.; Shaw, J.D. 2021. TreeMap 2016: A tree-level model of the forests of the conterminous United States circa 2016. Fort Collins, CO: Forest Service Research Data Archive. https://doi.org/10.2737/RDS-2021-0074

**LANDFIRE Fire Regime Groups.**  
LANDFIRE. Fire Regime Groups. https://www.landfire.gov/fire-regime/frg (LF 2016 Remap / FRG products as used in the study).

**Analysis code.**  
GitHub: https://github.com/WRI-Science/tree-traits-journal-of-ecology  
Zenodo (preferred citation): https://doi.org/10.5281/zenodo.23401676

### Draft Data Availability Statement (for the manuscript)

> Species trait data and analysis code supporting this study are available in Zenodo at https://doi.org/10.5281/zenodo.23401676. TreeMap 2016 data are available from the USDA Forest Service Research Data Archive (Riley et al. 2021; https://doi.org/10.2737/RDS-2021-0074). LANDFIRE Fire Regime Groups data are available from https://www.landfire.gov/fire-regime/frg. A working copy of the code is also available at https://github.com/WRI-Science/tree-traits-journal-of-ecology.

---

## Configure inputs (required for full pipeline)

Species traits default to `data/traits_95_west_conus.csv` in this repo. All other geospatial inputs must be set via environment variables (no machine-specific paths are hard-coded):

| Variable | What to point it at |
|----------|---------------------|
| `TREE_TRAITS_TREEMAP_TREE_CSV` | TreeMap 2016 tree table CSV |
| `TREE_TRAITS_STUDY_TREEMAP_TIF` | Study-region TreeMap raster (GeoTIFF) |
| `TREE_TRAITS_FRG_TIF` | FRG raster masked to the study region |
| `TREE_TRAITS_FRG_LOOKUP_CSV` | LANDFIRE FRG attribute lookup CSV (LF 2016) |
| `TREE_TRAITS_OUT_ROOT` | Output directory (default: `./output`) |
| `TREE_TRAITS_SPECIES_CSV` | Optional override of the traits CSV |

Also required on `PATH` for steps 03 and 05: **GDAL** (`gdalwarp`).

Example:

```bash
export TREE_TRAITS_TREEMAP_TREE_CSV=/path/to/TreeMap2016_tree_table.csv
export TREE_TRAITS_STUDY_TREEMAP_TIF=/path/to/study_treemap.tif
export TREE_TRAITS_FRG_TIF=/path/to/frg_masked.tif
export TREE_TRAITS_FRG_LOOKUP_CSV=/path/to/LF16_FRG_200.csv
export TREE_TRAITS_OUT_ROOT=./output
```

Fig 1 / Supplement S1 (species-level only) need only the traits CSV and can be run without TreeMap/FRG.

---

## How to run

From the repository root:

```bash
Rscript scripts/00_run_all.R
```

Or step by step: `Rscript scripts/01_….R` … `Rscript scripts/13_….R`.

### Run order

| # | Script | Role |
|---|--------|------|
| 00 | `00_run_all.R` | Runs 01–13 |
| 01 | `01_individual_trait_space.R` | Fig 1 — species groups + PCA |
| 02 | `02_supplement_individual_traits_by_group.R` | Supplement S1 |
| 03 | `03_rebuild_ba_cwmts.R` | BA CWMTs + majority FRG per `tm_id` |
| 04 | `04_paint_rasters.R` | Paint TreeMap with BA CWMTs |
| 05 | `05_pixel_sample.R` | Stratified pixel sample by FRG |
| 06–07 | Community PCA + CWMT boxplots | Fig 2 |
| 08–10 | Trait / resistance / resilience maps | Figs 3–5 |
| 11 | FRG boxplots | Fig 6 |
| 12 | PERMANOVA + ANOVA tables | Tables S1–S2 |
| 13 | FRG purity sensitivity | Supplement S2/S3 |

---

## Dependencies

See `DEPENDENCIES.md` and `packages.txt`.

---

## Zenodo (DOI for the journal)

Archived release: **https://doi.org/10.5281/zenodo.23401676** (`v1.0.1`).  
Record page: https://zenodo.org/records/23401676
