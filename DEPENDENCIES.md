# R dependencies (BA manuscript pipeline)

Pinned versions observed in the development environment (2026-07-12). No `renv` lockfile — install with the lines below if needed.

```r
install.packages(c(
  "cluster", "ggplot2", "dplyr", "tidyr", "data.table", "vegan",
  "terra", "sf", "rnaturalearth", "tidyterra", "viridis", "scales",
  "pals", "ggrepel", "patchwork", "cowplot"
))
```

System tools: **GDAL** (`gdalwarp`) for FRG→TreeMap alignment in `03` and `05`.

| Package | Version (session) | Used by |
|---------|-------------------|---------|
| cluster | 2.1.8.1 | 01 |
| ggplot2 | 4.0.2 | 01–02, 06–11, 13 |
| dplyr | 1.2.0 | 01–02 |
| tidyr | 1.3.1 | 01–02 |
| data.table | 1.18.4 | 03–07, 11–13 |
| vegan | 2.6.10 | 12 |
| terra | 1.9.11 | 03–05, 08–10 |
| sf | 1.0.16 | 08–10 |
| rnaturalearth | 1.0.1 | 08–10 |
| tidyterra | 0.6.1 | 08–10 |
| viridis | 0.6.5 | 08–10 |
| scales | 1.4.0 | 08–10 |
| pals | 1.10 | 09 |
| ggrepel | 0.9.5 | 01, 06, 13 (labeled PCA) |
| patchwork | (map composites) | 08 |
| cowplot | (map / Fig 6 composites) | 08–09, 11 |

See also `packages.txt` (one package per line).
