## ----setup, include = FALSE---------------------------------------------------
CAN_RUN <- requireNamespace("sf", quietly = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = CAN_RUN
)
library(hexify)

## ----h3-grid------------------------------------------------------------------
library(sf)
library(ggplot2)

# Create an H3 grid specification
grid_h3 <- hex_grid(resolution = 5, type = "h3")
grid_h3

## ----h3-hexify----------------------------------------------------------------
# Sample cities
cities <- data.frame(
  name = c("Vienna", "Paris", "Madrid", "Berlin", "Rome",
           "London", "Prague", "Warsaw", "Budapest", "Amsterdam"),
  lon = c(16.37, 2.35, -3.70, 13.40, 12.50,
          -0.12, 14.42, 21.01, 19.04, 4.90),
  lat = c(48.21, 48.86, 40.42, 52.52, 41.90,
          51.51, 50.08, 52.23, 47.50, 52.37)
)

result <- hexify(cities, lon = "lon", lat = "lat", grid = grid_h3)
result

## ----h3-cell-ids--------------------------------------------------------------
# Cell IDs are hexadecimal strings
result@cell_id

# All standard accessors work
cells(result)
n_cells(result)

## ----h3-by-area---------------------------------------------------------------
grid_area <- hex_grid(area_km2 = 500, type = "h3")
grid_area

## ----h3-grid-rect, fig.width=7, fig.height=5----------------------------------
# Generate H3 hexagons over Western Europe
grid_h3 <- hex_grid(resolution = 3, type = "h3")
europe_h3 <- grid_rect(c(-10, 35, 25, 60), grid_h3)

# Basemap
europe <- hexify_world[hexify_world$continent == "Europe", ]

ggplot() +
  geom_sf(data = europe, fill = "gray95", color = "gray60") +
  geom_sf(data = europe_h3, fill = NA, color = "#E6550D", linewidth = 0.4) +
  coord_sf(xlim = c(-10, 25), ylim = c(35, 60)) +
  labs(title = sprintf("H3 Resolution %d Grid (~%.0f km² avg cells)",
                       grid_h3@resolution, grid_h3@area_km2)) +
  theme_minimal()

## ----h3-grid-clip, fig.width=6, fig.height=6----------------------------------
# Clip H3 grid to France
france <- hexify_world[hexify_world$name == "France", ]
grid_h3 <- hex_grid(resolution = 4, type = "h3")
france_h3 <- grid_clip(france, grid_h3)

ggplot() +
  geom_sf(data = france, fill = "gray95", color = "gray40", linewidth = 0.5) +
  geom_sf(data = france_h3, fill = alpha("#E6550D", 0.3),
          color = "#E6550D", linewidth = 0.3) +
  coord_sf(xlim = c(-5, 10), ylim = c(41, 52)) +
  labs(title = sprintf("H3 Grid Clipped to France (res %d)", grid_h3@resolution)) +
  theme_minimal()

## ----h3-parent----------------------------------------------------------------
# Get parent cells (one resolution coarser)
grid_h3 <- hex_grid(resolution = 5, type = "h3")
child_ids <- lonlat_to_cell(
  lon = c(16.37, 2.35, 13.40),
  lat = c(48.21, 48.86, 52.52),
  grid = grid_h3
)

parent_ids <- get_parent(child_ids, grid_h3, levels = 1)
data.frame(child = child_ids, parent = parent_ids)

## ----h3-children--------------------------------------------------------------
# Get children of a single cell (one resolution finer)
grid_coarse <- hex_grid(resolution = 3, type = "h3")
coarse_id <- lonlat_to_cell(16.37, 48.21, grid_coarse)

children <- get_children(coarse_id, grid_coarse, levels = 1)
cat(length(children[[1]]), "children at resolution", grid_coarse@resolution + 1, "\n")
head(children[[1]])

## ----h3-hierarchy-plot, fig.width=7, fig.height=5-----------------------------
# Parent cell polygon
parent_poly <- cell_to_sf(coarse_id, grid_coarse)

# Children cell polygons
grid_fine <- hex_grid(resolution = 4, type = "h3")
children_poly <- cell_to_sf(children[[1]], grid_fine)

ggplot() +
  geom_sf(data = children_poly, fill = alpha("#E6550D", 0.3),
          color = "#E6550D", linewidth = 0.5) +
  geom_sf(data = parent_poly, fill = NA, color = "black", linewidth = 1.2) +
  labs(title = sprintf("H3 Hierarchy: 1 parent (res %d) → %d children (res %d)",
                       grid_coarse@resolution,
                       length(children[[1]]),
                       grid_fine@resolution)) +
  theme_minimal()

## ----h3-workflow, fig.width=8, fig.height=6, message=FALSE, warning=FALSE-----
set.seed(42)

# Simulate observations across Europe
obs <- data.frame(
  lon = c(rnorm(200, 10, 12), rnorm(100, 25, 8)),
  lat = c(rnorm(200, 48, 6), rnorm(100, 55, 4)),
  species = sample(c("Sp. A", "Sp. B", "Sp. C"), 300, replace = TRUE)
)
obs$lon <- pmax(-10, pmin(40, obs$lon))
obs$lat <- pmax(35, pmin(65, obs$lat))

# Hexify with H3
grid_h3 <- hex_grid(resolution = 3, type = "h3")
obs_hex <- hexify(obs, lon = "lon", lat = "lat", grid = grid_h3)

# Aggregate: species richness per cell
obs_df <- as.data.frame(obs_hex)
obs_df$cell_id <- obs_hex@cell_id

richness <- aggregate(species ~ cell_id, data = obs_df,
                      FUN = function(x) length(unique(x)))
names(richness)[2] <- "n_species"

# Map it
polys <- cell_to_sf(richness$cell_id, grid_h3)
polys <- merge(polys, richness, by = "cell_id")

europe <- hexify_world[hexify_world$continent == "Europe", ]

ggplot() +
  geom_sf(data = europe, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = polys, aes(fill = n_species), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "Species\nRichness") +
  coord_sf(xlim = c(-10, 40), ylim = c(35, 65)) +
  labs(title = "Species Richness on H3 Grid",
       subtitle = sprintf("H3 resolution %d (~%.0f km² avg cells)",
                          grid_h3@resolution, grid_h3@area_km2)) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## ----h3-crosswalk-------------------------------------------------------------
# Start with an ISEA grid and some cells
grid_isea <- hex_grid(resolution = 9, aperture = 3)
isea_ids <- lonlat_to_cell(
  lon = c(16.37, 2.35, 13.40, -3.70, 12.50),
  lat = c(48.21, 48.86, 52.52, 40.42, 41.90),
  grid = grid_isea
)

# Map ISEA cells to their closest H3 equivalents
xw <- h3_crosswalk(isea_ids, grid_isea)
xw[, c("isea_cell_id", "h3_cell_id", "isea_area_km2", "h3_area_km2")]

## ----h3-resolution-table------------------------------------------------------
h3_res <- hexify_compare_resolutions(type = "h3", res_range = 0:15)
h3_res$n_cells_fmt <- ifelse(
  h3_res$n_cells > 1e9,
  sprintf("%.1fB", h3_res$n_cells / 1e9),
  ifelse(h3_res$n_cells > 1e6,
         sprintf("%.1fM", h3_res$n_cells / 1e6),
         ifelse(h3_res$n_cells > 1e3,
                sprintf("%.1fK", h3_res$n_cells / 1e3),
                as.character(h3_res$n_cells)))
)
knitr::kable(
  h3_res[, c("resolution", "n_cells_fmt", "cell_area_km2", "cell_spacing_km")],
  col.names = c("Resolution", "# Cells", "Avg Area (km²)", "Spacing (km)"),
  digits = 1
)

