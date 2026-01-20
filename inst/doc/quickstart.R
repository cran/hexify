## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
library(hexify)

## ----cover-image, echo=FALSE, message=FALSE, warning=FALSE, fig.width=8, fig.height=4, fig.align='center'----
# Generate cover image: multi-resolution grids over different regions
library(sf)
library(ggplot2)

# Get world basemap
countries <- hexify_world

# Generate grids at different resolutions for specific regions (coarser for speed)
# Coarse grid: South America
grid_sa <- hex_grid(area_km2 = 300000)
sa_hexes <- grid_rect(c(-80, -40, -35, 10), grid_sa)

# Medium grid: Europe/Africa
grid_eu <- hex_grid(area_km2 = 80000)
eu_hexes <- grid_rect(c(-10, -10, 40, 50), grid_eu)

# Fine grid: East Asia
grid_asia <- hex_grid(area_km2 = 50000)
asia_hexes <- grid_rect(c(90, 10, 140, 50), grid_asia)

ggplot() +
  geom_sf(data = countries, fill = "gray95", color = "gray80", linewidth = 0.15) +
  geom_sf(data = sa_hexes, fill = NA, color = "#1B9E77", linewidth = 0.5) +
  geom_sf(data = eu_hexes, fill = NA, color = "#D95F02", linewidth = 0.4) +
  geom_sf(data = asia_hexes, fill = NA, color = "#7570B3", linewidth = 0.3) +
  coord_sf(xlim = c(-100, 150), ylim = c(-50, 60)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

## ----basic-usage--------------------------------------------------------------
# Sample data: European cities
cities <- data.frame(
  name = c("Vienna", "Paris", "Madrid", "Berlin", "Rome"),
  lon = c(16.37, 2.35, -3.70, 13.40, 12.50),
  lat = c(48.21, 48.86, 40.42, 52.52, 41.90)
)

# Create a grid specification
grid <- hex_grid(area_km2 = 10000)
grid

# Assign cities to hexagonal cells
result <- hexify(cities, lon = "lon", lat = "lat", grid = grid)
result

## ----hexdata-access-----------------------------------------------------------
# Get the grid specification
grid_info(result)

# Get unique cell IDs
cells(result)

# Count unique cells
n_cells(result)

# Access all cell IDs (one per row)
result@cell_id

# Access cell centers
head(result@cell_center)

# Extract original data as data.frame
head(as.data.frame(result))

## ----sf-usage, eval = requireNamespace("sf", quietly = TRUE)------------------
library(sf)

# Create sf object (any CRS works - hexify transforms automatically)
pts <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)

# hexify handles CRS transformation automatically
result_sf <- hexify(pts, area_km2 = 10000)
class(result_sf)

## ----bird-setup, message=FALSE, warning=FALSE---------------------------------
library(sf)
library(ggplot2)

# Simulate bird observation data (reduced for faster vignette build)
set.seed(123)
n_obs <- 500

# Generate observations with realistic spatial clustering
birds <- data.frame(
  lon = c(
    rnorm(150, mean = 10, sd = 15),    # Western Europe
    rnorm(100, mean = 25, sd = 10),    # Eastern Europe
    rnorm(100, mean = 20, sd = 20),    # Mediterranean
    rnorm(80, mean = 0, sd = 15),      # West Africa
    rnorm(70, mean = 35, sd = 10)      # East Africa
  ),
  lat = c(
    rnorm(150, mean = 50, sd = 8),     # Western Europe
    rnorm(100, mean = 55, sd = 6),     # Eastern Europe
    rnorm(100, mean = 42, sd = 5),     # Mediterranean
    rnorm(80, mean = 10, sd = 10),     # West Africa
    rnorm(70, mean = -5, sd = 12)      # East Africa
  ),
  species = sample(c("Passer domesticus", "Turdus merula", "Parus major",
                     "Columba palumbus", "Sturnus vulgaris"), n_obs, replace = TRUE)
)

# Clip to valid ranges
birds$lon <- pmax(-30, pmin(60, birds$lon))
birds$lat <- pmax(-35, pmin(70, birds$lat))

## ----bird-assign--------------------------------------------------------------
# Create grid and assign observations (coarser grid for faster build)
grid <- hex_grid(area_km2 = 50000)
birds_hex <- hexify(birds, lon = "lon", lat = "lat", grid = grid)

# Extract data with cell IDs
birds_gridded <- as.data.frame(birds_hex)
birds_gridded$cell_id <- birds_hex@cell_id

# Count observations per cell
obs_counts <- aggregate(
  species ~ cell_id,
  data = birds_gridded,
  FUN = length
)
names(obs_counts)[2] <- "n_observations"

# Species richness per cell
richness <- aggregate(
  species ~ cell_id,
  data = birds_gridded,
  FUN = function(x) length(unique(x))
)
names(richness)[2] <- "n_species"

obs_counts <- merge(obs_counts, richness, by = "cell_id")
head(obs_counts)

## ----bird-plot, fig.width=8, fig.height=6, message=FALSE, warning=FALSE-------
# Generate polygons for cells with data
cell_polys <- cell_to_sf(obs_counts$cell_id, grid)
cell_polys <- merge(cell_polys, obs_counts, by = "cell_id")

# Get relevant countries
region <- hexify_world[hexify_world$continent %in% c("Europe", "Africa"), ]

ggplot() +
  geom_sf(data = region, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = cell_polys, aes(fill = n_observations), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "Observations", trans = "sqrt") +
  coord_sf(xlim = c(-30, 60), ylim = c(-35, 70)) +
  labs(
    title = "Bird Observations in Equal-Area Hexagonal Cells",
    subtitle = sprintf("ISEA3H grid at resolution %d (~%.0f km² cells)",
                       grid@resolution, grid@area_km2)
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "gray90")
  )

## ----richness-plot, fig.width=8, fig.height=6, message=FALSE------------------
ggplot() +
  geom_sf(data = region, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = cell_polys, aes(fill = n_species), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "mako", name = "Species\nRichness", direction = -1) +
  coord_sf(xlim = c(-30, 60), ylim = c(-35, 70)) +
  labs(
    title = "Species Richness per Grid Cell",
    subtitle = "Number of unique species observed in each equal-area cell"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "gray90")
  )

