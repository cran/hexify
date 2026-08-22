## ----setup, include = FALSE---------------------------------------------------
source("_common.R")
library(hexify)

## ----lambert-geometry, echo=FALSE, fig.width=7, fig.height=6------------------
# Lambert Azimuthal Equal-Area Projection Geometry
# Reference: Snyder (1987), "Map Projections: A Working Manual", p. 182-185
#
# The projection formula ρ = 2R·sin(c/2) is derived analytically from
# the equal-area constraint, not from geometric chord distance.
# This diagram shows the geometric relationship that yields the same formula.

oldpar <- par(mar = c(1, 1, 2, 1), bg = "white")
plot(NULL, xlim = c(-1.5, 2.2), ylim = c(-1.5, 1.7), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "Lambert Azimuthal Equal-Area Projection Geometry")

R <- 1
theta_circle <- seq(0, 2*pi, length.out = 100)

# Draw sphere (great circle cross-section)
lines(R * cos(theta_circle), R * sin(theta_circle), lwd = 2, col = "gray30")

# Draw tangent plane at S (top of sphere)
lines(c(-1.4, 2.1), c(R, R), lwd = 2, col = "gray50")
text(-1.3, R + 0.1, "Tangent plane", adj = 0, cex = 0.9, col = "gray40")

# Mark center O
points(0, 0, pch = 19, cex = 1.2)
text(-0.12, -0.15, "O", cex = 1.1, font = 2)

# Mark tangent point S
points(0, R, pch = 19, cex = 1.2)
text(0.12, R + 0.12, "S", cex = 1.1, font = 2)

# Choose a point P on sphere at angle φ from S
phi <- 50 * pi/180  # Angular distance from S
P_x <- R * sin(phi)
P_y <- R * cos(phi)
points(P_x, P_y, pch = 19, cex = 1.2, col = "#E63946")
text(P_x + 0.1, P_y + 0.08, "P", cex = 1.1, font = 2, col = "#E63946")

# Draw radius to P
lines(c(0, P_x), c(0, P_y), lwd = 1.5, lty = 2, col = "gray50")

# Project P to tangent plane (perpendicular projection gives P')
# In Lambert: P' is at distance ρ = 2R·sin(φ/2) from S
rho <- 2 * R * sin(phi/2)
Pprime_x <- rho
Pprime_y <- R
points(Pprime_x, Pprime_y, pch = 19, cex = 1.2, col = "#457B9D")
text(Pprime_x + 0.1, Pprime_y + 0.12, "P'", cex = 1.1, font = 2, col = "#457B9D")

# Draw chord from S to P (this has length 2R·sin(φ/2))
lines(c(0, P_x), c(R, P_y), lwd = 2.5, col = "#E63946")

# Draw projected distance on tangent plane
lines(c(0, Pprime_x), c(R, R), lwd = 2.5, col = "#457B9D")

# Draw arc showing angle φ at center
arc_r <- 0.35
arc_theta <- seq(pi/2, pi/2 - phi, length.out = 30)
lines(arc_r * cos(arc_theta), arc_r * sin(arc_theta), lwd = 2, col = "#2A9D8F")
text(arc_r * 1.6, arc_r * 1.2, expression(italic(c)), cex = 1.2, col = "#2A9D8F")

# Label the chord d
mid_chord_x <- (0 + P_x)/2 - 0.12
mid_chord_y <- (R + P_y)/2 + 0.08
text(mid_chord_x, mid_chord_y, expression(italic(d)), cex = 1.2, col = "#E63946")

# Label the projected distance ρ
text(Pprime_x/2, R + 0.15, expression(rho), cex = 1.2, col = "#457B9D")

# Add formula box
rect(0.8, -1.3, 2.15, -0.7, col = "white", border = "gray70")
text(1.47, -0.85, expression(italic(d) == 2*italic(R)*sin(italic(c)/2)), cex = 1.1)
text(1.47, -1.15, expression(rho == italic(d)), cex = 1.1)

# Add reference note
text(1.47, -1.45, "Snyder (1987, p. 182)", cex = 0.9, col = "gray50")
par(oldpar)

## ----lambert-area-preservation, echo=FALSE, fig.width=7, fig.height=3.5-------
# Show equal-area property with concentric rings
oldpar <- par(mfrow = c(1, 2), mar = c(2, 1, 3, 1), bg = "white", cex = 1, cex.main = 1)

R <- 1

# Left: Sphere view (orthographic, looking from above at tangent point)
plot(NULL, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "Sphere (view from above S)")

# Draw outer circle (equator as seen from S)
theta <- seq(0, 2*pi, length.out = 100)
lines(cos(theta), sin(theta), lwd = 2)

# Draw concentric latitude circles and shade bands
n_bands <- 5
cols <- c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51")

for (i in n_bands:1) {
  phi_outer <- (i / n_bands) * (pi/2)
  phi_inner <- ((i-1) / n_bands) * (pi/2)

  r_outer <- sin(phi_outer)
  r_inner <- sin(phi_inner)

  theta_seq <- seq(0, 2*pi, length.out = 100)
  if (i == 1) {
    polygon(r_outer * cos(theta_seq), r_outer * sin(theta_seq),
            col = adjustcolor(cols[i], 0.5), border = "gray40")
  } else {
    polygon(c(r_outer * cos(theta_seq), rev(r_inner * cos(theta_seq))),
            c(r_outer * sin(theta_seq), rev(r_inner * sin(theta_seq))),
            col = adjustcolor(cols[i], 0.5), border = "gray40")
  }
}
points(0, 0, pch = 19, cex = 1.2)
text(0, -1.2, "Equal-area bands on sphere", cex = 0.9)

# Right: Lambert projection
plot(NULL, xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "Lambert Projection")

for (i in n_bands:1) {
  phi_outer <- (i / n_bands) * (pi/2)
  phi_inner <- ((i-1) / n_bands) * (pi/2)

  # Lambert projection: r = 2*sin(phi/2)
  r_outer <- 2 * sin(phi_outer / 2)
  r_inner <- 2 * sin(phi_inner / 2)

  theta_seq <- seq(0, 2*pi, length.out = 100)
  if (i == 1) {
    polygon(r_outer * cos(theta_seq), r_outer * sin(theta_seq),
            col = adjustcolor(cols[i], 0.5), border = "gray40")
  } else {
    polygon(c(r_outer * cos(theta_seq), rev(r_inner * cos(theta_seq))),
            c(r_outer * sin(theta_seq), rev(r_inner * sin(theta_seq))),
            col = adjustcolor(cols[i], 0.5), border = "gray40")
  }
}

r_eq <- 2 * sin(pi/4)
lines(r_eq * cos(theta), r_eq * sin(theta), lwd = 2, lty = 2, col = "gray50")

points(0, 0, pch = 19, cex = 1.2)
text(0, -1.35, "Same bands after projection\n(areas preserved)", cex = 0.9)
par(oldpar)

## ----icosahedron-projection, echo=FALSE, fig.width=7, fig.height=6------------
# Icosahedron Face Projection Geometry
# Reference: Snyder (1992), "An equal-area map projection for polyhedral globes", p. 10-12
# Reference: Coxeter (1973), "Regular Polytopes", p. 52-53
#
# Shows how points on the sphere project onto an icosahedral face tangent plane

oldpar <- par(mar = c(1, 1, 2, 1), bg = "white")
plot(NULL, xlim = c(-1.6, 1.6), ylim = c(-0.6, 1.8), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "Icosahedron Face Projection")

# Sphere cross-section (arc visible above the face)
R <- 1.2
theta_arc <- seq(20*pi/180, 160*pi/180, length.out = 50)
sphere_y_offset <- 0.3
lines(R * cos(theta_arc), R * sin(theta_arc) + sphere_y_offset - R,
      lwd = 2, col = "gray40")
text(-1.3, 1.1, "Sphere surface", cex = 0.9, col = "gray50", adj = 0)

# Draw triangular face (equilateral, base at bottom)
face_scale <- 1.3
v1 <- c(-face_scale * 0.866, -0.4)   # bottom-left vertex
v2 <- c(face_scale * 0.866, -0.4)    # bottom-right vertex
v3 <- c(0, face_scale * 1.0)          # top vertex

polygon(c(v1[1], v2[1], v3[1]), c(v1[2], v2[2], v3[2]),
        col = adjustcolor("gray90", 0.6), border = "gray30", lwd = 2)

# Label vertices
text(v1[1] - 0.12, v1[2] - 0.12, "V1", cex = 0.9, font = 2)
text(v2[1] + 0.12, v2[2] - 0.12, "V2", cex = 0.9, font = 2)
text(v3[1], v3[2] + 0.12, "V3", cex = 0.9, font = 2)

# Mark face center
face_center <- c((v1[1] + v2[1] + v3[1])/3, (v1[2] + v2[2] + v3[2])/3)
points(face_center[1], face_center[2], pch = 19, cex = 1)
text(face_center[1] - 0.2, face_center[2] - 0.08, "Face center", cex = 0.9, col = "gray50")

# A point P on sphere surface
P_sphere <- c(0.5, 0.95)
points(P_sphere[1], P_sphere[2], pch = 19, cex = 1.2, col = "#E63946")
text(P_sphere[1] + 0.12, P_sphere[2] + 0.08, "P", cex = 1.1, font = 2, col = "#E63946")

# Projected point P' on face
P_proj <- c(0.5, 0.45)
points(P_proj[1], P_proj[2], pch = 19, cex = 1.2, col = "#457B9D")
text(P_proj[1] + 0.12, P_proj[2] - 0.05, "P'", cex = 1.1, font = 2, col = "#457B9D")

# Draw projection line (dashed)
lines(c(P_sphere[1], P_proj[1]), c(P_sphere[2], P_proj[2]),
      lwd = 2, lty = 2, col = "#2A9D8F")

# Add text explanation
text(0, -0.55, "Each face covers ~1/20 of the sphere", cex = 0.9)
text(1.55, 0.7, "Tangent plane\n(icosa face)", cex = 0.9, col = "gray50", adj = 1)

# Reference note
text(0, -0.75, "Snyder (1992, p. 10-12)", cex = 0.9, col = "gray50")
par(oldpar)

## ----face-centers, echo=FALSE, fig.width=7, fig.height=3.5--------------------
library(sf)
library(ggplot2)

centers <- hexify_face_centers()

center_pts <- st_as_sf(
  data.frame(face = 0:19, lon = centers$lon * 180 / pi, lat = centers$lat * 180 / pi),
  coords = c("lon", "lat"), crs = 4326
)

ggplot() +
  geom_sf(data = hexify_world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = center_pts, color = "#E63946", size = 3) +
  geom_sf_text(data = center_pts, aes(label = face), nudge_y = 6, size = 5) +
  labs(title = "ISEA Icosahedron Face Centers",
       subtitle = "20 triangular faces, numbered 0-19") +
  theme_minimal(base_size = FIG_BASE_SIZE) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## ----cell-counts--------------------------------------------------------------
cat("Resolution  Aperture 3    Aperture 4    Aperture 7\n")
cat("---------  ----------    ----------    ----------\n")
for (res in 0:8) {
  cells_ap3 <- 10 * 3^res + 2
  cells_ap4 <- 10 * 4^res + 2
  cells_ap7 <- 10 * 7^res + 2
  cat(sprintf("    %d      %10s    %10s    %10s\n",
              res,
              format(cells_ap3, big.mark = ","),
              format(cells_ap4, big.mark = ","),
              format(cells_ap7, big.mark = ",")))
}

## ----orientation-classes, echo=FALSE, fig.width=7, fig.height=4---------------
oldpar <- par(mfrow = c(1, 3), mar = c(1, 0.2, 3, 0.2), bg = "white", cex = 1, cex.main = 1)

hex_v <- function(cx, cy, r, rotation = 0) {
  angles <- seq(rotation, 2*pi + rotation, length.out = 7)
  list(x = cx + r * cos(angles), y = cy + r * sin(angles))
}

# Class I: Flat-top
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class I\n(flat-top, 0°)")
h <- hex_v(0, 0, 1.2, rotation = 0)
polygon(h$x, h$y, col = adjustcolor("#457B9D", 0.4), border = "gray30", lwd = 2)
lines(h$x[1:2], h$y[1:2], col = "#E63946", lwd = 3)
text(0, -1.35, "Ap. 4 (all res)\nAp. 3 (even res)", cex = 0.9)

# Class II: Pointy-top (30° rotation)
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class II\n(pointy-top, 30°)")
h <- hex_v(0, 0, 1.2, rotation = pi/6)
polygon(h$x, h$y, col = adjustcolor("#2A9D8F", 0.4), border = "gray30", lwd = 2)
points(h$x[1], h$y[1], pch = 19, cex = 1.5, col = "#E63946")
text(0, -1.35, "Ap. 3 (odd res)", cex = 0.9)

# Class III: Skewed
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class III\n(skewed, ~19.1°)")
h_parent <- hex_v(0, 0, 1.2, rotation = 0)
polygon(h_parent$x, h_parent$y, col = NA, border = "gray50", lwd = 2, lty = 2)
# Correct rotation: arctan(sqrt(3/7))
rot_angle <- atan(sqrt(3/7))
h_child <- hex_v(0, 0, 0.8, rotation = rot_angle)
polygon(h_child$x, h_child$y, col = adjustcolor("#E9C46A", 0.4),
        border = "gray30", lwd = 2)
arc_r <- 0.5
arc_theta <- seq(0, rot_angle, length.out = 20)
lines(arc_r * cos(arc_theta), arc_r * sin(arc_theta), col = "#E63946", lwd = 2)
text(0.6, 0.15, "19.1°", cex = 0.9, col = "#E63946")
text(0, -1.35, "Ap. 7 (all res)\nRotation accumulates", cex = 0.9)
par(oldpar)

## ----pentagon-locations, echo=FALSE, fig.width=7, fig.height=3.5--------------
pentagon_coords <- data.frame(
  type = c("Pole", "Pole", rep("Upper ring", 5), rep("Lower ring", 5)),
  lon = c(0, 0, 0, 72, 144, 216, 288, 36, 108, 180, 252, 324),
  lat = c(90, -90, rep(26.57, 5), rep(-26.57, 5))
)

pentagon_pts <- st_as_sf(pentagon_coords, coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_sf(data = hexify_world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = pentagon_pts, aes(color = type), size = 4) +
  scale_color_manual(values = c("Pole" = "#E63946",
                                 "Upper ring" = "#457B9D",
                                 "Lower ring" = "#2A9D8F")) +
  labs(title = "Pentagon Cell Locations",
       subtitle = "12 pentagonal cells at icosahedron vertices\n(area = 5/6 of hexagons)",
       color = "Location") +
  theme_minimal(base_size = FIG_BASE_SIZE) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## ----z7-example---------------------------------------------------------------
# Z7 index encoding for aperture 7
g7 <- hex_grid(resolution = 4, aperture = 7)
cell <- lonlat_to_cell(16.37, 48.21, g7)
idx <- cell_to_index(cell, g7)
cat(sprintf("Cell %d -> Z7 index: %s\n", cell, idx))
cat(sprintf("  Leading field: %s (quad %d), Digits: %s\n",
            substr(idx, 1, 2), as.integer(substr(idx, 1, 2)) %% 12L,
            substr(idx, 3, nchar(idx))))

# Hierarchical property: parent is obtained by dropping the last digit
parent_idx <- substr(idx, 1, nchar(idx) - 1)
parent_info <- hexify_index_to_cell(parent_idx, 7, "z7")
cat(sprintf("  Parent index: %s (face %d, i=%d, j=%d)\n",
            parent_idx, parent_info$face,
            as.integer(parent_info$i), as.integer(parent_info$j)))

# A valid Z7 index decodes and re-encodes without changing
decoded <- hexify_index_to_cell(idx, 7, "z7")
stopifnot(identical(
  hexify_cell_to_index(decoded$face, decoded$i, decoded$j,
                       decoded$resolution, 7, "z7"),
  idx
))

## ----z3-example---------------------------------------------------------------
# Z3 index encoding for aperture 3
g3 <- hex_grid(resolution = 8, aperture = 3)
cell <- lonlat_to_cell(16.37, 48.21, g3)
idx <- cell_to_index(cell, g3)
cat(sprintf("Cell %d -> Z3 index: %s\n", cell, idx))
cat(sprintf("  Base cell: %s, Digits: %s (%d digit pairs)\n",
            substr(idx, 1, 2), substr(idx, 3, nchar(idx)),
            (nchar(idx) - 2) / 2))

## ----zorder-example-----------------------------------------------------------
# Z-order index for aperture 4
g4 <- hex_grid(resolution = 8, aperture = 4)
cell <- lonlat_to_cell(16.37, 48.21, g4)
idx <- cell_to_index(cell, g4)
cat(sprintf("Aperture 4: Cell %d -> Z-order index: %s\n", cell, idx))

## ----h3-area-variation, fig.width=7, fig.height=4.375-------------------------
# Compare ISEA (constant area) vs H3 (variable area) at similar resolutions
lats <- seq(-85, 85, by = 5)
lons <- rep(10, length(lats))

# ISEA: aperture 7, resolution 6 (~130 km² cells)
g_isea <- hex_grid(resolution = 6, aperture = 7)
isea_cells <- lonlat_to_cell(lons, lats, g_isea)
isea_areas <- cell_area(isea_cells, g_isea)

# H3: resolution 4 (~1,770 km² cells — different scale, but shows the pattern)
g_h3 <- hex_grid(resolution = 4, type = "h3")
h3_cells <- lonlat_to_cell(lons, lats, g_h3)
h3_areas <- cell_area(h3_cells, g_h3)

# Normalize to show relative variation
isea_rel <- isea_areas / mean(isea_areas)
h3_rel <- h3_areas / mean(h3_areas)

oldpar <- par(mar = c(4, 4, 2, 1), bg = "white")
plot(lats, h3_rel, type = "l", col = "#E63946", lwd = 2.5,
     xlab = "Latitude (degrees)", ylab = "Relative cell area",
     main = "Cell Area Variation by Latitude",
     ylim = range(c(isea_rel, h3_rel)))
lines(lats, isea_rel, col = "#457B9D", lwd = 2.5)
abline(h = 1, lty = 2, col = "gray50")
legend("topright", legend = c("H3 (gnomonic)", "ISEA (equal-area)"),
       col = c("#E63946", "#457B9D"), lwd = 2.5, bty = "n")
par(oldpar)

## ----resolution-mapping-------------------------------------------------------
# H3 resolution table: compare H3 and ISEA aperture-7 cell areas
cat("H3 Res  Avg Area (km²)  ISEA Ap7 Equivalent\n")
cat("------  --------------  -------------------\n")
for (h3_res in 0:8) {
  g_h3 <- hex_grid(resolution = h3_res, type = "h3")
  h3_area <- g_h3@area_km2

  # Find closest ISEA ap7 resolution by brute force
  best_res <- 0
  best_diff <- Inf
  for (r in 0:15) {
    g_test <- hex_grid(resolution = r, aperture = 7)
    d <- abs(log(g_test@area_km2) - log(h3_area))
    if (d < best_diff) { best_diff <- d; best_res <- r }
  }
  g_isea <- hex_grid(resolution = best_res, aperture = 7)
  cat(sprintf("   %2d   %14.1f  res %d (%.1f km²)\n",
              h3_res, h3_area, best_res, g_isea@area_km2))
}

## ----round-trip---------------------------------------------------------------
original_lon <- 16.37
original_lat <- 48.21

cat(sprintf("Original: (%.4f, %.4f)\n\n", original_lon, original_lat))

for (ap in c(3, 4, 7)) {
  res <- if (ap == 7) 6 else 10
  grid <- hex_grid(resolution = res, aperture = ap)
  cell_id <- lonlat_to_cell(original_lon, original_lat, grid)
  recovered <- cell_to_lonlat(cell_id, grid)
  error_km <- sqrt((recovered$lon - original_lon)^2 +
                   (recovered$lat - original_lat)^2) * 111
  cat(sprintf("Aperture %d (res %2d): cell %d -> (%.4f, %.4f), ~%.1f km from center\n",
              ap, res, cell_id,
              recovered$lon, recovered$lat, error_km))
}

