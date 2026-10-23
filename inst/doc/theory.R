## ----setup, include = FALSE---------------------------------------------------
source("_common.R")
library(hexify)

## ----eval = FALSE-------------------------------------------------------------
# # Install from CRAN
# install.packages("hexify")
# 
# # Or install development version from GitHub
# # install.packages("pak")
# pak::pak("gcol33/hexify")

## ----latlon-area-distortion, echo=FALSE, fig.width=7, fig.height=3.8----------
# A lat-lon cell spans a fixed lon x lat interval, so its area scales with
# cos(latitude); an ISEA cell has the same area everywhere.
oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), bg = "white")

# Left: orthographic hemisphere with a 15-degree graticule; the same
# lon x lat cell is shaded at the equator and at 60 degrees north
ortho_x <- function(lon, lat) cos(lat * pi/180) * sin(lon * pi/180)
ortho_y <- function(lon, lat) sin(lat * pi/180) + 0 * lon

par(mar = c(1, 1, 2, 1))
plot(NULL, xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Lat-lon grid")
theta <- seq(0, 2*pi, length.out = 200)
lines(cos(theta), sin(theta), lwd = 2, col = GREY)
for (lat in seq(-75, 75, by = 15)) {
  lon <- seq(-90, 90, length.out = 100)
  lines(ortho_x(lon, lat), ortho_y(lon, lat), col = "gray70")
}
for (lon in seq(-90, 90, by = 15)) {
  lat <- seq(-90, 90, length.out = 100)
  lines(ortho_x(lon, lat), ortho_y(lon, lat), col = "gray70")
}
shade_cell <- function(lon0, lon1, lat0, lat1) {
  edge <- rbind(
    cbind(seq(lon0, lon1, length.out = 25), lat0),
    cbind(lon1, seq(lat0, lat1, length.out = 25)),
    cbind(seq(lon1, lon0, length.out = 25), lat1),
    cbind(lon0, seq(lat1, lat0, length.out = 25)))
  polygon(ortho_x(edge[, 1], edge[, 2]), ortho_y(edge[, 1], edge[, 2]),
          col = adjustcolor("#E63946", 0.5), border = "#E63946", lwd = 1.5)
}
shade_cell(0, 15, 0, 15)
shade_cell(0, 15, 60, 75)

# Right: relative area of a 1-degree cell by latitude
par(mar = c(4, 4, 2, 1))
lat <- seq(-85, 85, by = 1)
plot(lat, cos(lat * pi/180), type = "l", col = "#E63946", lwd = 2.5,
     xlab = "Latitude (degrees)", ylab = "Relative cell area",
     main = "Cell area vs latitude", ylim = c(0, 1.05))
lines(lat, rep(1, length(lat)), col = "#457B9D", lwd = 2.5)
legend("bottom", legend = c("Lat-lon grid", "ISEA grid"),
       col = c("#E63946", "#457B9D"), lwd = 2.5, bty = "n")
par(oldpar)

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
lines(R * cos(theta_circle), R * sin(theta_circle), lwd = 2, col = GREY)

# Draw tangent plane at S (top of sphere)
lines(c(-1.4, 2.1), c(R, R), lwd = 2, col = GREY)
text(-1.3, R + 0.1, "Tangent plane", adj = 0, cex = 0.9, col = GREY)

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

# Draw the radii O-S and O-P that subtend the angle c
lines(c(0, 0), c(0, R), lwd = 1.5, lty = 2, col = GREY)
lines(c(0, P_x), c(0, P_y), lwd = 1.5, lty = 2, col = GREY)

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

# Rotate the chord about S into the tangent plane, so |SP| = |SP'| = rho
swing <- seq(-phi/2, 0, length.out = 60)
lines(rho * cos(swing), R + rho * sin(swing), lwd = 1.5, lty = 3, col = GREY)

# Draw arc showing angle c at center, labelled on the arc's bisector
arc_r <- 0.35
arc_theta <- seq(pi/2, pi/2 - phi, length.out = 30)
lines(arc_r * cos(arc_theta), arc_r * sin(arc_theta), lwd = 2, col = "#2A9D8F")
bisector <- pi/2 - phi/2
text(0.52 * cos(bisector), 0.52 * sin(bisector), expression(italic(c)),
     cex = 1.2, col = "#2A9D8F")

# Label the chord d on its inner side, between the chord and the radius O-P
chord_mid <- c(P_x, R + P_y) / 2
chord_normal <- c(R - P_y, P_x) / sqrt((R - P_y)^2 + P_x^2)
text(chord_mid[1] - 0.14 * chord_normal[1], chord_mid[2] - 0.14 * chord_normal[2],
     expression(italic(d)), cex = 1.2, col = "#E63946")

# Label the projected distance ρ
text(Pprime_x/2, R + 0.15, expression(rho), cex = 1.2, col = "#457B9D")

# Add formula box
rect(0.8, -1.3, 2.15, -0.7, col = "white", border = "gray70")
text(1.47, -0.85, expression(italic(d) == 2*italic(R)*sin(italic(c)/2)), cex = 1.1)
text(1.47, -1.15, expression(rho == italic(d)), cex = 1.1)

# Add reference note
text(1.47, -1.45, "Snyder (1987, p. 182)", cex = 0.9, col = GREY)
par(oldpar)

## ----lambert-area-preservation, echo=FALSE, fig.width=7, fig.height=3.5-------
# Five bands of equal area on the hemisphere around S. A spherical cap of
# angular radius phi has area 2*pi*R^2*(1 - cos(phi)), so equal-area bands
# have cos(phi_i) = 1 - i/n, and the Lambert radius 2*sin(phi_i/2) then has
# squared length 2*i/n: the projected annuli are equal-area as well.
oldpar <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1), bg = "white", cex = 1, cex.main = 1)

n_bands <- 5
phi_bands <- acos(1 - (0:n_bands) / n_bands)
cols <- c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51")
theta <- seq(0, 2*pi, length.out = 200)

draw_bands <- function(radius, main, caption, lim) {
  plot(NULL, xlim = c(-lim, lim), ylim = c(-lim - 0.45, lim), asp = 1,
       axes = FALSE, xlab = "", ylab = "", main = main)
  for (i in seq_len(n_bands)) {
    r_outer <- radius(phi_bands[i + 1])
    r_inner <- radius(phi_bands[i])
    polygon(c(r_outer * cos(theta), rev(r_inner * cos(theta))),
            c(r_outer * sin(theta), rev(r_inner * sin(theta))),
            col = adjustcolor(cols[i], 0.5), border = NA)
  }
  for (r in radius(phi_bands[-1])) lines(r * cos(theta), r * sin(theta), col = GREY)
  points(0, 0, pch = 19, cex = 1.2)
  text(0, -lim - 0.25, caption, cex = 0.9)
}

# Left: the hemisphere seen from above S (orthographic, r = R sin(phi))
draw_bands(function(phi) sin(phi), "Sphere (view from above S)",
           "Five equal-area bands", lim = 1.1)

# Right: the same bands after Lambert projection, r = 2R sin(phi/2)
draw_bands(function(phi) 2 * sin(phi / 2), "Lambert projection",
           "Same bands, areas preserved", lim = 1.5)
par(oldpar)

## ----icosahedron-projection, echo=FALSE, fig.width=7, fig.height=4------------
# Icosahedron face geometry, Snyder (1992, p. 10-14), unit sphere.
# The icosahedron has the same surface area as the sphere, so its face plane
# lies at the inradius R1 = 0.9103832815 from the centre and its vertices at
# angular distance g = 37.37736814 degrees from the face centre. Along the
# azimuth of a vertex the radial distance reduces to
# rho = 2 R1 f sin(z/2) with f = tan(g) / (2 sin(g/2)).
g  <- 37.37736814 * pi/180
R1 <- 0.9103832815
rv <- R1 * tan(g)
f  <- tan(g) / (2 * sin(g/2))
red <- "#E63946"; blue <- "#457B9D"; teal <- "#2A9D8F"
oldpar <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5), bg = "white",
              cex = 1, cex.main = 1)

# Left: cross-section through the centre O, the face centre C and a vertex V,
# showing the upper part of the sphere. In this plane the face runs from V on
# one side to the midpoint of the opposite edge (at rv/2 from C) on the other.
plot(NULL, xlim = c(-0.9, 0.9), ylim = c(-0.2, 1.25), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Cross-section through O, C, V")
theta <- seq(0, 2*pi, length.out = 200)
lines(cos(theta), sin(theta), lwd = 2, col = GREY)

# Face plane and the planar face itself
lines(c(-0.9, 0.9), c(R1, R1), lty = 3, col = "gray60")
lines(c(-rv/2, rv), c(R1, R1), lwd = 4)

# Ray from O through the sphere vertex to the planar vertex V
lines(c(0, rv), c(0, R1), lty = 2, col = GREY)
points(sin(g), cos(g), pch = 19, cex = 0.9, col = GREY)
points(rv, R1, pch = 19, cex = 1.1)
text(rv + 0.06, R1 + 0.11, "V", font = 2, cex = 1)

# Centre, face centre and the face inradius R1
points(0, 0, pch = 19, cex = 1.1)
text(0.09, -0.08, "O", font = 2, cex = 1, adj = 0)
lines(c(0, 0), c(0, 1), col = GREY)
points(0, 1, pch = 19, cex = 1.1)
text(0, 1.13, "C", font = 2, cex = 1)
text(-0.06, 0.72, expression(italic(R)[1]), adj = 1, cex = 0.9, col = GREY)

# Angle g from C to V, drawn at O
arc <- seq(pi/2, pi/2 - g, length.out = 30)
lines(0.3 * cos(arc), 0.3 * sin(arc), lwd = 2, col = GREY)
text(0.42 * cos(pi/2 - g/2), 0.42 * sin(pi/2 - g/2), expression(italic(g)), cex = 1)

# A point P at angular distance z from C on the azimuth of V, and its image P'
z <- 30 * pi/180
rho <- 2 * R1 * f * sin(z/2)
lines(c(0, sin(z)), c(0, cos(z)), lty = 2, col = red)
points(sin(z), cos(z), pch = 19, cex = 1.2, col = red)
text(sin(z) - 0.14, cos(z) - 0.02, "P", font = 2, cex = 1, col = red)
points(rho, R1, pch = 19, cex = 1.2, col = blue)
text(rho - 0.1, R1 + 0.11, "P'", font = 2, cex = 1, col = blue)
lines(c(sin(z), rho), c(cos(z), R1), lwd = 2, lty = 3, col = teal)
arc <- seq(pi/2, pi/2 - z, length.out = 20)
lines(0.62 * cos(arc), 0.62 * sin(arc), lwd = 2, col = red)
text(0.74 * cos(pi/2 - z/2), 0.74 * sin(pi/2 - z/2), expression(italic(z)),
     cex = 1, col = red)

# Right: the planar face seen from above, azimuth measured from the vertex ray
plot(NULL, xlim = c(-0.85, 0.85), ylim = c(-0.65, 0.9), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Planar face, view from above")
va <- pi/2 + c(0, 2, 4) * pi/3
vx <- rv * cos(va); vy <- rv * sin(va)
polygon(c(0, vx[1], vx[3]), c(0, vy[1], vy[3]),
        col = adjustcolor(teal, 0.15), border = NA)
polygon(vx, vy, col = NA, border = "black", lwd = 2)
for (k in 1:3) lines(c(0, vx[k]), c(0, vy[k]), lty = 2, col = GREY)
points(0, 0, pch = 19, cex = 1.1)
text(-0.08, -0.06, "C", font = 2, cex = 1)
text(vx[1], vy[1] + 0.09, "V", font = 2, cex = 1)

az <- 45 * pi/180
rq <- 0.33
points(rq * cos(pi/2 - az), rq * sin(pi/2 - az), pch = 19, cex = 1.2, col = blue)
text(rq * cos(pi/2 - az) + 0.1, rq * sin(pi/2 - az) - 0.02, "P'",
     font = 2, cex = 1, col = blue)
lines(c(0, rq * cos(pi/2 - az)), c(0, rq * sin(pi/2 - az)), col = blue)
arc <- seq(pi/2, pi/2 - az, length.out = 20)
lines(0.17 * cos(arc), 0.17 * sin(arc), lwd = 2, col = blue)
text(0.3 * cos(pi/2 - az/2), 0.3 * sin(pi/2 - az/2), "Az'", cex = 0.9, col = blue)
text(0.3, -0.5, "120° sector", cex = 0.8, col = teal)
par(oldpar)

## ----face-centers, echo=FALSE, message=FALSE, warning=FALSE, fig.width=7, fig.height=3.5----
library(sf)
library(ggplot2)

centers <- hexify_face_centers()

center_pts <- st_as_sf(
  data.frame(face = 0:19, lon = centers$lon * 180 / pi, lat = centers$lat * 180 / pi),
  coords = c("lon", "lat"), crs = 4326
)

ggplot() +
  geom_sf(data = hexify_world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = center_pts, color = "#E63946", size = 2.5) +
  geom_sf_text(data = center_pts, aes(label = face), fun.geometry = sf::st_centroid,
               nudge_y = 9, size = 4.5) +
  coord_sf(expand = FALSE) +
  labs(title = "ISEA Icosahedron Face Centers",
       subtitle = "20 triangular faces, numbered 0-19", x = NULL, y = NULL) +
  theme_minimal(base_size = FIG_BASE_SIZE) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## ----face-distortion, echo=FALSE, fig.width=7, fig.height=5.2-----------------
# Maximum angular deformation omega over one icosahedron face, computed from
# the forward projection's Jacobian (Snyder, 1992, eq. A17-A23). The face is
# drawn with the pole at its center and a 5-degree graticule.
G  <- 36 * pi/180
El <- 37.37736814 * pi/180
R1 <- 0.9103832815
cot30 <- 1 / tan(30 * pi/180)

fwd_sector <- function(Az, z) {
  H   <- acos(sin(Az) * sin(G) * cos(El) - cos(Az) * cos(G))
  AG  <- Az + G + H - pi
  Azp <- atan2(2 * AG, R1^2 * tan(El)^2 - 2 * AG * cot30)
  q   <- atan2(tan(El), cos(Az) + sin(Az) * cot30)
  dp  <- R1 * tan(El) / (cos(Azp) + sin(Azp) * cot30)
  f   <- dp / (2 * R1 * sin(q / 2))
  rho <- 2 * R1 * f * sin(z / 2)
  list(x = rho * sin(Azp), y = rho * cos(Azp))
}
fwd <- function(Az, z) {
  Az <- Az %% (2 * pi)
  k  <- floor(Az / (2 * pi / 3))
  p  <- fwd_sector(Az - k * 2 * pi / 3, z)
  a  <- k * 2 * pi / 3
  list(x = p$x * cos(a) + p$y * sin(a), y = p$y * cos(a) - p$x * sin(a))
}

# Tissot scale factors from a numerical Jacobian in the orthonormal frame;
# for an equal-area projection a * b = 1 and sin(omega/2) = (a - b)/(a + b)
distortion <- function(Az, z, h = 1e-6) {
  p1 <- fwd(Az, z + h); p0 <- fwd(Az, z - h)
  q1 <- fwd(Az + h, z); q0 <- fwd(Az - h, z)
  e11 <- (p1$x - p0$x) / (2 * h);  e21 <- (p1$y - p0$y) / (2 * h)
  e12 <- (q1$x - q0$x) / (2 * h) / sin(z)
  e22 <- (q1$y - q0$y) / (2 * h) / sin(z)
  tt <- e11^2 + e21^2 + e12^2 + e22^2
  dd <- abs(e11 * e22 - e12 * e21)
  s  <- sqrt(pmax(tt^2 - 4 * dd^2, 0))
  a  <- sqrt((tt + s) / 2); b <- sqrt(pmax((tt - s) / 2, 0))
  2 * asin((a - b) / (a + b)) * 180/pi
}

# Inverse of the azimuth adjustment: plane (x, y) -> (Az, z). Az' determines
# the target spherical area A_G directly; Newton iteration recovers Az.
inverse_az <- function(x, y) {
  azp_full <- atan2(x, y) %% (2 * pi)
  k <- floor(azp_full / (2 * pi / 3))
  azp <- azp_full - k * 2 * pi / 3
  rho <- sqrt(x^2 + y^2)
  AG_t <- R1^2 * tan(El)^2 * sin(azp) / (2 * cos(azp) + 2 * cot30 * sin(azp))
  Az <- azp
  for (it in 1:25) {
    u  <- sin(Az) * sin(G) * cos(El) - cos(Az) * cos(G)
    H  <- acos(pmin(pmax(u, -1), 1))
    Fv <- (Az + G + H - pi) - AG_t
    du <- cos(Az) * sin(G) * cos(El) + sin(Az) * cos(G)
    Fp <- 1 - du / sqrt(pmax(1 - u^2, 1e-12))
    Az <- Az - Fv / Fp
  }
  q  <- atan2(tan(El), cos(Az) + sin(Az) * cot30)
  dp <- R1 * tan(El) / (cos(azp) + sin(azp) * cot30)
  f  <- dp / (2 * R1 * sin(q / 2))
  z  <- 2 * asin(pmin(rho / (2 * R1 * f), 1))
  list(Az = Az + k * 2 * pi / 3, z = z)
}

oldpar <- par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
rv <- R1 * tan(El)
tri_x <- rv * sin(c(0, 2, 4) * pi/3); tri_y <- rv * cos(c(0, 2, 4) * pi/3)
plot(NULL, xlim = c(-0.66, 0.95), ylim = c(-0.39, 0.75), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "Angular deformation on one face")

in_tri <- function(x, y, tol = 1e-9) {
  s <- rep(TRUE, length(x))
  for (k in 1:3) {
    j <- k %% 3 + 1
    cr <- (tri_x[j] - tri_x[k]) * (y - tri_y[k]) -
          (tri_y[j] - tri_y[k]) * (x - tri_x[k])
    s <- s & (cr <= tol)
  }
  s
}

nx <- 420; ny <- 380
gx <- seq(-0.62, 0.62, length.out = nx)
gy <- seq(-0.36, 0.71, length.out = ny)
px <- rep(gx, times = ny); py <- rep(gy, each = nx)
inside <- in_tri(px, py) & (px^2 + py^2) > 1e-4
om <- rep(NA_real_, nx * ny)
iv <- inverse_az(px[inside], py[inside])
om[inside] <- distortion(iv$Az, iv$z)
ramp <- colorRampPalette(c("#FFFFFF", "#E9C46A", "#F4A261", "#E76F51",
                           "#E63946"))(64)
image(gx, gy, matrix(om, nrow = nx), zlim = c(0, 17.5), col = ramp,
      add = TRUE, useRaster = TRUE)

# 5-degree parallels and 10-degree meridians, pole at the face center
for (z_deg in seq(5, 35, by = 5)) {
  az <- seq(0, 2 * pi, length.out = 721)
  p <- fwd(az, z_deg * pi/180)
  ok <- in_tri(p$x, p$y)
  lines(ifelse(ok, p$x, NA), ifelse(ok, p$y, NA), col = "gray40", lwd = 0.7)
}
for (az_deg in seq(0, 355, by = 10)) {
  z <- seq(2.5, 37.37, length.out = 120) * pi/180
  p <- fwd(rep(az_deg * pi/180, length(z)), z)
  ok <- in_tri(p$x, p$y)
  lines(ifelse(ok, p$x, NA), ifelse(ok, p$y, NA), col = "gray40", lwd = 0.7)
}
polygon(tri_x, tri_y, border = "black", lwd = 2)

bx0 <- 0.78; bx1 <- 0.85; by0 <- -0.30; by1 <- 0.60
yy <- seq(by0, by1, length.out = 65)
rect(bx0, yy[-65], bx1, yy[-1], col = ramp, border = NA)
rect(bx0, by0, bx1, by1, border = GREY)
for (v in c(0, 5, 10, 15)) {
  yv <- by0 + (by1 - by0) * v / 17.5
  lines(c(bx1, bx1 + 0.02), c(yv, yv), col = GREY)
  text(bx1 + 0.04, yv, sprintf("%d°", v), adj = 0, cex = 0.8)
}
text((bx0 + bx1) / 2, by1 + 0.08, expression(omega), cex = 1.1)
par(oldpar)

## ----ellipsoid-ogc------------------------------------------------------------
ogc <- hex_grid(resolution = 8, aperture = 3, ellipsoid = "WGS84",
                orientation = "ogc")
ogc

## ----unfolded-projection, fig.width=7, fig.height=3.9-------------------------
draw_path <- function(lon, lat, col, lwd) {
  p <- hexify_lonlat_to_plane(lon, lat)
  jump <- c(FALSE, sqrt(diff(p$plane_x)^2 + diff(p$plane_y)^2) > 0.2)
  p$plane_x[jump] <- NA
  lines(p$plane_x, p$plane_y, col = col, lwd = lwd)
}

oldpar <- par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
plot(NULL, xlim = c(-0.1, 5.6), ylim = c(-0.05, 2.65), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Unfolded ISEA Projection")

# face outlines: side-1 triangles around each projected face centroid
centers <- hexify_face_centers()
pc <- hexify_lonlat_to_plane(centers$lon * 180/pi, centers$lat * 180/pi)
s3 <- sqrt(3)
for (f in 1:20) {
  base <- floor(pc$plane_y[f] / (s3/2)) * s3/2
  up <- (pc$plane_y[f] - base) < s3/4
  yy <- if (up) base + c(0, 0, s3/2) else base + c(s3/2, s3/2, 0)
  polygon(pc$plane_x[f] + c(-0.5, 0.5, 0), yy, border = "gray60", lwd = 0.8)
}

# 15-degree graticule
for (lon in seq(-180, 165, by = 15))
  draw_path(rep(lon, 181), seq(-90, 90, by = 1), "gray75", 0.6)
for (lat in seq(-75, 75, by = 15))
  draw_path(seq(-180, 180, by = 0.5), rep(lat, 721), "gray75", 0.6)

# coastlines
xy <- st_coordinates(st_cast(st_geometry(hexify_world), "MULTIPOLYGON"))
for (ring in split(seq_len(nrow(xy)),
                   interaction(xy[, "L1"], xy[, "L2"], xy[, "L3"], drop = TRUE)))
  draw_path(xy[ring, 1], xy[ring, 2], GREY, 0.8)
par(oldpar)

## ----aperture-subdivision, echo=FALSE, fig.width=7, fig.height=3--------------
# One parent cell (dashed, flat-top, circumradius 1) and the child cells of
# the next resolution for each aperture. Child centres sit on the child
# lattice, whose spacing is the parent spacing sqrt(3) divided by sqrt(a).
hex_xy <- function(cx, cy, r, rot = 0) {
  a <- rot + (0:5) * pi/3
  list(x = cx + r * cos(a), y = cy + r * sin(a))
}
draw_children <- function(main, child_r, child_rot, ring_r, ring_start, col) {
  plot(NULL, xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), asp = 1,
       axes = FALSE, xlab = "", ylab = "", main = main)
  ring_angles <- ring_start + (0:5) * pi/3
  cx <- c(0, ring_r * cos(ring_angles))
  cy <- c(0, ring_r * sin(ring_angles))
  for (k in seq_along(cx)) {
    h <- hex_xy(cx[k], cy[k], child_r, child_rot)
    polygon(h$x, h$y, col = adjustcolor(col, if (k == 1) 0.6 else 0.25),
            border = GREY, lwd = 1.5)
  }
  p <- hex_xy(0, 0, 1, 0)
  polygon(p$x, p$y, col = NA, border = "black", lwd = 2.5, lty = 2)
}

oldpar <- par(mfrow = c(1, 3), mar = c(0.2, 0.2, 2.5, 0.2), bg = "white",
              cex = 1, cex.main = 1)
# Aperture 3: children rotated 30 degrees, centred on the parent's vertices
draw_children("Aperture 3", child_r = 1/sqrt(3), child_rot = pi/6,
              ring_r = 1, ring_start = 0, col = "#457B9D")
# Aperture 4: same orientation, centred on the parent's edge midpoints
draw_children("Aperture 4", child_r = 1/2, child_rot = 0,
              ring_r = sqrt(3)/2, ring_start = pi/6, col = "#2A9D8F")
# Aperture 7: rotated by arctan(sqrt(3)/5), a rosette around the centre child
rot7 <- atan(sqrt(3)/5)
draw_children("Aperture 7", child_r = 1/sqrt(7), child_rot = rot7,
              ring_r = sqrt(3/7), ring_start = pi/6 + rot7, col = "#E9C46A")
par(oldpar)

## ----parent-shares------------------------------------------------------------
g7 <- hex_grid(resolution = 3, aperture = 7)
ids <- grid_global(g7)$cell_id
table(lengths(get_parent(ids, g7, overlapping = TRUE)))

# A uniform density adds up to each parent's area, the twelve pentagons included
up <- hex_aggregate(ids, cell_area(ids, g7), g7)
range(up$value / cell_area(up$cell_id, hex_grid(resolution = 2, aperture = 7)))

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
oldpar <- par(mfrow = c(1, 3), mar = c(0.2, 0.2, 3, 0.2), bg = "white", cex = 1, cex.main = 1)

hex_v <- function(cx, cy, r, rotation = 0) {
  angles <- rotation + (0:6) * pi/3
  list(x = cx + r * cos(angles), y = cy + r * sin(angles))
}
lim <- c(-1.55, 1.55)
ylim <- c(-1.7, 1.2)

# Class I: flat-top, vertices at 0, 60, ..., 300 degrees; the top edge is marked
plot(NULL, xlim = lim, ylim = ylim, asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class I\n(flat-top, 0°)")
h <- hex_v(0, 0, 1, rotation = 0)
polygon(h$x, h$y, col = adjustcolor("#457B9D", 0.4), border = GREY, lwd = 2)
lines(h$x[2:3], h$y[2:3], col = "#E63946", lwd = 4)
text(0, -1.4, "Ap. 4 (all res)\nAp. 3 (even res)", cex = 0.9)

# Class II: pointy-top, the Class I hexagon rotated by 30 degrees; the top vertex is marked
plot(NULL, xlim = lim, ylim = ylim, asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class II\n(pointy-top, 30°)")
h <- hex_v(0, 0, 1, rotation = pi/6)
polygon(h$x, h$y, col = adjustcolor("#2A9D8F", 0.4), border = GREY, lwd = 2)
points(h$x[2], h$y[2], pch = 19, cex = 1.5, col = "#E63946")
text(0, -1.4, "Ap. 3 (odd res)", cex = 0.9)

# Class III: the Class I hexagon rotated by arctan(sqrt(3)/5), with the
# rotation angle drawn outside the hexagons between the two vertex rays
plot(NULL, xlim = c(-1.3, 1.8), ylim = ylim, asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class III\n(skewed, 19.1°)")
rot_angle <- atan(sqrt(3)/5)
h_parent <- hex_v(0, 0, 1, rotation = 0)
polygon(h_parent$x, h_parent$y, col = NA, border = GREY, lwd = 2, lty = 2)
h_child <- hex_v(0, 0, 1, rotation = rot_angle)
polygon(h_child$x, h_child$y, col = adjustcolor("#E9C46A", 0.4),
        border = GREY, lwd = 2)
lines(c(0, 1.25), c(0, 0), col = GREY, lty = 2)
lines(c(0, 1.25 * cos(rot_angle)), c(0, 1.25 * sin(rot_angle)), col = "#E63946")
arc_theta <- seq(0, rot_angle, length.out = 20)
lines(1.15 * cos(arc_theta), 1.15 * sin(arc_theta), col = "#E63946", lwd = 2)
text(1.05, 0.62, "19.1°", cex = 0.9, col = "#E63946")
text(0, -1.4, "Ap. 7 (all res)\nRotation accumulates", cex = 0.9)
par(oldpar)

## ----pentagon-neighbors, echo=FALSE, fig.width=7, fig.height=3.8--------------
# Local topology around a pentagon cell and an ordinary hexagon cell.
# All cells share the same edge length; around a pentagon the neighbor
# ring leaves angular gaps that close up on the curved sphere.
oldpar <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5), bg = "white")

ngon <- function(cx, cy, r, n, start) {
  a <- start + (0:(n - 1)) * 2 * pi / n
  cbind(cx + r * cos(a), cy + r * sin(a))
}

# Pentagon: five hexagons, one on each edge (edge length 1 throughout)
plot(NULL, xlim = c(-2.85, 2.85), ylim = c(-2.85, 2.85), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Pentagon: 5 neighbors")
apothem_p <- cos(pi/5) / (2 * sin(pi/5))
apothem_h <- sqrt(3) / 2
for (k in 0:4) {
  mid_ang <- pi/2 + (2 * k + 1) * pi/5
  cx <- (apothem_p + apothem_h) * cos(mid_ang)
  cy <- (apothem_p + apothem_h) * sin(mid_ang)
  h <- ngon(cx, cy, 1, 6, mid_ang + pi/6)
  polygon(h[, 1], h[, 2], col = adjustcolor("#457B9D", 0.25),
          border = GREY, lwd = 1.5)
}
p <- ngon(0, 0, 1 / (2 * sin(pi/5)), 5, pi/2)
polygon(p[, 1], p[, 2], col = adjustcolor("#E63946", 0.5),
        border = GREY, lwd = 2)
text(0, 0, "5", font = 2, cex = 1.2)

# Hexagon: six hexagons tile around it exactly
plot(NULL, xlim = c(-2.85, 2.85), ylim = c(-2.85, 2.85), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Hexagon: 6 neighbors")
for (k in 0:5) {
  ang <- pi/6 + k * pi/3
  h <- ngon(sqrt(3) * cos(ang), sqrt(3) * sin(ang), 1, 6, 0)
  polygon(h[, 1], h[, 2], col = adjustcolor("#457B9D", 0.25),
          border = GREY, lwd = 1.5)
}
h0 <- ngon(0, 0, 1, 6, 0)
polygon(h0[, 1], h0[, 2], col = adjustcolor("#2A9D8F", 0.5),
        border = GREY, lwd = 2)
text(0, 0, "6", font = 2, cex = 1.2)
par(oldpar)

## ----pentagon-locations, echo=FALSE, message=FALSE, warning=FALSE, fig.width=7, fig.height=3.5----
pentagon_coords <- data.frame(
  type = factor(c("Pole", "Pole", rep("Upper ring", 5), rep("Lower ring", 5)),
                levels = c("Pole", "Upper ring", "Lower ring")),
  lon = c(0, 0, 0, 72, 144, 216, 288, 36, 108, 180, 252, 324),
  lat = c(90, -90, rep(26.57, 5), rep(-26.57, 5))
)
# Wrap longitudes above 180 into the map's -180..180 range
pentagon_coords$lon <- pentagon_coords$lon - 360 * (pentagon_coords$lon > 180)

pentagon_pts <- st_as_sf(pentagon_coords, coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_sf(data = hexify_world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = pentagon_pts, aes(color = type), size = 4) +
  scale_color_manual(values = c("Pole" = "#E63946",
                                 "Upper ring" = "#457B9D",
                                 "Lower ring" = "#2A9D8F")) +
  labs(title = "Pentagon Cell Locations",
       subtitle = "12 pentagonal cells at icosahedron vertices\n(area = 5/6 of hexagons)",
       color = "Location", x = NULL, y = NULL) +
  theme_minimal(base_size = FIG_BASE_SIZE) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## ----triangle-quad-pairing, echo=FALSE, fig.width=7, fig.height=4.7-----------
# Unfolded icosahedron as the 12 quads: each rhombus pairs an upper and a
# lower triangle; quads 0 and 11 hold only the polar pentagon cell.
oldpar <- par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
plot(NULL, xlim = c(-1.2, 10.2), ylim = c(-2.4, 5.2), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "20 Triangles Paired into 12 Quads")

hw <- 1     # rhombus half-width
hh <- 1.4   # rhombus half-height
draw_quad <- function(cx, cy, q, tri_top, tri_bot, fill) {
  polygon(c(cx, cx + hw, cx, cx - hw), c(cy + hh, cy, cy - hh, cy),
          col = fill, border = GREY, lwd = 1.5)
  lines(c(cx - hw, cx + hw), c(cy, cy), col = GREY)
  text(cx, cy + 0.62, tri_top, cex = 0.75, col = GREY)
  text(cx, cy - 0.62, tri_bot, cex = 0.75, col = GREY)
  text(cx, cy, paste0("Q", q), font = 2, cex = 0.9)
}
for (q in 1:5)
  draw_quad(2 * (q - 1), 2, q, q - 1, q + 4, adjustcolor("#457B9D", 0.15))
for (q in 6:10)
  draw_quad(2 * (q - 6) + 1, 0.6, q, q + 4, q + 9, adjustcolor("#2A9D8F", 0.15))

pentagon <- function(cx, cy, r) {
  a <- pi/2 + (0:4) * 2 * pi/5
  cbind(cx + r * cos(a), cy + r * sin(a))
}
p0 <- pentagon(4, 4.4, 0.45)
polygon(p0[, 1], p0[, 2], col = adjustcolor("#E63946", 0.4),
        border = GREY, lwd = 1.5)
text(4.65, 4.4, "Q0: north pole", adj = 0, cex = 0.85)
p11 <- pentagon(5, -1.75, 0.45)
polygon(p11[, 1], p11[, 2], col = adjustcolor("#E63946", 0.4),
        border = GREY, lwd = 1.5)
text(5.65, -1.75, "Q11: south pole", adj = 0, cex = 0.85)

text(-1.2, -2.2, "gray: triangle index (0-19)", adj = 0,
     cex = 0.8, col = GREY)
par(oldpar)

## ----z7-example---------------------------------------------------------------
# Z7 index encoding for aperture 7
g7 <- hex_grid(resolution = 4, aperture = 7)
cell <- lonlat_to_cell(16.37, 48.21, g7)
idx <- cell_to_index(cell, g7)
cat(sprintf("Cell %s -> Z7 index: %s\n", as.character(cell), idx))
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

## ----z7-integer---------------------------------------------------------------
g7 <- hex_grid(resolution = 4, aperture = 7)
cell <- index_to_cell("033000", g7)
cell_to_index(cell, g7, form = "hex")
cell_to_index(cell, g7, form = "int")
cell_to_index(cell, g7, form = "monotonic")

## ----ogc-zirs-----------------------------------------------------------------
cell <- index_to_cell("E6-317-A", ogc, form = "textZIRS")
cell_to_lonlat(cell, ogc)
cell_to_index(cell, ogc, form = "uint64ZIRS")

## ----z3-example---------------------------------------------------------------
# Z3 index encoding for aperture 3
g3 <- hex_grid(resolution = 8, aperture = 3)
cell <- lonlat_to_cell(16.37, 48.21, g3)
idx <- cell_to_index(cell, g3)
cat(sprintf("Cell %s -> Z3 index: %s\n", as.character(cell), idx))
cat(sprintf("  Base cell: %s, Digits: %s (%d digit pairs)\n",
            substr(idx, 1, 2), substr(idx, 3, nchar(idx)),
            (nchar(idx) - 2) / 2))

## ----zorder-example-----------------------------------------------------------
# Z-order index for aperture 4
g4 <- hex_grid(resolution = 8, aperture = 4)
cell <- lonlat_to_cell(16.37, 48.21, g4)
idx <- cell_to_index(cell, g4)
cat(sprintf("Aperture 4: Cell %s -> Z-order index: %s\n", as.character(cell), idx))

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
abline(h = 1, lty = 2, col = GREY)
legend("bottomleft", legend = c("H3 (gnomonic)", "ISEA (equal-area)"),
       col = c("#E63946", "#457B9D"), lwd = 2.5, bty = "n")
par(oldpar)

## ----resolution-mapping-------------------------------------------------------
# H3 resolution table: compare H3 and ISEA aperture-7 cell areas
cat("H3 Res  Avg Area (km²)  ISEA Ap7 Equivalent\n")
cat("------  --------------  -------------------\n")
for (h3_res in 0:8) {
  g_h3 <- hex_grid(resolution = h3_res, type = "h3")
  h3_area <- summary(g_h3)$area_km2

  # Find closest ISEA ap7 resolution by brute force
  best_res <- 0
  best_diff <- Inf
  for (r in 0:15) {
    g_test <- hex_grid(resolution = r, aperture = 7)
    d <- abs(log(summary(g_test)$area_km2) - log(h3_area))
    if (d < best_diff) { best_diff <- d; best_res <- r }
  }
  g_isea <- hex_grid(resolution = best_res, aperture = 7)
  cat(sprintf("   %2d   %14.1f  res %d (%.1f km²)\n",
              h3_res, h3_area, best_res, summary(g_isea)$area_km2))
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
  cat(sprintf("Aperture %d (res %2d): cell %s -> (%.4f, %.4f), ~%.1f km from center\n",
              ap, res, as.character(cell_id),
              recovered$lon, recovered$lat, error_km))
}

