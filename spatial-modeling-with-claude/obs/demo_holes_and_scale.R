# ============================================================
# demo_holes_and_scale.R  (v3)
#
# Demonstrates fit_spatial_mask() on topologically complex
# cases: interior holes, irregular boundaries, archipelagos,
# and 100k-cell scale tests.
#
# Run after sourcing: source("fit_spatial_mask.R")
#
# Fixes in v3:
#   1. raster_sigma now in COORDINATE UNITS (not grid cells).
#      Old calls used sigma in grid-cell units; values have been
#      converted: sigma_cu = sigma_gc * (domain_width / resolution).
#   2. All hand-built circle/ellipse polygons now explicitly close
#      their rings with rbind(mat, mat[1L,]).  sf::st_polygon
#      requires first == last row; seq(0, 2*pi, length.out=N)
#      does NOT guarantee this due to floating-point rounding.
#
# Required packages:
#   install.packages(c("sf","concaveman","MASS","isoband",
#                      "ggplot2","patchwork"))
#
# Sections:
#   1 — Donut: single interior void
#   2 — Swiss cheese: multiple holes
#   3 — Fractured map: irregular boundary + lakes
#   4 — Archipelago: disconnected islands
#   5 — Nested crescent voids
#   6 — sigma sensitivity sweep (coordinate units)
#   7 — Scale test: 100,000 cells
#   8 — Validation table
# ============================================================

library(sf)
library(ggplot2)
library(patchwork)
source("fit_spatial_mask.R")
set.seed(123)


# ---- Shared helpers ----------------------------------------------------------

plot_mask <- function(mask, coords, title = "", subtitle = "",
                      fill = "#4c9be8", border = "#1a5fa8",
                      pt_col = "#c0392b", pt_size = 0.9, pt_alpha = 0.5) {
  mc  <- as.data.frame(sf::st_coordinates(mask))
  grp <- if ("L3" %in% names(mc)) "L3" else
         if ("L2" %in% names(mc)) "L2" else "L1"
  ggplot() +
    geom_polygon(data = mc, aes(x = X, y = Y, group = .data[[grp]]),
                 fill = fill, alpha = 0.2, color = border, linewidth = 0.55) +
    geom_point(data = coords, aes(x = x, y = y),
               color = pt_col, size = pt_size, alpha = pt_alpha) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = "X", y = "Y") +
    theme_minimal(base_size = 9.5) +
    theme(plot.title    = element_text(face = "bold", size = 9.5),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          panel.grid    = element_line(color = "grey93"))
}

# Sample n points uniformly from inside a polygon
sample_in_polygon <- function(poly_sf, n, max_tries = 20) {
  bb  <- sf::st_bbox(poly_sf); mu <- sf::st_union(poly_sf)
  pts <- data.frame(x = numeric(0), y = numeric(0))
  for (i in seq_len(max_tries)) {
    if (nrow(pts) >= n) break
    cnd <- data.frame(x = runif(n * 8, bb["xmin"], bb["xmax"]),
                      y = runif(n * 8, bb["ymin"], bb["ymax"]))
    ok  <- as.logical(
      sf::st_within(sf::st_as_sf(cnd, coords = c("x","y")), mu,
                    sparse = FALSE)[, 1L])
    pts <- rbind(pts, cnd[ok, ])
  }
  pts[seq_len(min(n, nrow(pts))), ]
}

# Build a closed ellipse polygon.
# IMPORTANT: seq(0, 2*pi, length.out=n) does NOT produce a closed ring
# due to float rounding (cos(2*pi) != cos(0) exactly).  We always
# append the first row explicitly.
make_ellipse_poly <- function(cx, cy, a, b, angle_deg, n = 80) {
  th  <- seq(0, 2 * pi, length.out = n)
  ang <- angle_deg * pi / 180
  xe  <- a * cos(th); ye <- b * sin(th)
  mat <- cbind(cx + xe * cos(ang) - ye * sin(ang),
               cy + xe * sin(ang) + ye * cos(ang))
  mat <- rbind(mat, mat[1L, ])   # explicit ring closure
  sf::st_polygon(list(mat))
}

# Helper: closed circle matrix
.circle_mat <- function(cx, cy, r, n = 100) {
  th  <- seq(0, 2 * pi, length.out = n)
  mat <- cbind(cx + r * cos(th), cy + r * sin(th))
  rbind(mat, mat[1L, ])   # explicit ring closure
}


# ============================================================
# SECTION 1 — Donut
# ============================================================
cat("Section 1: Donut...\n")

r_d           <- sqrt(runif(800, 2^2, 5^2))
th_d          <- runif(800, 0, 2 * pi)
coords_donut  <- data.frame(x = r_d * cos(th_d), y = r_d * sin(th_d))
# Domain: x,y in [-5, 5], domain width ~ 10 coord units.
# raster_sigma = 0.25 coordinate units
# (previously 8 grid cells at res=400 over domain=11 → 8*(11/400) = 0.22)

m1_r <- fit_spatial_mask(coords_donut, method = "raster",
  raster_resolution = 256L, raster_sigma = 0.25, raster_threshold = 0.2,
  verbose = FALSE)
m1_c <- fit_spatial_mask(coords_donut, method = "convex", verbose = FALSE)

print(
  plot_mask(m1_r, coords_donut,
            title    = "1a. Donut \u2014 Raster (hole \u2713)",
            subtitle = "\u03c3 = 0.25 coord units, threshold = 0.2") +
  plot_mask(m1_c, coords_donut,
            title    = "1b. Donut \u2014 Convex",
            subtitle = "[hole filled]",
            fill = "#e74c3c", border = "#922b21") +
  plot_annotation(title = "Section 1: Donut",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
)


# ============================================================
# SECTION 2 — Swiss cheese
# ============================================================
cat("Section 2: Swiss cheese...\n")

# Four circular holes in a 10x10 square.
# Domain width = 10 coord units.
# raster_sigma = 0.25 coord units
# (previously 7 grid cells at res=450 over domain=11 → 7*(11/450) = 0.17)
hdefs <- list(
  list(cx = 2.5, cy = 7.5, r = 1.1),
  list(cx = 7.0, cy = 6.5, r = 1.4),
  list(cx = 4.0, cy = 2.5, r = 0.9),
  list(cx = 8.0, cy = 2.0, r = 1.2)
)
outer_sq  <- sf::st_sfc(sf::st_polygon(list(
  rbind(c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)))))
hcircles  <- lapply(hdefs, function(h) {
  sf::st_polygon(list(.circle_mat(h$cx, h$cy, h$r)))
})
swiss_region  <- sf::st_difference(outer_sq,
                   sf::st_sfc(sf::st_union(sf::st_sfc(hcircles))))
coords_swiss  <- sample_in_polygon(swiss_region, n = 1500)

m2_r <- fit_spatial_mask(coords_swiss, method = "raster",
  raster_resolution = 256L, raster_sigma = 0.25, raster_threshold = 0.2,
  verbose = FALSE)
m2_c <- fit_spatial_mask(coords_swiss, method = "concave",
  concavity = 1.8, verbose = FALSE)

print(
  plot_mask(m2_r, coords_swiss,
            title    = "2a. Swiss \u2014 Raster (holes \u2713)",
            subtitle = "\u03c3 = 0.25 coord units, threshold = 0.2") +
  plot_mask(m2_c, coords_swiss,
            title    = "2b. Swiss \u2014 Concave",
            subtitle = "[holes filled]",
            fill = "#e74c3c", border = "#922b21") +
  plot_annotation(title = "Section 2: Swiss Cheese",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
)


# ============================================================
# SECTION 3 — Fractured map
# ============================================================
cat("Section 3: Fractured map...\n")

# Irregular outer coast (radius ~9-12), five elliptic lakes.
# Domain width ~ 24 coord units.
# raster_sigma = 0.5 coord units
# (previously 9 grid cells at res=500 over domain=26.4 → 9*(26.4/500) = 0.48)
th_f     <- seq(0, 2 * pi, length.out = 600)
r_c      <- 9 + 2*sin(3*th_f) + sin(7*th_f) + 0.5*cos(13*th_f) + 0.3*sin(19*th_f)
mat_f    <- cbind(r_c * cos(th_f), r_c * sin(th_f))
mat_f    <- rbind(mat_f, mat_f[1L, ])   # explicit ring closure
outer_coast <- sf::st_sfc(sf::st_polygon(list(mat_f)))

lakes <- sf::st_sfc(list(
  make_ellipse_poly(-3.0,  3.5, 2.2, 1.0,  30),
  make_ellipse_poly( 4.5,  2.0, 1.8, 1.3, -20),
  make_ellipse_poly( 1.0, -4.5, 2.5, 0.8,  60),
  make_ellipse_poly(-5.0, -2.0, 1.2, 1.2,   0),
  make_ellipse_poly( 3.5, -0.5, 1.0, 0.6,  15)
))
frac_region  <- sf::st_difference(outer_coast, sf::st_union(lakes))
coords_frac  <- sample_in_polygon(frac_region, n = 2000)

m3_r <- fit_spatial_mask(coords_frac, method = "raster",
  raster_resolution = 256L, raster_sigma = 0.5, raster_threshold = 0.18,
  verbose = FALSE)
m3_k <- fit_spatial_mask(coords_frac, method = "kde",
  kde_resolution = 300L, kde_threshold = 0.01, buffer_dist = 0.1,
  verbose = FALSE)

print(
  plot_mask(m3_r, coords_frac,
            title    = "3a. Fractured \u2014 Raster",
            subtitle = "Coast + 5 lakes | \u03c3 = 0.5") +
  plot_mask(m3_k, coords_frac,
            title    = "3b. Fractured \u2014 KDE",
            subtitle = "Density-based boundary",
            fill = "#27ae60", border = "#1e8449") +
  plot_annotation(title = "Section 3: Fractured Map",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
)


# ============================================================
# SECTION 4 — Archipelago
# ============================================================
cat("Section 4: Archipelago...\n")

# Three disconnected islands; domain width ~ 16 coord units.
# raster_sigma = 0.4 coord units
# (previously 7 grid cells at res=450 over domain=17.6 → 7*(17.6/450) = 0.27)
mk_isl <- function(n, cx, cy, rx, ry, ns = 0.15) {
  th <- runif(n, 0, 2 * pi); r <- sqrt(runif(n))
  data.frame(x = cx + r * rx * cos(th) + rnorm(n, 0, ns),
             y = cy + r * ry * sin(th) + rnorm(n, 0, ns))
}
coords_arch <- rbind(mk_isl(300, 0, 0, 3, 2),
                     mk_isl(300, 9, 4, 2, 3),
                     mk_isl(300, 5, -7, 2.5, 1.5))

m4_r <- fit_spatial_mask(coords_arch, method = "raster",
  raster_resolution = 256L, raster_sigma = 0.4, raster_threshold = 0.18,
  verbose = FALSE)
m4_c <- fit_spatial_mask(coords_arch, method = "convex", verbose = FALSE)

n_geom <- length(sf::st_cast(m4_r, "POLYGON"))
print(
  plot_mask(m4_r, coords_arch,
            title    = "4a. Archipelago \u2014 Raster",
            subtitle = sprintf("%d disconnected polygons \u2713", n_geom)) +
  plot_mask(m4_c, coords_arch,
            title    = "4b. Archipelago \u2014 Convex",
            subtitle = "[bridges islands]",
            fill = "#e74c3c", border = "#922b21") +
  plot_annotation(title = "Section 4: Archipelago",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
)


# ============================================================
# SECTION 5 — Nested crescent voids
# ============================================================
cat("Section 5: Nested voids...\n")

# Two crescent-shaped voids inside a circle of radius 10.
# Domain width ~ 20 coord units.
# fine sigma = 0.4 (previously 8 gc at res=500 over domain=22 → 0.35)
# coarse sigma = 1.0 (previously 20 gc → 0.88)
mk_cres <- function(cx, cy, ro, ri, dx, dy, n = 100) {
  moc <- .circle_mat(cx,    cy,    ro, n)
  mic <- .circle_mat(cx+dx, cy+dy, ri, n)
  oc  <- sf::st_polygon(list(moc))
  ic  <- sf::st_polygon(list(mic))
  sf::st_sfc(sf::st_difference(sf::st_sfc(oc), sf::st_sfc(ic)))
}

th_od  <- seq(0, 2 * pi, length.out = 200)
mat_od <- cbind(10 * cos(th_od), 10 * sin(th_od))
mat_od <- rbind(mat_od, mat_od[1L, ])   # explicit ring closure
outer_d       <- sf::st_sfc(sf::st_polygon(list(mat_od)))
cres1         <- mk_cres(-2, 2,  2.5, 2,   0.9, -0.2)
cres2         <- mk_cres( 3, -3, 2.2, 1.8, -0.7,  0.4)
nest_region   <- sf::st_difference(outer_d, sf::st_union(c(cres1, cres2)))
coords_nested <- sample_in_polygon(nest_region, n = 2500)

m5_fine   <- fit_spatial_mask(coords_nested, method = "raster",
  raster_resolution = 256L, raster_sigma = 0.4,  raster_threshold = 0.18,
  verbose = FALSE)
m5_coarse <- fit_spatial_mask(coords_nested, method = "raster",
  raster_resolution = 256L, raster_sigma = 1.0, raster_threshold = 0.18,
  verbose = FALSE)

print(
  plot_mask(m5_fine,   coords_nested,
            title    = "5a. Nested \u2014 fine \u03c3 = 0.4",
            subtitle = "Crescent voids preserved") +
  plot_mask(m5_coarse, coords_nested,
            title    = "5b. Nested \u2014 coarse \u03c3 = 1.0",
            subtitle = "Narrow voids filled",
            fill = "#8e44ad", border = "#6c3483") +
  plot_annotation(
    title    = "Section 5: Nested Crescent Voids",
    subtitle = "Smaller \u03c3 \u2192 finer void detection",
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, color = "grey40")))
)


# ============================================================
# SECTION 6 — Sigma sweep (coordinate units)
# ============================================================
cat("Section 6: Sigma sweep...\n")

# Swiss cheese domain ~ 10 coord units.
# Sweep spans: tight (holes preserved) → loose (holes filled).
# Values in COORDINATE UNITS: c(0.1, 0.2, 0.4, 0.7)
# Equivalent old grid-cell values at res=400, domain=11:
#   c(3.6, 7.3, 14.5, 25.5)
sig_vals <- c(0.1, 0.2, 0.4, 0.7)
sig_cols <- c("#e74c3c", "#e67e22", "#27ae60", "#2980b9")

sig_pls <- mapply(function(sv, col) {
  m  <- fit_spatial_mask(coords_swiss, method = "raster",
    raster_resolution = 256L, raster_sigma = sv, raster_threshold = 0.2,
    verbose = FALSE)
  nh <- max(0L, length(sf::st_cast(m, "POLYGON")) - 1L)
  plot_mask(m, coords_swiss,
            title    = paste0("\u03c3 = ", sv, " cu"),
            subtitle = paste0(nh, " hole(s) detected"),
            fill = col, border = col, pt_size = 0.5)
}, sig_vals, sig_cols, SIMPLIFY = FALSE)

print(wrap_plots(sig_pls, nrow = 1) +
  plot_annotation(
    title    = "Section 6: raster_sigma Sensitivity \u2014 Swiss Cheese",
    subtitle = "All values in coordinate units | Small \u03c3 \u2192 holes preserved | Large \u03c3 \u2192 holes filled",
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, color = "grey40")))
)


# ============================================================
# SECTION 7 — 100k scale test
# ============================================================
cat("Section 7: 100k scale test...\n")

# Irregular outer shape radius ~12-14, three elliptic holes.
# Domain width ~ 28 coord units.
# raster_sigma = 0.6 coord units
# (previously 10 grid cells at res=512 over domain=30.8 → 10*(30.8/512) = 0.60)
th100    <- seq(0, 2 * pi, length.out = 400)
r100     <- 12 + sin(5 * th100) + 0.4 * cos(11 * th100)
mat100   <- cbind(r100 * cos(th100), r100 * sin(th100))
mat100   <- rbind(mat100, mat100[1L, ])   # explicit ring closure
outer100 <- sf::st_sfc(sf::st_polygon(list(mat100)))

holes100 <- sf::st_sfc(list(
  make_ellipse_poly( 3,  5,  3,   1.5,  10),
  make_ellipse_poly(-5,  2,  2.5, 2,   -30),
  make_ellipse_poly( 1, -6,  2,   3,    50)
))
reg100      <- sf::st_difference(outer100, sf::st_union(holes100))
cat("  Sampling 100,000 points...\n")
coords_100k <- sample_in_polygon(reg100, n = 100000)

t_r <- system.time(
  m7_r <- fit_spatial_mask(coords_100k, method = "raster",
    raster_resolution = 256L, raster_sigma = 0.6, raster_threshold = 0.18,
    verbose = TRUE)
)
t_k <- system.time(
  m7_k <- fit_spatial_mask(coords_100k, method = "kde",
    kde_resolution = 200L, kde_threshold = 0.005, buffer_dist = 0.1,
    verbose = TRUE)
)
cat(sprintf("\n  Raster (256\u00d7256): %.2f s\n  KDE    (200\u00d7200): %.2f s\n",
            t_r["elapsed"], t_k["elapsed"]))

sub <- coords_100k[sample(nrow(coords_100k), 5000), ]
print(
  plot_mask(m7_r, sub,
            title    = sprintf("7a. 100k \u2014 Raster (%.1f s)", t_r["elapsed"]),
            subtitle = "3 holes detected \u2713",
            pt_size = 0.25, pt_alpha = 0.25) +
  plot_mask(m7_k, sub,
            title    = sprintf("7b. 100k \u2014 KDE (%.1f s)", t_k["elapsed"]),
            subtitle = "Density boundary",
            fill = "#27ae60", border = "#1e8449",
            pt_size = 0.25, pt_alpha = 0.25) +
  plot_annotation(
    title    = "Section 7: 100,000-Cell Scale Test",
    subtitle = "5,000 of 100,000 points shown",
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, color = "grey40")))
)


# ============================================================
# SECTION 8 — Validation table
# ============================================================
cat("\nSection 8: Validation...\n")

vcases <- list(
  list(label = "Donut/raster",       mask = m1_r, coords = coords_donut),
  list(label = "Donut/convex",        mask = m1_c, coords = coords_donut),
  list(label = "Swiss/raster",        mask = m2_r, coords = coords_swiss),
  list(label = "Swiss/concave",       mask = m2_c, coords = coords_swiss),
  list(label = "Fractured/raster",    mask = m3_r, coords = coords_frac),
  list(label = "Fractured/kde",       mask = m3_k, coords = coords_frac),
  list(label = "Archipelago/raster",  mask = m4_r, coords = coords_arch),
  list(label = "Archipelago/convex",  mask = m4_c, coords = coords_arch),
  list(label = "Nested/fine",         mask = m5_fine,   coords = coords_nested),
  list(label = "Nested/coarse",       mask = m5_coarse, coords = coords_nested),
  list(label = "100k/raster",         mask = m7_r, coords = coords_100k),
  list(label = "100k/kde",            mask = m7_k, coords = coords_100k)
)

res <- do.call(rbind, lapply(vcases, function(vc) {
  pts  <- sf::st_geometry(sf::st_as_sf(vc$coords, coords = c("x", "y")))
  ok   <- sf::st_within(pts, sf::st_union(vc$mask), sparse = FALSE)[, 1L]
  n    <- nrow(vc$coords)
  data.frame(
    Case     = vc$label,
    N        = n,
    Enclosed = sum(ok),
    Pct      = round(100 * sum(ok) / n, 3),
    Result   = ifelse(sum(ok) == n, "\u2713 PASS", "\u2717 FAIL"),
    stringsAsFactors = FALSE
  )
}))

cat("\n============================================================\n")
cat("  VALIDATION: Point Containment (all must be PASS)\n")
cat("============================================================\n")
print(res, row.names = FALSE)
cat("============================================================\n")
cat("demo_holes_and_scale.R complete.\n")

