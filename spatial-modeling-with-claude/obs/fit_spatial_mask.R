# ============================================================
# fit_spatial_mask.R  (v3)
#
# Fits a spatial boundary mask to XY point coordinates.
# Methods: "convex" | "concave" | "kde" | "raster"
#
# RASTER METHOD (v3 rewrite)
# ---------------------------
# Previously: Gaussian smooth → isoband contour → manual hole assembly.
# Problem:    isoband winding is inconsistent; manual hole assembly
#             produced bowtie artefacts and diagonal fill lines.
#
# Now:        Gaussian smooth → binary threshold → union of "on" grid cells.
# Key insight: unioning axis-aligned rectangles via GEOS dissolve always
#             produces valid topology. Holes and islands emerge naturally —
#             no ring-winding logic needed at all.
#
# KEY PARAMETER CHANGE (v3):
#   raster_sigma is now in COORDINATE UNITS (same units as x/y),
#   not grid-cell units. This makes it directly interpretable:
#   "how wide is the Gaussian around each point?"
#   NULL → auto = 3 % of domain width.
#
# install.packages(c("sf", "concaveman", "MASS", "isoband", "ggplot2"))
# ============================================================


# ---- Internal helpers --------------------------------------------------------

.signed_area <- function(ring_mat) {
  n <- nrow(ring_mat); if (n < 3L) return(0)
  x <- ring_mat[, 1L]; y <- ring_mat[, 2L]
  0.5 * sum(x[seq_len(n - 1L)] * y[seq_len(n - 1L) + 1L] -
              x[seq_len(n - 1L) + 1L] * y[seq_len(n - 1L)])
}

.gaussian_kernel_1d <- function(sigma) {
  r <- ceiling(3 * sigma)
  k <- exp(-((-r):r)^2 / (2 * sigma^2))
  k / sum(k)
}

.smooth_rows <- function(mat, k) {
  t(apply(mat, 1L, function(row) {
    s <- stats::filter(row, k, method = "convolution", sides = 2L)
    s[is.na(s)] <- 0; as.numeric(s)
  }))
}

.smooth_cols <- function(mat, k) {
  apply(mat, 2L, function(col) {
    s <- stats::filter(col, k, method = "convolution", sides = 2L)
    s[is.na(s)] <- 0; as.numeric(s)
  })
}

# Separable 2-D Gaussian convolution.  sigma is in GRID-CELL units.
.gaussian_smooth_2d <- function(mat, sigma, n_cores = 1L) {
  k <- .gaussian_kernel_1d(sigma)
  if (n_cores > 1L && .Platform$OS.type == "unix") {
    row_list <- parallel::mclapply(seq_len(nrow(mat)), function(i) {
      s <- stats::filter(mat[i, ], k, method = "convolution", sides = 2L)
      s[is.na(s)] <- 0; as.numeric(s)
    }, mc.cores = n_cores)
    mat <- do.call(rbind, row_list)
  } else {
    mat <- .smooth_rows(mat, k)
  }
  .smooth_cols(mat, k)
}

# Used only by the "kde" method (isoband-based).
.build_sf_with_holes <- function(bands, crs) {
  rings <- list()
  for (band in bands) {
    if (length(band$x) == 0L) next
    ids <- if (!is.null(band$id)) band$id else rep(1L, length(band$x))
    for (uid in unique(ids[!is.na(ids)])) {
      sel <- !is.na(ids) & ids == uid
      rx  <- band$x[sel]; ry <- band$y[sel]
      if (length(rx) < 3L) next
      if (rx[1L] != rx[length(rx)] || ry[1L] != ry[length(ry)]) {
        rx <- c(rx, rx[1L]); ry <- c(ry, ry[1L])
      }
      mat  <- cbind(rx, ry)
      area <- abs(.signed_area(mat))
      rings <- c(rings, list(list(mat = mat, area = area)))
    }
  }
  if (length(rings) == 0L)
    stop("No rings extracted. Adjust kde_threshold / kde_bandwidth parameters.")
  ord   <- order(sapply(rings, `[[`, "area"), decreasing = TRUE)
  rings <- rings[ord]
  polys <- lapply(rings, function(r) {
    mat  <- rbind(r$mat, r$mat[1L, , drop = FALSE])
    dupl <- c(FALSE, rowSums(abs(diff(mat))) == 0)
    sf::st_sfc(sf::st_polygon(list(mat[!dupl, , drop = FALSE])), crs = crs)
  })
  result <- polys[[1L]]
  if (length(polys) > 1L) {
    for (i in seq(2L, length(polys))) {
      inside <- tryCatch(
        isTRUE(sf::st_within(polys[[i]], result, sparse = FALSE)[1L, 1L]),
        error = function(e) FALSE
      )
      result <- if (inside) sf::st_difference(result, polys[[i]])
                else         sf::st_union(result, polys[[i]])
    }
  }
  sf::st_make_valid(result)
}


# ---- Main function -----------------------------------------------------------

fit_spatial_mask <- function(
  coords,

  # Method
  method            = "raster",   # "convex" | "concave" | "kde" | "raster"

  # Concave-hull parameters
  concavity         = 2,
  length_threshold  = 0,

  # KDE parameters
  kde_bandwidth     = NULL,
  kde_threshold     = 0.05,
  kde_resolution    = 256L,

  # Raster parameters
  # raster_sigma: Gaussian spread in COORDINATE UNITS (same as x/y).
  #   NULL  → auto = 3 % of domain width.
  #   Larger sigma → holes fill in, islands merge, boundary smooths.
  #   Smaller sigma → holes and fine gaps are preserved.
  raster_resolution = 256L,
  raster_sigma      = NULL,
  raster_threshold  = 0.15,
  raster_min_pts    = 1L,

  # Post-processing
  buffer_dist       = 0,
  smooth_mask       = FALSE,
  smooth_tolerance  = NULL,

  # Misc
  n_cores           = 1L,
  crs               = NA,
  plot              = FALSE,
  verbose           = TRUE
) {

  # ---- 0. Package checks -----------------------------------------------------
  req  <- c("sf", "concaveman", "MASS", "isoband")
  miss <- req[!sapply(req, requireNamespace, quietly = TRUE)]
  if (length(miss) > 0L)
    stop("Install missing packages: install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))")

  # ---- 0b. Normalise CRS (sf::st_sfc rejects bare NA in sf >= 1.0) ----------
  crs <- if (length(crs) == 1L && is.na(crs)) sf::NA_crs_ else crs

  # ---- 1. Parse and validate coordinates ------------------------------------
  if (is.matrix(coords)) coords <- as.data.frame(coords)
  if (!is.data.frame(coords)) stop("`coords` must be a data.frame or matrix.")
  if (!all(c("x", "y") %in% names(coords))) {
    if (ncol(coords) >= 2L) {
      coords <- coords[, 1:2]; names(coords) <- c("x", "y")
    } else stop("Need columns 'x' and 'y'.")
  }
  coords <- coords[complete.cases(coords[, c("x", "y")]), c("x", "y")]
  n_pts  <- nrow(coords)
  if (n_pts < 3L) stop("Need at least 3 coordinate pairs.")
  if (n_pts > 50000L && method == "kde" && verbose)
    message("Large n = ", n_pts, " — consider method = 'raster'.")

  # ---- 2. Fit mask -----------------------------------------------------------
  mask_poly <- switch(method,

    # -- Convex hull ------------------------------------------------------------
    "convex" = {
      pts <- sf::st_as_sf(coords, coords = c("x", "y"), crs = crs)
      sf::st_convex_hull(sf::st_union(pts))
    },

    # -- Concave hull -----------------------------------------------------------
    "concave" = {
      hp <- concaveman::concaveman(as.matrix(coords),
                                   concavity = concavity,
                                   length_threshold = length_threshold)
      sf::st_sfc(sf::st_polygon(list(hp)), crs = crs)
    },

    # -- KDE contour (supports holes, slower for large n) ----------------------
    "kde" = {
      bw_x <- if (!is.null(kde_bandwidth)) kde_bandwidth[1L] else
               MASS::bandwidth.nrd(coords$x)
      bw_y <- if (!is.null(kde_bandwidth)) {
        if (length(kde_bandwidth) == 2L) kde_bandwidth[2L] else kde_bandwidth[1L]
      } else MASS::bandwidth.nrd(coords$y)
      xp  <- diff(range(coords$x)) * 0.1
      yp  <- diff(range(coords$y)) * 0.1
      kf  <- MASS::kde2d(coords$x, coords$y, h = c(bw_x, bw_y),
                         n = kde_resolution,
                         lims = c(range(coords$x) + c(-xp, xp),
                                  range(coords$y) + c(-yp, yp)))
      dv  <- as.vector(kf$z)
      bands <- isoband::isobands(kf$x, kf$y, kf$z,
                                  levels_low  = quantile(dv, kde_threshold),
                                  levels_high = max(dv) + 1)
      .build_sf_with_holes(bands, crs)
    },

    # -- Raster: Gaussian sum → threshold → cell union -------------------------
    #
    # Algorithm:
    #   1. Bin points onto a regular grid (raster_resolution x raster_resolution).
    #   2. Convolve the binary occupancy grid with a 2-D Gaussian of width
    #      raster_sigma (coordinate units).  Equivalent to placing an isotropic
    #      Gaussian at every point and summing them all.
    #   3. Threshold at raster_threshold x max(field).  Cells above threshold
    #      are "inside"; cells below are "outside".
    #   4. Convert each "inside" cell to an axis-aligned rectangle and dissolve
    #      via GEOS union.  Holes and islands emerge from the geometry naturally
    #      — no ring-winding logic required.
    #   5. Smooth the staircase boundary with a morphological close
    #      (buffer out then in by ~half the cell diagonal).
    #
    # Tuning guide:
    #   raster_sigma     up  → holes fill in, islands merge, edges smoother
    #   raster_sigma     down → holes preserved, fine structure retained
    #   raster_threshold up  → mask shrinks (requires denser coverage)
    #   raster_threshold down → mask grows (accepts sparse regions)
    "raster" = {

      # Grid cell edges
      xpad <- diff(range(coords$x)) * 0.05
      ypad <- diff(range(coords$y)) * 0.05
      xb   <- seq(min(coords$x) - xpad, max(coords$x) + xpad,
                  length.out = raster_resolution + 1L)
      yb   <- seq(min(coords$y) - ypad, max(coords$y) + ypad,
                  length.out = raster_resolution + 1L)
      hx   <- xb[2L] - xb[1L]   # cell width  (coordinate units)
      hy   <- yb[2L] - yb[1L]   # cell height (coordinate units)

      # Bin points: cm[yi, xi] = point count at cell (xi, yi)
      xi <- pmax(1L, pmin(raster_resolution,
                          findInterval(coords$x, xb, rightmost.closed = TRUE)))
      yi <- pmax(1L, pmin(raster_resolution,
                          findInterval(coords$y, yb, rightmost.closed = TRUE)))
      cv <- tabulate(yi + (xi - 1L) * raster_resolution,
                     nbins = raster_resolution^2L)
      cm <- matrix(cv, nrow = raster_resolution, ncol = raster_resolution)

      # Convert raster_sigma from coordinate units to grid-cell units
      domain_w <- diff(range(coords$x)) + 2 * xpad
      sigma_cu <- if (is.null(raster_sigma)) 0.03 * domain_w else raster_sigma
      sigma_gc <- sigma_cu / hx   # grid-cell units for the convolution kernel

      if (verbose)
        message("  raster_sigma = ", round(sigma_cu, 4),
                " coord units  (", round(sigma_gc, 1), " grid cells)")

      # Convolve binary occupancy with Gaussian — sum of point Gaussians
      indicator <- (cm >= raster_min_pts) * 1.0
      sm        <- .gaussian_smooth_2d(indicator, sigma = sigma_gc,
                                        n_cores = n_cores)

      # Threshold → binary
      thresh   <- raster_threshold * max(sm, na.rm = TRUE)
      binary   <- sm >= thresh
      on_cells <- which(binary, arr.ind = TRUE)   # [row = yi, col = xi]

      if (nrow(on_cells) == 0L)
        stop("No grid cells above threshold. Try: lower raster_threshold, ",
             "increase raster_sigma, or increase raster_resolution.")

      # Build one closed rectangle per on-cell
      ri <- on_cells[, 1L]   # y indices
      ci <- on_cells[, 2L]   # x indices
      rects <- lapply(seq_len(nrow(on_cells)), function(k) {
        sf::st_polygon(list(matrix(
          c(xb[ci[k]],      yb[ri[k]],
            xb[ci[k] + 1L], yb[ri[k]],
            xb[ci[k] + 1L], yb[ri[k] + 1L],
            xb[ci[k]],      yb[ri[k] + 1L],
            xb[ci[k]],      yb[ri[k]]),       # close ring
          ncol = 2L, byrow = TRUE
        )))
      })

      # GEOS dissolve — topology (holes, islands) emerges automatically
      raw <- sf::st_union(sf::st_sfc(rects, crs = crs))

      # Morphological close: smooth staircase boundary.
      # Buffer out then in by ~half the cell diagonal.
      cell_diag <- sqrt(hx^2 + hy^2)
      raw <- sf::st_buffer(raw, dist =  cell_diag * 0.6)
      raw <- sf::st_buffer(raw, dist = -cell_diag * 0.6)

      sf::st_make_valid(raw)
    },

    stop("Unknown method '", method, "'. Choose: convex, concave, kde, raster.")
  )

  if (!inherits(mask_poly, "sfc")) mask_poly <- sf::st_geometry(mask_poly)

  # ---- 3. Optional buffer ----------------------------------------------------
  if (buffer_dist > 0)
    mask_poly <- sf::st_buffer(mask_poly, dist = buffer_dist)

  # ---- 4. Containment guarantee ----------------------------------------------
  pts_sfc   <- sf::st_geometry(
    sf::st_as_sf(coords, coords = c("x", "y"), crs = crs))
  contained <- sf::st_within(pts_sfc, sf::st_union(mask_poly), sparse = FALSE)
  n_out     <- sum(!contained[, 1L])
  if (n_out > 0L) {
    warning(n_out, " point(s) outside mask — applying corrective buffer.")
    dists     <- sf::st_distance(pts_sfc[!contained[, 1L]],
                                  sf::st_union(mask_poly))
    mask_poly <- sf::st_buffer(mask_poly,
                                dist = max(as.numeric(dists)) * 1.05)
  }

  # ---- 5. Optional boundary smoothing ----------------------------------------
  if (smooth_mask) {
    if (is.null(smooth_tolerance)) {
      bb <- sf::st_bbox(mask_poly)
      smooth_tolerance <- sqrt((bb["xmax"] - bb["xmin"])^2 +
                               (bb["ymax"] - bb["ymin"])^2) * 0.01
    }
    mask_poly <- sf::st_simplify(mask_poly, dTolerance = smooth_tolerance)
  }
  sf::st_crs(mask_poly) <- crs

  # ---- 6. Verbose summary ----------------------------------------------------
  if (verbose) {
    bb   <- sf::st_bbox(mask_poly)
    area <- suppressWarnings(as.numeric(sf::st_area(mask_poly)))
    cat("--------------------------------------------------\n")
    cat("  Spatial mask fitted (v3)\n")
    cat("  Method         :", method, "\n")
    cat("  Points         :", n_pts, "\n")
    cat("  Sub-geometries :", length(mask_poly), "\n")
    cat("  Bounding box   : x [", round(bb["xmin"], 3), ",",
        round(bb["xmax"], 3), "]  y [",
        round(bb["ymin"], 3), ",", round(bb["ymax"], 3), "]\n")
    cat("  Area           :", round(area, 4), "\n")
    cat("  Buffer applied :", buffer_dist, "\n")
    cat("--------------------------------------------------\n")
  }

  # ---- 7. Optional quick plot ------------------------------------------------
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    md  <- as.data.frame(sf::st_coordinates(mask_poly))
    grp <- if ("L3" %in% names(md)) "L3" else
           if ("L2" %in% names(md)) "L2" else "L1"
    print(ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = md, ggplot2::aes(x = X, y = Y, group = .data[[grp]]),
        fill = "#4c9be8", alpha = 0.2, color = "#1a5fa8", linewidth = 0.55) +
      ggplot2::geom_point(
        data = coords, ggplot2::aes(x = x, y = y),
        color = "#c0392b", size = 0.7, alpha = 0.5) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = paste0("Mask | method='", method, "'"),
                    x = "X", y = "Y") +
      ggplot2::theme_minimal(base_size = 11))
  }

  return(mask_poly)
}
