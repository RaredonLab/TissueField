# Internal helpers for TissueField
# These functions are not exported and do not generate .Rd pages.

#' @noRd
.gauss_kernel_1d <- function(sigma, n_max = NULL) {
  r <- max(1L, ceiling(3 * sigma))
  if (!is.null(n_max)) r <- min(r, floor((n_max - 1L) / 2L))
  k <- exp(-(seq(-r, r))^2 / (2 * sigma^2))
  k / sum(k)
}

#' @noRd
#' @importFrom stats filter
.sep_gauss_smooth_2d <- function(mat, sigma_x, sigma_y) {
  nr <- nrow(mat); nc <- ncol(mat)
  kx <- .gauss_kernel_1d(sigma_x, n_max = nr)
  ky <- .gauss_kernel_1d(sigma_y, n_max = nc)
  out <- apply(mat, 2L, function(col) {
    s <- stats::filter(col, kx, method = "convolution", sides = 2L)
    s[is.na(s)] <- 0
    as.numeric(s)
  })
  t(apply(out, 1L, function(row) {
    s <- stats::filter(row, ky, method = "convolution", sides = 2L)
    s[is.na(s)] <- 0
    as.numeric(s)
  }))
}

#' @noRd
#' @importFrom sf st_within
.rasterize_mask <- function(gc_sf, mask_union, N, n_cores) {
  n_cells <- N * N
  if (n_cores > 1L && .Platform$OS.type == "unix") {
    chunk_size <- ceiling(n_cells / n_cores)
    chunks     <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / chunk_size))
    res_list   <- parallel::mclapply(
      chunks,
      function(idx)
        as.logical(sf::st_within(gc_sf[idx, ], mask_union, sparse = FALSE)[, 1L]),
      mc.cores = n_cores
    )
    in_msk <- unlist(res_list)
  } else {
    in_msk <- as.logical(sf::st_within(gc_sf, mask_union, sparse = FALSE)[, 1L])
  }
  matrix(in_msk, nrow = N, ncol = N)
}
