## Tests for estimate_concentration_field()
## Uses a simple square polygon so TissueMask is not required.

library(sf)

# ---- shared fixtures -------------------------------------------------------

sq   <- st_polygon(list(cbind(c(0, 10, 10, 0, 0), c(0, 0, 10, 10, 0))))
mask <- st_sfc(sq)

set.seed(42)
tc <- data.frame(x = runif(80, 1, 9), y = runif(80, 1, 9))

# run once for re-use
res_default <- estimate_concentration_field(
  mask, tc, D = 1, lambda = 0.2,
  grid_resolution = 64L, verbose = FALSE
)


# ---- return structure ------------------------------------------------------

test_that("returns a list with required elements", {
  expect_type(res_default, "list")
  expect_named(res_default, c("field", "x", "y", "hx", "hy",
                               "mask", "sources", "params", "diagnostics"))
})

test_that("field matrix has correct dimensions", {
  expect_equal(dim(res_default$field), c(64L, 64L))
})

test_that("x and y vectors have length N", {
  expect_length(res_default$x, 64L)
  expect_length(res_default$y, 64L)
})

test_that("mask matrix is logical with correct dimensions", {
  expect_type(res_default$mask, "logical")
  expect_equal(dim(res_default$mask), c(64L, 64L))
})

test_that("sources matrix is returned by default", {
  expect_equal(dim(res_default$sources), c(64L, 64L))
})

test_that("sources is NULL when return_sources = FALSE", {
  res <- estimate_concentration_field(
    mask, tc, grid_resolution = 32L, verbose = FALSE,
    return_sources = FALSE
  )
  expect_null(res$sources)
})


# ---- field values ----------------------------------------------------------

test_that("field is NA outside mask", {
  outside <- !res_default$mask
  expect_true(all(is.na(res_default$field[outside])))
})

test_that("field is non-negative inside mask (clip_negative = TRUE)", {
  inside_vals <- res_default$field[res_default$mask]
  expect_true(all(inside_vals >= 0, na.rm = TRUE))
})

test_that("field has positive peak concentration", {
  expect_gt(max(res_default$field, na.rm = TRUE), 0)
})


# ---- solvers ---------------------------------------------------------------

test_that("method = 'fd' direct solver works", {
  res <- estimate_concentration_field(
    mask, tc, method = "fd", fd_solver = "direct",
    grid_resolution = 32L, verbose = FALSE
  )
  expect_false(all(is.na(res$field)))
})

test_that("method = 'fd' iterative (SOR) solver works", {
  res <- estimate_concentration_field(
    mask, tc, method = "fd", fd_solver = "iterative",
    D = 1, lambda = 0.2, grid_resolution = 32L, verbose = FALSE
  )
  expect_false(all(is.na(res$field)))
})

test_that("method = 'green' works", {
  res <- estimate_concentration_field(
    mask, tc, method = "green", D = 1, lambda = 0.2,
    grid_resolution = 32L, verbose = FALSE
  )
  expect_false(all(is.na(res$field)))
})

test_that("method = 'kde' works", {
  res <- estimate_concentration_field(
    mask, tc, method = "kde", kde_bandwidth = 2,
    grid_resolution = 32L, verbose = FALSE
  )
  expect_false(all(is.na(res$field)))
})


# ---- diffusion_length override --------------------------------------------

test_that("diffusion_length override changes peak concentration", {
  res_sm <- estimate_concentration_field(
    mask, tc, diffusion_length = 0.5,
    grid_resolution = 32L, verbose = FALSE
  )
  res_lg <- estimate_concentration_field(
    mask, tc, diffusion_length = 8.0,
    grid_resolution = 32L, verbose = FALSE
  )
  pk_sm <- max(res_sm$field, na.rm = TRUE)
  pk_lg <- max(res_lg$field, na.rm = TRUE)
  # Smaller L -> tighter, higher peak
  expect_gt(pk_sm, pk_lg)
})

test_that("params$diffusion_length matches supplied value", {
  res <- estimate_concentration_field(
    mask, tc, diffusion_length = 3.5,
    grid_resolution = 32L, verbose = FALSE
  )
  expect_equal(res$params$diffusion_length, 3.5)
})


# ---- boundary conditions ---------------------------------------------------

test_that("Neumann BC produces higher peak than Dirichlet at same params", {
  rd <- estimate_concentration_field(
    mask, tc, method = "fd", fd_solver = "direct",
    boundary_condition = "dirichlet", D = 1, lambda = 0.2,
    grid_resolution = 32L, verbose = FALSE
  )
  rn <- estimate_concentration_field(
    mask, tc, method = "fd", fd_solver = "direct",
    boundary_condition = "neumann", D = 1, lambda = 0.2,
    grid_resolution = 32L, verbose = FALSE
  )
  expect_gt(max(rn$field, na.rm = TRUE), max(rd$field, na.rm = TRUE))
})


# ---- post-processing -------------------------------------------------------

test_that("normalize = TRUE gives field in [0, 1]", {
  res <- estimate_concentration_field(
    mask, tc, normalize = TRUE, grid_resolution = 32L, verbose = FALSE
  )
  inside <- res$field[res$mask]
  expect_lte(max(inside), 1 + 1e-10)
  expect_gte(min(inside), 0 - 1e-10)
})

test_that("log_transform = TRUE gives non-negative transformed field", {
  res <- estimate_concentration_field(
    mask, tc, log_transform = TRUE, grid_resolution = 32L, verbose = FALSE
  )
  inside <- res$field[res$mask]
  expect_true(all(inside >= 0))
})


# ---- weighting and external transcripts -----------------------------------

test_that("weight_col scales field proportionally", {
  tc2 <- tc
  tc2$umi <- 1
  tc2$umi[1:20] <- 10
  res_w <- estimate_concentration_field(
    mask, tc2, weight_col = "umi", production_rate = 1,
    grid_resolution = 32L, verbose = FALSE
  )
  # Simple unweighted version
  res_u <- estimate_concentration_field(
    mask, tc2, grid_resolution = 32L, verbose = FALSE
  )
  # Weighted peak should be higher (heavy transcripts have more influence)
  expect_gt(max(res_w$field, na.rm = TRUE), max(res_u$field, na.rm = TRUE))
})

test_that("include_external = TRUE does not error", {
  tc_ext <- data.frame(x = c(tc$x, 15, 16), y = c(tc$y, 15, 16))
  expect_no_error(
    estimate_concentration_field(
      mask, tc_ext, include_external = TRUE,
      grid_resolution = 32L, verbose = FALSE
    )
  )
})


# ---- matrix input ----------------------------------------------------------

test_that("matrix input for transcript_coords is accepted", {
  tc_mat <- as.matrix(tc[, c("x", "y")])
  expect_no_error(
    estimate_concentration_field(
      mask, tc_mat, grid_resolution = 32L, verbose = FALSE
    )
  )
})


# ---- NA rows in input ------------------------------------------------------

test_that("NA rows in transcript_coords are silently dropped", {
  tc_na <- tc
  tc_na[5, "x"] <- NA
  expect_no_error(
    estimate_concentration_field(
      mask, tc_na, grid_resolution = 32L, verbose = FALSE
    )
  )
})


# ---- verbose = FALSE -------------------------------------------------------

test_that("verbose = FALSE produces no messages", {
  expect_no_message(
    suppressWarnings(
      estimate_concentration_field(
        mask, tc, grid_resolution = 32L, verbose = FALSE
      )
    )
  )
})


# ---- error cases -----------------------------------------------------------

test_that("no transcripts in mask raises error", {
  tc_out <- data.frame(x = c(20, 21), y = c(20, 21))
  expect_error(
    estimate_concentration_field(mask, tc_out, verbose = FALSE),
    "No transcripts available"
  )
})

test_that("method = 'green' with lambda = 0 raises error", {
  expect_error(
    estimate_concentration_field(
      mask, tc, method = "green", lambda = 0,
      grid_resolution = 32L, verbose = FALSE
    ),
    "lambda"
  )
})

test_that("unknown method raises informative error", {
  expect_error(
    estimate_concentration_field(
      mask, tc, method = "banana",
      grid_resolution = 32L, verbose = FALSE
    ),
    "Unknown method"
  )
})


# ---- utility functions -----------------------------------------------------

test_that("field_to_df returns data.frame with required columns", {
  df <- field_to_df(res_default)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("x", "y", "field", "inside") %in% names(df)))
})

test_that("field_to_df inside_only = FALSE includes NA-field rows", {
  df_all  <- field_to_df(res_default, inside_only = FALSE)
  df_in   <- field_to_df(res_default, inside_only = TRUE)
  expect_gt(nrow(df_all), nrow(df_in))
})

test_that("sweep_diffusion_length returns data.frame with L column", {
  sw <- sweep_diffusion_length(
    c(1, 3), mask, tc, grid_resolution = 32L
  )
  expect_s3_class(sw, "data.frame")
  expect_true("L" %in% names(sw))
  expect_equal(sort(unique(sw$L)), c(1, 3))
})

test_that("plot_concentration_field returns ggplot when ggplot2 available", {
  skip_if_not_installed("ggplot2")
  p <- plot_concentration_field(res_default)
  expect_s3_class(p, "ggplot")
})
