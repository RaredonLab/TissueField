#' TissueField: Steady-State Molecular Concentration Fields for Spatial Transcriptomics
#'
#' Computes continuous steady-state molecular concentration fields from discrete
#' mRNA transcript or protein detection coordinates using a physically motivated
#' diffusion-clearance model. The screened Poisson PDE is solved numerically
#' using one of three methods: finite-difference sparse linear system (`"fd"`),
#' Green's function FFT convolution (`"green"`), or Gaussian kernel density
#' (`"kde"`).
#'
#' @keywords internal
"_PACKAGE"

## Suppress R CMD check NOTE for bare variable names used in ggplot2 aes()
utils::globalVariables(c("x", "y", "field", ".data"))
