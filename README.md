# TissueField

**TissueField** computes continuous steady-state molecular concentration fields from discrete mRNA
transcript or protein detection coordinates in spatial transcriptomics and spatial proteomics
experiments.

## Physical model

The package solves the **screened Poisson (diffusion-clearance) PDE**:

```
D * nabla^2 C(x) - lambda * C(x) + s(x) = 0
```

The key parameter is the **diffusion length** `L = sqrt(D / lambda)`, which sets the characteristic
decay distance from each source. Small `L` gives tight, cluster-specific fields; large `L` fills
the entire tissue.

## Solvers

| Method | Algorithm | Boundary conditions | Best for |
|--------|-----------|---------------------|----------|
| `"fd"` | Finite-difference sparse linear system | Exact Dirichlet or Neumann | Default; recommended |
| `"green"` | FFT convolution with analytic Green's function | None (infinite domain) | Parameter sweeps, large grids |
| `"kde"` | Gaussian kernel density | None | Quick visualisation baseline |

## Installation

```r
# Install from GitHub
remotes::install_github("RaredonLab/TissueField")

# Required: sf and Matrix
install.packages(c("sf", "Matrix"))

# Recommended companion package for mask fitting
remotes::install_github("RaredonLab/TissueMask")
```

## Quick start

```r
library(TissueField)
library(sf)

# 1. Create (or load) an sfc polygon mask
#    e.g. from TissueMask::fit_spatial_mask()
sq   <- st_polygon(list(cbind(c(0,10,10,0,0), c(0,0,10,10,0))))
mask <- st_sfc(sq)

# 2. Provide transcript coordinates
set.seed(1)
tc <- data.frame(x = runif(200, 1, 9), y = runif(200, 1, 9))

# 3. Estimate the concentration field
res <- estimate_concentration_field(
  mask              = mask,
  transcript_coords = tc,
  diffusion_length  = 2,     # L = 2 spatial units
  method            = "fd",  # finite-difference solver
  grid_resolution   = 128L,
  verbose           = TRUE
)

# 4. Visualise
plot_concentration_field(res, tc, show_contours = TRUE)

# 5. Get a tidy data frame for custom plotting
df <- field_to_df(res)
head(df)

# 6. Sweep over multiple diffusion lengths
sw <- sweep_diffusion_length(c(0.5, 1, 3, 8), mask, tc,
                              grid_resolution = 64L)
```

## Key parameters

| Parameter | Controls | Typical values |
|-----------|----------|----------------|
| `diffusion_length` | Spatial scale of field | 0.5 -- 20 (in coordinate units) |
| `D` | Amplitude (shape fixed by L) | 1.0 (leave at 1 unless you have physical units) |
| `lambda` | Clearance rate | Derived from L and D |
| `weight_col` | Per-transcript UMI weighting | e.g. `"umi"` |
| `boundary_condition` | Absorbing vs reflecting edge | `"dirichlet"` (default) |
| `grid_resolution` | Grid cells per axis | 128--512 |

## Vignettes

- **Getting started**: basic usage, `field_to_df()`, `plot_concentration_field()`, diffusion length intuition
- **Solver guide**: fd vs green vs kde, Dirichlet vs Neumann BCs, direct vs iterative
- **Parameter tuning**: L sweep, D vs lambda, UMI weighting, grid resolution, log-transform
- **Complete reference**: 12-section deep dive covering all features

## Citation

If you use TissueField in your research, please cite:

> Raredon Laboratory (2026). *TissueField: Steady-State Molecular Concentration Fields for
> Spatial Transcriptomics.* R package version 0.1.0.
> https://github.com/RaredonLab/TissueField

## License

MIT (c) 2026 Raredon Laboratory
