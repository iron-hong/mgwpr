# MGWPR

**MGWPR** is an R implementation of Multiscale Geographically Weighted Panel Regression. It provides pooled MGWR and multiscale geographically weighted fixed-effects regression (MGWFER), with location-specific coefficients and predictor-specific bandwidths.

The panel estimators support continuous outcomes with `family = "gaussian"` and binary outcomes with `family = "binomial"`. Binary models currently require `method = "pooled"`. Supported kernels are `bisquare`, `gaussian`, `exponential`, `tricube`, and `boxcar`.

Spatial data can be supplied as an `sf` object using `spatial`. Geometries are matched to panel units by ID. Distances between polygon locations are calculated from representative points, while the original geometries are retained for mapping. `mgwpr_sf()` joins location-level coefficients and available context effects to the `sf` object. Coefficient and significance maps can display polygon surfaces.

The package includes Moran's I, local F, local Hausman, and LM diagnostics through `diagnose_moran()`, `diagnose_local_f()`, `diagnose_hausman()`, `diagnose_lm()`, and the combined `diagnose_mgwpr()` interface. It also provides a cross-sectional MGWR interface through `GWmodel::gwr.multiscale()`.

## Installation

Install the released package from GitHub:

```r
install.packages("remotes")
remotes::install_github("iron-hong/mgwpr")
```

Or install a local source package:

```r
remotes::install_local("MGWPR_0.0.1.tar.gz")
```

## Simulate panel data and fit models

The example below generates a 15 by 15 spatial grid observed over three periods, then fits pooled MGWR and individual fixed-effects MGWFER models.

```r
library(MGWPR)

sim <- simulate_mgwpr_panel(
  n_grid = 15,
  n_time = 3,
  sigma_x = 0.5,
  spatial_coupling = 0.05,
  sigma_epsilon = 0.5,
  context_scale = 30,
  time_values = 2019:2021,
  seed = 42
)

panel <- sim$panel
units <- sim$units
model_formula <- y ~ x1 + x2 + x3 + x4

# pooled MGWR
fit_pmgwr <- mgwpr(
  formula = model_formula,
  data = panel,
  coords = c("grid_i", "grid_j"),
  id = "unit_id",
  time = "time_id",
  method = "pooled",
  family = "gaussian",
  adaptive = TRUE,
  kernel = "bisquare",
  scale = TRUE,
  max_iter = 30,
  tol = 1e-5
)

# multiscale geographically weighted fixed-effects regression 
fit_mgwfer <- mgwpr(
  formula = model_formula,
  data = panel,
  coords = c("grid_i", "grid_j"),
  id = "unit_id",
  time = "time_id",
  method = "fe",
  family = "gaussian",
  adaptive = TRUE,
  kernel = "bisquare",
  scale = TRUE,
  max_iter = 30,
  tol = 1e-5
)

summary(fit_pmgwr)
summary(fit_mgwfer)

fit_pmgwr$bandwidths
fit_mgwfer$bandwidths
```

Check the convergence status and iteration count in each summary. Increase `max_iter` if a model has not converged.

## Map coefficient surfaces

```r
plot_coefficients(
  fit_mgwfer,
  variables = c("x1", "x2", "x3", "x4"),
  tile_width = 1,
  tile_height = 1
)

plot_context_effects(
  fit_mgwfer,
  tile_width = 1,
  tile_height = 1
)
```

For polygon maps, provide one `sf` geometry per panel unit through `spatial`. Geometry matching uses the location ID, and the original polygon boundaries are preserved in the returned map data.

```r
# unit_sf must contain one geometry per unit and a column named unit_id
fit_sf <- mgwpr(
  y ~ x1 + x2 + x3 + x4,
  data = panel,
  spatial = unit_sf,
  id = "unit_id",
  time = "time_id",
  method = "fe",
  family = "gaussian"
)

map_data <- mgwpr_sf(fit_sf)
plot_coefficients(fit_sf, variables = c("x1", "x2"))
```

## Compare estimated coefficients with known surfaces

When the true coefficient surfaces are known, `compare_coefficient_surfaces()` returns location-level estimates and accuracy metrics.

```r
comparison <- compare_coefficient_surfaces(
  pooled = fit_pmgwr,
  fixed = fit_mgwfer,
  truth = units
)

head(comparison$coefficients)
comparison$metrics
```

To map a true surface alongside estimates:

```r
plot_mgwpr_comparison(
  pooled = fit_pmgwr,
  fixed = fit_mgwfer,
  variable = "x1",
  truth = setNames(units$beta1_true, units$unit_id),
  tile_width = 1,
  tile_height = 1
)
```

## Coefficient inference and significance maps

For Gaussian fits, `coefficient_inference()` returns local t-values, effective parameter counts, and degrees of freedom. `coefficient_significance()` applies the multiscale effective-number-of-parameters correction and returns results that can be mapped.

```r
inference <- coefficient_inference(fit_mgwfer)

sig <- coefficient_significance(
  object = fit_mgwfer,
  t_values = inference$t_values,
  enp = inference$enp,
  alpha = 0.05,
  df = inference$df
)

plot_coefficient_significance(sig, ncol = 2)
```

These calculations condition on the selected bandwidths and do not account for bandwidth-selection uncertainty.

## Spatial-panel diagnostics

`diagnose_mgwpr()` runs selected diagnostics and returns their results in a list. The following example requests Moran's I, local F, and LM tests.

```r
diagnostics <- diagnose_mgwpr(
  fit_mgwfer,
  tests = c("moran", "f_test", "lm_test"),
  k = 4,
  permutations = 999,
  seed = 123
)

diagnostics$moran[c("by_period", "statistic", "p_value")]
head(diagnostics$f_test)

lm_results <- diagnostics$lm_test
lm_results$p_BH <- p.adjust(lm_results$p_value, method = "BH")
head(lm_results)
```

## Cross-sectional MGWR with GWmodel

`GWmodel` is an optional dependency required only for the cross-sectional `mgwr()` interface:

```r
install.packages("GWmodel")
```

`mgwr()` wraps `GWmodel::gwr.multiscale()` for cross-sectional data. Supply one row per location; for a panel, select a period or aggregate before fitting. GWmodel selects a bandwidth for each model term and then optimizes the term-specific bandwidths during back-fitting, which can take time.

```r
cross_section <- panel[panel$time_id == 2019, ]

fit_mgwr <- mgwr(
  y ~ x1 + x2 + x3 + x4,
  data = cross_section,
  coords = c("grid_i", "grid_j"),
  id = "unit_id",
  adaptive = TRUE,
  kernel = "bisquare",
  approach = "AICc",
  criterion = "dCVR",
  max_iter = 2000,
  tol = 1e-5
)

summary(fit_mgwr)
fit_mgwr$bandwidths
plot_coefficients(fit_mgwr, variables = c("x1", "x2", "x3", "x4"))
```

## References

Li, Z. & Fotheringham, A. S. (2026). Spatial Context as a Time-Invariant Confounder: A Fixed-Effects Extension of MGWR. *Annals of the American Association of Geographers*, 116(9), 2294–2318. <https://doi.org/10.1080/24694452.2026.2654481>

Python reference implementation: <https://github.com/GeoZhipengLi/MGWPR>
