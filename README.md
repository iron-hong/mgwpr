# MGWPR

**MGWPR** is an R package for multiscale geographically weighted panel
regression. It provides pooled MGWR and individual fixed-effects MGWFER
estimators using spatially varying coefficients and predictor-specific
bandwidths.

## Installation

```r
install.packages(c("ggplot2", "sf"))
install.packages("remotes")
remotes::install_local("MGWPR_0.0.1.tar.gz")
```

Alternatively, from the unpacked source directory:

```sh
R CMD INSTALL MGWPR
```

## Cross-sectional MGWR with GWmodel

`mgwr()` provides a cross-sectional multiscale GWR interface through
`GWmodel::gwr.multiscale()`. Install `GWmodel` to use this optional estimator.
It selects a bandwidth for each model-matrix term, then optimizes the
term-specific bandwidths during back-fitting; bandwidth selection can take
time. Use one observation per spatial location (for panel data, select a
single period or aggregate first). This interface is distinct from the
panel-data `mgwpr()` estimators.

```r
install.packages("GWmodel")

# One cross-section: one row per unit
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
  tol = 1e-5,
  verbose = TRUE
)

summary(fit_mgwr)
fit_mgwr$bandwidths
plot_coefficients(fit_mgwr, variables = c("x1", "x2", "x3", "x4"))
```

For polygon maps, pass an `sf` object with one geometry per location through
`spatial`; `plot_coefficients()` then draws the coefficient surface over the
polygons. `fit_mgwr$gwmodel_fit` retains the original GWmodel result.

## Data requirements

Data should be in long format, with one row per unit and time period. The two
coordinate columns must be numeric and constant within a panel unit. Each
location needs at least four distinct units. Missing values are not accepted;
prepare a complete analysis sample before fitting.

Instead of coordinate columns, `coords` can be replaced by `spatial`, an `sf`
object with one geometry and the matching unit ID per panel unit. Polygon
representative points are used to calculate inter-location distances; the
original geometry is kept for maps and for `mgwpr_sf()` output.

## Simulate panel data

`simulate_mgwpr_panel()` generates a regular grid with known spatial
coefficient surfaces, a time-invariant context effect, four correlated
predictors, and a continuous response. It returns a list containing `units`
(the location-level truth table), `panel` (long-format observations),
simulation diagnostics, and settings. The defaults reproduce the 15 by 15,
three-period example.

```r
sim <- simulate_mgwpr_panel(seed = 42)
units <- sim$units
panel <- sim$panel
sim$diagnostics

fit_mgwfer <- mgwpr(
  y ~ x1 + x2 + x3 + x4,
  data = panel,
  coords = c("grid_i", "grid_j"),
  id = "unit_id",
  time = "time_id",
  method = "fe",
  family = "gaussian",
  adaptive = TRUE,
  kernel = "bisquare",
  scale = TRUE
)
```

## Example

```r
library(MGWPR)

# panel_data contains outcome, time-varying predictors, id, period, east, north
fit <- mgwpr(
  outcome ~ investment + population + service_access,
  data = panel_data,
  coords = c("east", "north"),
  id = "id",
  time = "period",
  method = "fe",
  adaptive = TRUE,
  kernel = "bisquare",
  max_iter = 30,
  verbose = TRUE
)

summary(fit)
head(coef(fit))
head(fitted(fit))

# Coefficient-surface map
plot_coefficients(fit, variables = c("investment", "population"))

# Context-effect map for an FE fit
plot_context_effects(fit)
```

Use `family = "binomial"` with a pooled model for binary outcomes. The
response must be coded 0/1 or supplied as a two-level factor. The positive
factor class is the second factor level. Binomial fixed-effects estimation is
not currently supported. Available kernels are `bisquare`, `gaussian`,
`exponential`, `tricube`, and `boxcar`.

```r
fit_binary <- mgwpr(
  adopted ~ income + age + access,
  data = panel_data,
  spatial = unit_polygons,
  id = "id",
  time = "period",
  method = "pooled",
  family = "binomial",
  kernel = "tricube"
)

# Coefficients and context effects joined to sf geometry
map_data <- mgwpr_sf(fit_binary)
plot_coefficients(fit_binary, variables = c("income", "age"))
```

### Coefficient significance maps

For Gaussian fits, `coefficient_inference()` calculates the local t-value matrix, panel-unit
cluster sandwich standard errors, and `ENP_j`. `coefficient_significance()`
then applies the multiscale effective-number-of-parameters correction of da
Silva and Fotheringham (2016). The adjusted level is `alpha / ENP_j`; the
result includes the critical t-value, significance classification, and
filtered t-value (set to zero when not significant). The result can be mapped
with `plot_coefficient_significance()`.

```r
inference <- coefficient_inference(fit)
t_matrix <- inference$t_values
enp_by_term <- inference$enp

sig <- coefficient_significance(
  fit,
  t_values = t_matrix,
  enp = enp_by_term,
  alpha = 0.05,
  df = inference$df
)

head(sig)
plot_coefficient_significance(sig, ncol = 2)
```

The inference calculations use fixed-bandwidth linear smoother operators and
a panel-unit cluster sandwich variance. They condition on the selected
bandwidths and do not include bandwidth-selection uncertainty. The reported
degrees of freedom are the number of fitted locations minus one, so interpret
the results as conditional local inference.

Use `method = "pooled"` to fit a pooled model without unit fixed effects.
For fixed-effects models, time-invariant predictors are removed by the within
transformation. Bandwidths are selected separately for each model-matrix
coefficient using leave-one-location-out cross-validation. Set `bw` to a
positive scalar or coefficient-length vector to supply bandwidths directly.
Adaptive bandwidth values count nearest distinct panel locations; fixed
bandwidth values use the coordinate units.

## Estimation notes

The estimator applies iterative univariate backfitting to spatial coefficient
surfaces with one of five spatial kernels. Fixed effects are
removed by subtracting each unit's time mean from the outcome and regressors.
The pooled and fixed-effects functions return location-level coefficients,
recovered location-specific context effects for fixed-effects fits,
observation-level fitted values and residuals, bandwidths, and convergence
status.

The panel estimators are a new R implementation of the documented MGWPR
approach, not a port of the repository's Python source. They use CV bandwidth
selection; the separate cross-sectional `mgwr()` interface delegates MGWR
estimation to GWmodel. `coefficient_inference()` derives conditional local inference
but does not account for bandwidth-selection uncertainty or provide
prediction intervals. The reported R-squared and RMSE are descriptive fit
summaries. Compare results with published reference implementations before
using them for substantive inference.

## Coefficient-surface accuracy

When simulated or otherwise known coefficient surfaces are available,
compare_coefficient_surfaces() aligns the fitted coefficients and truth by
location ID. It returns a coefficients table and a metrics table containing
PMGWR and MGWFER RMSE and correlation by coefficient. By default, terms named
x1, x2, etc. are matched to truth columns beta1_true, beta2_true, etc.

```r
comparison <- compare_coefficient_surfaces(
  pooled = fit_pmgwr,
  fixed = fit_mgwfer,
  truth = units
)

comparison$coefficients
comparison$metrics

# Return only the metrics table
coef_metrics(
  pooled = fit_pmgwr,
  fixed = fit_mgwfer,
  truth = units
)
```

## Spatial plots

The package provides three map functions. plot_coefficients() maps one or more
coefficient surfaces from a fitted model. plot_context_effects() maps the
recovered unit-level context effects from an FE model and can optionally show
a truth vector. plot_mgwpr_comparison() compares one coefficient surface from
pooled and FE fits, optionally with a true surface. The ordinary plot() method
also maps coefficients by default; use plot(fit, type = "context") for context
effects. With `spatial = sf_object`, maps use the supplied polygon geometries.
With coordinate columns, maps use filled tiles centered on fitted coordinates;
tile dimensions are inferred from coordinate spacing by default. Set
tile_width and tile_height explicitly for regular grids with known cell sizes.
The plotting functions use ggplot2.

```r
plot(fit)
plot(fit, type = "context")

# All three panels share a scale for this coefficient
plot_mgwpr_comparison(
  fit_pmgwr,
  fit_mgwfer,
  variable = "x1",
  truth = setNames(units$beta1_true, units$unit_id)
)
```

## Diagnostics

`diagnose_mgwpr()` provides residual Moran's I, local F tests comparing pooled
and unit fixed-effects fits, local Hausman comparisons of within and random
slopes, and local Breusch-Pagan LM tests for random intercepts. Moran weights
default to inverse-distance k-nearest neighbours, and permutation replicates
provide a reference p-value. For binomial fits, Moran's I uses Pearson
residuals and is exploratory; local F, Hausman, and LM are restricted to
Gaussian outcomes. The local Hausman test uses a Swamy-Arora-style local
variance-components approximation. LM calculations use balanced local panels;
unbalanced neighbourhoods are explicitly marked as not testable.

```r
diagnostics <- diagnose_mgwpr(
  fit,
  tests = c("moran", "f_test", "hausman", "lm_test"),
  k = 4,
  permutations = 999,
  seed = 123
)

diagnostics$moran$by_period
diagnostics$f_test
diagnostics$hausman
diagnostics$lm_test
```

## Method reference

Li, Z. & Fotheringham, A. S. (2026). Spatial Context as a Time-Invariant
Confounder: A Fixed-Effects Extension of MGWR. *Annals of the American
Association of Geographers*, 116(9), 2294–2318.
https://doi.org/10.1080/24694452.2026.2654481

Reference implementation: https://github.com/GeoZhipengLi/MGWPR
