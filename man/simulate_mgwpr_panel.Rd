\name{simulate_mgwpr_panel}
\alias{simulate_mgwpr_panel}
\title{Simulate a spatial panel with known coefficient surfaces}
\usage{
simulate_mgwpr_panel(n_grid = 15L, n_time = 3L, sigma_x = 0.5,
  spatial_coupling = 0.05, sigma_epsilon = 0.5, context_scale = 30,
  time_values = NULL, seed = NULL)
}
\arguments{
  \item{n_grid}{Number of rows and columns in the square spatial grid.}
  \item{n_time}{Number of periods per spatial unit.}
  \item{sigma_x}{Standard deviation of the independent component of each predictor.}
  \item{spatial_coupling}{Strength of the time-invariant context effect in each predictor.}
  \item{sigma_epsilon}{Standard deviation of the outcome disturbance.}
  \item{context_scale}{Multiplier controlling the scale of the true spatial context effect.}
  \item{time_values}{Optional vector of distinct period labels. Defaults to 0 through n_time - 1.}
  \item{seed}{Optional non-negative integer random seed. When supplied, the caller's existing random-number state is restored.}
}
\description{
Generates a regular spatial grid, known location-specific coefficients, a
time-invariant spatial context effect, four spatially confounded predictors,
and a continuous outcome in long panel format.
}
\value{
A list of class mgwpr_simulation with units (the location-level coefficient
truth table), panel (long-format observations), diagnostics, and settings.
}
\examples{
sim <- simulate_mgwpr_panel(seed = 42)
head(sim$panel)
head(sim$units)
}
\seealso{mgwpr, compare_coefficient_surfaces}
\keyword{datasets}
\keyword{spatial}
