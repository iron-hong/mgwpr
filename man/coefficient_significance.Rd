\name{coefficient_significance}
\alias{coefficient_significance}
\alias{plot_coefficient_significance}
\title{Correct and Map Local Coefficient Significance}
\description{
Applies the multiscale effective-number-of-parameters correction of da Silva
and Fotheringham (2016) to supplied local t-values, or to t-values computed
from supplied standard errors. The plotting function maps positive,
negative, and non-significant coefficient estimates.
}
\usage{
coefficient_significance(object, t_values = NULL, standard_errors = NULL,
  enp, alpha = 0.05, df = Inf)

plot_coefficient_significance(x, ncol = 2, tile_width = NULL,
  tile_height = NULL, xlab = "X coordinate", ylab = "Y coordinate")
}
\arguments{
  \item{object}{An object returned by \code{mgwpr()}.}
  \item{t_values}{Optional numeric matrix or data frame of local t-values,
  with fitted locations in rows and coefficients in columns. Supply either
  this or \code{standard_errors}.}
  \item{standard_errors}{Optional numeric matrix or data frame of coefficient
  standard errors with the same dimensions and dimnames as the fitted
  coefficient matrix.}
  \item{enp}{Positive finite effective number of parameters, one per
  coefficient. Names must match coefficient names, unless the vector is in
  coefficient order.}
  \item{alpha}{Family-wise target significance level between zero and one.}
  \item{df}{Positive degrees of freedom used for the two-sided critical
  t-value. The default \code{Inf} uses the normal critical value. May be a
  named vector by coefficient.}
  \item{x}{An object returned by \code{coefficient_significance()}.}
  \item{ncol}{Number of facet columns.}
  \item{tile_width, tile_height}{Optional tile dimensions in coordinate
  units.}
  \item{xlab, ylab}{Axis labels.}
}
\value{
\code{coefficient_significance()} returns a long data frame with estimates,
t-values, adjusted significance levels, critical values, significance flags,
directions, and filtered t-values. \code{plot_coefficient_significance()}
returns a ggplot object.
}
\details{
The per-coefficient adjusted significance level is \code{alpha / ENP_j}.
This function does not estimate standard errors, t-values, or effective
numbers of parameters. Those quantities must be obtained from an inferential
procedure appropriate to the fitted panel model.
}
\references{
da Silva, A. R. and Fotheringham, A. S. (2016). The multiple testing issue in
geographically weighted regression. \emph{Geographical Analysis}, 48(3),
233-247. doi:10.1111/gean.12084.
}
\seealso{\code{mgwpr()}, \code{plot_coefficients()}}
