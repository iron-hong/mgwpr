\name{coefficient_inference}
\alias{coefficient_inference}
\alias{mgwpr_tvalues}
\alias{mgwpr_enp}
\title{Calculate Conditional Local Inference for MGWPR}
\description{
Calculates local coefficient t-values, panel-unit cluster sandwich standard
errors, and coefficient-specific effective parameter counts from a fitted
MGWPR model while holding the selected bandwidths fixed.
}
\usage{
coefficient_inference(object, tol = 1e-8, max_iter = 500L,
  finite_sample = TRUE)

mgwpr_tvalues(object, ...)

mgwpr_enp(object, ...)
}
\arguments{
  \item{object}{An object returned by \code{mgwpr()}.}
  \item{tol}{Convergence tolerance for linear smoother backfitting.}
  \item{max_iter}{Maximum number of smoother backfitting iterations.}
  \item{finite_sample}{Whether to apply a finite-sample correction to the
  panel-unit cluster sandwich variance.}
  \item{...}{Additional arguments passed to \code{coefficient_inference()}.}
}
\value{
\code{coefficient_inference()} returns a list containing the location by
coefficient t-value matrix, standard-error matrix, named ENP vector, degrees
of freedom, iteration count, and convergence status. \code{mgwpr_tvalues()}
and \code{mgwpr_enp()} return the corresponding matrix and vector.
}
\details{
The coefficient-specific effective parameter count is calculated as
\code{tr(H_j)} from the term-specific linear smoother operator. Standard
errors use a panel-unit cluster sandwich variance. Both
calculations condition on the selected bandwidths and do not account for
bandwidth-selection uncertainty.
}
\seealso{\code{coefficient_significance()}, \code{mgwpr()}}
