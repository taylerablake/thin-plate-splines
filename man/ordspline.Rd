\name{ordspline}
\alias{ordspline}
\title{
Fits Ordinal Smoothing Spline
}
\description{
Given a real-valued response vector \eqn{\mathbf{y}=\{y_{i}\}_{n\times1}} and an ordinal predictor vector \eqn{\mathbf{x}=\{x_{i}\}_{n\times 1}} with \eqn{x_{i} \in \{1,\ldots,K\} \ \forall i}, an ordinal smoothing spline model has the form \deqn{y_{i}=\eta(x_{i})+e_{i}} where \eqn{y_{i}} is the \eqn{i}-th observation's respone, \eqn{x_{i}} is the \eqn{i}-th observation's predictor, \eqn{\eta} is an unknown function relating the response and predictor, and \eqn{e_{i}\sim\mathrm{N}(0,\sigma^{2})} is iid Gaussian error.
}
\usage{
ordspline(x, y, knots, weights, lambda, monotone=FALSE)
}
\arguments{
  \item{x}{
	Predictor vector.
}
  \item{y}{
  Response vector. Must be same length as \code{x}.
}
  \item{knots}{
	Either a scalar giving the number of equidistant knots to use, or a vector of values to use as the spline knots. If left blank, the number of knots is \code{min(50, nu)} where \code{nu = length(unique(x)).}
}
  \item{weights}{
	Weights vector (for weighted penalized least squares). Must be same length as \code{x} and contain non-negative values.
}
  \item{lambda}{
  Smoothing parameter. If left blank, \code{lambda} is tuned via Generalized Cross-Validation.
}
  \item{monotone}{
  If \code{TRUE}, the relationship between \code{x} and \code{y} is constrained to be monotonic increasing.
}
}
\details{
  To estimate \eqn{\eta} I minimize the penalized least-squares functional \deqn{\frac{1}{n}\sum_{i=1}^{n}(y_{i}-\eta(x_{i}))^{2}+\lambda \sum_{x=2}^K [\eta(x)-\eta(x-1)]^2 dx} where \eqn{\lambda\geq0} is a smoothing parameter that controls the trade-off between fitting and smoothing the data. 
  
  Default use of the function estimates \eqn{\lambda} by minimizing the GCV score: \deqn{\mbox{GCV}(\lambda) = \frac{n\|(\mathbf{I}_{n}-\mathbf{S}_{\lambda})\mathbf{y}\|^{2}}{[n-\mathrm{tr}(\mathbf{S}_{\lambda})]^2}} where \eqn{\mathbf{I}_{n}} is the identity matrix and \eqn{\mathbf{S}_{\lambda}} is the smoothing matrix.
}
\value{
	\item{fitted.values}{Vector of fitted values.}
	\item{se.fit}{Vector of standard errors of \code{fitted.values}.}
  \item{sigma}{Estimated error standard deviation, i.e., \eqn{\hat{\sigma}}.}
  \item{lambda}{Chosen smoothing parameter.}
  \item{info}{Model fit information: vector containing the GCV, R-squared, AIC, and BIC of fit model (assuming Gaussian error).}
  \item{coef}{Spline basis function coefficients.}
  \item{coef.csqrt}{Matrix square-root of covariace matrix of \code{coef}. Use \code{tcrossprod(coef.csqrt)} to get covariance matrix of \code{coef}.}
  \item{n}{Number of data points, i.e., \code{length(x)}.}
  \item{df}{Effective degrees of freedom (trace of smoothing matrix).}
  \item{xunique}{Unique elements of \code{x}.}
  \item{x}{Predictor vector (same as input).}
	\item{y}{Response vector (same as input).}
	\item{residuals}{Residual vector, i.e., \code{y - fitted.values}.}
	\item{knots}{Spline knots used for fit.}
	\item{monotone}{Logical (same as input).}
}
\references{
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.

  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.
	  
	Helwig, N. E. and Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics, 24}, 715-732.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
	The spline is estimated using penalized least-squares, which does not require the Gaussian error assumption. However, the spline inference information (e.g., standard errors and fit information) requires the Gaussian error assumption.
}
\section{Warnings }{
	When inputting user-specified \code{knots}, all values in \code{knots} must match a corresponding value in \code{x}.
}

\examples{

##########   EXAMPLE   ##########

# generate some data
n <- 100
nk <- 50
x <- seq(-3,3,length.out=n)
eta <- (sin(2*x/pi) + 0.25*x^3 + 0.05*x^5)/15
set.seed(1)
y <- eta + rnorm(n, sd=0.5)

# plot data and true eta
plot(x, y)
lines(x, eta, col="blue", lwd=2)

# fit ordinal smoothing spline
ossmod <- ordspline(x, y, knots=nk)
lines(ossmod$x, ossmod$fit, col="red", lwd=2)

# fit monotonic smoothing spline
mssmod <- ordspline(x, y, knots=nk, monotone=TRUE)
lines(mssmod$x, mssmod$fit, col="purple", lwd=2)

}

