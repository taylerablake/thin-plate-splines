\name{ssBasis}
\alias{ssBasis}
\title{
Smoothing Spline Basis for Polynomial Splines
}
\description{
Generate the smoothing spline basis matrix for a polynomial spline.
}
\usage{
ssBasis(x, knots, m=2, d=0, xmin=min(x), xmax=max(x), periodic=FALSE, intercept=FALSE)
}
\arguments{
  \item{x}{
Predictor variable.
}
  \item{knots}{
Spline knots.
}
  \item{m}{
Penalty order. 'm=1' for linear smoothing spline, 'm=2' for cubic, and 'm=3' for quintic.
}
  \item{d}{
Derivative order. 'd=0' for smoothing spline basis, 'd=1' for 1st derivative of basis, and 'd=2' for 2nd derivative of basis.
}
\item{xmin}{
Minimum value of 'x'.
}
\item{xmax}{
Maximum value of 'x'.
}
\item{periodic}{
If \code{TRUE}, the smoothing spline basis is periodic w.r.t. the interval [\code{xmin}, \code{xmax}].
}
\item{intercept}{
If \code{TRUE}, the first column of the basis will be a column of ones.
}
}
\value{
\item{X}{Spline Basis.}
\item{knots}{Spline knots.}
\item{m}{Penalty order.}
\item{d}{Derivative order.}
\item{xlim}{Inputs \code{xmin} and \code{xmax}.}
\item{periodic}{Same as input.}
\item{intercept}{Same as input.}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Inputs \code{x} and \code{knots} should be within the interval [\code{xmin}, \code{xmax}].
}

\examples{
##########   EXAMPLE   ##########

# define function and its derivatives
n <- 500
x <- seq(0, 1, length.out=n)
knots <- seq(0, 1, length=20)
y <- sin(4 * pi * x)
d1y <- 4 * pi * cos(4 * pi * x)
d2y <- - (4 * pi)^2 * sin(4 * pi * x)

# linear smoothing spline
linmat0 <- ssBasis(x, knots, m=1)
lincoef <- pinvsm(crossprod(linmat0$X)) \%*\% crossprod(linmat0$X, y)
linyhat <- linmat0$X \%*\% lincoef
linmat1 <- ssBasis(x, knots, m=1, d=1)
linyd1 <- linmat1$X \%*\% lincoef

# plot linear smoothing spline results
par(mfrow=c(1,2))
plot(x, y, type="l", main="Function")
lines(x, linyhat, lty=2, col="red")
plot(x, d1y, type="l", main="First Derivative")
lines(x, linyd1, lty=2, col="red")

# cubic smoothing spline
cubmat0 <- ssBasis(x, knots)
cubcoef <- pinvsm(crossprod(cubmat0$X)) \%*\% crossprod(cubmat0$X, y)
cubyhat <- cubmat0$X \%*\% cubcoef
cubmat1 <- ssBasis(x, knots, d=1)
cubyd1 <- cubmat1$X \%*\% cubcoef
cubmat2 <- ssBasis(x, knots, d=2)
cubyd2 <- cubmat2$X \%*\% cubcoef

# plot cubic smoothing spline results
par(mfrow=c(1,3))
plot(x, y, type="l", main="Function")
lines(x, cubyhat, lty=2, col="red")
plot(x, d1y, type="l", main="First Derivative")
lines(x, cubyd1, lty=2, col="red")
plot(x, d2y, type="l", main="Second Derivative")
lines(x, cubyd2, lty=2, col="red")

# quintic smoothing spline
quimat0 <- ssBasis(x, knots, m=3)
quicoef <- pinvsm(crossprod(quimat0$X)) \%*\% crossprod(quimat0$X, y)
quiyhat <- quimat0$X \%*\% quicoef
quimat1 <- ssBasis(x, knots, m=3, d=1)
quiyd1 <- quimat1$X \%*\% quicoef
quimat2 <- ssBasis(x, knots, m=3, d=2)
quiyd2 <- quimat2$X \%*\% quicoef

# plot quintic smoothing spline results
par(mfrow=c(1,3))
plot(x, y, type="l", main="Function")
lines(x, quiyhat, lty=2, col="red")
plot(x, d1y, type="l", main="First Derivative")
lines(x, quiyd1, lty=2, col="red")
plot(x, d2y, type="l", main="Second Derivative")
lines(x, quiyd2, lty=2, col="red")

}
