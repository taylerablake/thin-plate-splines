\name{predict.ordspline}
\alias{predict.ordspline}
\title{
Predicts for "ordspline" Objects
}
\description{
Get fitted values and standard error estimates for ordinal smoothing splines.
}
\usage{
\method{predict}{ordspline}(object,newdata=NULL,se.fit=FALSE,...)
}
\arguments{
  \item{object}{
	Object of class "ordspline", which is output from \code{\link{ordspline}}.
}
  \item{newdata}{
	Vector containing new data points for prediction. See Details and Example. Default of \code{newdata=NULL} uses original data in \code{object} input.
}
  \item{se.fit}{
	Logical indicating whether the standard errors of the fitted values should be estimated. Default is \code{se.fit=FALSE}.
}
  \item{\dots}{
	Ignored.
}
}
\value{
  If \code{se.fit=FALSE}, returns vector of fitted values.

  Otherwise returns list with elements:
	\item{fit}{Vector of fitted values}
  \item{se.fit}{Vector of standard errors of fitted values}
}
\references{
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.

  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.

	Helwig, N. E. and Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics, 24}, 715-732.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
	Uses the coefficient and smoothing parameter estimates from a fit ordinal smoothing spline (estimated by \code{\link{ordspline}}) to predict for new data.
}

\examples{

##########   EXAMPLE   ##########

# define univariate function and data
set.seed(773)
myfun <- function(x){ 2 + x/2 + sin(x) }
x <- sample(1:20, size=500, replace=TRUE)
y <- myfun(x) + rnorm(500)

# fit ordinal spline model
ordmod <- ordspline(x, y)
monmod <- ordspline(x, y, monotone=TRUE)
crossprod( predict(ordmod) - myfun(x) ) / 500
crossprod( predict(monmod) - myfun(x) ) / 500

# plot truth and predictions
ordfit <- predict(ordmod, 1:20, se.fit=TRUE)
monfit <- predict(monmod, 1:20, se.fit=TRUE)
plotci(1:20, ordfit$fit, ordfit$se.fit, ylab="f(x)")
plotci(1:20, monfit$fit, monfit$se.fit, col="red", col.ci="pink", add=TRUE)
points(1:20, myfun(1:20))

}
