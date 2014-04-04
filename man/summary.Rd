\name{summary}
\alias{summary.css}
\alias{summary.ssa}
\alias{summary.ssp}
\alias{summary.tps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes fit information for big splines model
}
\description{
This function summarizes basic model fit information for a fit \code{bigsplines} model.
}
\usage{
\method{summary}{css}(object,fitresid=TRUE,chunksize=10000,...)
\method{summary}{ssa}(object,fitresid=TRUE,chunksize=10000,...)
\method{summary}{ssp}(object,fitresid=TRUE,chunksize=10000,...)
\method{summary}{tps}(object,fitresid=TRUE,chunksize=10000,...)
}
\arguments{
  \item{object}{
  	Object of class "css" (output from \code{\link{bigspline}}), class "ssa" (output from \code{\link{bigssa}}), class "ssp" (output from \code{\link{bigssp}}), or class "tps" (output from \code{\link{bigtps}}).
}
  \item{fitresid}{
  Logical indicating whether the fitted values and residuals should be calculated for all data points in input \code{object}. 
}
  \item{chunksize}{
  If \code{fitresid=TRUE}, fitted values are calculated in chunks of size \code{chunksize}.
}
  \item{\ldots}{
  Ignored.
}
}
\value{
	\item{call}{Called model in input \code{formula}.}
  \item{type}{Type of smoothing spline that was used for each predictor.}
  \item{fitted.values}{Vector of fitted values (if \code{fitresid=TRUE}).}
  \item{residuals}{Vector of residuals (if \code{fitresid=TRUE}).}
  \item{sigma}{Estimated error standard deviation.}
  \item{n}{Total sample size.}
  \item{df}{Effective degrees of freedom of the model.}
  \item{info}{Model fit information: vector containing the GCV, multiple R-squared, AIC, and BIC of fit model.}
  \item{converged}{Convergence status: \code{converged=TRUE} if the iterative theta update converged, \code{converged=FALSE} if the iterative theta update failed to converge, and \code{converged=NA} if option \code{skip.iter=TRUE} was used.}
  \item{iter}{Number of iterative updates (\code{iter=NA} if option \code{skip.iter=TRUE} was used).}
  \item{rparm}{Rounding parameters used for model fitting.}
  \item{lambda}{Global smoothing parameter used for model fitting.}
  \item{gammas}{Vector of additional smoothing parameters (only for class "ssa").}
  \item{thetas}{Vector of additional smoothing parameters (only for class "ssp").}
}
\details{
  See \code{\link{bigspline}}, \code{\link{bigssa}}, \code{\link{bigssp}}, and \code{\link{bigtps}} for more details.
}
\author{
Nathaniel E. Helwig <nhelwig2@illinois.edu>
}
\note{
For "css" and "tps" objects, the outputs \code{call}, \code{converged}, and \code{iter} are NA.
}

\examples{

##########   EXAMPLE 1   ##########

# define relatively smooth function
set.seed(773)
myfun<-function(x){sin(2*pi*x)}
x<-runif(10^4)
y<-myfun(x)+rnorm(10^4)

# cubic spline
cubmod<-bigspline(x,y)
summary(cubmod)


##########   EXAMPLE 2   ##########

# function with two continuous predictors
set.seed(773)
myfun<-function(x1v,x2v){sin(2*pi*x1v)+log(x2v+.1)+cos(pi*(x1v-x2v))}
x1v<-runif(10^4);     x2v<-runif(10^4)
y<-myfun(x1v,x2v)+rnorm(10^4)

# cubic splines with 100 randomly selected knots (efficient parameterization)
cubmod<-bigssa(y~x1v*x2v,type=list(x1v="cub",x2v="cub"),nknots=100)
summary(cubmod)


##########   EXAMPLE 3   ##########

# function with two continuous predictors
set.seed(773)
myfun<-function(x1v,x2v){sin(2*pi*x1v)+log(x2v+.1)+cos(pi*(x1v-x2v))}
x1v<-runif(10^4);     x2v<-runif(10^4)
y<-myfun(x1v,x2v)+rnorm(10^4)

# cubic splines with 100 randomly selected knots (classic parameterization)
cubmod<-bigssp(y~x1v*x2v,type=list(x1v="cub",x2v="cub"),nknots=100)
summary(cubmod)


##########   EXAMPLE 4   ##########

# define relatively smooth function
set.seed(773)
myfun<-function(x){sin(2*pi*x)}
x<-runif(10^4)
y<-myfun(x)+rnorm(10^4)

# thin-plate with default (30 knots)
tpsmod<-bigtps(x,y)
summary(tpsmod)

}

