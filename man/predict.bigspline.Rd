\name{predict.bigspline}
\alias{predict.bigspline}
\title{
Predicts for "bigspline" Objects
}
\description{
Get fitted values and standard error estimates for cubic smoothing splines.
}
\usage{
\method{predict}{bigspline}(object,newdata=NULL,se.fit=FALSE,
        effect=c("all","0","lin","non"),
        design=FALSE,smoothMatrix=FALSE,...)
}
\arguments{
  \item{object}{
	Object of class "bigspline", which is output from \code{\link{bigspline}}.
}
  \item{newdata}{
	Vector containing new data points for prediction. See Details and Example. Default of \code{newdata=NULL} uses original data in \code{object} input.
}
  \item{se.fit}{
	Logical indicating whether the standard errors of the fitted values should be estimated. Default is \code{se.fit=FALSE}.
}
\item{effect}{
  Which effect to estimate: \code{effect="all"} gives full \eqn{\hat{y}}, \code{effect="0"} gives the intercept (constant) portion of \eqn{\hat{y}},  \code{effect="lin"} gives linear portion of \eqn{\hat{y}}, and \code{effect="non"} gives nonlinear portion of \eqn{\hat{y}}. 
}
  \item{design}{
  Logical indicating whether the design matrix should be returned.
}
  \item{smoothMatrix}{
  Logical indicating whether the smoothing matrix should be returned.
}
  \item{\dots}{
	Ignored.
}
}
\value{
  If \code{se.fit=FALSE}, \code{design=FALSE}, and \code{smoothMatrix=FALSE}, returns vector of fitted values.

  Otherwise returns list with elements:
	\item{fit}{Vector of fitted values}
  \item{se.fit}{Vector of standard errors of fitted values (if \code{se.fit=TRUE})}
  \item{X}{Design matrix used to create fitted values (if \code{design=TRUE})}
  \item{ix}{Index vector such that \code{fit=X\%*\%object$coef[ix]} (if \code{design=TRUE})}
  \item{S}{Smoothing matrix corresponding to fitted values (if \code{smoothMatrix=TRUE})}
  
}
\references{
  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.

	Helwig, N. E. and Ma, P. (in press). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
  
  Helwig, N. E. and Ma, P. (in revision). Smoothing spline ANOVA for super-large samples: Scalable computation via rounding parameters. \emph{Statistics and Its Interface}.
  
  Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
  
  Helwig, N. E. and Ma, P. (in preparation). Optimal eigenbasis truncation for thin-plate regression splines. 
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
	Uses the coefficient and smoothing parameter estimates from a fit cubic smoothing spline (estimated by \code{\link{bigspline}}) to predict for new data.
}

\examples{

##########   EXAMPLE 1   ##########

# define univariate function and data
set.seed(773)
myfun <- function(x){ 2 + x + sin(2*pi*x) }
x <- runif(10^4)
y <- myfun(x) + rnorm(10^4)

# fit cubic spline model
cubmod <- bigspline(x,y)
crossprod( predict(cubmod) - myfun(x) )/10^4

# define new data for prediction
newdata <- data.frame(x=seq(0,1,length.out=100))

# get fitted values and standard errors for new data
yc <- predict(cubmod,newdata,se.fit=TRUE)

# plot results with 95\% Bayesian confidence interval
plot(newdata$x,yc$fit,type="l")
lines(newdata$x,yc$fit+qnorm(.975)*yc$se.fit,lty=3)
lines(newdata$x,yc$fit-qnorm(.975)*yc$se.fit,lty=3)

# predict constant, linear, and nonlinear effects
yc0 <- predict(cubmod,newdata,se.fit=TRUE,effect="0")
ycl <- predict(cubmod,newdata,se.fit=TRUE,effect="lin")
ycn <- predict(cubmod,newdata,se.fit=TRUE,effect="non")
crossprod( yc$fit - (yc0$fit + ycl$fit + ycn$fit) )

# plot results with 95\% Bayesian confidence intervals
par(mfrow=c(1,2))
plot(newdata$x,ycl$fit,type="l",main="Linear effect")
lines(newdata$x,ycl$fit+qnorm(.975)*ycl$se.fit,lty=3)
lines(newdata$x,ycl$fit-qnorm(.975)*ycl$se.fit,lty=3)
plot(newdata$x,ycn$fit,type="l",main="Nonlinear effect")
lines(newdata$x,ycn$fit+qnorm(.975)*ycn$se.fit,lty=3)
lines(newdata$x,ycn$fit-qnorm(.975)*ycn$se.fit,lty=3)


##########   EXAMPLE 2   ##########

# define (same) univariate function and data
set.seed(773)
myfun <- function(x){ 2 + x + sin(2*pi*x) }
x <- runif(10^4)
y <- myfun(x) + rnorm(10^4)

# fit a different cubic spline model
cubamod <- bigspline(x,y,type="cub0")
crossprod( predict(cubamod) - myfun(x) )/10^4

# define (same) new data for prediction
newdata <- data.frame(x=seq(0,1,length.out=100))

# get fitted values and standard errors for new data
ya <- predict(cubamod,newdata,se.fit=TRUE)

# plot results with 95\% Bayesian confidence interval
plot(newdata$x,ya$fit,type="l")
lines(newdata$x,ya$fit+qnorm(.975)*ya$se.fit,lty=3)
lines(newdata$x,ya$fit-qnorm(.975)*ya$se.fit,lty=3)

# predict constant, linear, and nonlinear effects
ya0 <- predict(cubamod,newdata,se.fit=TRUE,effect="0")
yal <- predict(cubamod,newdata,se.fit=TRUE,effect="lin")
yan <- predict(cubamod,newdata,se.fit=TRUE,effect="non")
crossprod( ya$fit - (ya0$fit + yal$fit + yan$fit) )

# plot results with 95\% Bayesian confidence intervals
par(mfrow=c(1,2))
plot(newdata$x,yal$fit,type="l",main="Linear effect")
lines(newdata$x,yal$fit+qnorm(.975)*yal$se.fit,lty=3)
lines(newdata$x,yal$fit-qnorm(.975)*yal$se.fit,lty=3)
plot(newdata$x,yan$fit,type="l",main="Nonlinear effect")
lines(newdata$x,yan$fit+qnorm(.975)*yan$se.fit,lty=3)
lines(newdata$x,yan$fit-qnorm(.975)*yan$se.fit,lty=3)

}
