\name{predict.bigtps}
\alias{predict.bigtps}
\title{
Predicts for "bigtps" Objects
}
\description{
Get fitted values and standard error estimates for thin-plate splines.
}
\usage{
\method{predict}{bigtps}(object,newdata=NULL,se.fit=FALSE,
        effect=c("all","0","lin","non"),
        design=FALSE,smoothMatrix=FALSE,...)
}
\arguments{
  \item{object}{
	Object of class "bigtps", which is output from \code{\link{bigtps}}.
}
  \item{newdata}{
	Vector or matrix containing new data points for prediction. See Details and Example. Default of \code{newdata=NULL} uses original data in \code{object} input.
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
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.

	Helwig, N. E. and Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics, 24}, 715-732.
  
  Helwig, N. E. and Ma, P. (2016). Smoothing spline ANOVA for super-large samples: Scalable computation via rounding parameters. \emph{Statistics and Its Interface, 9}, 433-444.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
	Uses the coefficient and smoothing parameter estimates from a fit thin-plate spline (estimated by \code{\link{bigtps}}) to predict for new data.
}

\examples{

##########   EXAMPLE 1   ##########

# define univariate function and data
set.seed(773)
myfun <- function(x){ 2 + x + sin(2*pi*x) }
x <- runif(10^4)
y <- myfun(x) + rnorm(10^4)

# fit thin-plate spline (default 1 dim: 30 knots)
tpsmod <- bigtps(x,y)
crossprod( predict(tpsmod) - myfun(x) )/10^4

# define new data for prediction
newdata <- data.frame(x=seq(0,1,length.out=100))

# get fitted values and standard errors for new data
yc <- predict(tpsmod,newdata,se.fit=TRUE)

# plot results with 95\% Bayesian confidence interval
plot(newdata$x,yc$fit,type="l")
lines(newdata$x,yc$fit+qnorm(.975)*yc$se.fit,lty=3)
lines(newdata$x,yc$fit-qnorm(.975)*yc$se.fit,lty=3)

# predict constant, linear, and nonlinear effects
yc0 <- predict(tpsmod,newdata,se.fit=TRUE,effect="0")
ycl <- predict(tpsmod,newdata,se.fit=TRUE,effect="lin")
ycn <- predict(tpsmod,newdata,se.fit=TRUE,effect="non")
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

# function with two continuous predictors
set.seed(773)
myfun <- function(x1v,x2v){
  sin(2*pi*x1v) + log(x2v+.1) + cos(pi*(x1v-x2v))
}
x <- cbind(runif(10^4),runif(10^4))
y <- myfun(x[,1],x[,2]) + rnorm(10^4)

# fit thin-plate spline (default 2 dim: 100 knots)
tpsmod <- bigtps(x,y)

# define new data
newdata <- as.matrix(expand.grid(seq(0,1,length=50),seq(0,1,length=50)))

# get fitted values for new data
yp <- predict(tpsmod,newdata)

# plot results
imagebar(seq(0,1,length=50),seq(0,1,length=50),matrix(yp,50,50),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))))

# predict linear and nonlinear effects
yl <- predict(tpsmod,newdata,effect="lin")
yn <- predict(tpsmod,newdata,effect="non")

# plot results
par(mfrow=c(1,2))
imagebar(seq(0,1,length=50),seq(0,1,length=50),matrix(yl,50,50),
         main="Linear effect",xlab=expression(italic(x)[1]),
         ylab=expression(italic(x)[2]),zlab=expression(hat(italic(y))))
imagebar(seq(0,1,length=50),seq(0,1,length=50),matrix(yn,50,50),
         main="Nonlinear effect",xlab=expression(italic(x)[1]),
         ylab=expression(italic(x)[2]),zlab=expression(hat(italic(y))))

}
