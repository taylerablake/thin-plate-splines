\name{predict.bigssg}
\alias{predict.bigssg}
\title{
Predicts for "bigssg" Objects
}
\description{
Get fitted values and standard error estimates for generalized smoothing spline anova models. 
}
\usage{
\method{predict}{bigssg}(object,newdata=NULL,se.lp=FALSE,include=object$tnames,
        effect=c("all","0","lin","non"),includeint=FALSE,
        design=FALSE,smoothMatrix=FALSE,intercept=NULL,...)
}
\arguments{
  \item{object}{
	Object of class "bigssg", which is output from \code{\link{bigssg}}.
}
  \item{newdata}{
	Data frame or list containing the new data points for prediction. Variable names must match those used in the \code{formula} input of \code{\link{bigssg}}. See Details and Example. Default of \code{newdata=NULL} uses original data in \code{object} input.
}
  \item{se.lp}{
	Logical indicating if the standard errors of the linear predictors (\eqn{\eta}) should be estimated. Default is \code{se.lp=FALSE}.
}
  \item{include}{
	Which terms to include in the estimate. You can get fitted values for any combination of terms in the \code{tnames} element of an "bigssg" object.
}
  \item{effect}{
  Which effect to estimate: \code{effect="all"} gives \eqn{\hat{y}} for given terms in \code{include}, \code{effect="lin"} gives linear portion of \eqn{\hat{y}} for given terms in \code{include}, and \code{effect="non"} gives nonlinear portion of \eqn{\hat{y}} for given terms in \code{include}. Use \code{effect="0"} to return the intercept.
}
  \item{includeint}{
  Logical indicating whether the intercept should be included in the prediction. If \code{include=object$tnames} and \code{effect="all"} (default), then this input is ignored and the intercept is automatically included in the prediction.
}
  \item{design}{
  Logical indicating whether the design matrix should be returned.
}
  \item{smoothMatrix}{
  Logical indicating whether the smoothing matrix should be returned.
}
  \item{intercept}{
  Logical indicating whether the intercept should be included in the prediction. When used, this input overrides the \code{includeint} input.
}
  \item{\dots}{
	Ignored.
}
}
\value{
  Returns list with elements:
  \item{fitted.values}{Vector of fitted values (on data scale)}
  \item{linear.predictors}{Vector of fitted values (on link scale)}
  \item{se.lp}{Vector of standard errors of linear predictors (if \code{se.lp=TRUE})}
  \item{X}{Design matrix used to create linear predictors (if \code{design=TRUE})}
  \item{ix}{Index vector such that \code{linear.predictors=X\%*\%object$modelspec$coef[ix]} (if \code{design=TRUE})}
  \item{S}{Smoothing matrix corresponding to fitted values (if \code{smoothMatrix=TRUE})}
}
\references{
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.

  Gu, C. and Xiang, D. (2001). Cross-validating non-Gaussian data: Generalized approximate cross-validation revisited. \emph{Journal of Computational and Graphical Statistics, 10}, 581-591.

  Helwig, N. E. (in preparation). Nonparametric exponential family regression for ultra large samples: Scalable computation via rounding parameters. 

	Helwig, N. E. and Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics, 24}, 715-732.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
	Uses the coefficient and smoothing parameter estimates from a fit generalized smoothing spline anova (estimated by \code{\link{bigssg}}) to predict for new data.
}

\examples{

##########   EXAMPLE 1   ##########

# define univariate function and data
set.seed(1)
myfun <- function(x){ sin(2*pi*x) }
ndpts <- 1000
x <- runif(ndpts)

# negative binomial response (unknown dispersion)
set.seed(773)
lp <- myfun(x)
mu <- exp(lp)
y <- rnbinom(n=ndpts,size=2,mu=mu)

# fit cubic spline model
cubmod <- bigssg(y~x,family="negbin",type="cub",nknots=20)
1/cubmod$dispersion   ## dispersion = 1/size
crossprod( lp - cubmod$linear.predictor )/length(lp)

# define new data for prediction
newdata <- data.frame(x=seq(0,1,length.out=100))

# get fitted values and standard errors for new data
yc <- predict(cubmod,newdata,se.lp=TRUE)

# plot results with 95\% Bayesian confidence interval (link scale)
plot(newdata$x,yc$linear.predictor,type="l")
lines(newdata$x,yc$linear.predictor+qnorm(.975)*yc$se.lp,lty=3)
lines(newdata$x,yc$linear.predictor-qnorm(.975)*yc$se.lp,lty=3)

# plot results with 95\% Bayesian confidence interval (data scale)
plot(newdata$x,yc$fitted,type="l")
lines(newdata$x,exp(yc$linear.predictor+qnorm(.975)*yc$se.lp),lty=3)
lines(newdata$x,exp(yc$linear.predictor-qnorm(.975)*yc$se.lp),lty=3)

# predict constant, linear, and nonlinear effects
yc0 <- predict(cubmod,newdata,se.lp=TRUE,effect="0")
ycl <- predict(cubmod,newdata,se.lp=TRUE,effect="lin")
ycn <- predict(cubmod,newdata,se.lp=TRUE,effect="non")
crossprod( yc$linear - (yc0$linear + ycl$linear + ycn$linear) )

# plot results with 95\% Bayesian confidence intervals (link scale)
par(mfrow=c(1,2))
plot(newdata$x,ycl$linear,type="l",main="Linear effect")
lines(newdata$x,ycl$linear+qnorm(.975)*ycl$se.lp,lty=3)
lines(newdata$x,ycl$linear-qnorm(.975)*ycl$se.lp,lty=3)
plot(newdata$x,ycn$linear,type="l",main="Nonlinear effect")
lines(newdata$x,ycn$linear+qnorm(.975)*ycn$se.lp,lty=3)
lines(newdata$x,ycn$linear-qnorm(.975)*ycn$se.lp,lty=3)
         
# plot results with 95\% Bayesian confidence intervals (data scale)
par(mfrow=c(1,2))
plot(newdata$x,ycl$fitted,type="l",main="Linear effect")
lines(newdata$x,exp(ycl$linear+qnorm(.975)*ycl$se.lp),lty=3)
lines(newdata$x,exp(ycl$linear-qnorm(.975)*ycl$se.lp),lty=3)
plot(newdata$x,ycn$fitted,type="l",main="Nonlinear effect")
lines(newdata$x,exp(ycn$linear+qnorm(.975)*ycn$se.lp),lty=3)
lines(newdata$x,exp(ycn$linear-qnorm(.975)*ycn$se.lp),lty=3)

         
         
##########   EXAMPLE 2   ##########

# define bivariate function and data
set.seed(1)
myfun <- function(x1v,x2v){
  sin(2*pi*x1v) + log(x2v+.1) + cos(pi*(x1v-x2v))
}
ndpts <- 1000
x1v <- runif(ndpts)
x2v <- runif(ndpts)

# binomial response (with weights)
set.seed(773)
lp <- myfun(x1v,x2v)
p <- 1/(1+exp(-lp))
w <- sample(c(10,20,30,40,50),length(p),replace=TRUE)
y <- rbinom(n=ndpts,size=w,p=p)/w   ## y is proportion correct
cubmod <- bigssg(y~x1v*x2v,family="binomial",type=list(x1v="cub",x2v="cub"),nknots=100,weights=w)
crossprod( lp - cubmod$linear.predictor )/length(lp)

# define new data for prediction
xnew <- as.matrix(expand.grid(seq(0,1,length=50),seq(0,1,length=50)))
newdata <- list(x1v=xnew[,1],x2v=xnew[,2])

# get fitted values for new data
yp <- predict(cubmod,newdata)

# plot linear predictor and fitted values
par(mfrow=c(2,2))
imagebar(seq(0,1,l=50),seq(0,1,l=50),matrix(myfun(newdata$x1v,newdata$x2v),50,50),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))),zlim=c(-4.5,1.5),main="True Linear Predictor")
imagebar(seq(0,1,l=50),seq(0,1,l=50),matrix(yp$linear.predictor,50,50),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))),zlim=c(-4.5,1.5),main="Estimated Linear Predictor")
newprob <- 1/(1+exp(-myfun(newdata$x1v,newdata$x2v)))
imagebar(seq(0,1,l=50),seq(0,1,l=50),matrix(newprob,50,50),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))),zlim=c(0,0.8),main="True Probabilities")
imagebar(seq(0,1,l=50),seq(0,1,l=50),matrix(yp$fitted.values,50,50),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))),zlim=c(0,0.8),main="Estimated Probabilities")         

# predict linear and nonlinear effects for x1v (link scale)
newdata <- list(x1v=seq(0,1,length.out=100))
yl <- predict(cubmod,newdata,include="x1v",effect="lin",se.lp=TRUE)
yn <- predict(cubmod,newdata,include="x1v",effect="non",se.lp=TRUE)

# plot results with 95\% Bayesian confidence intervals (link scale)
par(mfrow=c(1,2))
plot(newdata$x1v,yl$linear,type="l",main="Linear effect")
lines(newdata$x1v,yl$linear+qnorm(.975)*yl$se.lp,lty=3)
lines(newdata$x1v,yl$linear-qnorm(.975)*yl$se.lp,lty=3)
plot(newdata$x1v,yn$linear,type="l",main="Nonlinear effect")
lines(newdata$x1v,yn$linear+qnorm(.975)*yn$se.lp,lty=3)
lines(newdata$x1v,yn$linear-qnorm(.975)*yn$se.lp,lty=3)

}
