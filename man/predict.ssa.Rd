\name{predict.ssa}
\alias{predict.ssa}
\title{
Predicts for "ssa" objects
}
\description{
Get fitted values and standard error estimates for smoothing spline anova models. 
}
\usage{
\method{predict}{ssa}(object,newdata=NULL,se.fit=FALSE,
        include=c(object$tnames),
        effect=c("all","0","lin","non"),...)
}
\arguments{
  \item{object}{
	Object of class "ssa", which is output from \code{\link{bigssa}}.
}
  \item{newdata}{
	Data frame or list containing the new data points for prediction. Variable names must match those used in the \code{formula} input of \code{\link{bigssa}}. See Details and Example. Default of \code{newdata=NULL} uses original data in \code{object} input.
}
  \item{se.fit}{
	Logical indicating whether the standard errors of the fitted values should be estimated. Default is \code{se.fit=FALSE}.
}
  \item{include}{
	Which terms to include in the estimate. You can get fitted values for any combination of terms in the \code{tnames} element of an "ssa" object.
}
\item{effect}{
  Which effect to estimate: \code{effect="all"} gives \eqn{\hat{y}} for given terms in \code{include}, \code{effect="lin"} gives linear portion of \eqn{\hat{y}} for given terms in \code{include}, and \code{effect="non"} gives nonlinear portion of \eqn{\hat{y}} for given terms in \code{include}. Use \code{effect="0"} to return the intercept.
}
  \item{\dots}{
	Ignored.
}
}
\value{
  For \code{se.fit=FALSE}, returns vector of fitted values.

  For \code{se.fit=TRUE}, returns list with two components:
	\item{fit}{
  Vector of fitted values
}
\item{se.fit}{
  Vector of standard errors of fitted values
}
}
\references{
  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.

	Helwig, N. E. and Ma, P. (under review). Fast thin-plate splines via bin-sampled knots. \emph{Computational Statistics & Data Analysis}.
	
	Helwig, N. E. and Ma, P. (under review). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
	
	Helwig, N. E. and Ma, P. (in preparation). Nonparametric regression for ultra large samples: Scalable computation via rounding parameters.
}
\author{
Nathaniel E. Helwig <nhelwig2@illinois.edu>
}
\details{
	Uses the coefficient and smoothing parameter estimates from a fit smoothing spline anova (estimated by \code{\link{bigssa}}) to predict for new data.
}

\examples{

##########   EXAMPLE 1   ##########

# define univariate function and data
set.seed(773)
myfun<-function(x){2+x+sin(2*pi*x)}
x<-runif(500)
y<-myfun(x)+rnorm(500)

# fit cubic spline model
cubmod<-bigssa(y~x,type="cub",nknots=30)
crossprod(predict(cubmod)-myfun(x))/500

# define new data for prediction
newdata<-data.frame(x=seq(0,1,length.out=100))

# get fitted values and standard errors for new data
yc<-predict(cubmod,newdata,se.fit=TRUE)

# plot results with 95\% Bayesian confidence interval
plot(newdata$x,yc$fit,type="l")
lines(newdata$x,yc$fit+qnorm(.975)*yc$se.fit,lty=3)
lines(newdata$x,yc$fit-qnorm(.975)*yc$se.fit,lty=3)

# predict constant, linear, and nonlinear effects
yc0<-predict(cubmod,newdata,se.fit=TRUE,effect="0")
ycl<-predict(cubmod,newdata,se.fit=TRUE,effect="lin")
ycn<-predict(cubmod,newdata,se.fit=TRUE,effect="non")
crossprod(yc$fit-(yc0$fit+ycl$fit+ycn$fit))

# plot results with 95\% Bayesian confidence intervals
par(mfrow=c(1,2))
plot(newdata$x,ycl$fit,type="l",main="Linear effect")
lines(newdata$x,ycl$fit+qnorm(.975)*ycl$se.fit,lty=3)
lines(newdata$x,ycl$fit-qnorm(.975)*ycl$se.fit,lty=3)
plot(newdata$x,ycn$fit,type="l",main="Nonlinear effect")
lines(newdata$x,ycn$fit+qnorm(.975)*ycn$se.fit,lty=3)
lines(newdata$x,ycn$fit-qnorm(.975)*ycn$se.fit,lty=3)
         
         
##########   EXAMPLE 2   ##########

# define bivariate function and data
set.seed(773)
myfun<-function(x){2+x[,1]/10-x[,2]/5+2*sin(sqrt(x[,1]^2+x[,2]^2+.1))/sqrt(x[,1]^2+x[,2]^2+.1)}
x1v<-runif(500)*16-8;   x2v<-runif(500)*16-8
y<-myfun(cbind(x1v,x2v))+rnorm(500)

# tensor product cubic splines with 50 knots
cubmod<-bigssa(y~x1v*x2v,type=list(x1v="cub",x2v="cub"),nknots=50)
crossprod(predict(cubmod)-myfun(cbind(x1v,x2v)))/500

# define new data for prediction
xnew<-as.matrix(expand.grid(seq(-8,8,l=50),seq(-8,8,l=50)))
newdata<-list(x1v=xnew[,1],x2v=xnew[,2])

# get fitted values for new data
yp<-predict(cubmod,newdata)

# plot results
imagebar(seq(-8,8,l=50),seq(-8,8,l=50),matrix(yp,50,50),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))))

# predict linear and nonlinear effects for x1v
newdata<-list(x1v=seq(-8,8,length.out=100))
yl<-predict(cubmod,newdata,include="x1v",effect="lin",se.fit=TRUE)
yn<-predict(cubmod,newdata,include="x1v",effect="non",se.fit=TRUE)

# plot results with 95\% Bayesian confidence intervals
par(mfrow=c(1,2))
plot(newdata$x1v,yl$fit,type="l",main="Linear effect")
lines(newdata$x1v,yl$fit+qnorm(.975)*yl$se.fit,lty=3)
lines(newdata$x1v,yl$fit-qnorm(.975)*yl$se.fit,lty=3)
plot(newdata$x1v,yn$fit,type="l",main="Nonlinear effect",ylim=c(-.3,.4))
lines(newdata$x1v,yn$fit+qnorm(.975)*yn$se.fit,lty=3)
lines(newdata$x1v,yn$fit-qnorm(.975)*yn$se.fit,lty=3)

}
