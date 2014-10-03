\name{predict.ssp}
\alias{predict.ssp}
\title{
Predicts for "ssp" Objects
}
\description{
Get fitted values and standard error estimates for smoothing splines with parametric effects.
}
\usage{
\method{predict}{ssp}(object,newdata=NULL,se.fit=FALSE,include=object$tnames,
        effect=c("all","0","lin","non"),includeint=FALSE,...)
}
\arguments{
  \item{object}{
	Object of class "ssp", which is output from \code{\link{bigssp}}.
}
  \item{newdata}{
	Data frame or list containing the new data points for prediction. Variable names must match those used in the \code{formula} input of \code{\link{bigssp}}. See Details and Example. Default of \code{newdata=NULL} uses original data in \code{object} input.
}
  \item{se.fit}{
	Logical indicating whether the standard errors of the fitted values should be estimated. Default is \code{se.fit=FALSE}.
}
  \item{include}{
	Which terms to include in the estimate. You can get fitted values for any combination of terms in the \code{tnames} element of an "ssp" object.
}
  \item{effect}{
  Which effect to estimate: \code{effect="all"} gives \eqn{\hat{y}} for given terms in \code{include}, \code{effect="lin"} gives linear portion of \eqn{\hat{y}} for given terms in \code{include}, and \code{effect="non"} gives nonlinear portion of \eqn{\hat{y}} for given terms in \code{include}. Use \code{effect="0"} to return the intercept.
}
  \item{includeint}{
  Logical indicating whether the intercept should be included in the prediction. If \code{include=object$tnames} and \code{effect="all"} (default), then this input is ignored and the intercept is automatically included in the prediction.
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

	Helwig, N. E. and Ma, P. (in press). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
  
  Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
  
  Helwig, N. E. and Ma, P. (in preparation). Optimal eigenbasis truncation for thin-plate regression splines. 
	
	Helwig, N. E. and Ma, P. (in preparation). Smoothing spline ANOVA for super-large samples: Scalable computation via rounding parameters.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
	Uses the coefficient and smoothing parameter estimates from a fit smoothing spline with parametric effects (estimated by \code{\link{bigssp}}) to predict for new data.
}

\examples{

##########   EXAMPLE 1   ##########

# define univariate function and data
set.seed(773)
myfun<-function(x){2+x+sin(2*pi*x)}
x<-runif(500)
y<-myfun(x)+rnorm(500)

# fit cubic spline model
cubmod<-bigssp(y~x,type="cub",nknots=30)
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
sum(yc$fit-(yc0$fit+ycl$fit+ycn$fit))

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
x<-cbind(runif(500),runif(500))*16-8
y<-myfun(x)+rnorm(500)

# bidimensional thin-plate spline with 50 knots
tpsmod<-bigssp(y~x,type="tps",nknots=50)
crossprod(predict(tpsmod)-myfun(x))/500

# define new data for prediction
xnew<-as.matrix(expand.grid(seq(-8,8,l=50),seq(-8,8,l=50)))
newdata<-list(x=xnew)

# get fitted values for new data
yp<-predict(tpsmod,newdata)

# plot results
imagebar(seq(-8,8,l=50),seq(-8,8,l=50),matrix(yp,50,50),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))))

# predict linear and nonlinear effects
yl<-predict(tpsmod,newdata,effect="lin")
yn<-predict(tpsmod,newdata,effect="non")

# plot results
par(mfrow=c(1,2))
imagebar(seq(-8,8,l=50),seq(-8,8,l=50),matrix(yl,50,50),
         main="Linear effect",xlab=expression(italic(x)[1]),
         ylab=expression(italic(x)[2]),zlab=expression(hat(italic(y))))
imagebar(seq(-8,8,l=50),seq(-8,8,l=50),matrix(yn,50,50),
         main="Nonlinear effect",xlab=expression(italic(x)[1]),
         ylab=expression(italic(x)[2]),zlab=expression(hat(italic(y))))
         

}
