\name{bigtps}
\alias{bigtps}
\title{
Fits cubic thin-plate splines
}
\description{
Given a real-valued response vector \eqn{\mathbf{y}=\{y_{i}\}_{n\times1}}, a thin-plate spline model has the form \deqn{y_{i}=\eta(\mathbf{x}_{i})+e_{i}} where \eqn{y_{i}} is the \eqn{i}-th observation's respone, \eqn{\mathbf{x}_{i}=(x_{i1},\ldots,x_{id})} is the \eqn{i}-th observation's nonparametric predictor vector, \eqn{\eta} is an unknown smooth function relating the response and predictor, and \eqn{e_{i}\sim\mathrm{N}(0,\sigma^{2})} is iid Gaussian error. Function only fits interaction models.
}
\usage{
bigtps(x,y,nknots=NULL,nvec=NULL,rparm=NA,
       alpha=1,lambdas=NULL,se.fit=FALSE,
       rseed=1234,knotcheck=TRUE)
}
\arguments{
  \item{x}{
	Predictor vector or matrix with three or less columns.
}
  \item{y}{
  Response vector. Must be same length as \code{x} has rows.
}
  \item{nknots}{
	Two possible options: (a) scalar giving total number of random knots to sample, or (b) vector indexing which rows of \code{x} to use as knots.
} 
  \item{nvec}{
  Number of eigenvectors (and eigenvalues) to use in approximation. Must be less than or equal to the number of knots and greater than or equal to \code{ncol(x)+2}. Default sets \code{nvec<-nknots}. Can also input \code{0<nvec<1} to retain \code{nvec} percentage of eigenbasis variation.
}
  \item{rparm}{
	Rounding parameter(s) for \code{x}. Use \code{rparm=NA} to fit unrounded solution. Can provide one (positive) rounding parameter for each column of \code{x}.
}
  \item{alpha}{
  Manual tuning parameter for GCV score. Using \code{alpha=1} gives unbaised esitmate. Using a larger alpha enforces a smoother estimate.
}
  \item{lambdas}{
	Vector of global smoothing parameters to try. Default estimates smoothing parameter that minimizes GCV score.
}
  \item{se.fit}{
	Logical indicating if the standard errors of fitted values should be estimated. 
}
  \item{rseed}{
  Random seed for knot sampling. Input is ignored if \code{nknots} is an input vector of knot indices. Set \code{rseed=NULL} to obtain a different knot sample each time, or set \code{rseed} to any positive integer to use a different seed than the default.
}
  \item{knotcheck}{
If \code{TRUE}, only unique knots are used (for stability).  
}
}
\details{
  To estimate \eqn{\eta} I minimize the penalized least-squares functional \deqn{\frac{1}{n}\sum_{i=1}^{n}(y_{i}-\eta(\mathbf{x}_{i}))^{2}+\lambda J(\eta)} where \eqn{J(\eta)} is the thin-plate penalty (see Helwig and Ma) and \eqn{\lambda\geq0} is a smoothing parameter that controls the trade-off between fitting and smoothing the data. Default use of the function estimates \eqn{\lambda} by minimizing the GCV score (see \code{\link{bigspline}}).
  
  Using the rounding parameter input \code{rparm} can greatly speed-up and stabilize the fitting for large samples. When \code{rparm} is used, the spline is fit to a set of unique data points after rounding; the unique points are determined using the efficient algorithm described in Helwig (2013). Rounding parameter should be on the raw data scale. 
}
\value{
	\item{fitted.values}{Vector of fitted values corresponding to the original data points in \code{x} (if \code{rparm=NA}) or the rounded data points in \code{xunique} (if \code{rparm} is used).}
	\item{se.fit}{Vector of standard errors of \code{fitted.values} (if input \code{se.fit=TRUE)}.}
  \item{x}{Predictor vector (same as input).}
	\item{y}{Response vector (same as input).}
  \item{xunique}{Unique elements of \code{x} after rounding (if \code{rparm} is used).}
  \item{yunique}{Mean of \code{y} for unique elements of \code{x} after rounding (if \code{rparm} is used).}
  \item{funique}{Vector giving frequency of each element of \code{xunique} (if \code{rparm} is used).}
	\item{sigma}{Estimated error standard deviation, i.e., \eqn{\hat{\sigma}}.}
  \item{ndf}{Data frame with two elements: \code{n} is total sample size, and \code{df} is effective degrees of freedom of fit model (trace of smoothing matrix).}
  \item{info}{Model fit information: vector containing the GCV, multiple R-squared, AIC, and BIC of fit model (assuming Gaussian error).}
  \item{myknots}{Spline knots used for fit.}
  \item{nvec}{Number of eigenvectors used for solution.}
  \item{rparm}{Rounding parameter for \code{x} (same as input).}
  \item{lambda}{Optimal smoothing parameter.}
  \item{coef}{Spline basis function coefficients.}
  \item{coef.csqrt}{Matrix square-root of covariace matrix of \code{coef}. Use \code{tcrossprod(coef.csqrt)} to get covariance matrix of \code{coef}.}
}
\references{
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.
	
	Helwig, N. E. and Ma, P. (in press). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
  
  Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
  
  Helwig, N. E. and Ma, P. (in preparation). Optimal eigenbasis truncation for thin-plate regression splines. 
	
	Helwig, N. E. and Ma, P. (in preparation). Nonparametric regression for ultra large samples: Scalable computation via rounding parameters. 
}
\author{
Nathaniel E. Helwig <nhelwig2@illinois.edu>
}
\note{
	The spline is estimated using penalized least-squares, which does not require the Gaussian error assumption. However, the spline inference information (e.g., standard errors and fit information) requires the Gaussian error assumption.
}
\section{Warnings }{
  Input \code{nvec} must be greater than \code{ncol(x)+1}.

  When using rounding parameters, output \code{fitted.values} corresponds to unique rounded predictor scores in output \code{xunique}. Use \code{\link{predict.tps}} function to get fitted values for full \code{y} vector.
}

\section{Computational Details }{
  According to thin-plate spline theory, the function \eqn{\eta} can be approximated as \deqn{\eta(x) = \sum_{k=1}^{M}d_{k}\phi_{k}(\mathbf{x}) + \sum_{h=1}^{q}c_{h}\xi(\mathbf{x},\mathbf{x}_{h}^{*})} where the \eqn{\{\phi_{k}\}_{k=1}^{M}} are linear functions, \eqn{\xi} is the thin-plate spline semi-kernel, \eqn{\{\mathbf{x}_{h}^{*}\}_{h=1}^{q}} are the knots, and the \eqn{c_{h}} coefficients are constrained to be orthongonal to the \eqn{\{\phi_{k}\}_{k=1}^{M}} functions.
  
  This implies that the penalized least-squares functional can be rewritten as \deqn{ \|\mathbf{y} - \mathbf{K}\mathbf{d} - \mathbf{J}\mathbf{c}\|^{2} + n\lambda\mathbf{c}'\mathbf{Q}\mathbf{c} }
where \eqn{\mathbf{K}=\{\phi(\mathbf{x}_{i})\}_{n \times M}} is the null space basis function matrix, \eqn{\mathbf{J}=\{\xi(\mathbf{x}_{i},\mathbf{x}_{h}^{*})\}_{n \times q}} is the contrast space basis funciton matrix, \eqn{\mathbf{Q}=\{\xi(\mathbf{x}_{g}^{*},\mathbf{x}_{h}^{*})\}_{q \times q}} is the penalty matrix, and \eqn{\mathbf{d}=(d_{0},\ldots,d_{M})'} and \eqn{\mathbf{c}=(c_{1},\ldots,c_{q})'} are the unknown basis function coefficients, where \eqn{\mathbf{c}} are constrained to be orthongonal to the \eqn{\{\phi_{k}\}_{k=1}^{M}} functions.

See Helwig and Ma for specifics about how the constrained estimation is handled.

}



\examples{

##########   EXAMPLE 1   ##########

# define relatively smooth function
set.seed(773)
myfun<-function(x){sin(2*pi*x)}
x<-runif(500)
y<-myfun(x)+rnorm(500)

# fit thin-plate spline (default 1 dim: 30 knots)
tpsmod<-bigtps(x,y)
tpsmod


##########   EXAMPLE 2   ##########

# define more jagged function
set.seed(773)
myfun<-function(x){2*x+cos(2*pi*x)}
x<-runif(500)*4
y<-myfun(x)+rnorm(500)

# try different numbers of knots
r1mod<-bigtps(x,y,nknots=20,rparm=0.01)
crossprod(myfun(r1mod$xunique)-r1mod$fitted)/length(r1mod$fitted)
r2mod<-bigtps(x,y,nknots=35,rparm=0.01)
crossprod(myfun(r2mod$xunique)-r2mod$fitted)/length(r2mod$fitted)
r3mod<-bigtps(x,y,nknots=50,rparm=0.01)
crossprod(myfun(r3mod$xunique)-r3mod$fitted)/length(r3mod$fitted)


##########   EXAMPLE 3   ##########

# function with two continuous predictors
set.seed(773)
myfun<-function(x1v,x2v){sin(2*pi*x1v)+log(x2v+.1)+cos(pi*(x1v-x2v))}
x<-cbind(runif(500),runif(500))
y<-myfun(x[,1],x[,2])+rnorm(500)

# fit thin-plate spline with 50 knots (default 2 dim: 100 knots)
tpsmod<-bigtps(x,y,nknots=50)
tpsmod
crossprod(myfun(x[,1],x[,2])-tpsmod$fitted.values)/500


##########   EXAMPLE 4   ##########

# function with three continuous predictors
set.seed(773)
myfun<-function(x1v,x2v,x3v){
  sin(2*pi*x1v)+log(x2v+.1)+cos(pi*x3v)
  }
x=cbind(runif(500),runif(500),runif(500))
y<-myfun(x[,1],x[,2],x[,3])+rnorm(500)

# fit thin-plate spline with 50 knots (default 3 dim: 200 knots)
tpsmod<-bigtps(x,y,nknots=50)
tpsmod
crossprod(myfun(x[,1],x[,2],x[,3])-tpsmod$fitted.values)/500


}

