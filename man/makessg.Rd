\name{makessg}
\alias{makessg}
\title{
Makes Objects to Fit Generalized Smoothing Spline ANOVA Models
}
\description{
This function creates a list containing the necessary information to fit a generalized smoothing spline anova model (see \code{\link{bigssg}}). 
}
\usage{
makessg(formula,family,data,type=NULL,nknots=NULL,rparm=NA,
        lambdas=NULL,skip.iter=TRUE,se.lp=FALSE,rseed=1234,
        gcvopts=NULL,knotcheck=TRUE,gammas=NULL,weights=NULL,
        gcvtype=c("acv","gacv","gacv.old"))
}
\arguments{
  \item{formula}{
  An object of class "\code{formula}": a symbolic description of the model to be fitted (see Details and Examples for more information).
}
  \item{family}{
  Distribution for response. One of five options: \code{"binomial"}, \code{"poisson"}, \code{"Gamma"}, \code{"inverse.gaussian"}, or \code{"negbin"}. See \code{\link{bigssg}}.
}
  \item{data}{
  Optional data frame, list, or environment containing the variables in \code{formula}.
}
  \item{type}{
	List of smoothing spline types for predictors in \code{formula} (see Details). Options include \code{type="cub"} for cubic, \code{type="acub"} for another cubic, \code{type="per"} for cubic periodic, \code{type="tps"} for cubic thin-plate, and \code{type="nom"} for nominal.
}
  \item{nknots}{
	Two possible options: (a) scalar giving total number of random knots to sample, or (b) vector indexing which rows of \code{data} to use as knots.
}
  \item{rparm}{
	List of rounding parameters for each predictor. See Details. 
}
  \item{lambdas}{
	Vector of global smoothing parameters to try. Default uses \code{lambdas=10^-c(9:0)}
}
  \item{skip.iter}{
	Logical indicating whether to skip the iterative smoothing parameter update. Using \code{skip.iter=FALSE} should provide a more optimal solution, but the fitting time may be substantially longer. See Computational Details.
}
  \item{se.lp}{
  Logical indicating if the standard errors of the linear predictors (\eqn{\eta}) should be estimated.
}
  \item{rseed}{
  Random seed for knot sampling. Input is ignored if \code{nknots} is an input vector of knot indices. Set \code{rseed=NULL} to obtain a different knot sample each time, or set \code{rseed} to any positive integer to use a different seed than the default.
}
  \item{gcvopts}{
  Control parameters for optimization. List with 6 elements: (i) \code{maxit}: maximum number of outer iterations, (ii) \code{gcvtol}: covergence tolerance for iterative GACV update, (iii) \code{alpha}: tuning parameter for GACV minimization, (iv) \code{inmaxit}: maximum number of inner iterations for iteratively reweighted fitting, (v) \code{intol}: inner convergence tolerance for iteratively reweighted fitting, and (vi) \code{insub}: number of data points to subsample when checking inner convergence. \code{gcvopts=list(maxit=5,gcvtol=10^-5,alpha=1,inmaxit=100,intol=10^-5,insub=10^4)}
}
  \item{knotcheck}{
If \code{TRUE}, only unique knots are used (for stability).  
}
  \item{gammas}{
List of initial smoothing parameters for each predictor. See Details. 
}
  \item{weights}{
Vector of positive weights for fitting (default is vector of ones).  
}
  \item{gcvtype}{
Cross-validation criterion for selecting smoothing parameters (see Details).
}
}
\details{
 See \code{\link{bigssg}} and below example for more details.
}
\value{
	An object of class "makessg", which can be input to \code{\link{bigssg}}.
}
\references{
  Gu, C. and Xiang, D. (2001). Cross-validating non-Gaussian data: Generalized approximate cross-validation revisited. \emph{Journal of Computational and Graphical Statistics, 10}, 581-591.
  
  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.

  Helwig, N. E. (in preparation). Nonparametric exponential family regression for ultra large samples: Scalable computation via rounding parameters. 

	Helwig, N. E. and Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
  
  Helwig, N. E. and Ma, P. (in press). Smoothing spline ANOVA for super-large samples: Scalable computation via rounding parameters. \emph{Statistics and Its Interface}.
	
  Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
  
  Helwig, N. E. and Ma, P. (in preparation). Optimal eigenbasis truncation for thin-plate regression splines. 
  
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\section{Warning }{
  When inputting a "makessg" class object into \code{\link{bigssg}}, the formula input to \code{bigssg} must be a nested version of the original formula input to \code{makessg}. In other words, you cannot add any new effects after a "makessg" object has been created, but you can drop (remove) effects from the model.
}


\examples{

##########   EXAMPLE  ##########

# function with two continuous predictors
set.seed(1)
myfun <- function(x1v,x2v){
  sin(2*pi*x1v) + log(x2v+.1) + cos(pi*(x1v-x2v))
}
ndpts <- 1000
x1v <- runif(ndpts)
x2v <- runif(ndpts)

# binomial response (no weights)
set.seed(773)
lp <- myfun(x1v,x2v)
p <- 1/(1+exp(-lp))
y <- rbinom(n=ndpts,size=1,p=p)

# fit 2 possible models (create information 2 separate times)
system.time({
  intmod <- bigssg(y~x1v*x2v,family="binomial",type=list(x1v="cub",x2v="cub"),nknots=50)
  addmod <- bigssg(y~x1v+x2v,family="binomial",type=list(x1v="cub",x2v="cub"),nknots=50)
})

# fit 2 possible models (create information 1 time)
system.time({
  makemod <- makessg(y~x1v*x2v,family="binomial",type=list(x1v="cub",x2v="cub"),nknots=50)
  int2mod <- bigssg(y~x1v*x2v,data=makemod)
  add2mod <- bigssg(y~x1v+x2v,data=makemod)
})

# check difference (no difference)
crossprod( intmod$fitted.values - int2mod$fitted.values )
crossprod( addmod$fitted.values - add2mod$fitted.values )

}

