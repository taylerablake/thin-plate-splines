\name{makessa}
\alias{makessa}
\title{
Makes Objects to Fit Smoothing Spline ANOVA Models
}
\description{
This function creates a list containing the necessary information to fit a smoothing spline anova model (see \code{\link{bigssa}}). 
}
\usage{
makessa(formula,data=NULL,type=NULL,nknots=NULL,rparm=NA,
        lambdas=NULL,skip.iter=TRUE,se.fit=FALSE,rseed=1234,
        gcvopts=NULL,knotcheck=TRUE,gammas=NULL,weights=NULL,
        random=NULL,remlalg=c("FS","NR","EM","none"),remliter=500,
       remltol=10^-4,remltau=NULL)
}
\arguments{
  \item{formula}{
  An object of class "\code{formula}": a symbolic description of the model to be fitted (see Details and Examples for more information).
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
  \item{se.fit}{
	Logical indicating if the standard errors of the fitted values should be estimated.
}
  \item{rseed}{
  Random seed for knot sampling. Input is ignored if \code{nknots} is an input vector of knot indices. Set \code{rseed=NULL} to obtain a different knot sample each time, or set \code{rseed} to any positive integer to use a different seed than the default.
}
  \item{gcvopts}{
	Control parameters for optimization. List with 3 elements: (a) \code{maxit}: maximum number of algorithm iterations, (b) \code{gcvtol}: covergence tolerance for iterative GCV update, and (c) \code{alpha}: tuning parameter for GCV minimization. Default: \code{gcvopts=list(maxit=5,gcvtol=10^-5,alpha=1)}
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
  \item{random}{
Adds random effects to model (see Random Effects section).
}
  \item{remlalg}{
REML algorithm for estimating variance components (see Random Effects section). Input is ignored if \code{is.null(random)}.  
}
  \item{remliter}{
Maximum number of iterations for REML estimation of variance components. Input is ignored if \code{random=NULL}.  
}
  \item{remltol}{
Convergence tolerance for REML estimation of variance components. Input is ignored if \code{random=NULL}.  
}
  \item{remltau}{
Initial estimate of variance parameters for REML estimation of variance components. Input is ignored if \code{random=NULL}.  
}
}
\details{
 See \code{\link{bigssa}} and below example for more details.
}
\value{
	An object of class "makessa", which can be input to \code{\link{bigssa}}.
}
\references{
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.
  
  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.

  Helwig, N. E. (2016). Efficient estimation of variance components in nonparametric mixed-effects models with large samples. \emph{Statistics and Computing, 26}, 1319-1336.

	Helwig, N. E. and Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics, 24}, 715-732.
  
  Helwig, N. E. and Ma, P. (2016). Smoothing spline ANOVA for super-large samples: Scalable computation via rounding parameters. \emph{Statistics and Its Interface, 9}, 433-444.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\section{Warning }{
  When inputting a "makessa" class object into \code{\link{bigssa}}, the formula input to \code{bigssa} must be a nested version of the original formula input to \code{makessa}. In other words, you cannot add any new effects after a "makessa" object has been created, but you can drop (remove) effects from the model.
}


\examples{

##########   EXAMPLE  ##########

# function with two continuous predictors
set.seed(773)
myfun <- function(x1v,x2v){
  sin(2*pi*x1v) + log(x2v+.1) + cos(pi*(x1v-x2v))
}
x1v <- runif(500)
x2v <- runif(500)
y <- myfun(x1v,x2v) + rnorm(500)

# fit 2 possible models (create information 2 separate times)
system.time({
  intmod <- bigssa(y~x1v*x2v,type=list(x1v="cub",x2v="cub"),nknots=50)
  addmod <- bigssa(y~x1v+x2v,type=list(x1v="cub",x2v="cub"),nknots=50)
})

# fit 2 possible models (create information 1 time)
system.time({
  makemod <- makessa(y~x1v*x2v,type=list(x1v="cub",x2v="cub"),nknots=50)
  int2mod <- bigssa(y~x1v*x2v,makemod)
  add2mod <- bigssa(y~x1v+x2v,makemod)
})

# check difference (no difference)
crossprod( intmod$fitted.values - int2mod$fitted.values )
crossprod( addmod$fitted.values - add2mod$fitted.values )

}

