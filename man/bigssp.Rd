\name{bigssp}
\alias{bigssp}
\title{
Fits smoothing splines with parametric effects
}
\description{
Given a real-valued response vector \eqn{\mathbf{y}=\{y_{i}\}_{n\times1}}, a semiparametric regression model has the form \deqn{y_{i}= \eta(\mathbf{x}_{i}) + \sum_{j=1}^{t}b_{j}z_{ij} + e_{i}} where \eqn{y_{i}} is the \eqn{i}-th observation's respone, \eqn{\mathbf{x}_{i}=(x_{i1},\ldots,x_{ip})} is the \eqn{i}-th observation's nonparametric predictor vector, \eqn{\eta} is an unknown smooth function relating the response and nonparametric predictors, \eqn{\mathbf{z}_{i}=(z_{i1},\ldots,z_{it})} is the \eqn{i}-th observation's parametric predictor vector, and \eqn{e_{i}\sim\mathrm{N}(0,\sigma^{2})} is iid Gaussian error. Function can fit both additive and interactive non/parametric effects, and allows for 2-way and 3-way interactions between nonparametric and parametric effects (see Details and Examples).
}
\usage{
bigssp(formula,data=NULL,type=NULL,nknots=NULL,rparm=NA,
       lambdas=NULL,skip.iter=TRUE,se.fit=FALSE,rseed=1234,
       gcvopts=NULL,knotcheck=TRUE,thetas=NULL)
}
\arguments{
  \item{formula}{
	An object of class "\code{formula}": a symbolic description of the model to be fitted (see Details and Examples for more information).
}
  \item{data}{
  Optional data frame, list, or environment containing the variables in \code{formula}. Or an object of class "makessp", which is output from \code{\link{makessp}}.
}
  \item{type}{
	List of smoothing spline types for predictors in \code{formula} (see Details). Options include \code{type="cub"} for cubic, \code{type="cub0"} for another cubic, \code{type="per"} for cubic periodic, \code{type="tps"} for cubic thin-plate, and \code{type="nom"} for nominal. Use \code{type="prm"} for parametric effect.
}
  \item{nknots}{
	Two possible options: (a) scalar giving total number of random knots to sample, or (b) vector indexing which rows of \code{data} to use as knots.
}
  \item{rparm}{
	List of rounding parameters for each predictor. See Details. 
}
  \item{lambdas}{
	Vector of global smoothing parameters to try. Default \code{lambdas=10^-c(9:0)}.
}
  \item{skip.iter}{
	Logical indicating whether to skip the iterative smoothing parameter update. Using \code{skip.iter=FALSE} should provide a more optimal solution, but the fitting time may be substantially longer. See Skip Iteration section.
}
  \item{se.fit}{
	Logical indicating if the standard errors of the fitted values should be estimated.
}
  \item{rseed}{
  Random seed for knot sampling. Input is ignored if \code{nknots} is an input vector of knot indices. Set \code{rseed=NULL} to obtain a different knot sample each time, or set \code{rseed} to any positive integer to use a different seed than the default.
}
  \item{gcvopts}{
	Control parameters for optimization. List with 3 elements: (a) \code{maxit}: maximum number of algorithm iterations, (b) \code{gcvtol}: covergence tolerange for iterative GCV update, and (c) \code{alpha}: tuning parameter for GCV minimization. Default: \code{gcvopts=list(maxit=5,gcvtol=10^-5,alpha=1)}
}
  \item{knotcheck}{
If \code{TRUE}, only unique knots are used (for stability).  
}
  \item{thetas}{
List of initial smoothing parameters for each predictor subspace. See Details. 
}
}
\details{
 The \code{formula} syntax is similar to that used in \code{\link{lm}} and many other R regression functions. Use \code{y~x} to predict the response \code{y} from the predictor \code{x}. Use \code{y~x1+x2} to fit an additive model of the predictors \code{x1} and \code{x2}, and use \code{y~x1*x2} to fit an interaction model. The syntax \code{y~x1*x2} includes the interaction and main effects, whereas the syntax \code{y~x1:x2} only includes the interaction. See Computational Details for specifics about how non/parametric effects are estimated.
 
 See \code{\link{bigspline}} for definitions of \code{type="cub"}, \code{type="cub0"}, and \code{type="per"} splines, which can handle one-dimensional predictors. See Appendix of Helwig and Ma (under review) for information about \code{type="tps"} and \code{type="nom"} splines. Note that \code{type="tps"} can handle one-, two-, or three-dimensional predictors. I recommend using \code{type="cub"} if the predictor scores have no extreme outliers; when outliers are present, \code{type="tps"} may produce a better result. 
  
  Using the rounding parameter input \code{rparm} can greatly speed-up and stabilize the fitting for large samples. For typical cases, I recommend using \code{rparm=0.01} for cubic and periodic splines, but smaller rounding parameters may be needed for particularly jagged functions. For thin-plate splines, the data are NOT transformed to the interval [0,1] before fitting, so rounding parameter should be on raw data scale. Also, for \code{type="tps"} you can enter one rounding parameter for each predictor dimension. Use \code{rparm=1} for nominal splines.
}
\value{
	\item{fitted.values}{Vector of fitted values corresponding to the original data points in \code{xvars} (if \code{rparm=NA}) or the rounded data points in \code{xunique} (if \code{rparm} is used).}
	\item{se.fit}{Vector of standard errors of \code{fitted.values} (if input \code{se.fit=TRUE)}.}
  \item{yvar}{Response vector.}
	\item{xvars}{List of predictors.}
  \item{type}{Type of smoothing spline that was used for each predictor.}
  \item{yunique}{Mean of \code{yvar} for unique points after rounding (if \code{rparm} is used).}
  \item{xunique}{Unique rows of \code{xvars} after rounding (if \code{rparm} is used).}
	\item{sigma}{Estimated error standard deviation, i.e., \eqn{\hat{\sigma}}.}
  \item{ndf}{Data frame with two elements: \code{n} is total sample size, and \code{df} is effective degrees of freedom of fit model (trace of smoothing matrix).}
  \item{info}{Model fit information: vector containing the GCV, multiple R-squared, AIC, and BIC of fit model (assuming Gaussian error).}
  \item{modelspec}{List containing specifics of fit model (needed for prediction).}
  \item{converged}{Convergence status: \code{converged=TRUE} if iterative update converged, \code{converged=FALSE} if iterative update failed to converge, and \code{converged=NA} if option \code{skip.iter=TRUE} was used.}
  \item{tnames}{Names of the terms in model.}
  \item{call}{Called model in input \code{formula}.}
}
\references{
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.

  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.
	
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
	Cubic and cubic periodic splines transform the predictor to the interval [0,1] before fitting.
  
  When using rounding parameters, output \code{fitted.values} corresponds to unique rounded predictor scores in output \code{xunique}. Use \code{\link{predict.ssp}} function to get fitted values for full \code{yvar} vector.
}

\section{Computational Details }{
	To estimate \eqn{\eta} I minimize the penalized least-squares functional \deqn{\frac{1}{n}\sum_{i=1}^{n}\left(y_{i} - \eta(\mathbf{x}_{i}) - \textstyle\sum_{j=1}^{t}b_{j}z_{ij} \right)^{2} + \lambda J(\eta)} where \eqn{J(\cdot)} is a nonnegative penalty functional quantifying the roughness of \eqn{\eta} and \eqn{\lambda>0} is a smoothing parameter controlling the trade-off between fitting and smoothing the data. Note that for \eqn{p>1} nonparametric predictors, there are additional \eqn{\theta_{k}} smoothing parameters embedded in \eqn{J}. 
  
  The penalized least squares functioncal can be rewritten as \deqn{ \|\mathbf{y} - \mathbf{K}\mathbf{d} - \mathbf{J}_{\theta}\mathbf{c}\|^{2} + n\lambda\mathbf{c}'\mathbf{Q}_{\theta}\mathbf{c} }
where \eqn{\mathbf{K}=\{\phi(x_{i}),\mathbf{z}_{i}\}_{n \times m}} is the parametric space basis function matrix, \eqn{\mathbf{J}_{\theta}=\sum_{k=1}^{s}\theta_{k}\mathbf{J}_{k}} with \eqn{\mathbf{J}_{k}=\{\rho_{k}(\mathbf{x}_{i},\mathbf{x}_{h}^{*})\}_{n \times q}} denoting the \eqn{k}-th contrast space basis funciton matrix, \eqn{\mathbf{Q}_{\theta}=\sum_{k=1}^{s}\theta_{k}\mathbf{Q}_{k}} with \eqn{\mathbf{Q}_{k}=\{\rho_{k}(\mathbf{x}_{g}^{*},\mathbf{x}_{h}^{*})\}_{q \times q}} denoting the \eqn{k}-th penalty matrix, and \eqn{\mathbf{d}=(d_{0},\ldots,d_{m})'} and \eqn{\mathbf{c}=(c_{1},\ldots,c_{q})'} are the unknown basis function coefficients. The optimal smoothing parameters are chosen by minimizing the GCV score (see \code{\link{bigspline}}). 

Note that this function uses the classic smoothing spline parameterization (see Gu, 2013), so there is more than one smoothing parameter per predictor (if interactions are included in the model). To evaluate the GCV score, this function uses the improved (scalable) SSA algorithm discussed in Helwig (2013) and Helwig and Ma (under review).
}

\section{Skip Iteration }{
For \eqn{p>1} predictors, initial values for the \eqn{\theta_{k}} parameters are estimated using Algorithm 3.2 described in Gu and Wahba (1991). 

Default use of this function (\code{skip.iter=TRUE}) fixes the \eqn{\theta_{k}} parameters afer the smart start, and then finds the global smoothing parameter \eqn{\lambda} (among the input \code{lambdas}) that minimizes the GCV score. This approach typically produces a solution very similar to the more optimal solution using \code{skip.iter=FALSE}.

Setting \code{skip.iter=FALSE} uses the same smart starting algorithm as setting \code{skip.iter=TRUE}. However, instead of fixing the \eqn{\theta_{k}} parameters afer the smart start, using \code{skip.iter=FALSE} iterates between estimating the optimal \eqn{\lambda} and the optimal \eqn{\theta_{k}} parameters. The R function \code{\link{nlm}} is used to minimize the GCV score with respect to the \eqn{\theta_{k}} parameters, which can be time consuming for models with many predictors.
}

\examples{

##########   EXAMPLE   ##########

# function with four continuous predictors
set.seed(773)
myfun<-function(x1v,x2v,x3v,x4v){
  sin(2*pi*x1v)+log(x2v+.1)+x3v*cos(pi*(x4v))
  }
x1v<-runif(500);     x2v<-runif(500)
x3v<-runif(500);     x4v<-runif(500)
y<-myfun(x1v,x2v,x3v,x4v)+rnorm(500)

# fit cubic spline model with x3v*x4v interaction and x3v as "cub" 
# (includes x3v and x4v main effects)
cubmod<-bigssp(y~x1v+x2v+x3v*x4v,type=list(x1v="cub",x2v="cub",x3v="cub",x4v="cub"),nknots=50)
crossprod(myfun(x1v,x2v,x3v,x4v)-cubmod$fitted.values)/500

# fit cubic spline model with x3v*x4v interaction and x3v as "cub0"
# (includes x3v and x4v main effects)
cubmod<-bigssp(y~x1v+x2v+x3v*x4v,type=list(x1v="cub",x2v="cub",x3v="cub0",x4v="cub"),nknots=50)
crossprod(myfun(x1v,x2v,x3v,x4v)-cubmod$fitted.values)/500

# fit model with x3v*x4v interaction treating x3v as parametric effect
# (includes x3v and x4v main effects)
cubmod<-bigssp(y~x1v+x2v+x3v*x4v,type=list(x1v="cub",x2v="cub",x3v="prm",x4v="cub"),nknots=50)
crossprod(myfun(x1v,x2v,x3v,x4v)-cubmod$fitted.values)/500

# fit cubic spline model with x3v:x4v interaction and x3v as "cub"
# (excludes x3v and x4v main effects)
cubmod<-bigssp(y~x1v+x2v+x3v:x4v,type=list(x1v="cub",x2v="cub",x3v="cub",x4v="cub"),nknots=50)
crossprod(myfun(x1v,x2v,x3v,x4v)-cubmod$fitted.values)/500

# fit cubic spline model with x3v:x4v interaction and x3v as "cub0"
# (excludes x3v and x4v main effects)
cubmod<-bigssp(y~x1v+x2v+x3v:x4v,type=list(x1v="cub",x2v="cub",x3v="cub0",x4v="cub"),nknots=50)
crossprod(myfun(x1v,x2v,x3v,x4v)-cubmod$fitted.values)/500

# fit model with x3v:x4v interaction treating x3v as parametric effect
# (excludes x3v and x4v main effects)
cubmod<-bigssp(y~x1v+x2v+x3v:x4v,type=list(x1v="cub",x2v="cub",x3v="prm",x4v="cub"),nknots=50)
crossprod(myfun(x1v,x2v,x3v,x4v)-cubmod$fitted.values)/500


}

