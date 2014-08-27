\name{bigssg}
\alias{bigssg}
\title{
Fits generalized smoothing spline anova models
}
\description{
Given an exponential family response vector \eqn{\mathbf{y}=\{y_{i}\}_{n\times1}}, a Generalized Smoothing Spline Anova (GSSA) has the form \deqn{g(\mu_{i}) = \eta(\mathbf{x}_{i})} where \eqn{\mu_{i}} is the expected value of the \eqn{i}-th observation's respone, \eqn{g(\cdot)} is some invertible link function, \eqn{\mathbf{x}_{i}=(x_{i1},\ldots,x_{ip})} is the \eqn{i}-th observation's nonparametric predictor vector, and \eqn{\eta} is an unknown smooth function relating the response and nonparametric predictors. Function can fit additive models, and also allows for 2-way and 3-way interactions between any number of predictors. Response can be one of five non-Gaussian distributions: Binomial, Poisson, Gamma, Inverse Gaussian, or Negative Binomial (see Details and Examples).
}
\usage{
bigssg(formula,family,data=NULL,type=NULL,nknots=NULL,rparm=NA,
       lambdas=NULL,skip.iter=TRUE,se.lp=FALSE,rseed=1234,
       gcvopts=NULL,knotcheck=TRUE,gammas=NULL,weights=NULL)
}
\arguments{
  \item{formula}{
	An object of class "\code{formula}": a symbolic description of the model to be fitted (see Details and Examples for more information).
}
  \item{family}{
  Distribution for response. One of five options: \code{"binomial"}, \code{"poisson"}, \code{"Gamma"}, \code{"inverse.gaussian"}, or \code{"negbin"}. See Response section.
}
  \item{data}{
  Optional data frame, list, or environment containing the variables in \code{formula}. Or an object of class "makessg", which is output from \code{\link{makessg}}.
}
  \item{type}{
	List of smoothing spline types for predictors in \code{formula} (see Details). Options include \code{type="cub"} for cubic, \code{type="cub0"} for another cubic, \code{type="per"} for cubic periodic, \code{type="tps"} for cubic thin-plate, and \code{type="nom"} for nominal.
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
}
\details{
 The \code{formula} syntax is similar to that used in \code{\link{lm}} and many other R regression functions. Use \code{y~x} to predict the response \code{y} from the predictor \code{x}. Use \code{y~x1+x2} to fit an additive model of the predictors \code{x1} and \code{x2}, and use \code{y~x1*x2} to fit an interaction model. The syntax \code{y~x1*x2} includes the interaction and main effects, whereas the syntax \code{y~x1:x2} is not supported. See Computational Details for specifics about how nonparametric effects are estimated.
 
 See \code{\link{bigspline}} for definitions of \code{type="cub"}, \code{type="cub0"}, and \code{type="per"} splines, which can handle one-dimensional predictors. See Appendix of Helwig and Ma (in press) for information about \code{type="tps"} and \code{type="nom"} splines. Note that \code{type="tps"} can handle one-, two-, or three-dimensional predictors. I recommend using \code{type="cub"} if the predictor scores have no extreme outliers; when outliers are present, \code{type="tps"} may produce a better result. 
  
 Using the rounding parameter input \code{rparm} can greatly speed-up and stabilize the fitting for large samples. For typical cases, I recommend using \code{rparm=0.01} for cubic and periodic splines, but smaller rounding parameters may be needed for particularly jagged functions. For thin-plate splines, the data are NOT transformed to the interval [0,1] before fitting, so rounding parameter should be on raw data scale. Also, for \code{type="tps"} you can enter one rounding parameter for each predictor dimension. Use \code{rparm=1} for nominal splines.
}
\value{
	\item{fitted.values}{Vector of fitted values (data scale) corresponding to the original data points in \code{xvars} (if \code{rparm=NA}) or the rounded data points in \code{xunique} (if \code{rparm} is used).}
  \item{linear.predictors}{Vector of fitted values (link scale) corresponding to the original data points in \code{xvars} (if \code{rparm=NA}) or the rounded data points in \code{xunique} (if \code{rparm} is used).}
	\item{se.lp}{Vector of standard errors of \code{linear.predictors} (if input \code{se.lp=TRUE)}.}
  \item{yvar}{Response vector.}
	\item{xvars}{List of predictors.}
  \item{type}{Type of smoothing spline that was used for each predictor.}
  \item{yunique}{Mean of \code{yvar} for unique points after rounding (if \code{rparm} is used).}
  \item{xunique}{Unique rows of \code{xvars} after rounding (if \code{rparm} is used).}
	\item{dispersion}{Estimated dispersion parameter (see Response section).}
  \item{ndf}{Data frame with two elements: \code{n} is total sample size, and \code{df} is effective degrees of freedom of fit model (trace of smoothing matrix).}
  \item{info}{Model fit information: vector containing the GCV, multiple R-squared, AIC, and BIC of fit model.}
  \item{modelspec}{List containing specifics of fit model (needed for prediction).}
  \item{converged}{Convergence status: \code{converged=TRUE} if iterative update converged, \code{converged=FALSE} if iterative update failed to converge, and \code{converged=NA} if option \code{skip.iter=TRUE} was used.}
  \item{tnames}{Names of the terms in model.}
  \item{family}{Distribution family (same as input).}
  \item{call}{Called model in input \code{formula}.}
}
\references{
  Gu, C. and Xiang, D. (2001). Cross-validating non-Gaussian data: Generalized approximate cross-validation revisited. \emph{Journal of Computational and Graphical Statistics, 10}, 581-591.
  
  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.
  
  Helwig, N. E. (in preparation). Nonparametric exponential family regression for ultra large samples: Scalable computation via rounding parameters. 
	
	Helwig, N. E. and Ma, P. (in press). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
	
  Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
  
  Helwig, N. E. and Ma, P. (in preparation). Optimal eigenbasis truncation for thin-plate regression splines. 
  
	Helwig, N. E. and Ma, P. (in preparation). Nonparametric regression for ultra large samples: Scalable computation via rounding parameters. 
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
	The spline is estimated using penalized likelihood estimation. Standard errors of the linear predictors are formed using Bayesian confidence intervals.
}
\section{Warnings }{
	Cubic and cubic periodic splines transform the predictor to the interval [0,1] before fitting.
  
  When using rounding parameters, output \code{fitted.values} corresponds to unique rounded predictor scores in output \code{xunique}. Use \code{\link{predict.ssg}} function to get fitted values for full \code{yvar} vector.
}

\section{Response}{
Only one link is permitted for each family:

  \code{family="binomial"}
    Logit link. Response should be vector of proportions in the interval [0,1]. If response is a sample proportion, the total count should be input through \code{weights} argument.
  
  \code{family="poisson"}
    Log link. Response should be vector of counts (non-negative integers).
  
  \code{family="Gamma"}
    Inverse link. Response should be vector of positive real-valued data. Estimated \code{dispersion} parameter is the inverse of the \code{shape} parameter, so that the variance of the response increases as \code{dispersion} increases.
  
  \code{family="inverse.gaussian"}
    Inverse-square link. Response should be vector of positive real-valued data. Estimated \code{dispersion} parameter is the inverse of the \code{shape} parameter, so that the variance of the response increases as \code{dispersion} increases.
  
  \code{family="negbin"}
    Log link. Response should be vector of counts (non-negative integers). Estimated \code{dispersion} parameter is the inverse of the \code{size} parameter, so that the variance of the response increases as \code{dispersion} increases. 
  
  \code{family=list("negbin",2)}
    Log link. Response should be vector of counts (non-negative integers). Second element is the known (common) \code{dispersion} parameter (2 in this case). The input \code{dispersion} parameter should be the inverse of the \code{size} parameter, so that the variance of the response increases as \code{dispersion} increases.
  
}

\section{Computational Details }{
	To estimate \eqn{\eta} I minimize the (negative of the) penalized log likelihood \deqn{-\frac{1}{n}\sum_{i=1}^{n}\left\{y_{i}\eta(\mathbf{x}_{i}) - b(\eta(\mathbf{x}_{i})) \right\} + \frac{\lambda}{2} J(\eta)} where \eqn{J(\cdot)} is a nonnegative penalty functional quantifying the roughness of \eqn{\eta} and \eqn{\lambda>0} is a smoothing parameter controlling the trade-off between fitting and smoothing the data. Note that for \eqn{p>1} nonparametric predictors, there are additional \eqn{\theta_{k}} smoothing parameters embedded in \eqn{J}. 
  
  Following standard exponential family theory, \eqn{\mu_{i} = \dot{b}(\eta(\mathbf{x}_{i}))} and \eqn{v_{i} = \ddot{b}(\eta(\mathbf{x}_{i}))a(\xi)}, where \eqn{\dot{b}(\cdot)} and \eqn{\ddot{b}(\cdot)} denote the first and second derivatives of \eqn{b(\cdot)}, \eqn{v_{i}} is the variance of \eqn{y_{i}},and \eqn{\xi} is the dispersion parameter. Given fixed smoothing parameters, the optimal \eqn{\eta} can be estimated by iteratively minimizing the penalized reweighted least-squares functional \deqn{\frac{1}{n}\sum_{i=1}^{n}v_{i}^{*}\left(y_{i}^{*} - \eta(\mathbf{x}_{i}) \right)^{2} + \lambda J(\eta)} where \eqn{v_{i}^{*}=v_{i}/a(\xi)} is the weight, \eqn{y_{i}^{*}=\hat{\eta}(\mathbf{x}_{i})+(y_{i}-\hat{\mu}_{i})/v_{i}^{*}} is the adjusted dependent variable, and \eqn{\hat{\eta}(\mathbf{x}_{i})} is the current estimate of \eqn{\eta}.
  
  The optimal smoothing parameters are chosen by minimizing (an approximation of) the alternative Generalized Approximate Cross-Validation (GACV) score of Gu and Xiang (2001): \deqn{-\frac{1}{n}\sum_{i=1}^{n}\{y_{i}\hat{\eta}(\mathbf{x}_{i}) - b(\hat{\eta}(\mathbf{x}_{i}))\} + \frac{\mathrm{tr}(\mathbf{S}_{\boldsymbol\lambda}\mathbf{V}^{-1})}{n-\mathrm{tr}(\mathbf{S}_{\boldsymbol\lambda})}\frac{1}{n}\sum_{i=1}^{n}y_{i}(y_{i}-\hat{\mu}_{i}) } where \eqn{\mathbf{S}_{\boldsymbol\lambda}} is the smoothing matrix, and \eqn{\mathbf{V}=\mathrm{diag}(v_{1}^{*},\ldots,v_{n}^{*})}. We approximate \eqn{\frac{1}{n}\mathrm{tr}(\mathbf{S}_{\boldsymbol\lambda}\mathbf{V}^{-1})} using \eqn{\frac{1}{n^2}\mathrm{tr}(\mathbf{S}_{\boldsymbol\lambda})\mathrm{tr}(\mathbf{V}^{-1})} for faster computation.

Note that this function uses the efficient SSA reparameterization described in Helwig (2013) and Helwig and Ma (in press); using is parameterization, there is one unique smoothing parameter per predictor (\eqn{\gamma_{j}}), and these \eqn{\gamma_{j}} parameters determine the structure of the \eqn{\theta_{k}} parameters in the tensor product space. To evaluate the GACV score, this function uses the improved (scalable) GSSA algorithm discussed in Helwig (in preparation).
}

\section{Skip Iteration }{
For \eqn{p>1} predictors, initial values for the \eqn{\gamma_{j}} parameters (that determine the structure of the \eqn{\theta_{k}} parameters) are estimated using an extension of the smart starting algorithm described in Helwig (2013) and Helwig and Ma (in press). 

Default use of this function (\code{skip.iter=TRUE}) fixes the \eqn{\gamma_{j}} parameters afer the smart start, and then finds the global smoothing parameter \eqn{\lambda} (among the input \code{lambdas}) that minimizes the GCV score. This approach typically produces a solution very similar to the more optimal solution using \code{skip.iter=FALSE}.

Setting \code{skip.iter=FALSE} uses the same smart starting algorithm as setting \code{skip.iter=TRUE}. However, instead of fixing the \eqn{\gamma_{j}} parameters afer the smart start, using \code{skip.iter=FALSE} iterates between estimating the optimal \eqn{\lambda} and the optimal \eqn{\gamma_{j}} parameters. The R function \code{\link{nlm}} is used to minimize the approximate GACV score with respect to the \eqn{\gamma_{j}} parameters, which can be time consuming for models with many predictors and/or a large number of knots.
}

\examples{

##########   EXAMPLE 1 (1-way GSSA)   ##########

# define univariate function and data
set.seed(1)
myfun<-function(x){sin(2*pi*x)}
ndpts<-1000
x<-runif(ndpts)

# binomial response (no weights)
set.seed(773)
lp<-myfun(x)
p<-1/(1+exp(-lp))
y<-rbinom(n=ndpts,size=1,p=p)     ## y is binary data
gmod<-bigssg(y~x,family="binomial",type="cub",nknots=20)
crossprod(lp-gmod$linear.predictor)/length(lp)

# binomial response (with weights)
set.seed(773)
lp<-myfun(x)
p<-1/(1+exp(-lp))
w<-sample(c(10,20,30,40,50),length(p),replace=TRUE)
y<-rbinom(n=ndpts,size=w,p=p)/w   ## y is proportion correct
gmod<-bigssg(y~x,family="binomial",type="cub",nknots=20,weights=w)
crossprod(lp-gmod$linear.predictor)/length(lp)

# poisson response
set.seed(773)
lp<-myfun(x)
mu<-exp(lp)
y<-rpois(n=ndpts,lambda=mu)
gmod<-bigssg(y~x,family="poisson",type="cub",nknots=20)
crossprod(lp-gmod$linear.predictor)/length(lp)

# Gamma response
set.seed(773)
lp<-myfun(x)+2
mu<-1/lp
y<-rgamma(n=ndpts,shape=4,scale=mu/4)
gmod<-bigssg(y~x,family="Gamma",type="cub",nknots=20)
1/gmod$dispersion   ## dispersion = 1/shape
crossprod(lp-gmod$linear.predictor)/length(lp)

# inverse gaussian response (not run: requires statmod package)
# require(statmod)
# set.seed(773)
# lp<-myfun(x)+2
# mu<-sqrt(1/lp)
# y<-rinvgauss(n=ndpts,mean=mu,shape=2)
# gmod<-bigssg(y~x,family="inverse.gaussian",type="cub",nknots=20)
# 1/gmod$dispersion   ## dispersion = 1/shape
# crossprod(lp-gmod$linear.predictor)/length(lp)

# negative binomial response (known dispersion)
set.seed(773)
lp<-myfun(x)
mu<-exp(lp)
y<-rnbinom(n=ndpts,size=.5,mu=mu)
gmod<-bigssg(y~x,family=list("negbin",2),type="cub",nknots=20)
1/gmod$dispersion   ## dispersion = 1/size
crossprod(lp-gmod$linear.predictor)/length(lp)

# negative binomial response (unknown dispersion)
set.seed(773)
lp<-myfun(x)
mu<-exp(lp)
y<-rnbinom(n=ndpts,size=.5,mu=mu)
gmod<-bigssg(y~x,family="negbin",type="cub",nknots=20)
1/gmod$dispersion   ## dispersion = 1/size
crossprod(lp-gmod$linear.predictor)/length(lp)


##########   EXAMPLE 2 (2-way GSSA)   ##########

# function with two continuous predictors
set.seed(1)
myfun<-function(x1v,x2v){sin(2*pi*x1v)+log(x2v+.1)+cos(pi*(x1v-x2v))}
ndpts<-1000;    x1v<-runif(ndpts);    x2v<-runif(ndpts)

# binomial response (no weights)
set.seed(773)
lp<-myfun(x1v,x2v)
p<-1/(1+exp(-lp))
y<-rbinom(n=ndpts,size=1,p=p)     ## y is binary data
gmod<-bigssg(y~x1v*x2v,family="binomial",type=list(x1v="cub",x2v="cub"),nknots=50)
crossprod(lp-gmod$linear.predictor)/length(lp)

# binomial response (with weights)
set.seed(773)
lp<-myfun(x1v,x2v)
p<-1/(1+exp(-lp))
w<-sample(c(10,20,30,40,50),length(p),replace=TRUE)
y<-rbinom(n=ndpts,size=w,p=p)/w   ## y is proportion correct
gmod<-bigssg(y~x1v*x2v,family="binomial",type=list(x1v="cub",x2v="cub"),nknots=50,weights=w)
crossprod(lp-gmod$linear.predictor)/length(lp)

# poisson response
set.seed(773)
lp<-myfun(x1v,x2v)
mu<-exp(lp)
y<-rpois(n=ndpts,lambda=mu)
gmod<-bigssg(y~x1v*x2v,family="poisson",type=list(x1v="cub",x2v="cub"),nknots=50)
crossprod(lp-gmod$linear.predictor)/length(lp)

# Gamma response
set.seed(773)
lp<-myfun(x1v,x2v)+6
mu<-1/lp
y<-rgamma(n=ndpts,shape=4,scale=mu/4)
gmod<-bigssg(y~x1v*x2v,family="Gamma",type=list(x1v="cub",x2v="cub"),nknots=50)
1/gmod$dispersion   ## dispersion = 1/shape
crossprod(lp-gmod$linear.predictor)/length(lp)

# inverse gaussian response (not run: requires 'statmod' package)
# require(statmod)
# set.seed(773)
# lp<-myfun(x1v,x2v)+6
# mu<-sqrt(1/lp)
# y<-rinvgauss(n=ndpts,mean=mu,shape=2)
# gmod<-bigssg(y~x1v*x2v,family="inverse.gaussian",type=list(x1v="cub",x2v="cub"),nknots=50)
# 1/gmod$dispersion   ## dispersion = 1/shape
# crossprod(lp-gmod$linear.predictor)/length(lp)

# negative binomial response (known dispersion)
set.seed(773)
lp<-myfun(x1v,x2v)
mu<-exp(lp)
y<-rnbinom(n=ndpts,size=.5,mu=mu)
gmod<-bigssg(y~x1v*x2v,family=list("negbin",2),type=list(x1v="cub",x2v="cub"),nknots=50)
1/gmod$dispersion   ## dispersion = 1/size
crossprod(lp-gmod$linear.predictor)/length(lp)

# negative binomial response (unknown dispersion)
set.seed(773)
lp<-myfun(x1v,x2v)
mu<-exp(lp)
y<-rnbinom(n=ndpts,size=.5,mu=mu)
gmod<-bigssg(y~x1v*x2v,family="negbin",type=list(x1v="cub",x2v="cub"),nknots=50)
1/gmod$dispersion   ## dispersion = 1/size
crossprod(lp-gmod$linear.predictor)/length(lp)


}

