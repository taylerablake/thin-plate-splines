\name{bigspline}
\alias{bigspline}
\title{
Fits Cubic Smoothing Splines
}
\description{
Given a real-valued response vector \eqn{\mathbf{y}=\{y_{i}\}_{n\times1}} and a real-valued predictor vector \eqn{\mathbf{x}=\{x_{i}\}_{n\times 1}} with \eqn{a \leq x_{i} \leq b \ \forall i}, a smoothing spline model has the form \deqn{y_{i}=\eta(x_{i})+e_{i}} where \eqn{y_{i}} is the \eqn{i}-th observation's respone, \eqn{x_{i}} is the \eqn{i}-th observation's predictor, \eqn{\eta} is an unknown smooth function relating the response and predictor, and \eqn{e_{i}\sim\mathrm{N}(0,\sigma^{2})} is iid Gaussian error.
}
\usage{
bigspline(x,y,type="cub",nknots=30,rparm=0.01,xmin=min(x),
          xmax=max(x),alpha=1,lambdas=NULL,se.fit=FALSE,
          rseed=1234,knotcheck=TRUE)
}
\arguments{
  \item{x}{
	Predictor vector.
}
  \item{y}{
  Response vector. Must be same length as \code{x}.
}
  \item{type}{
	Type of spline for \code{x}. Options include \code{type="cub"} for cubic, \code{type="cub0"} for different cubic, and \code{type="per"} for cubic periodic. See Spline Types section.
}
  \item{nknots}{
	Scalar giving maximum number of knots to bin-sample. Use more knots for more jagged functions.
}
  \item{rparm}{
	Rounding parameter for \code{x}. Use \code{rparm=NA} to fit unrounded solution. Rounding parameter must be in interval (0,1].
}
  \item{xmin}{
  Minimum \code{x} value (i.e., \eqn{a}). Used to transform data to interval [0,1].
}
  \item{xmax}{
  Maximum \code{x} value (i.e., \eqn{b}). Used to transform data to interval [0,1].
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
  Random seed. Input to \code{\link{set.seed}} to reproduce same knots when refitting same model. Use \code{rseed=NULL} to generate a different sample of knots each time.
}
  \item{knotcheck}{
If \code{TRUE}, only unique knots are used (for stability).  
}
}
\details{
  To estimate \eqn{\eta} I minimize the penalized least-squares functional \deqn{\frac{1}{n}\sum_{i=1}^{n}(y_{i}-\eta(x_{i}))^{2}+\lambda \int [\ddot{\eta}(x)]^2 dx} where \eqn{\ddot{\eta}} denotes the second derivative of \eqn{\eta} and \eqn{\lambda\geq0} is a smoothing parameter that controls the trade-off between fitting and smoothing the data. 
  
  Default use of the function estimates \eqn{\lambda} by minimizing the GCV score: \deqn{\mbox{GCV}(\lambda) = \frac{n\|(\mathbf{I}_{n}-\mathbf{S}_{\lambda})\mathbf{y}\|^{2}}{[n-\mathrm{tr}(\mathbf{S}_{\lambda})]^2}} where \eqn{\mathbf{I}_{n}} is the identity matrix and \eqn{\mathbf{S}_{\lambda}} is the smoothing matrix (see Computational Details).
  
  Using the rounding parameter input \code{rparm} can greatly speed-up and stabilize the fitting for large samples. When \code{rparm} is used, the spline is fit to a set of unique data points after rounding; the unique points are determined using the efficient algorithm described in Helwig (2013). For typical cases, I recommend using \code{rparm=0.01}, but smaller rounding parameters (e,g., \code{rparm=0.001}) may be needed for particularly jagged functions (or when \code{x} has outliers). 
}
\value{
	\item{fitted.values}{Vector of fitted values corresponding to the original data points in \code{x} (if \code{rparm=NA}) or the rounded data points in \code{xunique} (if \code{rparm} is used).}
	\item{se.fit}{Vector of standard errors of \code{fitted.values} (if input \code{se.fit=TRUE)}.}
  \item{x}{Predictor vector (same as input).}
	\item{y}{Response vector (same as input).}
  \item{type}{Type of spline that was used.}
  \item{xunique}{Unique elements of \code{x} after rounding (if \code{rparm} is used).}
  \item{yunique}{Mean of \code{y} for unique elements of \code{x} after rounding (if \code{rparm} is used).}
  \item{funique}{Vector giving frequency of each element of \code{xunique} (if \code{rparm} is used).}
	\item{sigma}{Estimated error standard deviation, i.e., \eqn{\hat{\sigma}}.}
  \item{ndf}{Data frame with two elements: \code{n} is total sample size, and \code{df} is effective degrees of freedom of fit model (trace of smoothing matrix).}
  \item{info}{Model fit information: vector containing the GCV, multiple R-squared, AIC, and BIC of fit model (assuming Gaussian error).}
  \item{xrng}{Predictor range: \code{xrng=c(xmin,xmax)}. }
  \item{myknots}{Bin-sampled spline knots used for fit.}
  \item{rparm}{Rounding parameter for \code{x} (same as input).}
  \item{lambda}{Optimal smoothing parameter.}
  \item{coef}{Spline basis function coefficients.}
  \item{coef.csqrt}{Matrix square-root of covariace matrix of \code{coef}. Use \code{tcrossprod(coef.csqrt)} to get covariance matrix of \code{coef}.}
}
\references{
  Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.

  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.
	
	Helwig, N. E. and Ma, P. (in press). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
	
  Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
  
  Helwig, N. E. and Ma, P. (in preparation). Optimal eigenbasis truncation for thin-plate regression splines. 
  
	Helwig, N. E. and Ma, P. (in preparation). Smoothing spline ANOVA for super-large samples: Scalable computation via rounding parameters. 
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
	The spline is estimated using penalized least-squares, which does not require the Gaussian error assumption. However, the spline inference information (e.g., standard errors and fit information) requires the Gaussian error assumption.
}
\section{Warnings }{
	Cubic and cubic periodic splines transform the predictor to the interval [0,1] before fitting. So input \code{xmin} must be less than or equal to \code{min(x)}, and input \code{xmax} must be greater than or equal to \code{max(x)}.
    
  When using rounding parameters, output \code{fitted.values} corresponds to unique rounded predictor scores in output \code{xunique}. Use \code{\link{predict.bigspline}} function to get fitted values for full \code{y} vector.
}

\section{Computational Details }{
  According to smoothing spline theory, the function \eqn{\eta} can be approximated as \deqn{\eta(x) = d_{0} + d_{1}\phi_{1}(x) + \sum_{h=1}^{q}c_{h}\rho(x,x_{h}^{*})} where the \eqn{\phi_{1}}  is a linear function, \eqn{\rho} is the reproducing kernel of the contrast (nonlinear) space, and \eqn{\{x_{h}^{*}\}_{h=1}^{q}} are the selected spline knots.
  
  This implies that the penalized least-squares functional can be rewritten as \deqn{ \|\mathbf{y} - \mathbf{K}\mathbf{d} - \mathbf{J}\mathbf{c}\|^{2} + n\lambda\mathbf{c}'\mathbf{Q}\mathbf{c} }
where \eqn{\mathbf{K}=\{\phi(x_{i})\}_{n \times 2}} is the null space basis function matrix, \eqn{\mathbf{J}=\{\rho(x_{i},x_{h}^{*})\}_{n \times q}} is the contrast space basis funciton matrix, \eqn{\mathbf{Q}=\{\rho(x_{g}^{*},x_{h}^{*})\}_{q \times q}} is the penalty matrix, and \eqn{\mathbf{d}=(d_{0},d_{1})'} and \eqn{\mathbf{c}=(c_{1},\ldots,c_{q})'} are the unknown basis function coefficients.

Given the smoothing parameter \eqn{\lambda}, the optimal basis function coefficients have the form 
\deqn{ \left(\begin{array}{cc} \hat{\mathbf{d}} \\ \hat{\mathbf{c}} \end{array}\right) =
\left(\begin{array}{cc} \mathbf{K'K} & \mathbf{K}'\mathbf{J} \\
 \mathbf{J}'\mathbf{K} & \mathbf{J}'\mathbf{J} + n\lambda\mathbf{Q} \end{array}\right)^{\dagger} \left(\begin{array}{c} \mathbf{K}' \\ \mathbf{J}' \end{array}\right)\mathbf{y} }
where \eqn{(\cdot)^{\dagger}} denotes the pseudoinverse of the input matrix.

Given the optimal coefficients, the fitted values are given by \eqn{\hat{\mathbf{y}} = \mathbf{K}\hat{\mathbf{d}}+\mathbf{J}\hat{\mathbf{c}} = \mathbf{S}_{\lambda}\mathbf{y}}, where \deqn{ \mathbf{S}_{\lambda} = \left(\begin{array}{cc} \mathbf{K} & \mathbf{J} \end{array}\right)
\left(\begin{array}{cc} \mathbf{K'K} & \mathbf{K}'\mathbf{J} \\
 \mathbf{J}'\mathbf{K} & \mathbf{J}'\mathbf{J} + n\lambda\mathbf{Q} \end{array}\right)^{\dagger} \left(\begin{array}{c} \mathbf{K}' \\ \mathbf{J}' \end{array}\right) }
 is the smoothing matrix, which depends on \eqn{\lambda}.

}
\section{Spline Types }{
  For a cubic spline (\code{type="cub"}) with \eqn{x \in [0,1]}, the needed functions are 
\deqn{ \phi_{1}(x) = k_{1}(x) \qquad \mbox{and} \qquad \rho(x,z) =  k_{2}(x)k_{2}(z)-k_{4}(|x-z|)} where \eqn{k_{1}(x)=x-0.5}, \eqn{k_{2}(x)=\frac{1}{2}\left(k_{1}^{2}(x) - \frac{1}{12} \right)}, and \eqn{k_{4}(x)=\frac{1}{24}\left(k_{1}^{4}(x) - \frac{k_{1}^{2}(x)}{2} + \frac{7}{240} \right)}. 
  
  For a different cubic spline (\code{type="cub0"}) with \eqn{x \in [0,1]}, the needed functions are \deqn{ \phi_{1}(x) = x \qquad \mbox{and} \qquad \rho(x,z) = (x \wedge z)^2[3(x \vee z) - (x \wedge z)]/6}
where \eqn{(x \wedge z) = \min(x,z)} and \eqn{(x \vee z) = \max(x,z)}. 

Note that \code{type="cub"} and \code{type="cub0"} use different definitions of the averaging operator in the null space. The overall spline estimates should be the same (up to approximation accuracy), but the null and constrast space effect functions will differ (see \code{\link{predict.bigspline}}). See Helwig (2013) and Gu (2013) for a further discussion of polynomial splines.  

For a periodic cubic spline (\code{type="per"}) with \eqn{x \in [0,1]}, the needed functions are 
\deqn{ \phi_{1}(x) = 0 \qquad \mbox{and} \qquad \rho(x,z) = -k_{4}(|x-z|)} where \eqn{k_{4}(x)} is defined as it was for \code{type="cub"}; in this case \eqn{\mathbf{K}=\mathbf{1}_{n}} and \eqn{\mathbf{d}=d_{0}}.

}




\examples{

##########   EXAMPLE 1   ##########

# define relatively smooth function
set.seed(773)
myfun<-function(x){sin(2*pi*x)}
x<-runif(10^6)
y<-myfun(x)+rnorm(10^6)

# cubic, different cubic, and periodic splines
cubmod<-bigspline(x,y)
cubmod
cub0mod<-bigspline(x,y,type="cub0")
cub0mod
permod<-bigspline(x,y,type="per")
permod


##########   EXAMPLE 2   ##########

# define more jagged function
set.seed(773)
myfun<-function(x){2*x+cos(4*pi*x)}
x<-runif(10^6)*4
y<-myfun(x)+rnorm(10^6)

# try different numbers of knots
r1mod<-bigspline(x,y,nknots=20)
crossprod(myfun(r1mod$xunique)-r1mod$fitted)/length(r1mod$fitted)
r2mod<-bigspline(x,y,nknots=30)
crossprod(myfun(r2mod$xunique)-r2mod$fitted)/length(r2mod$fitted)
r3mod<-bigspline(x,y,nknots=40)
crossprod(myfun(r3mod$xunique)-r3mod$fitted)/length(r3mod$fitted)

}

