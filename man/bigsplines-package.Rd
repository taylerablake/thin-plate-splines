\name{bigsplines-package}
\alias{bigsplines-package}
\alias{bigsplines}
\docType{package}
\title{Smoothing Splines for Large Samples}
\description{
Fits smoothing spline regression models using scalable algorithms designed for large samples. Six marginal spline types are supported: cubic, different cubic, cubic periodic, cubic thin-plate, ordinal, and nominal. Random effects and parametric predictors are also supported. Response can be Gaussian or non-Gaussian: Binomial, Poisson, Gamma, Inverse Gaussian, or Negative Binomial.
}
\details{
%\tabular{ll}{
%Package: \tab bigsplines\cr
%Type: \tab Package\cr
%Version: \tab 1.0\cr
%Date: \tab 2013-08-31\cr
%License: \tab GPL (>=2)\cr
%}
The function \code{\link{bigspline}} fits one-dimensional cubic smoothing splines (unconstrained or periodic). The function \code{\link{bigssa}} fits Smoothing Spline Anova (SSA) models (Gaussian data). The function \code{\link{bigssg}} fits Generalized Smoothing Spline Anova (GSSA) models (non-Gaussian data). The function \code{\link{bigssp}} is for fitting Smoothing Splines with Parametric effects (semi-parametric regression). The function \code{\link{bigtps}} fits one-, two-, and three-dimensional cubic thin-plate splines. There are corresponding predict, print, and summary functions for these methods.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>

Maintainer: Nathaniel E. Helwig <helwig@umn.edu>
}
\references{
	Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.
  
  Gu, C. and Wahba, G. (1991). Minimizing GCV/GML scores with multiple smoothing parameters via the Newton method. \emph{SIAM Journal on Scientific and Statistical Computing, 12}, 383-398.
	
  Gu, C. and Xiang, D. (2001). Cross-validating non-Gaussian data: Generalized approximate cross-validation revisited. \emph{Journal of Computational and Graphical Statistics, 10}, 581-591.
  
  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.
  
  Helwig, N. E. (2016). Efficient estimation of variance components in nonparametric mixed-effects models with large samples. \emph{Statistics and Computing, 26}, 1319-1336.
  
	Helwig, N. E. and Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics, 24}, 715-732.
	
  Helwig, N. E. and Ma, P. (2016). Smoothing spline ANOVA for super-large samples: Scalable computation via rounding parameters. \emph{Statistics and Its Interface, 9}, 433-444.
}
\keyword{ package }
\examples{
# See examples for bigspline, bigssa, bigssg, bigssp, and bigtps
}
