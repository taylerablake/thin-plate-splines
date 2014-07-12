\name{bigsplines-package}
\alias{bigsplines-package}
\alias{bigsplines}
\docType{package}
\title{Smoothing splines for large samples}
\description{
Fits smoothing spline regression models using scalable algorithms designed for large samples. Options include: (a) spline knots can be user-selected, randomly sampled, or bin-sampled throughout the covariate domain, (b) rounding parameters can be provided, and (c) parametric effects can be included. Supports 2-way and 3-way interactions between any combination of nonparametric and parametric predictors.
}
\details{
%\tabular{ll}{
%Package: \tab bigsplines\cr
%Type: \tab Package\cr
%Version: \tab 1.0\cr
%Date: \tab 2013-08-31\cr
%License: \tab GPL (>=2)\cr
%}
The function \code{\link{bigspline}} fits one-dimensional cubic smoothing splines (unconstrained or periodic). The function \code{\link{bigssa}} fits Smoothing Spline Anova (SSA) models. The function \code{\link{bigssp}} is a more general function for fitting Smoothing Splines with Parametric effects (semi-parametric regression). The function \code{\link{bigtps}} fits one-, two-, and three-dimensional cubic thin-plate splines. There are corresponding predict, print, and summary functions for these methods.
}
\author{
Nathaniel E. Helwig <nhelwig2@illinois.edu>

Maintainer: Nathaniel E. Helwig <nhelwig2@illinois.edu>
}
\references{
	Gu, C. (2013). \emph{Smoothing spline ANOVA models, 2nd edition}. New York: Springer.
	
  Helwig, N. E. (2013). \emph{Fast and stable smoothing spline analysis of variance models for large samples with applications to electroencephalography data analysis}. Unpublished doctoral dissertation. University of Illinois at Urbana-Champaign.
  
	Helwig, N. E. and Ma, P. (in press). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics}.
	
  Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
  
  Helwig, N. E. and Ma, P. (in preparation). Optimal eigenbasis truncation for thin-plate regression splines. 
  
	Helwig, N. E. and Ma, P. (in preparation). Nonparametric regression for ultra large samples: Scalable computation via rounding parameters. 
}
\keyword{ package }
\examples{
# See examples for bigspline, bigssa, bigssp, and bigtps
}
