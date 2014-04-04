\name{print}
\alias{print.css}
\alias{print.ssa}
\alias{print.ssp}
\alias{print.tps}
\alias{print.summary.css}
\alias{print.summary.ssa}
\alias{print.summary.ssp}
\alias{print.summary.tps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints fit information for big splines model
}
\description{
This function prints basic model fit information for a fit \code{bigsplines} model.
}
\usage{
\method{print}{css}(x,...)
\method{print}{ssa}(x,...)
\method{print}{ssp}(x,...)
\method{print}{tps}(x,...)
\method{print}{summary.css}(x,digits=4,...)
\method{print}{summary.ssa}(x,digits=4,...)
\method{print}{summary.ssp}(x,digits=4,...)
\method{print}{summary.tps}(x,digits=4,...)
}
\arguments{
  \item{x}{
  	Object of class "css" (output from \code{\link{bigspline}}), class "summary.css" (output from \code{\link{summary.css}}), class "ssa" (output from \code{\link{bigssa}}), class "summary.ssa" (output from \code{\link{summary.ssa}}), class "ssp" (output from \code{\link{bigssp}}), class "summary.ssp" (output from \code{\link{summary.ssp}}), class "tps" (output from \code{\link{bigtps}}), class "summary.tps" (output from \code{\link{summary.tps}}).
}
  \item{digits}{
  Number of decimal places to print.
} 
  \item{\ldots}{
  Ignored.
}
}
\value{
  "css" objects: prints Spline Type, Fit Statistic information, and Smoothing Parameter.
  
  "summary.css" objects: prints Spline Type, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameter.
  
  "ssa" objects: prints Spline Types, Fit Statistic information, and Algorithm Convergence status.
  
  "summary.ssa" objects: prints the formula Call, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameters.
  
  "ssp" objects: prints Predictor Types, Fit Statistic information, and Algorithm Convergence status.
  
  "summary.ssp" objects: prints formula Call, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameters.
  
  "tps" objects: prints Spline Type, Fit Statistic information, and Smoothing Parameter.
  
  "summary.tps" objects: prints Spline Type, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameter.
}
\details{
	See \code{\link{bigspline}}, \code{\link{bigssa}}, \code{\link{bigssp}}, and \code{\link{bigtps}} for more details.
}
\author{
Nathaniel E. Helwig <nhelwig2@illinois.edu>
}
\examples{

### see examples for bigspline, bigssa, bigssp, and bigtps

}

