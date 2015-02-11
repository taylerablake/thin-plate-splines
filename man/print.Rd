\name{print}
\alias{print.bigspline}
\alias{print.bigssa}
\alias{print.bigssg}
\alias{print.bigssp}
\alias{print.bigtps}
\alias{print.summary.bigspline}
\alias{print.summary.bigssa}
\alias{print.summary.bigssg}
\alias{print.summary.bigssp}
\alias{print.summary.bigtps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints Fit Information for bigsplines Model
}
\description{
This function prints basic model fit information for a fit \code{bigsplines} model.
}
\usage{
\method{print}{bigspline}(x,...)
\method{print}{bigssa}(x,...)
\method{print}{bigssg}(x,...)
\method{print}{bigssp}(x,...)
\method{print}{bigtps}(x,...)
\method{print}{summary.bigspline}(x,digits=4,...)
\method{print}{summary.bigssa}(x,digits=4,...)
\method{print}{summary.bigssg}(x,digits=4,...)
\method{print}{summary.bigssp}(x,digits=4,...)
\method{print}{summary.bigtps}(x,digits=4,...)
}
\arguments{
  \item{x}{
  	Object of class "bigspline" (output from \code{\link{bigspline}}), class "summary.bigspline" (output from \code{\link{summary.bigspline}}), class "bigssa" (output from \code{\link{bigssa}}), class "summary.bigssa" (output from \code{\link{summary.bigssa}}), class "bigssg" (output from \code{\link{bigssg}}), class "summary.bigssg" (output from \code{\link{summary.bigssg}}), class "bigssp" (output from \code{\link{bigssp}}), class "summary.bigssp" (output from \code{\link{summary.bigssp}}), class "bigtps" (output from \code{\link{bigtps}}), class "summary.bigtps" (output from \code{\link{summary.bigtps}}).
}
  \item{digits}{
  Number of decimal places to print.
} 
  \item{\ldots}{
  Ignored.
}
}
\value{
  "bigspline" objects: prints Spline Type, Fit Statistic information, and Smoothing Parameter.
  
  "summary.bigspline" objects: prints Spline Type, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameter.
  
  "bigssa" objects: prints Spline Types, Fit Statistic information, and Algorithm Convergence status.
  
  "summary.bigssa" objects: prints the formula Call, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameters.
  
  "bigssg" objects: prints Family, Spline Types, Fit Statistic information, and Algorithm Convergence status.
  
  "summary.bigssg" objects: prints the Family, formula Call, five number summary of Residuals, Dispersion Estimate, Fit Statistics, and Smoothing Parameters (with selection criterion).
  
  "bigssp" objects: prints Predictor Types, Fit Statistic information, and Algorithm Convergence status.
  
  "summary.bigssp" objects: prints formula Call, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameters.
  
  "bigtps" objects: prints Spline Type, Fit Statistic information, and Smoothing Parameter.
  
  "summary.bigtps" objects: prints Spline Type, five number summary of Residuals, Error Standard Deviation Estimate, Fit Statistics, and Smoothing Parameter.
}
\details{
	See \code{\link{bigspline}}, \code{\link{bigssa}}, \code{\link{bigssg}}, \code{\link{bigssp}}, and \code{\link{bigtps}} for more details.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\examples{

### see examples for bigspline, bigssa, bigssg, bigssp, and bigtps

}

