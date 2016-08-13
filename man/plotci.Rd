\name{plotci}
\alias{plotci}
\title{
Generic X-Y Plotting with Confidence Intervals
}
\description{
This is a modification to the R function \code{\link{plot}} that adds confidence intervals to the plot.
}
\usage{
plotci(x, y, se, level=0.95, cval=NULL, col="blue",
       col.ci="cyan", alpha=0.65, add=FALSE, 
       type="l", link=function(y){y}, axes=TRUE, 
       bars=FALSE, barlty=1, barlwd=2, bw=0.2, ...)
}
\arguments{
  \item{x, y}{
	The x and y coordinates of the points to plot.
}
  \item{se}{
  Numeric vector the same length as \code{x} and \code{y} containing the standard errors of the \code{y} values.
}
  \item{level}{
  	Significance level for the confidence interval. Default forms 95\% interval.
}
  \item{cval}{
  	Critical value for the confidence interval. Default uses \code{cval=qnorm(1-(1-level)/2)}.
}
  \item{col}{
  	Color for plotting the relationship between \code{x} and \code{y}.
}
  \item{col.ci}{
  	Color for plotting the confidence interval.
}
  \item{alpha}{
  Transparency used for plotting confidence polygons. Only used when \code{bars=FALSE}.
}
  \item{add}{
  Logical indicating whether lines should be added to current plot.
}
  \item{type}{
  	Type of plot to create (defaults to "l" for lines).
}
  \item{link}{
  	Link function to apply. See Details.
}
  \item{axes}{
    Logical indicating if the axes should be drawn.
}
  \item{bars}{
    Logical indicating if confidence bars should be plotted instead of polygons.
}
  \item{barlty, barlwd}{
    Line type and width for confidence bars. Only used when \code{bars=TRUE}.
}
  \item{bw}{
    Positive scalar giving the width of the confidence bars. Only used when \code{bars=TRUE}.
}
  \item{\dots}{
  	Additional arguments to be passed to \code{plot} (e.g., \code{xlab}, \code{ylab}, \code{main}, \code{cex}, \code{cex.axis}, \code{cex.lab}, etc.)
}
}

\value{
	Produces a \code{plot} with a colorbar.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\details{
The plotted confidence interval is \code{c(link(y-cval*se), link(y+cval*se))} where \code{link} is the user-specified link function and \code{cval} is the user-sepcified critival value, which defaults to \code{cval = qnorm(1-(1-level)/2)}.
}

\examples{
##########   EXAMPLE   ##########

# define relatively smooth function
set.seed(773)
myfun <- function(x){ sin(2*pi*x) }
x <- runif(10^4)
y <- myfun(x) + rnorm(10^4)

# fit cubic smoothing spline
cubmod <- bigspline(x,y)
newdata <- data.frame(x=seq(0,1,length=20))
ypred <- predict(cubmod, newdata, se.fit=TRUE)

# plot predictions with CIs in two ways
plotci(newdata$x, ypred$fit, ypred$se.fit)
plotci(newdata$x, ypred$fit, ypred$se.fit, type="p", bars=TRUE, bw=0.02)

}