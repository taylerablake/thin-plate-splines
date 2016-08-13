\name{plotbar}
\alias{plotbar}
\title{
Generic X-Y Plotting with Colorbar
}
\description{
This is a modification to the R function \code{\link{plot}} that adds a colorbar to the margin.
}
\usage{
plotbar(x,y,z,xlim=NULL,ylim=NULL,zlim=NULL,
        zlab=NULL,zcex.axis=NULL,zcex.lab=NULL,
        zaxis.at = NULL, zaxis.labels = TRUE,
        col=NULL,ncolor=21,drawbar=TRUE,zline=2,
        pltimage=c(.2,.8,.2,.8),pltbar=c(.82,.85,.2,.8),...)
}
\arguments{
  \item{x, y}{
	The x and y coordinates of the points to plot.
}
  \item{z}{
  Numeric vector the same length as \code{x} and \code{y} containing the values to be plotted in color.
}
  \item{xlim, ylim}{
  	Ranges for the plotted \code{x} and \code{y} values, defaulting to the ranges of \code{x} and \code{y}.
}
  \item{zlim}{
  	The minimum and maximum \code{z} values for which colors should be plotted, defaulting to the range of the finite values of \code{z}.
}
  \item{zlab}{
  	Label for the colorbar.
}
  \item{zcex.axis}{
  The magnification to be used for the z-axis annotation (colorbar scale).
}
  \item{zcex.lab}{
  The magnification to be used for the z-axis label (\code{zlab}).
}
  \item{zaxis.at}{
  The points at which tick-marks are to be drawn for the colorbar. Points outside of the range of \code{zlim} will not be plotted.
}
  \item{zaxis.labels}{
  This can either be a logical value specifying whether (numerical) annotations are to be made at the tickmarks, or a character or expression vector of labels to be placed at the tickpoints.
}  
  \item{col}{
  	Color scheme to use. Default is from \code{blueviolet} (low) to \code{red} (high).
}
  \item{ncolor}{
  	The number of colors to use in the color scheme.
}
  \item{drawbar}{
    Logical indicating if the colorbar should be drawn.
}
  \item{zline}{
    Number of lines into the margin at which the axis line will be drawn (see \code{\link{axis}}).
}
  \item{pltimage}{
    A vector of the form c(x1, x2, y1, y2) giving the coordinates of the image region as fractions of the current figure region (see \code{\link{par}}).
}
  \item{pltbar}{
    A vector of the form c(x1, x2, y1, y2) giving the coordinates of the colorbar region as fractions of the current figure region (see \code{\link{par}}).
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

\examples{
##########   EXAMPLE 1   ##########

myfun <- function(x){
  2*sin(sqrt(x[,1]^2+x[,2]^2+.1))/sqrt(x[,1]^2+x[,2]^2+.1)
}
x <- expand.grid(seq(-8,8,l=100),seq(-8,8,l=100))
plotbar(x[,1],x[,2],myfun(x),
        xlab=expression(italic(x)[1]),
        ylab=expression(italic(x)[2]),
        zlab=expression(hat(italic(y))))


##########   EXAMPLE 2   ##########

myfun <- function(x1v,x2v){
  sin(2*pi*x1v) + 2*sin(sqrt(x2v^2+.1))/sqrt(x2v^2+.1)
}
x <- expand.grid(x1v=seq(0,1,l=100),x2v=seq(-8,8,l=100))
plotbar(x[,1],x[,2],myfun(x$x1v,x$x2v),
         col=c("red","orange","yellow","white"),
         xlab="x1v",ylab="x2v",zlab=expression(hat(italic(y))))


##########   EXAMPLE 3   ##########

myfun <- function(x1v,x2v){
  sin(3*pi*x1v) + sin(2*pi*x2v) + 3*cos(pi*(x1v-x2v))
}
x <- expand.grid(x1v=seq(-1,1,l=100),x2v=seq(-1,1,l=100))
plotbar(x[,1],x[,2],myfun(x$x1v,x$x2v),
         col=c("blue","green","light green","yellow"),
         xlab="x1v",ylab="x2v",zlab=expression(hat(italic(y))))

}

