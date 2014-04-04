\name{imagebar}
\alias{imagebar}
\title{
Displays a color image with colorbar
}
\description{
This is a modification to the R function \code{\link{image}} that adds a colorbar to the right margin.
}
\usage{
imagebar(x,y,z,xlim=NULL,ylim=NULL,zlim=NULL,
         zlab=NULL,zcex.axis=NULL,zcex.lab=NULL,
         col=NULL,ncolor=100,drawbar=TRUE,...)
}
\arguments{
  \item{x, y}{
	Locations of grid lines at which the values in \code{z} are measured. These must be finite, non-missing and in (strictly) ascending order.
}
  \item{z}{
  	A matrix containing the values to be plotted (\code{NA}s are allowed).
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
  \item{col}{
  	Color scheme to use. Default is reverse rainbow from violet (low) to red (high).
}
  \item{ncolor}{
  	The number of colors to use in the color scheme.
}
  \item{drawbar}{
    Logical indicating if the colorbar should be drawn.
}
\item{...}{
  	Additional arguments to be passed to \code{image} (e.g., \code{xlab}, \code{ylab}, \code{main}, \code{cex}, \code{cex.axis}, \code{cex.lab}, etc.)
}
}

\value{
	Produces an \code{image} plot with a colorbar.
}
\author{
Nathaniel E. Helwig <nhelwig2@illinois.edu>
}

\examples{
##########   EXAMPLE 1   ##########

myfun<-function(x){2*sin(sqrt(x[,1]^2+x[,2]^2+.1))/sqrt(x[,1]^2+x[,2]^2+.1)}
x<-expand.grid(seq(-8,8,l=100),seq(-8,8,l=100))
imagebar(seq(-8,8,l=100),seq(-8,8,l=100),matrix(myfun(x),100,100),
         xlab=expression(italic(x)[1]),ylab=expression(italic(x)[2]),
         zlab=expression(hat(italic(y))))


##########   EXAMPLE 2   ##########

myfun<-function(x1v,x2v){sin(2*pi*x1v)+2*sin(sqrt(x2v^2+.1))/sqrt(x2v^2+.1)}
x<-expand.grid(x1v=seq(0,1,l=100),x2v=seq(-8,8,l=100))
imagebar(seq(0,1,l=100),seq(-8,8,l=100),matrix(myfun(x$x1v,x$x2v),100,100),
         col=c("red","orange","yellow","white"),xlab="x1v",ylab="x2v",
         zlab=expression(hat(italic(y))))


##########   EXAMPLE 3   ##########

myfun<-function(x1v,x2v){sin(3*pi*x1v)+sin(2*pi*x2v)+3*cos(pi*(x1v-x2v))}
x<-expand.grid(x1v=seq(-1,1,l=100),x2v=seq(-1,1,l=100))
imagebar(seq(-1,1,l=100),seq(-1,1,l=100),matrix(myfun(x$x1v,x$x2v),100,100),
         col=c("blue","green","light green","yellow"),xlab="x1v",ylab="x2v",
         zlab=expression(hat(italic(y))))

}

