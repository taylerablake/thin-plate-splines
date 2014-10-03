\name{binsamp}
\alias{binsamp}
\title{
Bin-Samples Strategic Knot Indices
}
\description{
Breaks the predictor domain into a user-specified number of disjoint subregions, and randomly samples a user-specified number of observations from each (nonempty) subregion. 
}
\usage{
binsamp(x,xrng=NULL,nmbin=11,nsamp=1,alg=c("new","old"))
}
\arguments{
  \item{x}{
Matrix of predictors \eqn{\mathbf{X}=\{x_{ij}\}_{n \times p}} where \eqn{n} is the number of observations, and \eqn{p} is the number of predictors.
}
  \item{xrng}{
Optional matrix of predictor ranges: \eqn{\mathbf{R}=\{r_{kj}\}_{2 \times p}} where \eqn{r_{1j}=\min_{i}x_{ij}} and \eqn{r_{2j}=\max_{i}x_{ij}}.
}
  \item{nmbin}{
Vector \eqn{\mathbf{b}=(b_{1},\ldots,b_{p})'}, where \eqn{b_{j}\geq1} is the number of marginal bins to use for the \eqn{j}-th predictor. If \code{length(nmbin)<ncol(x)}, then \code{nmbin[1]} is used for all columns. Default is \code{nmbin=11} marginal bins for each dimension.
}
  \item{nsamp}{
Scalar \eqn{s\geq1} giving the number of observations to sample from each bin. Default is sample \code{nsamp=1} observation from each bin.
}
\item{alg}{
Bin-sampling algorithm. New algorithm forms equidistant grid, whereas old algorithm forms approximately equidistant grid. New algorithm is default for versions 1.0-1 and later.
}
}
\value{
Returns an index vector indicating the rows of \code{x} that were bin-sampled. 
}
\references{
Helwig, N. E. and Ma, P. (in preparation). Stable smoothing spline approximation via bin-sampled knots.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The number of returned knots will depend on the distribution of the covariate scores. The maximum number of possible bin-sampled knots is \eqn{s\prod_{j=1}^{p}b_{j}}, but fewer knots will be returned if one (or more) of the bins is empty (i.e., if there is no data in one or more bins).
}
\section{Warnings }{
If \eqn{x_{ij}} is nominal with \eqn{g} levels, the function requires \eqn{b_{j}=g} and \eqn{x_{ij}\in\{1,\ldots,g\}} for \eqn{i\in\{1,\ldots,n\}}.
}

\examples{
##########   EXAMPLE 1   ##########

# create 2-dimensional predictor (both continuous)
set.seed(123)
xmat<-cbind(runif(10^6),runif(10^6))

# Default use:
#   10 marginal bins for each predictor
#   sample 1 observation from each subregion
xind<-binsamp(xmat)

# get the corresponding knots
bknots<-xmat[xind,]

# compare to randomly-sampled knots
rknots<-xmat[sample(1:(10^6),100),]
par(mfrow=c(1,2))
plot(bknots,main="bin-sampled")
plot(rknots,main="randomly sampled")



##########   EXAMPLE 2   ##########

# create 2-dimensional predictor (continuous and nominal)
set.seed(123)
xmat<-cbind(runif(10^6),sample(1:3,10^6,replace=TRUE))

# use 10 marginal bins for x1 and 3 marginal bins for x2 
# and sample one observation from each subregion
xind<-binsamp(xmat,nmbin=c(10,3))

# get the corresponding knots
bknots<-xmat[xind,]

# compare to randomly-sampled knots
rknots<-xmat[sample(1:(10^6),30),]
par(mfrow=c(1,2))
plot(bknots,main="bin-sampled")
plot(rknots,main="randomly sampled")



##########   EXAMPLE 3   ##########

# create 3-dimensional predictor (continuous, continuous, nominal)
set.seed(123)
xmat<-cbind(runif(10^6),runif(10^6),sample(1:2,10^6,replace=TRUE))

# use 10 marginal bins for x1 and x2, and 2 marginal bins for x3 
# and sample one observation from each subregion
xind<-binsamp(xmat,nmbin=c(10,10,2))

# get the corresponding knots
bknots<-xmat[xind,]

# compare to randomly-sampled knots
rknots<-xmat[sample(1:(10^6),200),]
par(mfrow=c(2,2))
plot(bknots[1:100,1:2],main="bin-sampled, x3=1")
plot(bknots[101:200,1:2],main="bin-sampled, x3=2")
plot(rknots[rknots[,3]==1,1:2],main="randomly sampled, x3=1")
plot(rknots[rknots[,3]==2,1:2],main="randomly sampled, x3=2")

}
