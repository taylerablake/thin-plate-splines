ssBasis <- 
  function(x, knots, m=2, d=0, xmin=min(x), xmax=max(x), periodic=FALSE, intercept=FALSE){
    # smoothing spline basis for polynomial splines
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: January 30, 2017
    
    ## define Bernoulli polynomials
    k1fun <- function(x) { 
      x - 1/2 
    }
    k2fun <- function(x) { 
      (k1fun(x)^2 - (1/12))/2 
    }
    k3fun <- function(x){
      (4*(k1fun(x)^3) - k1fun(x))/24
    }
    k4fun <- function(x) { 
      (k1fun(x)^4 - ((k1fun(x)^2)/2) + 7/240)/24 
    }
    k5fun <- function(x){
      ((k1fun(x)^5)/5 - (k1fun(x)^3)/6 + 7*k1fun(x)/240)/24
    }
    k6fun <- function(x) { 
      ( (k1fun(x)^6)/30 - (k1fun(x)^4)/24 + 7*(k1fun(x)^2)/480 - 31/40320 ) / 24
    }
    
    ## check penalty order 'm'
    m <- as.integer(m[1])
    if(m < 1L) stop("Input 'm' must be a positive integer between 1 and 3.")
    if(m > 3L) stop("Input 'm' must be a positive integer between 1 and 3.")
    
    ## check derivative 'd'
    d <- as.integer(d[1])
    if(d < 0L) stop("Input 'd' must be a positive integer between 0 and 2.")
    if(d > 2L) stop("Input 'd' must be a positive integer between 0 and 2.")
    
    ## check 'x' and 'knots'
    x <- as.matrix(x)
    knots <- as.matrix(knots)
    nx <- nrow(x)
    nknots <- nrow(knots)
    
    ## transform 'x' and 'knots'
    x <- (x - xmin) / (xmax - xmin)
    knots <- (knots - xmin) / (xmax - xmin)
    
    ## check 'periodic' input
    periodic <- periodic[1]
    if(!is.logical(periodic)) stop("Input 'periodic' should be a logical (TRUE/FALSE) variable.")
    
    ## linear smoothing spline
    if(m == 1L){
      
      # check derivative
      if(d == 0L){
        Xn <- matrix(1, nx, 1)
        colnames(Xn) <- "null.0"
        Xc <- matrix(k1fun(x),nx,nknots) * matrix(k1fun(knots),nx,nknots,byrow=TRUE) + k2fun( abs( matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE) ) )
        colnames(Xc) <- paste0("knot.",1:nknots)
      } else if(d == 1L){
        Xn <- matrix(0, nx, 1)
        colnames(Xn) <- "null.0"
        dmat <- matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE)
        smat <- sign(dmat)
        smat[smat==0L] <- 1L
        Xc <- matrix(k1fun(knots),nx,nknots,byrow=TRUE) + smat * k1fun(abs(dmat))
        colnames(Xc) <- paste0("knot.",1:nknots)
      } else {
        stop("Cannot set 'd=2' when 'm=1' (need d <= m).")
      } # end if(d == 0L)
      
    } # end if(m == 1L)
    
    ## cubic smoothing spline
    if(m == 2L){
      
      # check derivative
      if(d == 0L){
        Xn <- cbind(1, k1fun(x))
        colnames(Xn) <- paste0("null.",0:1)
        Xc <- matrix(k2fun(x),nx,nknots) * matrix(k2fun(knots),nx,nknots,byrow=TRUE) - k4fun( abs( matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE) ) )
        colnames(Xc) <- paste0("knot.",1:nknots)
      } else if(d == 1L){
        Xn <- matrix(c(0,1), nx, 2, byrow=TRUE)
        colnames(Xn) <- paste0("null.",0:1)
        dmat <- matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE)
        smat <- sign(dmat)
        smat[smat==0L] <- 1L
        Xc <- matrix(k1fun(x),nx,nknots) * matrix(k2fun(knots),nx,nknots,byrow=TRUE) - smat * k3fun(abs(dmat))
        colnames(Xc) <- paste0("knot.",1:nknots)
      } else {
        Xn <- matrix(0, nx, 2)
        colnames(Xn) <- paste0("null.",0:1)
        Xc <- matrix(k2fun(knots),nx,nknots,byrow=TRUE) - k2fun( abs( matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE) ) )
        colnames(Xc) <- paste0("knot.",1:nknots)
      } # end if(d == 0L)
      
    } # end if(m == 2L)
    
    ## quintic smoothing spline
    if(m == 3L) {
      
      # check derivative
      if(d == 0L){
        Xn <- cbind(1, k1fun(x), k2fun(x))
        colnames(Xn) <- paste0("null.",0:2)
        Xc <- matrix(k3fun(x),nx,nknots) * matrix(k3fun(knots),nx,nknots,byrow=TRUE) + k6fun( abs( matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE) ) )
        colnames(Xc) <- paste0("knot.",1:nknots)
      } else if(d == 1L){
        Xn <- cbind(0, 1, k1fun(x))
        colnames(Xn) <- paste0("null.",0:2)
        dmat <- matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE)
        smat <- sign(dmat)
        smat[smat==0L] <- 1L
        Xc <- matrix(k2fun(x),nx,nknots) * matrix(k3fun(knots),nx,nknots,byrow=TRUE) + smat * k5fun(abs(dmat))
        colnames(Xc) <- paste0("knot.",1:nknots)
      } else {
        Xn <- matrix(c(0, 0, 1), nx, 3, byrow=TRUE)
        colnames(Xn) <- paste0("null.",0:2)
        Xc <- matrix(k1fun(x),nx,nknots) * matrix(k3fun(knots),nx,nknots,byrow=TRUE) + k4fun( abs( matrix(x,nx,nknots) - matrix(knots,nx,nknots,byrow=TRUE) ) )
        colnames(Xc) <- paste0("knot.",1:nknots)
      } # end if(d == 0L)
      
    } # end if(m == 3L)
    
    if(intercept){
      return(list(X=cbind(Xn,Xc), knots=knots, m=m, d=d, xlim=c(xmin,xmax), 
                  periodic=periodic, intercept=intercept))
    } else {
      return(list(X=cbind(Xn,Xc)[,-1], knots=knots, m=m, d=d, xlim=c(xmin,xmax),
                  periodic=periodic, intercept=intercept))
    }
    
  
}
