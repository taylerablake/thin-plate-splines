ordspline <-
  function(x, y, knots, weights, lambda, monotone=FALSE){
    ###### Ordinal Smoothing Spline
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: January 31, 2017
    
    ### check 'x' and 'y'
    nx <- length(x)
    if(length(y) != nx) stop("Lengths of 'x' and 'y' must match.")
    y <- as.numeric(y)
    
    ### recode 'x' as integers {1, ..., length(unique(x))}
    xuni <- sort(unique(x))
    nlev <- length(xuni)
    xo <- x
    x <- match(x, xuni)
    
    ### check 'knots'
    if(missing(knots)){
      nk <- min(nlev, 50L)
      knots <- unique(as.integer(seq(1, nlev, length.out=nk)))
    } else if(length(knots)==1L){
      knots <- as.integer(knots)
      if(knots < 2L) stop("Need at least 2 knots.")
      nk <- knots
      knots <- unique(as.integer(seq(1, nlev, length.out=nk)))
    } else {
      nk <- length(knots)
      if(nk > nx) stop("Input too many knots: length(knots) > length(x).")
      iknots <- match(knots, xuni)
      if(any(is.na(iknots))) stop("Input 'knots' contains values not found in 'x'.")
      knots <- iknots
    }
    
    ### check 'weights'
    if(missing(weights)){
      weights <- 1
    } else {
      if(length(weights) != nx) stop("Lengths of 'x' and 'weights' must match.")
      if(any(weights<0)) stop("Input 'weights' must be non-negative.")
    }
    wsqrt <- sqrt(weights)
    
    ### check 'lambda'
    if(missing(lambda)){
      fixlambda <- FALSE
      lambda <- 1e-5
    } else {
      lambda <- as.numeric(lambda[1])
      if(lambda < 0) stop("Input 'lambda' must be a non-negative scalar.")
      fixlambda <- TRUE
    }
    
    ### check 'monotone'
    monotone <- monotone[1]
    if(!is.logical(monotone)) stop("Input 'monotone' must be logical (TRUE/FALSE).")
    
    ### make RK matrix
    if(monotone){
      X <- cbind(1, ((.Fortran("ordkermon", x, knots, nx, nk, matrix(0, nx, nk), PACKAGE = "bigsplines"))[[5]])[,-nk])
      Q <- diag(nk)
      Q[1,1] <- 0
    } else {
      X <- cbind(1, (.Fortran("ordker", x, knots, nx, nk, nlev, matrix(0, nx, nk), PACKAGE = "bigsplines"))[[6]])
      Q <- cbind(0, rbind(0, (.Fortran("ordkersym", knots, nk, nlev, matrix(0, nk, nk), PACKAGE = "bigsplines"))[[4]]))
    }
    
    ### make crossproduct matrix
    XtX <- crossprod(wsqrt * X)
    Xty <- crossprod(X, weights * y)
    
    ### find optimal lambda
    if(!fixlambda){
      gcvopt <- optimize(f=gcvoss, interval=c(0,1), yty=sum(y^2),
                         xtx=XtX, xty=Xty, qmat=Q, ndpts=nx, monotone=monotone)
      lambda <- gcvopt$minimum
    }
    
    ### get covariance matrix
    cpmat <- XtX + nx * lambda * Q
    ceig <- eigen(cpmat, symmetric=TRUE)
    nze <- sum(ceig$val > (.Machine$double.eps * max(ceig$val)))
    covsqrt <- ceig$vec[,1:nze] %*% diag(1/sqrt(ceig$val[1:nze]))
    cvar <- tcrossprod(covsqrt)
    
    ### check monotone input
    if(monotone){
      
      # define inequality constraint matrix
      Gmat <- diag(nk-1)
      Gmat <- cbind(0, Gmat)
      
      # get inequality-constrained coefficients
      qpfit <- solve.QP(Dmat=cpmat, dvec=Xty, Amat=t(Gmat))
      coefs <- qpfit$solution
      
      # notes on organization of solve.QP output: 
      #  (a) qpfit$solution == qpfit$unconstrained.solution + solve(cpmat) %*% crossprod(Gmat, qpfit$Lagrangian)
      #  (b) qpfit$solution == (diag(ncol(cpmat)) + solve(cpmat) %*% crossprod(Gmat[lix,], M2) %*% Gmat) %*% qpfit$unconstrained.solution
      #  (c) qpfit$Lagrangian == MPinv(t(Gmat)) %*% cpmat %*% (qpfit$solution - qpfit$unconstrained.solution)
      #  (d) sse == sum(y^2) + 2*qpfit$value - nx*lambda*sum(qpfit$solution[-1]^2)
      
      # get covariance matrix information (C.K. Liew, 1976 - JASA)
      Wmat <- Gmat %*% cvar %*% t(Gmat)
      lam <- qpfit$Lagrangian
      lix <- which(lam > 0L)
      lenlix <- length(lix)
      if(lenlix == 0L){
        CM <- diag(nk)
      } else if(lenlix == 1L){
        M1M2 <- MPinv(cbind(diag(nk-1)[,-lix],-Wmat[,lix]))
        M2 <- M1M2[(nk-lenlix):(nk-1),]
        CM <- diag(nk) +  cvar %*% outer(Gmat[lix,],M2) %*% Gmat
      } else {
        M1M2 <- MPinv(cbind(diag(nk-1)[,-lix],-Wmat[,lix]))
        M2 <- M1M2[(nk-lenlix):(nk-1),]
        CM <- diag(nk) +  cvar %*% crossprod(Gmat[lix,], M2) %*% Gmat
      }
      covsqrt <- CM %*% covsqrt
      edf <- sum(diag(CM %*% cvar %*% XtX))
      
    } else {
      
      coefs <- cvar %*% Xty
      edf <- sum(diag(XtX %*% cvar))
      
    }
    
    # get fit information
    fit <- X %*% coefs
    sse <- sum( (y - fit)^2 )
    GCV <- (sse / nx) / (1 - edf/nx)^2
    aic <- nx * (1 + log(2 * pi)) + nx * log(sse/nx) + edf * 2
    bic <- nx * (1 + log(2 * pi)) + nx * log(sse/nx) + edf * log(nx)
    Rsq <- 1 - sse / sum((y-mean(y))^2)
    sigma <- sqrt(sse / (nx - edf))
    covsqrt <- sigma * covsqrt
    se.fit <- sqrt(rowSums((X %*% covsqrt )^2))
    
    # collect results
    ossfit <- list(fitted.values=fit, se.fit=se.fit, sigma=sigma,
                   lambda=lambda, info=c(gcv=GCV,rsq=Rsq,aic=aic,bic=bic), 
                   coef=coefs, coef.csqrt=covsqrt, n=nx, df=edf, xunique=xuni,  
                   x=xo, y=y, residuals=y-fit, knots=knots, monotone=monotone)
    class(ossfit) <- "ordspline"
    ossfit
    
  }