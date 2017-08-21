fittps <-  function(xgrid,y,D2,nknots=NULL,nvec=NULL,rparm=NA,
           alpha=1,lambdas=NULL,se.fit=FALSE,
           rseed=1234,knotcheck=TRUE){
      
    ### initial info
    if(is.null(rseed)==FALSE){set.seed(rseed)}
    x <- as.matrix(xgrid+0.0)
    n <- nrow(x)
    nx <- ncol(x)
    N <- nrow(y)
    M <- ncol(y)
    regressors <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
    no.skip <- 0
    for (t in 2:M){
          regressors[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
          no.skip <- no.skip + t - 1
    }
    y <- matrix(as.vector(t(y[,-1])),ncol=1,nrow=N*(M-1))
    
    yty <- sum((y*rep(diag(sqrt(D2))[-1],N))^2)
    ysm <- sum(y)
    if(nrow(y)!=n){stop("Lengths of 'x' and 'y' must match.")}
    if(ncol(y)>1){stop("Response must be unidimensional (vector).")}
    if(nx>3){stop("Too many predictors. Use another function (bigssa or bigssp).")}

    ### check knots
    
    if(is.null(nknots)){
      if(nx==1L){nknots <- 30L} else if(nx==2L){nknots <- 100L} else {nknots <- 200L}
    }
    if(!is.null(rseed) && length(nknots)==1L){set.seed(rseed)}
    kidx <- NA
    theknots <- x 
    nknots <- nrow(theknots)

    nunewr <- n
    w <- diag(rep(diag(D2)[-1],N))
    if(nknots<1){stop("Input 'nknots' must be positive integer.")}
    
    if(knotcheck){
      theknots <- unique(theknots)
      nknots <- nrow(theknots)
    }
    if(nknots<nx+2){stop("Input 'nknots' is too small. Need nknots-2 >= ncol(x).")}

    ### create penalty matrix
    if(nx<3){gconst <- 1} else{gconst <- (-1)}
    mqed <- gconst*(.Fortran("tpskersym",theknots,nknots,nx,matrix(0,nknots,nknots)))[[4]]
    sqed <- gconst*(.Fortran("tpsker",x,theknots,nunewr,nx,nknots,matrix(0,nunewr,nknots)))[[6]]

      # tranforms problem to unconstrained
      Cqrd <- qr(cbind(1,theknots),LAPACK=TRUE)
      CqrQ <- qr.Q(Cqrd,complete=TRUE)
      CZ <- CqrQ[,(nx+2):nknots]
      rm(Cqrd,CqrQ)
      Qmat <- crossprod(CZ,mqed%*%CZ)
      # make design matrix and crossproduct matrix
      wsqrt <- sqrt(w)
      Kmat <- cbind(1,x)
      Jmat <- sqed%*%CZ
      rm(sqed)
      KtK <- crossprod(wsqrt %*% regressors %*% Kmat)
      KtJ <- crossprod(wsqrt %*% regressors %*% Kmat,
                       wsqrt %*% regressors %*% Jmat)
      JtJ <- crossprod(wsqrt %*% regressors %*% Jmat)
      Kty <- crossprod(wsqrt %*% regressors %*% Kmat,wsqrt %*% y)
      Jty <- crossprod(wsqrt %*% regressors %*% Jmat,wsqrt %*% y)
     
    # # ### find optimal smoothing parameter 
    nbf <- nx+1L
    ncdim <- nrow(Qmat)
    if(is.null(lambdas)){
      nqmat <- matrix(0,nbf+ncdim,nbf+ncdim)
      nqmat[(nbf+1):(ncdim+nbf),(nbf+1):(ncdim+nbf)] <- n*Qmat
      newlam <- lamloop(10^-seq(9,3,length.out=30),1,Kty,Jty,KtK,KtJ,JtJ,
                        Qmat,ncdim,n,alpha,yty,nbf)
      gcvopt <- nlm(f=gcvcss,p=log(newlam),yty=yty,
                    xtx=rbind(cbind(KtK,KtJ),cbind(t(KtJ),JtJ)),
                    xty=rbind(Kty,Jty),nqmat=nqmat,ndpts=n,alpha=alpha)
      lambda <- exp(gcvopt$est)
    } else if(length(lambdas)>1){
      if(any(lambdas<0)){stop("Input 'lambdas' must be nonnegative.")}
      lambda <- lamloop(lambdas,1,Kty,Jty,KtK,KtJ,JtJ,
                        Qmat,ncdim,n,alpha,yty,nbf)
    } else {
      lambda <- lambdas[1]
      if(lambda<0){stop("Input 'lambdas' must be nonnegative.")}
    }

    # # # ### get final estimates
    fxhat <- lamcoef(lambda,1,Kty,Jty,KtK,KtJ,JtJ,Qmat,ncdim,n,alpha,yty,nbf)
    fhat <- fxhat[[1]]
    dchat <- fhat[1:(nbf+ncdim)]
    yhat <- regressors %*% cbind(Kmat,Jmat) %*% dchat
    sseval <- fhat[nbf+ncdim+1]
    effdf <- fhat[nbf+ncdim+2]
    mevar <- sseval/(n-effdf)
    gcv <- n*sseval/((n-alpha*effdf)^2)
    aic <- n*(1+log(2*pi)) + n*log(sseval/n) + effdf*2
    bic <- n*(1+log(2*pi)) + n*log(sseval/n) + effdf*log(n)
    csqrt <- sqrt(mevar)*fxhat[[2]]
    
    # # # ### posterior variance
    pse <- NA
    if(se.fit){pse <- sqrt(postvar(Kmat,Jmat,csqrt))}

    #  dchat <- c(dchat[1:(nx+1)],CZ%*%dchat[(nx+2):(ncdim+nbf)])
    #  csqrt <- rbind(csqrt[1:(nx+1),],CZ%*%csqrt[(nx+2):(ncdim+nbf),])

    # # ### calculate vaf
    mval <- ysm/n
    vaf <- 1 - sseval/(yty-n*(mval^2))

    ndf <- data.frame(n=n,df=effdf,row.names="")
    tpsfit <- list(fitted.values=yhat,
                   se.fit=pse,x=x,y=y,
                   funique=w,sigma=sqrt(mevar),
                   ndf=ndf,
                   info=c(gcv=gcv,rsq=vaf,aic=aic,bic=bic),
                   myknots=theknots,
                   nvec=nvec,rparm=rparm,
                   lambda=lambda,
                   coef=dchat,
                   Kmat=Kmat,
                   Jmat=Jmat,
                   fxhat=fxhat)
    tpsfit

}