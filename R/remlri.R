remlri <-
  function(yty,Xty,Zty,XtX,ZtZ,XtZ,ndf,tau=1,imx=100,tol=10^-5,alg=c("FS","EM")){
    ###### REML Estimation of Random Intercept (via Fisher Scoring or EM)
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: March 10, 2015
    
    ### Inputs:
    # yty: crossprod(y)
    # Xty: crossprod(X,y)
    # Zty: crossprod(Z,y)
    # XtX: crossprod(X)
    # ZtZ: diag(crossprod(Z))
    # XtZ: crossprod(X,Z)
    # ndf: nrow(X) - ncol(X)
    # tau: initial estimate of variance component
    # imx: maximum number of iterations
    # tol: convergence tolerance
    # alg: estimation algorithm
    
    ### Notes:
    #   y = X%*%alpha + Z%*%beta + e   where 
    #
    #   y is n by 1 response vector
    #   X is n by p fixed effect design matrix
    #     note: df is typically ncol(X)
    #   alpha is p by 1 fixed effect coefficients
    #   Z is n by q random effect design matrix 
    #   beta is q by 1 random effect coefficients
    #     note: beta ~ N(0,sig*tau*diag(q))
    #   e is n by 1 residual vector
    #     note: e ~ N(0,sig*diag(n))
    
    ### initial info
    XtXi <- pinvsm(XtX)
    XtXiZ <- XtXi%*%XtZ
    ZtZX <- (-1)*crossprod(XtZ,XtXiZ)
    diag(ZtZX) <- diag(ZtZX)+ZtZ
    ZtyX <- Zty-crossprod(XtZ,XtXi)%*%Xty
    nz <- length(Zty)
    Deig <- eigen(ZtZX,symmetric=TRUE)
    
    ### iterative update
    alg <- alg[1]
    vtol <- 1
    iter <- 0
    tauold <- tau
    if(alg=="FS"){
      # Fisher scoring
      while(vtol>tol && iter<imx) {
        # update fixed and random effects
        newval <- 1/(Deig$val+1/tau)
        Dmat <- Deig$vec%*%tcrossprod(diag(newval),Deig$vec)
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty-Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
        # update score and information
        trval <- sum(newval)
        gg <- (crossprod(beta)/((tau^2)*sig)) - ((nz/tau)-(trval/(tau^2)))
        hh <- (nz/(tau^2)) - 2*(trval/(tau^3)) + sum(newval^2)/(tau^4)
        # update parameter and check for convergence
        tau <- tau+gg/hh
        vtol <- abs((tau-tauold)/tau)
        tauold <- tau
        iter <- iter+1
      }
    } else {
      # Expectation Maximization
      if(is.matrix(ZtZ)){
        WtW <- rbind(cbind(XtX,XtZ),cbind(t(XtZ),ZtZ))
      } else {
        WtW <- rbind(cbind(XtX,XtZ),cbind(t(XtZ),diag(ZtZ)))
      }
      nx <- nrow(XtX)
      sig <- 0
      while(vtol>tol && iter<imx) {
        # update fixed and random effects
        newval <- 1/(Deig$val+1/tau)
        Dmat <- Deig$vec%*%tcrossprod(diag(newval),Deig$vec)
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty-Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        abhat <- c(alpha,beta)
        sig <- (yty - 2*t(c(Xty,Zty))%*%abhat + crossprod(abhat,WtW%*%abhat) + sig*(nx+sum(diag(Dmat%*%ZtZX))))/ndf
        # update parameters and check for convergence
        tau <- (mean(beta^2)/sig) + sum(newval)/nz
        vtol <- abs((tau-tauold)/tau)
        tauold <- tau
        iter <- iter+1
      }
    } # end if(alg=="FS")
    
    list(tau=as.numeric(tau),sig=as.numeric(sig),
         iter=iter,cnvg=as.logical(ifelse(vtol>tol,FALSE,TRUE)))
    
  }