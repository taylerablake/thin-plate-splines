lamloopg <-
  function(lambdas,gamvec,family,Kmat,Jmats,yvar,Qmats,
           nknots,ndpts,alpha,yty,nbf,fweights,weights,
           maxit,intol,subsamp,dispersion) {
    
    # initialize things
    lamgcv=vector("numeric",length(lambdas))
    nqmat=matrix(0,nbf+nknots,nbf+nknots)
    thmat=kronecker(gamvec,diag(nknots))
    nqmat[(nbf+1):(nknots+nbf),(nbf+1):(nknots+nbf)]=ndpts*(Qmats%*%thmat)
    Jmat=Jmats%*%thmat
    rm(Jmats,Qmats)
    nunewr=length(yvar)
    if(nunewr>subsamp){idx=sample.int(nunewr,subsamp)} else {idx=1:nunewr}
    if(family=="binomial" & length(weights)>1){gcvndpts=sum(weights*fweights)} else {gcvndpts=ndpts}
    
    # loop through different lambdas
    for(jj in 1:length(lambdas)){
      lamgcv[jj]=tryCatch({
        
        # get starting values
        if(family=="binomial"){
          mu0=(yvar+0.5)/(fweights*weights+1)
          eta0=log(mu0/(1-mu0))
          vwts=as.numeric(mu0*(1-mu0)*weights)
          zvar=eta0*fweights+(yvar-fweights*mu0*weights)/vwts
        } else if(family=="poisson"){
          mu0=(yvar+0.1)/fweights
          eta0=log(mu0)
          vwts=as.numeric(mu0*weights)
          zvar=eta0*fweights+(yvar-fweights*mu0)/vwts
        } else if(family=="Gamma"){
          mu0=(yvar+0.1)/fweights
          eta0=1/mu0
          vwts=as.numeric((mu0^2)*weights)
          zvar=eta0*fweights-(yvar-fweights*mu0)/vwts
        } else if(family=="inverse.gaussian"){
          mu0=(yvar+0.1)/fweights
          eta0=1/(2*(mu0^2))
          vwts=as.numeric((mu0^3)*weights)
          zvar=eta0*fweights-(yvar-fweights*mu0)/vwts
        } else if(family=="negbin"){
          mu0=(yvar+0.1)/fweights
          size=1/dispersion
          eta0=log(mu0)
          vwts=as.numeric(((mu0*size)/(mu0+size))*weights)
          zvar=eta0*fweights-1+yvar/(fweights*mu0)
        }
        
        # iterative reweighted least squares udpate
        fitchng=1;   iters=0L;   yold=yty
        while(fitchng>intol & iters<maxit){
          
          # get crossproduct matrices
          KtJ=crossprod(Kmat*(vwts*fweights),Jmat)
          KtK=crossprod(Kmat*(vwts*fweights),Kmat)
          JtJ=crossprod(Jmat*(vwts*fweights),Jmat)
          Kty=crossprod(Kmat*vwts,zvar)
          Jty=crossprod(Jmat*vwts,zvar)
          xty=c(Kty,Jty)
          xtx=rbind(cbind(KtK,KtJ),cbind(t(KtJ),JtJ))
          
          # update coefficients
          ceig=eigen(xtx+lambdas[jj]*nqmat,symmetric=TRUE)
          nze=sum(ceig$val>ceig$val[1]*.Machine$double.eps)
          isqrts=ceig$vec[,1:nze]%*%diag(ceig$val[1:nze]^-0.5)
          chi=tcrossprod(isqrts)
          bhat=chi%*%xty
          yhat=cbind(Kmat,Jmat)%*%bhat
          
          # check solution and update iteration
          if(family=="binomial"){
            mu0=exp(yhat)/(1+exp(yhat))
            mu0[mu0<=0]=10^-4
            mu0[mu0>=1]=1-10^-4
            vwts=as.numeric(mu0*(1-mu0)*weights)
            zvar=yhat*fweights+(yvar-fweights*mu0*weights)/vwts
          } else if(family=="poisson"){
            mu0=exp(yhat)
            mu0[mu0<=0]=10^-4
            vwts=as.numeric(mu0*weights)
            zvar=yhat*fweights+(yvar-fweights*mu0)/vwts
          } else if (family=="Gamma"){
            yhat[yhat<=0]=10^-4
            mu0=1/yhat
            vwts=as.numeric((mu0^2)*weights)
            zvar=yhat*fweights-(yvar-fweights*mu0)/vwts
          } else if(family=="inverse.gaussian"){
            yhat[yhat<=0]=10^-4
            mu0=sqrt(1/(2*yhat))
            vwts=as.numeric((mu0^3)*weights)
            zvar=yhat*fweights-(yvar-fweights*mu0)/vwts
          } else if(family=="negbin"){
            mu0=exp(yhat)
            mu0[mu0<=0]=10^-4
            vwts=as.numeric(((mu0*size)/(mu0+size))*weights)
            zvar=yhat*fweights-1+yvar/(fweights*mu0)
          }
          fitchng=crossprod(yold[idx]-yhat[idx])/crossprod(yold[idx])
          iters=iters+1L;   yold=yhat
          
        }
                
        # update solution (using direct GCV -- Gu and Xiang, 2001)
        trval=sum(diag(chi%*%xtx))
        if(family=="binomial"){
          cbeta=sum(log(1+exp(yhat))*fweights*weights)
          p1gcv=crossprod(yvar,yhat)-cbeta
          p2gcv=(trval/(gcvndpts-trval))*sum((weights^2/vwts)*fweights)/gcvndpts
          p3gcv=crossprod(yvar,1-mu0)
        } else if(family=="poisson"){
          cbeta=sum(exp(yhat)*fweights)
          p1gcv=crossprod(yvar,yhat)-cbeta
          p2gcv=(trval/(ndpts-trval))*sum((1/vwts)*fweights)/ndpts
          p3gcv=sum(yty)-crossprod(yvar,mu0)
        } else if(family=="Gamma"){
          cbeta=sum(log(mu0)*fweights)
          p1gcv=crossprod(yvar,-yhat)-cbeta
          p2gcv=(trval/(ndpts-trval))*sum((1/vwts)*fweights)/ndpts
          p3gcv=sum(yty)-crossprod(yvar,mu0)
        } else if(family=="inverse.gaussian"){
          cbeta=-sum((1/mu0)*fweights)
          p1gcv=crossprod(yvar,-yhat)-cbeta
          p2gcv=(trval/(ndpts-trval))*sum((1/vwts)*fweights)/ndpts
          p3gcv=sum(yty)-crossprod(yvar,mu0)
        } else if(family=="negbin"){
          cbeta=-sum(fweights*size*log(size/mu0))
          p0=size/(size+mu0)
          p1gcv=sum((yvar+fweights*size)*log(1-p0))-cbeta
          p2gcv=(trval/(ndpts-trval))*sum((1/vwts)*fweights)/ndpts
          p3gcv=sum((yty+yvar*size)*(p0^2))-sum(yvar*size*p0)
        }
        
        (1/gcvndpts)*(alpha*p2gcv*p3gcv-p1gcv)
        
      }, error = function(e) sum(yty))
    }
    newlam=lambdas[which.min(lamgcv)]
    
  }
