binsamp <-
  function(x,xrng=NULL,nmbin=10,nsamp=1){
    # bin-sampled knots
    
    x=as.matrix(x)
    xdim=dim(x)
    if(is.null(xrng)){
      if(xdim[2]>1){xrng=apply(x,2,range)} else{xrng=matrix(range(x),2,1)}
    }
    mysamp<-function(z){sample(z,size=min(nsamp,length(z)))}
    nmbin=as.integer(nmbin)
    if(length(nmbin)!=xdim[2]){nmbin=rep(nmbin[1],xdim[2])}
    gvec=matrix(1,xdim[1],1)
    kconst=1
    for(kk in 1:xdim[2]){
      gvec = gvec + kconst*pmin(floor(nmbin[kk]*((x[,kk]-xrng[1,kk])/(xrng[2,kk]-xrng[1,kk]))),nmbin[kk]-1L)
      kconst = kconst*nmbin[kk]
    }
    gvec=as.integer(gvec)
    kidx=unlist(tapply(1:xdim[1],gvec,mysamp))
    
  }
