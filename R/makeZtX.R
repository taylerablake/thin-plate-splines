makeZtX <-
  function(lev1var,lev2var,Ku,Ju,uidx,l1names=NULL){
    
    v2lev <- levels(lev2var)
    v2nlev <- length(v2lev)
    nKcol <- ncol(Ku)
    nJcol <- ncol(Ju)
    if(length(lev1var)==1L){
      # random intercept model
      ZtK <- matrix(0,v2nlev,ncol(Ku))
      ZtJ <- matrix(0,v2nlev,ncol(Ju))
      for(k in 1:v2nlev){
        idx <- which(lev2var==v2lev[k])
        lidx <- length(idx)
        ZtK[k,] <- colSums(matrix(Ku[uidx[idx],],nrow=lidx,ncol=nKcol))
        ZtJ[k,] <- colSums(matrix(Ju[uidx[idx],],nrow=lidx,ncol=nJcol))
      }
    } else {
      # subjects nested in groups
      v1lev <- l1names
      v1nlev <- length(v1lev)
      ZtK <- matrix(0,v1nlev,ncol(Ku))
      ZtJ <- matrix(0,v1nlev,ncol(Ju))
      for(k in 1:v1nlev){
        idx <- which(lev1var==v1lev[k])
        lidx <- length(idx)
        ZtK[k,] <- colSums(matrix(Ku[uidx[idx],],nrow=lidx,ncol=nKcol))
        ZtJ[k,] <- colSums(matrix(Ju[uidx[idx],],nrow=lidx,ncol=nJcol))
      }
    }
    
    return(cbind(ZtK,ZtJ))
    
  }