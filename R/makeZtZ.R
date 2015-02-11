makeZtZ <-
  function(rs,lev1var,lev2var,yvar){
    
    v2lev=levels(lev2var)
    v2nlev=length(v2lev)
    if(rs[1]=="1"){
      # random intercept model
      ysp=split(yvar,lev2var)
      ZtZ=unlist(lapply(ysp,length))
      Zty=unlist(lapply(ysp,sum))
      nspg=length(Zty)
    } else {
      # subjects nested in groups
      v1lev=levels(lev1var)
      v1nlev=length(v1lev)
      ysp=split(yvar,list(lev1var,lev2var))
      zt=matrix(unlist(lapply(ysp,length)),v1nlev,v2nlev)
      zt[zt==0]=NA
      nspg=colSums(!is.na(zt))
      rszt=rowSums(zt)
      if(any(is.na(rszt)==FALSE)){stop(paste("Group 2 variable",rs[1],"must be nested within Group 1 variable",rs[2],"."))}
      ZtZ=zt[!is.na(c(zt))]
      zs=matrix(unlist(lapply(ysp,sum)),v1nlev,v2nlev)
      zname=rep(v1lev,v2nlev)
      idx=which(zs==0)
      Zty=zs[-idx]
      names(Zty)=zname[-idx]
    }
    
    return(list(Zty,ZtZ,nspg))
    
  }