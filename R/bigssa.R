bigssa <-
  function(formula,data=NULL,type=NULL,nknots=NULL,rparm=NA,
           lambdas=NULL,skip.iter=TRUE,se.fit=FALSE,rseed=1234,
           gcvopts=NULL,knotcheck=TRUE,gammas=NULL) {
    ###### Fits Smoothing Spline ANOVA models
    ###### Nathaniel E. Helwig (nhelwig2@illinois.edu)
    ###### Last modified: July 11, 2014
    
    if(class(data)=="makessa"){
      ssafit=ssawork(formula,data)
    } else{
      ssamk=makessa(formula,data,type,nknots,rparm,
                    lambdas,skip.iter,se.fit,rseed,
                    gcvopts,knotcheck,gammas)
      ssafit=ssawork(formula,ssamk)
    }
    ssafit=c(ssafit,list(call=formula))
    class(ssafit) <- "ssa"
    return(ssafit)
    
  }