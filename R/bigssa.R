bigssa <-
  function(formula,data=NULL,type=NULL,nknots=NULL,rparm=NA,
           lambdas=NULL,skip.iter=TRUE,se.fit=FALSE,rseed=1234,
           gcvopts=NULL,knotcheck=TRUE,gammas=NULL,weights=NULL,
           random=NULL,remlalg=c("FS","EM","none"),remlopts=NULL) {
    ###### Fits Smoothing Spline ANOVA models
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: December 5, 2014
    
    if(class(data)=="makessa"){
      ssafit=ssawork(formula,data)
    } else{
      ssamk=makessa(formula,data,type,nknots,rparm,
                    lambdas,skip.iter,se.fit,rseed,
                    gcvopts,knotcheck,gammas,weights,
                    random,remlalg,remlopts)
      ssafit=ssawork(formula,ssamk)
    }
    ssafit=c(ssafit,list(call=formula))
    class(ssafit) <- "bigssa"
    return(ssafit)
    
  }