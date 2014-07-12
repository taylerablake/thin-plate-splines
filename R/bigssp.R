bigssp <-
  function(formula,data=NULL,type=NULL,nknots=NULL,rparm=NA,
           lambdas=NULL,skip.iter=TRUE,se.fit=FALSE,rseed=1234,
           gcvopts=NULL,knotcheck=TRUE,thetas=NULL) {
    ###### Fits Smoothing Splines with Parametric effects
    ###### Nathaniel E. Helwig (nhelwig2@illinois.edu)
    ###### Last modified: July 11, 2014
    
    if(class(data)=="makessp"){
      sspfit=sspwork(formula,data)
    } else{
      sspmk=makessp(formula,data,type,nknots,rparm,
                    lambdas,skip.iter,se.fit,rseed,
                    gcvopts,knotcheck,thetas)
      sspfit=sspwork(formula,sspmk)
    }
    sspfit=c(sspfit,list(call=formula))
    class(sspfit) <- "ssp"
    return(sspfit)
    
  }