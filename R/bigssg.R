bigssg <-
  function(formula,family,data=NULL,type=NULL,nknots=NULL,rparm=NA,
           lambdas=NULL,skip.iter=TRUE,se.lp=FALSE,rseed=1234,
           gcvopts=NULL,knotcheck=TRUE,gammas=NULL,weights=NULL) {
    ###### Fits Generalized Smoothing Spline ANOVA models
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: August 26, 2014
    
    if(class(data)=="makessg"){
      ssgfit=ssgwork(formula,data)
    } else{
      ssgmk=makessg(formula,family,data,type,nknots,rparm,
                    lambdas,skip.iter,se.lp,rseed,gcvopts,
                    knotcheck,gammas,weights)
      ssgfit=ssgwork(formula,ssgmk)
    }
    ssgfit=c(ssgfit,list(call=formula))
    class(ssgfit) <- "ssg"
    return(ssgfit)
    
  }