print.summary.tps <- 
  function(x,digits=4,...){
    
    ehat=round(quantile(x$residuals),digits)
    names(ehat)=c("Min","1Q","Median","3Q","Max")
    x$info=round(x$info,digits)
    cat("\nSpline Type:\n")
    cat(x$type,"\n")
    cat("\nResiduals:\n")
    print(ehat)
    cat("\nError Std Dev Estimate:\n")
    cat(x$sigma," on ",as.numeric(x$n-x$df)," degrees of freedom","\n")
    cat("\nFit Statistics:")
    cat("\nGCV:  ",x$info[1])
    cat("\nR^2:  ",x$info[2])
    cat("\nAIC:  ",x$info[3])
    cat("\nBIC:  ",x$info[4],"\n ")
    cat("\nSmoothing Parameter:\n")
    cat(x$lambda,"\n ")
    
  }