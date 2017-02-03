print.ordspline <- 
  function(x,...){
    
    cat("\nMonotonic:\n")
    cat(x$monotone,"\n")
    cat("\nFit Statistics:\n")
    print(x$info)
    cat("\nSmoothing Parameter:\n")
    cat(x$lambda,"\n ")
    
  }