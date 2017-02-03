predict.ordspline <-
  function(object, newdata=NULL, se.fit=FALSE, ...) {
    
    ### check if 'newdata' is provided
    if(is.null(newdata)) {
      if(se.fit){
        return(list(fit=object$fitted.values, se.fit=object$se.fit))
      } else {
        return(object$fitted.values)
      }
    }
    
    ### recode 'x' as integers {1, ..., length(unique(x))}
    newdata <- match(newdata, object$xunique)
    if(any(is.na(newdata))) stop("Input 'newdata' contains values not found in 'x'.")
    
    ### get number of new data points and knots
    nx <- length(newdata)
    nk <- length(object$knots)
    nlev <- length(object$xunique)
    
    ### make design matrix
    if(object$monotone){
      X <- cbind(1, ((.Fortran("ordkermon", newdata, object$knots, nx, nk, matrix(0, nx, nk), PACKAGE = "bigsplines"))[[5]])[,-nk])
    } else {
      X <- cbind(1, (.Fortran("ordker", newdata, object$knots, nx, nk, nlev, matrix(0, nx, nk), PACKAGE = "bigsplines"))[[6]])
    }
    
    ### output predictions
    fit <- X %*% object$coef
    if(se.fit){
      se.fit <- sqrt(rowSums((X %*% object$coef.csqrt)^2))
      return(list(fit=fit, se.fit=se.fit))
    } else {
      return(fit)
    }
    
}
