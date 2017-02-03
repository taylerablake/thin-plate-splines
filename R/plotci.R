plotci <- 
  function(x, y, se, level=0.95, cval=NULL, col="blue",
           col.ci="cyan", alpha=0.65, add=FALSE, 
           type="l", link=function(y){y}, axes=TRUE, 
           bars=FALSE, barlty=1, barlwd=2, bw=0.2, ...){
  
  if (any(se < 0)) {
    stop("Input 'se' must contain nonnegative standard errors.")
  }
  
  if (level <= 0 || level >= 1) {
    stop("Input 'level' must statisfy 0 < level < 1.")
  }
  if (is.null(cval)) cval <- qnorm(1 - (1 - level)/2)
  
  if (add) {
    lines(x, link(y), col=col, type=type, ...)
  }
  else {
    plot(x, link(y), col=col, type=type, axes=axes, ...)
  }
  
  if(bars){
    if(length(barlty)!=length(x)) barlty <- rep(barlty[1],length(x))
    if(length(barlwd)!=length(x)) barlwd <- rep(barlwd[1],length(x))
    if(length(col.ci)!=length(x)) col.ci <- rep(col.ci[1],length(x))
    for(k in 1:length(x)){
      lines(rep(x[k],2), link(c(y[k]-cval*se[k], y[k]+cval*se[k])), lty=barlty[k], lwd=barlwd[k], col=col.ci[k])
      lines(c(x[k]-bw, x[k]+bw), rep(link(y[k]-cval*se[k]),2), lty=barlty[k], lwd=barlwd[k], col=col.ci[k])
      lines(c(x[k]-bw, x[k]+bw), rep(link(y[k]+cval*se[k]),2), lty=barlty[k], lwd=barlwd[k], col=col.ci[k])
    }
  } else {
    myrgb <- col2rgb(col.ci)/255
    polycoords <- rbind(cbind(x, link(y + cval*se)), 
                        cbind(rev(x), link(rev(y - cval*se))))
    polygon(polycoords[, 1], polycoords[, 2], 
            col = rgb(myrgb[1], myrgb[2], myrgb[3], alpha), border = NA)
    lines(x, link(y), col=col, ...)
  }
  
}