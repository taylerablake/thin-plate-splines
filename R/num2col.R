num2col <- 
  function(x, col=heat.colors(12), xlim=range(x)){
    ncol <- length(col)
    q <- seq(xlim[1], xlim[2], length.out=(ncol+1L))
    xq <- outer(x,q,FUN=">=")
    nx <- rowSums(xq)
    if(any(is.na(nx))) nx[which(is.na(nx))] <- ncol
    if(any(nx>ncol)) nx[nx>ncol] <- ncol
    col[nx]
}