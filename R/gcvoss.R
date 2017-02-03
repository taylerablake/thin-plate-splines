gcvoss <-
  function(lambda,yty,xtx,xty,qmat,ndpts,monotone) {
  
  gcv <- tryCatch({
    cpmat <- xtx + ndpts * lambda * qmat
    chi <- pinvsm(cpmat)
    if(monotone){
      nk <- length(xty)
      Gmat <- diag(nk-1)
      Gmat <- cbind(0, Gmat)
      qpfit <- solve.QP(Dmat = cpmat, dvec = xty, Amat = t(Gmat))
      Wmat <- Gmat %*% chi %*% t(Gmat)
      lam <- qpfit$Lagrangian
      lix <- which(lam > 0L)
      lenlix <- length(lix)
      if(lenlix == 0L){
        CM <- diag(nk)
      } else if(lenlix == 1L){
        M1M2 <- MPinv(cbind(diag(nk-1)[,-lix],-Wmat[,lix]))
        M2 <- M1M2[(nk-lenlix):(nk-1),]
        CM <- diag(nk) +  chi %*% outer(Gmat[lix,],M2) %*% Gmat
      } else {
        M1M2 <- MPinv(cbind(diag(nk-1)[,-lix],-Wmat[,lix]))
        M2 <- M1M2[(nk-lenlix):(nk-1),]
        CM <- diag(nk) +  chi %*% crossprod(Gmat[lix,], M2) %*% Gmat
      }
      mse <- (yty + 2*qpfit$value - ndpts*lambda*sum(qpfit$solution[-1]^2)) / ndpts
      gcv <- mse / ( (1 - sum(diag(CM%*%chi%*%xtx))/ndpts )^2 )
    } else {
      parta <- chi %*% xty
      mse <- (yty - 2*crossprod(xty,parta) + crossprod(parta,xtx%*%parta)) / ndpts
      gcv <- mse / ( (1 - sum(diag(chi%*%xtx))/ndpts )^2 )
    }
    gcv
  }, error = function(e) yty)
  
}
