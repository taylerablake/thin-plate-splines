## Fit Multiple Smoothing Parameter Varying Coefficient (Gaussian) REGression


mspvcreg1 <- function(s,r,id.basis,y,wt,method,alpha,varht,random,skip.iter,w,d2)
{
  qr.trace <- FALSE
  
  random <- NULL
  if ((alpha<0)&(method%in%c("u","v"))) qr.trace <- TRUE
  alpha <- abs(alpha)
  ## get dimensions
  nobs <- nrow(r)
  nxi <- ncol(r)
  if (!is.null(s)) {
    if (is.vector(s)) nnull <- 1
    else nnull <- ncol(s)
  }
  else nnull <- 0
  if (!is.null(random)) nz <-ncol(as.matrix(random$z))
  else nz <- 0
  nxiz <- nxi + nz
  nn <- nxiz + nnull
  nq <- dim(r)[3]
  ## cv function
  cv <- function(theta) {
    ind.wk <- theta[1:nq]!=theta.old
    if (sum(ind.wk)==nq) {
      r.wk0 <- 0
      for (i in 1:nq) {
        r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
      }
      assign("r.wk",r.wk0+0,inherits=TRUE)
      assign("theta.old",theta[1:nq]+0,inherits=TRUE)
    }
    else {
      r.wk0 <- r.wk
      for (i in (1:nq)[ind.wk]) {
        theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
        r.wk0 <- r.wk0 + theta.wk*r[,,i]
      }
    }
    qq.wk <- r.wk0[id.basis,]
    q.wk <- 10^nlambda*qq.wk

    if (!is.null(wt)) {
      y.wk <- wt*y
      s.wk <- wt*s
      r.wk0 <- wt*r.wk0
    }
    if (qr.trace) {
      qq.wk <- chol(q.wk,pivot=TRUE)
      sr <- cbind(s.wk,r.wk0[,attr(qq.wk,"pivot")])
      sr <- rbind(sr,cbind(matrix(0,nxiz,nnull),qq.wk))
      sr <- qr(sr,tol=0)
      rss <- mean(qr.resid(sr,c(y.wk,rep(0,nxiz)))[1:nobs]^2)
      trc <- sum(qr.Q(sr)[1:nobs,]^2)/nobs
      if (method=="u") score <- rss + alpha*2*varht*trc
      if (method=="v") score <- rss/(1-alpha*trc)^2
      alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
      alpha.wk <- min(alpha.wk,3)
      if (alpha.wk>alpha) {
        if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*trc
        if (method=="v") score <- rss/(1-alpha.wk*trc)^2
      }
      if (return.fit) {
        z <- .Fortran("reg",
                      as.double(cbind(s.wk,r.wk0)), as.integer(nobs), as.integer(nnull),
                      as.double(q.wk), as.integer(nxiz), as.double(y.wk),
                      as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                      as.double(alpha), varht=as.double(varht),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                      wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
        z$score <- score
        assign("fit",z[c(1:5,7)],inherits=TRUE)
      }
    }
    else {
      z <- .Fortran("reg",
                    as.double(cbind(s.wk,r.wk0)), as.integer(nobs), as.integer(nnull),
                    as.double(q.wk), as.integer(nxiz), as.double(y.wk),
                    as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                    as.double(alpha), varht=as.double(varht),
                    score=double(1), dc=double(nn),
                    as.double(.Machine$double.eps),
                    chol=double(nn*nn), double(nn),
                    jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                    wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                    PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
      if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
      assign("fit",z[c(1:5,7)],inherits=TRUE)
      score <- z$score
      alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
      alpha.wk <- min(alpha.wk,3)
      if (alpha.wk>alpha) {
        if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*z$wk[2]
        if (method=="v") score <- z$wk[1]/(1-alpha.wk*z$wk[2])^2
      }
    }
    score
  }
  cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
  ## initialization
  theta <- -log10(apply(r[id.basis,,],3,function(x)sum(diag(x))))
  r.wk <- 0
  for (i in 1:nq) {
    r.wk <- r.wk + 10^theta[i]*r[,,i]
  }
  ## theta adjustment
  return.fit <- FALSE
  z <- sspreg1(s,r.wk,r.wk[id.basis,],y,wt,method,alpha,varht,random)
  theta <- theta + z$theta
  r.wk <- 0
  for (i in 1:nq) {
    theta[i] <- 2*theta[i] + log10(t(z$c)%*%r[id.basis,,i]%*%z$c)
    r.wk <- r.wk + 10^theta[i]*r[,,i]
  }
  if (!is.null(wt)) q.wk <- wt*r.wk
  else q.wk <- r.wk
  log.la0 <- log10(sum(q.wk^2)/sum(diag(r.wk[id.basis,])))
  log.th0 <- theta-log.la0
  ## lambda search
  z <- sspreg1(s,r.wk,r.wk[id.basis,],y,wt,method,alpha,varht,random)
  nlambda <- z$nlambda
  log.th0 <- log.th0 + z$nlambda
  theta <- theta + z$theta
  if (!is.null(random)) ran.scal <- z$ran.scal
  ## early return
  if (skip.iter) {
    z$theta <- theta
    return(z)
  }
  ## theta search
  fit <- NULL
  counter <- 0
  y.wk <- y
  s.wk <- s
  r.wk <- 0
  for (i in 1:nq) {
    r.wk <- r.wk + 10^theta[i]*r[,,i]
  }
  theta.old <- theta
  if (!is.null(random)) theta <- c(theta,z$zeta)
  ## scale and shift cv
  tmp <- abs(cv(theta))
  cv.scale <- 1
  cv.shift <- 0
  if (tmp<1&tmp>10^(-4)) {
    cv.scale <- 10/tmp
    cv.shift <- 0
  }
  if (tmp<10^(-4)) {
    cv.scale <- 10^2
    cv.shift <- 10
  }
  repeat {
    zz <- nlm(cv.wk,theta,stepmax=1,ndigit=7)
    if (zz$code<=3)  break
    theta <- zz$est        
    counter <- counter + 1
    if (counter>=5) {
      warning("gss warning in ssanova: iteration for model selection fails to converge")
      break
    }
  }
  ## return
  return.fit <- TRUE
  jk1 <- cv(zz$est)
  r.wk <- 0
  for (i in 1:nq) {
    r.wk <- r.wk + 10^zz$est[i]*r[,,i]
  }
  qq.wk <- r.wk[id.basis,]
  if (is.null(random)) q.wk <- qq.wk
  else {
    r.wk <- cbind(r.wk,10^(ran.scal)*random$z)
    q.wk <- matrix(0,nxiz,nxiz)
    q.wk[1:nxi,1:nxi] <- qq.wk
    q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
      10^(2*ran.scal-nlambda)*random$sigma$fun(zz$est[-(1:nq)],random$sigma$env)
  }
  if (!is.null(wt)) {
    s <- wt*s
    r.wk <- wt*r.wk
  }
  se.aux <- regaux(s,r.wk,q.wk,nlambda,fit)
  c <- fit$dc[nnull+(1:nxi)]
  if (nnull) d <- fit$dc[1:nnull]
  else d <- NULL
  if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
  else b <- NULL
  c(list(method=method,theta=zz$est[1:nq],c=c,d=d,b=b,nlambda=nlambda,
         zeta=zz$est[-(1:nq)]),fit[-3],list(se.aux=se.aux))
}
