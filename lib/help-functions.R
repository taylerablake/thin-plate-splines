

# Compute B-spline base matrix
bspline <- function(X., XL., XR., NDX., BDEG.) {
      require(splines)
      dx <- (XR. - XL.)/NDX.
      knots <- seq(XL. - BDEG. * dx, XR. + BDEG. * dx, by = dx)
      B <- spline.des(knots, X., BDEG. + 1, 0 * X.)$design
      B
}

# row-tensor product (or Box product of two matrices)
rowtens <- function(X, B) {
      one.1 <- matrix(1, 1, ncol(X))
      one.2 <- matrix(1, 1, ncol(B))
      kronecker(X, one.2) * kronecker(one.1, B)
}


# Mixed Model Basis
MM.basis <- function(x, xl, xr, ndx, bdeg, pord) {
      
      B <- bspline(x, xl, xr, ndx, bdeg)
      m <- ncol(B)
      
      if (pord > 0) {
            D <- diff(diag(m), differences = pord)
      } else {
            D <- diag(m)
      }
      
      P.svd <- svd(t(D) %*% D)
      U <- (P.svd$u)[, 1:(m - pord)]
      d <- (P.svd$d)[1:(m - pord)]
      
      Z <- B %*% U
      
      X <- NULL
      if (pord > 0) {
            for (i in 0:(pord - 1)) {
                  X <- cbind(X, x^i)
            }
      } else if (pord == 0) {
            X <- cbind(X,x^0)
      }
      
      list(X = X, Z = Z, d = d, B = B, m = m, D = D)
}

# Construct 2 x 2 block symmetric matrices:
construct.block2 <- function(A1, A2, A4) {
      block <- rbind(cbind(A1, A2), cbind(t(A2), A4))
      return(block)
}


# construct the mixed model matrices from marginal B-spline bases
mmodel.frame2d <- function(x1, x2,
                   nseg1, nseg2,
                   div, bdeg,
                   pord1, pord2) {
      
      pord1n <- pord1
      pord2n <- pord2
      
      MM1 <- MM.basis(x1, min(x1) - 0.01, max(x1) + 0.01, nseg1, 
                      bdeg, pord1)  # bases for x1
      X1 <- MM1$X
      G1 <- MM1$Z
      d1 <- MM1$d
      B1 <- MM1$B
      
      MM2 <- MM.basis(x2, min(x2) - 0.01, max(x2) + 0.01, nseg2, 
                      bdeg, pord2)  # bases for x2
      X2 <- MM2$X
      G2 <- MM2$Z
      d2 <- MM2$d
      B2 <- MM2$B
      
      MM1n <- MM.basis(x1, min(x1) - 0.01, max(x1) + 0.01, nseg1/div, 
                       bdeg, pord1n)  # Nested bases for x1
      G1n <- MM1n$Z
      d1n <- MM1n$d
      B1n <- MM1n$B
      
      MM2n <- MM.basis(x2, min(x2) - 0.01, max(x2) + 0.01, nseg2/div, 
                       bdeg, pord2n)  # Nested bases for x2
      G2n <- MM2n$Z
      d2n <- MM2n$d
      B2n <- MM2n$B
      
      c1 <- ncol(B1)
      c2 <- ncol(B2)
      c1n <- ncol(B1n)
      c2n <- ncol(B2n)
      
      X <- rowtens(X2, X1)  # -> Fixed effects
      
      #####################
      
      d3 <- c(rep(1,
                  c2n - pord2) %x% d1n + d2n %x% rep(1,
                                                     c1n - pord1))
      Delta1 <- diag(1/sqrt(d1))
      Delta2 <- diag(1/sqrt(d2))
      Delta3 <- diag(1/sqrt(d3))

      # random effects matrices
      Z1 <- G1 %*% Delta1  # smooth random comp. fx1
      Z2 <- G2 %*% Delta2  # smooth random comp. fx2
      Z12 <- rowtens(G2n, G1n) %*% Delta3  # smooth interaction  fx1:fx2

      Z2ListNames <- c("Z2x1","Z2x1sq","Z2x1cubed")
      Z1ListNames <- c("Z1x2","Z1x2sq","Z1x2cubed")
      
      Z2squareX1 <- NULL
      X2squareZ1 <- NULL
      if ( pord1>0 ) {
            Z2squareX1 <- list()
            for(power in 2:ncol(X1)) {
                  Z2squareX1 <- list.append(Z2squareX1,
                                            rowtens(G2n, matrix(X1[,power],ncol=1)))     # linear:smooth comp. fx2 x [x1, x1^2,...] 
            }
            names(Z2squareX1) <- Z2ListNames[1:length(Z2squareX1)]
      }
      if ( pord2>0 ) {
            X2squareZ1 <- list()
            for ( power in 2:ncol(X2) ) {
                  X2squareZ1 <- list.append(X2squareZ1,
                                            rowtens(matrix(X2[,power],ncol=1), G1n) ) # linear:smooth comp. fx1 x [x2, x2^2,...] 
            }
            names(X2squareZ1) <- Z1ListNames[1:length(X2squareZ1)]
      }

      
      
      
      
      ## BUILD RANDOM EFFECTS MATRIX
      
      Zlist <- list(Z1=Z1,Z2=Z2)
      if (pord1 > 0){
            for (this_col in 1:length(Z2squareX1)) {
                  Zlist <- list.append(Zlist,Z2squareX1[[this_col]])
                  names(Zlist)[2+this_col] <- names(Z2squareX1)[this_col]
            }
      }
      lngth <- length(Zlist)
      if (pord2 > 0){
            for (this_col in 1:length(X2squareZ1)) {
                  Zlist <- list.append(Zlist,X2squareZ1[[this_col]])
                  names(Zlist)[lngth+this_col] <- names(X2squareZ1)[this_col]
            }
      }
      Zlist <- list.append(Zlist,Z12=Z12)
      Z <- list.cbind(Zlist)
      M <- cbind(X, Z)

      list(M=M,
           X=X,
           X1=X1,
           X2=X2,
           Zlist=Zlist,
           B1=B1,
           B2=B2,
           Z=Z)
}









# construct the mixed model matrices from marginal B-spline bases
mmodel.pred.frame2d <- function(x1, x2,
                           nseg1, nseg2,
                           div, bdeg,
                           pord1, pord2,
                           ngrid) {

      g1 <- rep(seq(min(x1), max(x1), length = ngrid), ngrid)
      g2 <- rep(seq(min(x2), max(x2), length = ngrid), rep(ngrid, ngrid))
      
      mmodel.frame2d(g1, g2,
                     nseg1, nseg2,
                     div, bdeg,
                     pord1, pord2)
}




