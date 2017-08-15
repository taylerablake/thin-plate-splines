


##############################################################################
##
##    Arguments: N - sample size (number of independent
##                   vectors to be generated)
##               M - dimension of the random vectors
##               rho - off-diagonal elements  
##               tausq - term to be added to the covariace
##                       to be added to the common variance  
##
##    Return valuess: y - N x M matrix of simulated vectors
##                    Sigma - M x M covariance matrix                    
##                    Omega - inverse covariance matrix
##                    T_mat - Cholesky factor of inverse covariance
##                    D - diagonal matrix of innovation variances
##                    phi - vector of the generalized varying coefficient
##                          function evaluated at the observed design points
##                    grid - dataframe containing the unscaled and scaled 
##                           observation coordinates
##
##############################################################################

generate_compound_symmetry <- function(N, M, rho, tausq) { 
  
  Grid <- build_grid(M)
  
  ## scale the predictors to lie within (0,1)
  
  Grid <- Grid %>%
    transform(.,l_s=l/(max(Grid$l)+min(Grid$l)),
                                 m_s=m/(max(Grid$m)+min(Grid$m)))
  
  ## define the covariace and precision matrices
  Sigma <- Sigma <- matrix(rho,nrow=M,ncol=M) + diag(rep(tausq,M))
  Omega <- solve(Sigma)
  
  ## construct the cholesky decomposition
  C <- t(chol(Sigma))
  D <- diag(diag(C))
  Dsq <- diag(diag(C)^2)
  L <- C%*%solve(D)
  T_mat <- solve(L)
  phi <- -as.vector(T_mat[lower.tri(T_mat)])    
  
  y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)

  list(grid=Grid,
       y=y,
       Sigma=Sigma,
       Omega=Omega,
       D=D,
       Dsq=Dsq,
       T_mat=T_mat,
       phi_vec=phi)
  
  }


