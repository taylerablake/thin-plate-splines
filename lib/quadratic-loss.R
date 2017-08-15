
####################################################
#      
#      file: quadratic-loss.R
#      author: Tayler Blake
#
#      inputs: omegaHat - an estimated inverse 
#                         covariance matrix
#              Sig    - the true covariance matrix
#
#      output: quadratic_loss - non-negative loss 
#                             measure
#
####################################################


quadratic_loss <- function(omegaHat, Sig) {
      I_hat <- omegaHat%*%Sig
      sum(diag(I_hat)^2)
}