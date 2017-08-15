

####################################################
#      
#      file: entropy-loss.R
#      author: Tayler Blake
#
#      inputs: omegaHat - an estimated inverse 
#                         covariance matrix
#              Sig    - the true covariance matrix
#
#      output: entropy_loss - non-negative loss 
#                             measure
#
####################################################


entropy_loss <- function(omegaHat, Sig) {
            I_hat <- omegaHat%*%Sig
            sum(diag(I_hat)) - log(det(I_hat)) - ncol(Sig)
}