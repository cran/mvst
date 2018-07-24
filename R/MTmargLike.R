MTmargLike = function(y, X=NULL, LOG=FALSE){
 # Checks
 if(is.data.frame(y)) y = as.matrix(y)
 if(is.vector(y)) y = matrix(y, ncol=1)
 dimnames(y) = NULL
 n = nrow(y)
 p = ncol(y)
 #
 XFlag = !is.null(X)
 if(XFlag){
  if(!is.matrix(X)) stop('(MTmargLike.R) the object X is not a matrix.\n')
  if(nrow(X) != n) stop('(MTmargLike.R) matrices y and X should have the same number of rows.\n')
  k = ncol(X)
 }

 #
 if(!XFlag){
  nS = (n-1) * var(y)
  log.psi = 0
  for(j in 1:p){
   log.psi = log.psi + lgamma(0.5*(n-j))
  }
  log.m0 = -0.5*p*log(n) + log.psi - 0.25*p*(2*n-p-1)*log(pi) - (0.5*(n-1)) * as.numeric(determinant(nS, logarithm=TRUE)$modulus)
  # log.m0 = ((3-2*n)/2)*log(pi)+log.psi-n*log(n)-((n-1)/2)*log(det(S))
  m0 = exp(log.m0)
 } else {
  OmegaX = solve(t(X) %*% X)
  logdet.OmegaX = as.numeric(determinant(OmegaX, logarithm=TRUE)$modulus)
  nS = t(y) %*% y - t(y) %*% X %*% OmegaX %*% t(X) %*% y
  logdet.nS = as.numeric(determinant(nS, logarithm=TRUE)$modulus)
  log.psi = 0
  for(j in 1:p){
   log.psi = log.psi + lgamma(0.5*(n-k+1-j))
  }
  log.m0 = 0.5*p*logdet.OmegaX + log.psi - 0.25*p*(2*(n-k)-p+1)*log(pi) - (0.5*(n-k)) * logdet.nS
  m0 = exp(log.m0)
 }
 if(LOG) return(log.m0) else return(m0)
}
