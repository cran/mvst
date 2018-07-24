cmlN = function(y, X, latentVars=NULL){
# Complete Maximum Likelihood estimates for a p-variate Normal model.
# Arguments:
# y: data matrix
# X: a design matrix
# latentVars: an empty list
 n = nrow(y)
 p = ncol(y)
 #
# pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
# pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
# n.pmat.indices.w = p * (p+1) / 2
# n.pmat.indices.wo = p * (p-1) / 2
 # Complete Maximum Likelihood estimates
 if(is.null(X)){
  # xi & G
  xiCML = apply(y, 2, mean)
  GCML = (n-1) * var(y) / n
 # SigmaCML = GCML
 #	omegaCML = diag(sqrt(diag(SigmaCML)))
 # OmegaCML = cov2cor(SigmaCML)
 # rho.ijCML = OmegaCML[pmat.indices.wo]
  thetaCML = list(xi=xiCML, G=GCML)
 } else {
  # B & G
  k = ncol(X)
  BCML = matrix(0, k, p)
  GCML = matrix(0, p, p)
  constFlag = identical(as.integer(X[,1]), as.integer(rep(1,n)))
  if(constFlag){
   lmFit = stats::lm(y ~ X - 1)
  } else {
   lmFit = stats::lm(y ~ X)
  }
  lmCoeff = lmFit$coefficients
  if(any(is.na(lmCoeff))) stop('cmlN.R: X is misspecified.\n')
  BCML = matrix(lmCoeff, k, p)
  BX = lmFit$fitted.values   # Conditional mean of y given X
  #
  e = y - BX
  GCML = (n-1) * var(e) / n
  thetaCML = list(B=BCML, G=GCML)
 }
 return(thetaCML)
}
