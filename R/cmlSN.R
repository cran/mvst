cmlSN = function(y, X, latentVars){
# Complete Maximum Likelihood estimates for a p-variate Skew-Normal model.
# Arguments:
# y: data matrix
# X: a design matrix
# latentVars: a list containing (at least) the
#  element 'z', the value of the latent
#  variables z
 n = nrow(y)
 p = ncol(y)
 #
 if(is.null(latentVars$z)){
  stop('cmlSN.R: values of z are missing.\n')
 } else {
  z = latentVars$z
 }
 #
# pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
# pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
# n.pmat.indices.w = p * (p+1) / 2
# n.pmat.indices.wo = p * (p-1) / 2
 #
 abs.z = abs(z)
 mean.absz = mean(abs.z)
 mean.y = apply(y, 2, mean)
 # Complete Maximum Likelihood estimates
 # psi
 if(all(z == 0)){
  psiCML = rep(0, p)
 } else {
  psiCML = numeric(p)
 }
 for(icol in 1:p){
  psiNUM = sum((abs.z-mean.absz) * (y[,icol] - mean.y[icol]))
  psiDEN = sum((abs.z-mean.absz)^2)
  if(all(z == 0)){
   psiCML[icol] = 0
  } else {
   psiCML[icol] = psiNUM / psiDEN
  }
 }
 if(is.null(X)){
  # xi & G
  xiCML = numeric(p)
  for(icol in 1:p){
   xiCML = mean.y - psiCML[icol] * mean.absz
  }
  e = y - matrix(rep(xiCML,each=n), n, p, byrow=F) - matrix(rep(psiCML,each=n),n,p,byrow=F) * matrix(rep(abs(z), each=p),n,p,byrow=T)
  GCML = t(e) %*% e / n
 # SigmaCML = GCML + psiCML %*% t(psiCML)
 #	omegaCML = diag(sqrt(diag(SigmaCML)))
 # OmegaCML = cov2cor(SigmaCML)
 # rho.ijCML = OmegaCML[pmat.indices.wo]
 #	deltaCML = psiCML/sqrt(diag(GCML)+psiCML^2)
 # alphaCML = as.numeric((1-t(deltaCML) %*% solve(OmegaCML) %*% deltaCML)^(-0.5)) * solve(OmegaCML) %*% deltaCML
  thetaCML = list(xi=xiCML, G=GCML, psi=psiCML)
 } else {
  # B & G
  k = ncol(X)
  BCML = matrix(0, k, p)
  GCML = matrix(0, p, p)
  constFlag = identical(as.integer(X[,1]), as.integer(rep(1,n)))
  yStar = matrix(0, n, p)
  for(iCol in 1:p){
   yStar[,iCol] = y[,iCol] - rep(psiCML[iCol], n) * abs.z
  }
  if(constFlag){
   lmFit = stats::lm(yStar ~ X - 1)
  } else {
   lmFit = stats::lm(yStar ~ X)
  }
  lmCoeff = lmFit$coefficients
  if(any(is.na(lmCoeff))) stop('cmlSN.R: X is misspecified.\n')
  BCML = matrix(lmCoeff, k, p)
  BX = lmFit$fitted.values   # Conditional mean of y given X
  #
  e = yStar - BX
  GCML = (n-1) * var(e) / n
  thetaCML = list(B=BCML, G=GCML, psi=psiCML)
 }
 return(thetaCML)
}
