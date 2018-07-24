cmlST = function(y, X, latentVars){
# Complete Maximum Likelihood estimates for a p-variate Skew-t model.
# Arguments:
# y: data matrix
# X: a design matrix
# latentVars: list of (at least) two elements
#  'z': the value of the latent variable z
#  'v': the value of the latent variable v
 n = nrow(y)
 p = ncol(y)
 
 # Latent variables
 if(is.null(latentVars$v)){
  stop('cmlST.R: values of v are missing.\n')
 } else {
  v = latentVars$v
 }
 if(is.null(latentVars$z)){
  stop('cmlST.R: values of z are missing.\n')
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
 sum.v = sum(v)
 sqrt.v = sqrt(v)
 sum.sqrtv.absz = sum(sqrt.v * abs.z)
 sum.sqrtv.absz.y = apply(matrix(rep(sqrt.v*abs.z, each=p), n, p, byrow=T) * y, 2, sum)
 sum.v.y = apply(matrix(rep(v, each=p), n, p, byrow=T) * y, 2, sum)
 sum.z2 = sum(z^2)

 # Complete Maximum Likelihood estimates
 # nu
 Lnu = function(x){
  L = n * log(x/2) - n * digamma(x/2) - sum(v) + sum(log(v)) + n
  return(L)
 }
 if(all(v == 1)){
  nuCML = Inf
 } else {
  nuCML = uniroot(Lnu, c(0.0001,1000))$root
 }
 # psi
 psiNUM = (sum.v * sum.sqrtv.absz.y - sum.sqrtv.absz * sum.v.y)
 psiDEN = sum.z2 * sum.v - sum.sqrtv.absz^2
 if(all(z == 0)){
  psiCML = rep(0, p)
 } else {
  psiCML = psiNUM / psiDEN
 }
 if(is.null(X)){
  # xi & G
  xiCML = (sum.v.y - psiCML * sum.sqrtv.absz) / sum.v
  e = y - matrix(rep(xiCML,each=n), n, p, byrow=F) - matrix(rep(psiCML,each=n),n,p,byrow=F) * matrix(rep(abs.z, each=p),n,p,byrow=T) / matrix(rep(sqrt.v, each=p),n,p,byrow=T)
  sqrtv.e = matrix(rep(sqrt.v, each=p),n,p,byrow=T) * e
  GCML = (n-1) * var(sqrtv.e) / n # GCML = t(sqrtv.e) %*% sqrtv.e / n
  SigmaCML = GCML + psiCML %*% t(psiCML)
  omegaCML = diag(sqrt(diag(SigmaCML))) # diag(SigmaCML)^(0.5)
  thetaCML = list(xi=xiCML, G=GCML, psi=psiCML, nu= nuCML)
# OmegaCML = cov2cor(SigmaCML)
# rho.ijCML = OmegaCML[pmat.indices.wo]
#	deltaCML = psiCML/sqrt(diag(GCML)+psiCML^2) # psiCML / diag(omegaCML)
# alphaCML = as.numeric((1-t(deltaCML) %*% solve(OmegaCML) %*% deltaCML)^(-0.5)) * solve(OmegaCML) %*% deltaCML
 } else {
  # B & G
  k = ncol(X)
  BCML = matrix(0, k, p)
  GCML = matrix(0, p, p)
  constFlag = identical(as.integer(X[,1]), as.integer(rep(1,n)))
  yStar = matrix(0, n, p)
  for(iCol in 1:p){
   yStar[,iCol] = y[,iCol] - rep(psiCML[iCol], n) * abs.z / sqrt.v
  }
  if(constFlag){
   lmFit = stats::lm(yStar ~ X - 1)
  } else {
   lmFit = stats::lm(yStar ~ X)
  }
  lmCoeff = lmFit$coefficients
  if(any(is.na(lmCoeff))) stop('cmlST.R: X is misspecified.\n')
  BCML = matrix(lmCoeff, k, p)
  BX = lmFit$fitted.values   # Conditional mean of y given X
  #
  e = yStar - BX
  sqrtv.e = matrix(rep(sqrt.v, each=p), ncol=p, byrow=T) * e
  GCML = (n-1) * var(sqrtv.e) / n
  thetaCML = list(B=BCML, G=GCML, psi=psiCML, nu=nuCML)
 }
 return(thetaCML)
}
