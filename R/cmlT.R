cmlT = function(y, X, latentVars){
# Complete Maximum Likelihood estimates for a p-variate Student-t model.
# Arguments:
# y: data matrix
# X: a design matrix
# latentVars: list of (at least) the element 'v',
#  the value of the latent variables v
 n = nrow(y)
 p = ncol(y)
 
 # Latent variables
 if(is.null(latentVars$v)){
  stop('cmlT.R: values of v are missing.\n')
 } else {
  v = latentVars$v
 }
 #
# pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
# pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
# n.pmat.indices.w = p * (p+1) / 2
# n.pmat.indices.wo = p * (p-1) / 2
 #
 sum.v = sum(v)
 sqrt.v = sqrt(v)
 sum.v.y = apply(matrix(rep(v, each=p), n, p, byrow=T) * y, 2, sum)
 # Complete Maximum Likelihood estimates
 # nu
 Lnu = function(x){
  L = n * log(x/2) - n * digamma(x/2) - sum(v) + sum(log(v)) + n
  return(L)
 }
 nuCML = uniroot(Lnu, c(0.0001,1000))$root
 if(is.null(X)){
  # xi & G
  xiCML = sum.v.y / sum.v
  e = y - matrix(rep(xiCML,each=n), n, p, byrow=F)
  sqrtv.e = matrix(rep(sqrt.v, each=p),n,p,byrow=T) * e
  GCML = (n-1) * var(sqrtv.e) / n # t(sqrtv.e) %*% sqrtv.e / n
 # SigmaCML = GCML
 #	omegaCML = diag(sqrt(diag(SigmaCML))) # diag(SigmaCML)^(0.5)
 # OmegaCML = cov2cor(SigmaCML)
 # rho.ijCML = OmegaCML[pmat.indices.wo]
  thetaCML = list(xi=xiCML, G=GCML, nu=nuCML)
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
  if(any(is.na(lmCoeff))) stop('cmlT.R: X is misspecified.\n')
  BCML = matrix(lmCoeff, k, p)
  BX = lmFit$fitted.values   # Conditional mean of y given X
  #
  e = y - BX
  sqrtv.e = matrix(rep(sqrt.v, each=p), ncol=p, byrow=T) * e
  GCML = (n-1) * var(sqrtv.e) / n
  thetaCML = list(B=BCML, G=GCML, nu=nuCML)
 }
 return(thetaCML)
}
