cmlST = function(y, latentVars){
# Complete Maximum Likelihood estimates for a p-variate Skew-Elliptical model.
# Arguments:
# y: data matrix
# latentVars: list of two elements:
#  'z': real value of the vector z
#  'v': real value of the vector v
 n = nrow(y)
 p = ncol(y)
 
 # Latent variables
 if(is.null(latentVars$v)){
  v = rep(1, n)
 } else {
  v = latentVars$v
 }
 if(is.null(latentVars$z)){
  z = rep(0, n)
 } else {
  z = latentVars$z
 }

 #
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
 n.pmat.indices.wo = p * (p-1) / 2
 #
 sum.v = sum(v)
 sum.sqrtv.absz = sum(sqrt(v) * abs(z))
 sum.sqrtv.absz.y = apply(matrix(rep(sqrt(v)*abs(z), each=p), n, p, byrow=T) * y, 2, sum)
 sum.v.y = apply(matrix(rep(v, each=p), n, p, byrow=T) * y, 2, sum)
 sum.z2 = sum(z^2)
 # Complete Maximum Likelihood estimates
 psiNUM = (sum.v * sum.sqrtv.absz.y - sum.sqrtv.absz * sum.v.y)
 psiDEN = sum.z2 * sum.v - sum.sqrtv.absz^2
 if(all(z == 0)){
  psiCML = rep(0, p)
 } else {
  psiCML = psiNUM / psiDEN
 }
 xiCML = (sum.v.y - psiCML * sum.sqrtv.absz) / sum.v
 e = y - matrix(rep(xiCML,each=n), n, p, byrow=F) - matrix(rep(psiCML,each=n),n,p,byrow=F) * matrix(rep(abs(z), each=p),n,p,byrow=T) / matrix(rep(sqrt(v), each=p),n,p,byrow=T)
 sqrtv.e = matrix(rep(sqrt(v), each=p),n,p,byrow=T) * e
 GCML = t(sqrtv.e) %*% sqrtv.e / n
 SigmaCML = GCML + psiCML %*% t(psiCML)
	omegaCML = diag(sqrt(diag(SigmaCML))) # diag(SigmaCML)^(0.5)
# rhoCML = SigmaCML[1,2] / prod(omegaCML)
# OmegaCML = matrix(c(1, rhoCML, rhoCML, 1), ncol=2)
 OmegaCML = cov2cor(SigmaCML)
 rho.ijCML = OmegaCML[pmat.indices.wo]
	deltaCML = psiCML/sqrt(diag(GCML)+psiCML^2) # psiCML / diag(omegaCML)
 alphaCML = as.numeric((1-t(deltaCML) %*% solve(OmegaCML) %*% deltaCML)^(-0.5)) * solve(OmegaCML) %*% deltaCML
# Lnu = function(x, n=n, v=v){
 Lnu = function(x){
  L = n * log(x/2) - n * digamma(x/2) - sum(v) + sum(log(v)) + n
  return(L)
 }
 if(all(v == 1)){
  nuCML = Inf
 } else {
  nuCML = uniroot(Lnu, c(0.0001,1000))$root
 }
 thetaCML = c(xiCML, psiCML, GCML[pmat.indices.w], nuCML, deltaCML, alphaCML, diag(omegaCML), rho.ijCML)
 names(thetaCML) = c(paste(rep(c('xi', 'psi'),each=p), 1:p, sep=''), paste(rep('G', n.pmat.indices.w), paste(pmat.indices.w[,1], pmat.indices.w[,2], sep='.'), sep=''), 'nu', paste(rep(c('delta', 'alpha', 'omega'),each=p), 1:p, sep=''), paste(rep('rho', n.pmat.indices.wo), paste(pmat.indices.wo[,1], pmat.indices.wo[,2], sep='.'), sep=''))

 return(thetaCML)
}
