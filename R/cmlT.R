cmlT = function(y, latentVars){
 n = nrow(y)
 p = ncol(y)
 v = latentVars[['v']]
 #
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
 n.pmat.indices.wo = p * (p-1) / 2
 #
 sum.v = sum(v)
# sum.sqrtv.absz = sum(sqrt(v) * abs(z))
# sum.sqrtv.absz.y = apply(matrix(rep(sqrt(v)*abs(z), each=p), n, p, byrow=T) * y, 2, sum)
 sum.v.y = apply(matrix(rep(v, each=p), n, p, byrow=T) * y, 2, sum)
# sum.z2 = sum(z^2)
 # Complete Maximum Likelihood estimates
# psiCML = (sum.v * sum.sqrtv.absz.y - sum.sqrtv.absz * sum.v.y) / (sum.z2 * sum.v - sum.sqrtv.absz^2)
 xiCML = sum.v.y / sum.v
 e = y - matrix(rep(xiCML,each=n), n, p, byrow=F)
 sqrtv.e = matrix(rep(sqrt(v), each=p),n,p,byrow=T) * e
 GCML = t(sqrtv.e) %*% sqrtv.e / n
 SigmaCML = GCML
	omegaCML = diag(sqrt(diag(SigmaCML))) # diag(SigmaCML)^(0.5)
# rhoCML = SigmaCML[1,2] / prod(omegaCML)
# OmegaCML = matrix(c(1, rhoCML, rhoCML, 1), ncol=2)
 OmegaCML = cov2cor(SigmaCML)
 rho.ijCML = OmegaCML[pmat.indices.wo]
 Lnu = function(x){
  L = n * log(x/2) - n * digamma(x/2) - sum(v) + sum(log(v)) + n
  return(L)
 }
 nuCML = uniroot(Lnu, c(0.0001,1000))$root
 thetaCML = c(xiCML, GCML[pmat.indices.w], nuCML, diag(omegaCML), rho.ijCML)
 names(thetaCML) = c(paste(rep('xi',each=p), 1:p, sep=''), paste(rep('G', n.pmat.indices.w), paste(pmat.indices.w[,1], pmat.indices.w[,2], sep='.'), sep=''), 'nu', paste(rep('omega',each=p), 1:p, sep=''), paste(rep('rho', n.pmat.indices.wo), paste(pmat.indices.wo[,1], pmat.indices.wo[,2], sep='.'), sep=''))

 return(thetaCML)
}
