cmlSN = function(y, latentVars){
 n = nrow(y)
 p = ncol(y)
 
 #
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
 absz = abs(z)
 mean.absz = mean(absz)
 mean.y = apply(y, 2, mean)
 # Complete Maximum Likelihood estimates
 xiCML = numeric(p)
 psiCML = numeric(p)
 for(icol in 1:p){
  psiNUM = sum((absz-mean.absz) * (y[,icol] - mean.y[icol]))
  psiDEN = sum((absz-mean.absz)^2)
  if(all(z == 0)){
   psiCML[icol] = 0
  } else {
   psiCML[icol] = psiNUM / psiDEN
  }
  xiCML = mean.y - psiCML[icol] * mean.absz
 }
 e = y - matrix(rep(xiCML,each=n), n, p, byrow=F) - matrix(rep(psiCML,each=n),n,p,byrow=F) * matrix(rep(abs(z), each=p),n,p,byrow=T)
 GCML = t(e) %*% e / n
 SigmaCML = GCML + psiCML %*% t(psiCML)
	omegaCML = diag(sqrt(diag(SigmaCML)))
 OmegaCML = cov2cor(SigmaCML)
 rho.ijCML = OmegaCML[pmat.indices.wo]
	deltaCML = psiCML/sqrt(diag(GCML)+psiCML^2)
 alphaCML = as.numeric((1-t(deltaCML) %*% solve(OmegaCML) %*% deltaCML)^(-0.5)) * solve(OmegaCML) %*% deltaCML
 #
 thetaCML = c(xiCML, psiCML, GCML[pmat.indices.w], deltaCML, alphaCML, diag(omegaCML), rho.ijCML)
 names(thetaCML) = c(paste(rep(c('xi', 'psi'),each=p), 1:p, sep=''), paste(rep('G', n.pmat.indices.w), paste(pmat.indices.w[,1], pmat.indices.w[,2], sep='.'), sep=''), paste(rep(c('delta', 'alpha', 'omega'),each=p), 1:p, sep=''), paste(rep('rho', n.pmat.indices.wo), paste(pmat.indices.wo[,1], pmat.indices.wo[,2], sep='.'), sep=''))

 return(thetaCML)
}
