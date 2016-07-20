cmlN = function(y, latentVars=NULL){
 n = nrow(y)
 p = ncol(y)
 z = latentVars[['z']]
 #
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
 n.pmat.indices.wo = p * (p-1) / 2
 #
#  absz = abs(z)
#  mean.absz = mean(absz)
#  mean.y = apply(y, 2, mean)
 # Complete Maximum Likelihood estimates
 xiCML = apply(y, 2, mean)
 GCML = (n-1) * var(y) / n
 SigmaCML = GCML
	omegaCML = diag(sqrt(diag(SigmaCML)))
 OmegaCML = cov2cor(SigmaCML)
 rho.ijCML = OmegaCML[pmat.indices.wo]
 #
 thetaCML = c(xiCML, GCML[pmat.indices.w], diag(omegaCML), rho.ijCML)
 names(thetaCML) = c(paste(rep('xi', each=p), 1:p, sep=''), paste(rep('G', n.pmat.indices.w), paste(pmat.indices.w[,1], pmat.indices.w[,2], sep='.'), sep=''), paste(rep('omega', each=p), 1:p, sep=''), paste(rep('rho', n.pmat.indices.wo), paste(pmat.indices.wo[,1], pmat.indices.wo[,2], sep='.'), sep=''))

 return(thetaCML)
}
