reparamT = function(N, p, particles){
 # Reparameterization function
 n.pmat.indices.wo = p * (p-1) / 2
 #
# alpha = matrix(0, N, p)
# delta = matrix(0, N, p)
 omega = matrix(0, N, p)
 rho = matrix(0, N, n.pmat.indices.wo)
# h = matrix(0, N, p)
# out = numeric(N)
 log.detSigma = numeric(N)
 log.detOmega = numeric(N)
 detG = numeric(N)
 diagIndices = intersect(triangleIndices(p, side='u', dgn=T, dataframe=F), triangleIndices(p, side='l', dgn=T, dataframe=F))
 pmat.indices.woDF = triangleIndices(p, side='u', dgn=F, dataframe=T)
  #
 G = particles$G
# psi = particles$psi
 for(icol in 1:p){
#  h[,icol] = (G[,diagIndices[icol]] + psi[,icol]^2)^(0.5)
#  delta[,icol] = psi[,icol] / h[,icol] # (h[,icol])^(-0.5) * psi[,icol]
  omega[,icol] = sqrt(G[,diagIndices[icol]])       # G = Sigma in the symmetric cases
  # omega[,icol] = sqrt(G[,diagIndices[icol]]/(1-delta[,icol]^2))	# metodo alternativo per ottenere omega
 }

 for(iN in 1:N){
  GiN = matrix(G[iN,], p, p)
  Sigma.iN = GiN # + psi[iN,] %*% t(psi[iN,])
#  if(!all(eigen(Sigma.iN)$values>0)) out[iN] = 1
  log.detSigma[iN] = log(det(Sigma.iN))
  Omega.iN = cov2cor(Sigma.iN)
  rho[iN,] = as.numeric(Omega.iN[pmat.indices.woDF])
  log.detOmega[iN] = log(det(Omega.iN))
#  alpha[iN,] = (1 - as.numeric(t(delta[iN,]) %*% solve(Omega.iN, delta[iN,]))) * solve(Omega.iN, delta[iN,])
#  detG[iN] = det(GiN)
 }
# pos = 1 - out
# if(sum(out, na.rm=T)>0) browser()
 return(list(omega=omega, rho=rho, log.detOmega=log.detOmega, log.detSigma=log.detSigma)) # alpha=alpha, delta=delta, h=h, pos=pos
}
