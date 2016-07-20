samplepsiSN = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable Psi, in the p-variate skew-N model.
 n = nrow(y)
 p = ncol(y)
#
# pmat.indices.wVEC = triangleIndices(p, side='u', dgn=T, dataframe=F)
 Ry = matrix(0, n*N, p)
 Rxi = matrix(0, n*N, p)
 for(icol in 1:p){
  Ry[,icol] = rep(y[,icol], N)
  Rxi[,icol] = rep(particles$xi[,icol], each=n)
 }
 abs.z = abs(as.numeric(t(particles$z)))
 sums.z2 = apply((particles$z)^2, 1, sum)
# sqrt.v = sqrt(as.numeric(t(particles$v)))
 G = particles$G
 rm(particles)
 #
 mpsi = matrix(0,N,p)
 for(icol in 1:p){
  mpsi[,icol] = ((apply(matrix(Ry[,icol] * abs.z, N, n, byrow=T),1,sum) - apply(matrix(Rxi[,icol] * abs.z, N, n, byrow=T), 1, sum))) / sums.z2
 }

 psi = matrix(0, N, p)
 log.dpsi = numeric(N)
 for(iN in 1:N){
  G.iN = matrix(G[iN,], p, p)
  varPsi.iN = G.iN / sums.z2[iN]
		varPsi.iN = (varPsi.iN + t(varPsi.iN))/2	# force symmetry
  psiValue = as.numeric(rmnorm(1, mpsi[iN,], varPsi.iN))
  psi[iN,] = psiValue
  log.dpsi[iN] = dmnorm(psiValue, mpsi[iN,], varPsi.iN, log=TRUE)
 }
 return(list(values=psi, log.dq=log.dpsi))
}
