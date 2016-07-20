samplepsiST = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable Psi, in the p-variate skew-t model.
 n = nrow(y)
 p = ncol(y)
#
 pmat.indices.wVEC = triangleIndices(p, side='u', dgn=T, dataframe=F)
 Ry = matrix(0, n*N, p)
 Rxi = matrix(0, n*N, p)
 for(icol in 1:p){
  Ry[,icol] = rep(y[,icol], N)
  Rxi[,icol] = rep(particles$xi[,icol], each=n)
 }
 abs.z = abs(as.numeric(t(particles$z)))
 sums.z2 = apply((particles$z)^2, 1, sum)
 sqrt.v = sqrt(as.numeric(t(particles$v)))
 G = particles$G
 rm(particles)
 #
 mpsi = matrix(0,N,p)
 for(icol in 1:p){
  mpsi[,icol] = ((apply(matrix(Ry[,icol] * abs.z * sqrt.v,N,n,byrow=T),1,sum) - apply(matrix(Rxi[,icol] * abs.z * sqrt.v,N,n,byrow=T),1,sum)))/sums.z2
 }

 psi = matrix(0, N, p)
 log.dpsi = numeric(N)
# out = numeric(N)
 for(iN in 1:N){
  G.iN = matrix(G[iN,], p, p)
  varPsi.iN = G.iN / sums.z2[iN]
		varPsi.iN = (varPsi.iN + t(varPsi.iN))/2	# force symmetry

#   if(!all(eigen(varPsi.iN)$values>0)) stop(c('varPsi not PD.\n varPsi.iN = ', as.numeric(varPsi.iN)))
#   if(!isSymmetric(varPsi.iN)){
#    cat('avviso 2: varPsi =', as.numeric(varPsi.iN))
#    browser()
#   }
#   if (max(abs(varPsi.iN - t(varPsi.iN))) > .Machine$double.eps){
#    cat('avviso 2: varPsi =', as.numeric(varPsi.iN))
#    browser()
#   }
  
  psiValue = as.numeric(rmnorm(1, mpsi[iN,], varPsi.iN))
  psi[iN,] = psiValue
  log.dpsi[iN] = dmnorm(psiValue, mpsi[iN,], varPsi.iN, log=TRUE)
#  out[iN] = as.numeric(!all(eigen(G.iN + psiValue %*% t(psiValue))$values > 0))
 }
 # PSI = rmnorm.mia2(mpsi,G[,pmat.indices.wVEC]/sums.z2)		# Attenzione a non confondere con l'oggetto 'Psi'
 # psi = PSI$x
 # log.dpsi = PSI$log.dens

 # Indices of the columns of G related to diagonal elements of the matrices
# diagIndices = intersect(triangleIndices(p, side='u', dgn=T, dataframe=F), triangleIndices(p, side='l', dgn=T, dataframe=F))
# out = sort(unique(which(psi^2 > G[,diagIndices],arr.ind=T)[,1]))

#  if(sum(out) > 0) psi[out,] = matrix(0,sum(out),p)
#  pos = rep(1,N)				# Indicatore: 1 se Sigma (altra param) e' pos-def
#  pos[out] = rep(0,length(out))

 return(list(values=psi, log.dq=log.dpsi)) #,reppsi=Rpsi,detGsign=pos))
}
