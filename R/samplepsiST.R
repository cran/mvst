samplepsiST = function(y, X, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable Psi, in the p-variate skew-t model.
 n = nrow(y)
 p = ncol(y)
 #
 G = particles$G
 z = particles$z
 v = particles$v
 Ry = y[rep(1:n, times=N),,drop=F]
 XFlag = !is.null(X)
 if(XFlag){
  k = ncol(X)
  B = particles$B
  RX = X[rep(1:n, times=N),,drop=F]
 } else {
  xi = particles$xi
 }
 RM = matrix(0, n*N, p) # xi or B'X
 for(icol in 1:p){
  if(XFlag){
   Bcol = matrix(B[,icol,], ncol=k, byrow=T)
   RB = Bcol[rep(1:N, each=n),]
   RM[,icol] = apply(RB*RX,1,sum)
  } else {
   RM[,icol] = xi[rep(1:N, each=n),icol]
  }
 }
 abs.z = abs(as.numeric(t(z)))
 sums.z2 = apply(z^2, 1, sum)
 sqrt.v = sqrt(as.numeric(t(v)))
 rm(particles)
 #
 mpsi = matrix(0,N,p)
 for(icol in 1:p){
  mpsi[,icol] = ((apply(matrix(Ry[,icol] * abs.z * sqrt.v,N,n,byrow=T),1,sum) - apply(matrix(RM[,icol] * abs.z * sqrt.v,N,n,byrow=T),1,sum))) / sums.z2
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
 # PSI = rmnorm.mia2(mpsi,G[,pmat.indices.wVEC]/sums.z2)		# Attenzione a non confondere con l'oggetto 'Psi'
 # psi = PSI$x
 # log.dpsi = PSI$log.dens

 return(list(values=psi, log.dq=log.dpsi))
}
