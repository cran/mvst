samplezST = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the latent variable Z, in the p-variate skew-t model.
 n = nrow(y)
 p = ncol(y)
 #
 yVec = as.numeric(t(y))
 xiVec = as.numeric(t(particles$xi))
 psiVec = as.numeric(t(particles$psi))
 GVec = as.numeric(t(particles$G))
 # zVec = as.numeric(t(particles$z))
 values = as.double(rep(0.0, n*N))
 propdens = as.double(rep(0.0, N))
 vVec = as.numeric(t(particles$v))
 # nu = particles$nu

 # pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T) 
 #
 M = as.double(rep(0.0,n*N))
 RV = as.double(rep(0.0,n*N))
 absz = as.double(rep(0.0, n*N))
 zPars = .C('rzST', as.double(absz), as.double(M), as.double(RV), as.integer(n), as.integer(p), as.integer(N), as.double(yVec), as.double(xiVec), as.double(psiVec), as.double(GVec), as.double(vVec), NAOK=TRUE)
 absz = zPars[[1]]
 M = zPars[[2]]
 RV = zPars[[3]]
 truncProbs = pnorm(0, M, sqrt(RV))
 zSigns = sample(c(1,-1), size=n*N, replace=T)
 truncProbs = pnorm(0, M, sqrt(RV))
 # norm.densities = dnorm(zi.ass,moda,sqrt(var.z),log=T)-log(2)-log(1-p.tronca)
 log.dZI = dnorm(absz, M, sqrt(RV), log=T)		 # density of normal rv's with mean M and variance RV
 norm.dens = apply(matrix(log.dZI, N, n, byrow=T), 1, sum)
 truncs = apply(matrix(log(1-truncProbs), N, n, byrow=T), 1, sum)
 log.dz = norm.dens - rep(n*log(2), N) - truncs	# pdf of Z
 Rz = zSigns * absz
 z = matrix(Rz, N, n, byrow=T)
# mean.absz = apply(abs(z), 1, mean) # TROVA E SOSTITUISCI mean.absz
# sums.z2 = apply(z^2, 1, sum)
 return(list(values=z, log.dq=log.dz)) #, repz=abs.z, mean.abs=mean.abs.z, sumsq=sums.z2))
}
