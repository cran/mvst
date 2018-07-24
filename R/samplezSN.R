samplezSN = function(y, X, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the latent variable Z, in the p-variate skew-N model.
 n = nrow(y)
 p = ncol(y)
 #
 yVec = as.numeric(t(y))
 psiVec = as.numeric(t(particles$psi))
 GVec = as.numeric(t(particles$G))
 # vVec = as.numeric(t(particles$v))
 values = as.double(rep(0.0, n*N))
 propdens = as.double(rep(0.0, N))
 #
 M = as.double(rep(0.0,n*N))
 RV = as.double(rep(0.0,n*N))
 absz = as.double(rep(0.0, n*N))
 #
 XFlag = !is.null(X)
 if(XFlag){
  k = ncol(X)
  XVec = as.numeric(t(X))
  BVec = as.numeric(particles$B)
  zPars = .C(C_rzSNX, as.double(absz), as.double(M), as.double(RV), as.integer(n), as.integer(p), as.integer(k), as.integer(N), as.double(yVec), as.double(XVec), as.double(BVec), as.double(psiVec), as.double(GVec), NAOK=TRUE)
 } else {
  xiVec = as.numeric(t(particles$xi))
  zPars = .C(C_rzSN, as.double(absz), as.double(M), as.double(RV), as.integer(n), as.integer(p), as.integer(N), as.double(yVec), as.double(xiVec), as.double(psiVec), as.double(GVec), NAOK=TRUE)
 }
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

 return(list(values=z, log.dq=log.dz))
}
