samplevST = function(y, X, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the latent variable V, in the p-variate skew-t model.
 n = nrow(y)
 p = ncol(y)
 #
 yVec = as.numeric(t(y))
 psiVec = as.numeric(t(particles$psi))
 GVec = as.numeric(t(particles$G))
 nu = as.numeric(particles$nu)
 zVec = as.numeric(t(particles$z))
 values = rep(0, n*N)
 propdens = rep(0, N)
 #
 XFlag = !is.null(X)
 if(XFlag){
  k = ncol(X)
  XVec = as.numeric(t(X))
  BVec = as.numeric(particles$B)
  vList = .C(C_rvSTX, as.double(values), as.double(propdens), as.integer(n), as.integer(p), as.integer(k), as.integer(N), as.double(yVec), as.double(XVec), as.double(nu), as.double(BVec), as.double(psiVec), as.double(GVec), as.double(zVec), NAOK=TRUE)
 } else {
  xiVec = as.numeric(t(particles$xi))
  vList = .C(C_rvST, as.double(values), as.double(propdens), as.integer(n), as.integer(p), as.integer(N), as.double(yVec), as.double(nu), as.double(xiVec), as.double(psiVec), as.double(GVec), as.double(zVec), NAOK=TRUE)
 }
 values = matrix(vList[[1]], ncol=n, byrow=T)
 propdens = vList[[2]]

 return(list(values=values, log.dq=propdens))
}
