samplevST = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the latent variable V, in the p-variate skew-t model.
 n = nrow(y)
 p = ncol(y)
 #
 yVec = as.numeric(t(y))
 xiVec = as.numeric(t(particles$xi))
 psiVec = as.numeric(t(particles$psi))
 GVec = as.numeric(t(particles$G))
# invGVec = as.numeric(t(invG))
 nu = as.numeric(particles$nu)
 zVec = as.numeric(t(particles$z))
 values = rep(0, n*N)
 propdens = rep(0, N)
 	# TEMP
 #	lista = list(y=y, N=N, particles=particles, priorList=priorList, n=as.integer(n), p=as.integer(p), yVec=as.double(yVec), nu=as.double(nu), xiVec=as.double(xiVec), psiVec=as.double(psiVec), GVec=as.double(GVec), zVec=as.double(zVec))
 #	saveRDS(lista, 'listaFromsamplevST.RDS')
 #
 vList = .C('rvST', as.double(values), as.double(propdens), as.integer(n), as.integer(p), as.integer(N), as.double(yVec), as.double(nu), as.double(xiVec), as.double(psiVec), as.double(GVec), as.double(zVec), NAOK=TRUE)
 values = matrix(vList[[1]], ncol=n, byrow=T)
 propdens = vList[[2]]
# mean.v = apply(values, 1, mean)
# sqrt.v = as.numeric(t(sqrt(values)))		# sqrt(v), as a vector
# mean.yv = (values %*% y) / n
 #
 return(list(values=values, log.dq=propdens)) #, mv=mean.v, myv=mean.yv, sqroot=sqrt.v))
}
