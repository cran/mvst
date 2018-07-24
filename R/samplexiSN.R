samplexiSN = function(y, X, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable Xi, in the p-variate skew-N model.
 n = nrow(y)
 p = ncol(y)
 #
 pmat.indices.wVEC = triangleIndices(p, side='u', dgn=T, dataframe=F)
 pmat.indices.wDF = triangleIndices(p, side='u', dgn=T, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2

 G = particles$G
 psi = particles$psi
 z = particles$z
 #
 mean.y = apply(y, 2, mean)
 mean.absz = apply(abs(z), 1, mean)
 mxi = matrix(0, N, p)
 for(icol in 1:p){
  mxi[,icol] = mean.y[icol] - psi[,icol] * mean.absz
 }
 vxi = matrix(0, N, n.pmat.indices.w)
 for(icol in 1:n.pmat.indices.w){
  vxi[,icol] = G[,pmat.indices.wVEC[icol]] / n
 }
 xi = matrix(0, N, p)
 log.dxi = numeric(N)
 for(iN in 1:N){
  vxi.iN = matrix(0, p, p)
  vxi.iN[pmat.indices.wDF] = vxi[iN,]
  vxi.iN[pmat.indices.wDF[,2:1]] = vxi[iN,]
  xi.iN = rmnorm(1, mxi[iN,], vxi.iN)
  xi[iN,] = xi.iN
  log.dxi[iN] = dmnorm(as.numeric(xi.iN), mxi[iN,], vxi.iN, log=TRUE)
 }
 #
 return(list(values=xi, log.dq=log.dxi))
}
