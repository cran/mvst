samplexiT = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable Xi, in the p-variate skew-t model.
 n = nrow(y)
 p = ncol(y)
#
 pmat.indices.wVEC = triangleIndices(p, side='u', dgn=T, dataframe=F)
 pmat.indices.wDF = triangleIndices(p, side='u', dgn=T, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2

 G = particles$G
# psi = particles$psi
 v = particles$v
 mean.v = apply(v, 1, mean)
 mean.yv = (v %*% y) / n
 
# mean.absz.sqrtv = apply(abs(particles$z) * sqrt(v), 1, mean)
 mxi = matrix(0, N, p)
 for(icol in 1:p){
  mxi[,icol] = mean.yv[,icol] / mean.v
 }
 vxi = matrix(0, N, n.pmat.indices.w)
 for(icol in 1:n.pmat.indices.w){
  vxi[,icol] = G[,pmat.indices.wVEC[icol]] / (n * mean.v)
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
# Rxi = matrix(0, n*N, p)
# for(icol in 1:p){
#  Rxi[,icol] = rep(xi[,icol], each=n)
# }
 return(list(values=xi, log.dq=log.dxi))# , repxi=Rxi))
}
