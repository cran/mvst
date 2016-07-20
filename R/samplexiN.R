samplexiN = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable Xi, in the p-variate Normal model.
 n = nrow(y)
 p = ncol(y)
 G = particles$G
 #
 ybar = apply(y, 2, mean)
 xi = matrix(0, N, p)
 log.dxi = numeric(N)
 for(iN in 1:N){
  temp = matrix(G[iN,], p, p)/n
  vxi.iN = 0.5 * (temp + t(temp))
  xi.iN = rmnorm(1, ybar, vxi.iN)
  xi[iN,] = xi.iN
  log.dxi[iN] = dmnorm(as.numeric(xi.iN), ybar, vxi.iN, log=TRUE)
 }
 return(list(values=xi, log.dq=log.dxi))
}
