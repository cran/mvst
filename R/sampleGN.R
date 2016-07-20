sampleGN = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable G, in the p-variate Normal model.
 n = nrow(y)
 p = ncol(y)
 xi = particles$xi
#
 ybar = apply(y, 2, mean)
 G = matrix(0, N, p^2)
 log.dG = numeric(N)
 for(iN in 1:N){
  Lambda = n * ((n-1) * cov(y)/n + (ybar - xi[iN,]) %*% t(ybar - xi[iN,]))
  G.iN = riwish(v=(n+priorList$m), S=(Lambda+priorList$W))
  G[iN,] = as.numeric(G.iN)
  log.dG[iN] = diwish.mia3(G.iN, (n+priorList$m), (Lambda+priorList$W), LOG=TRUE)
 }
 return(list(values=G, log.dq=log.dG))
}
