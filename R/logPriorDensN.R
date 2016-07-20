logPriorDensN = function(N, particles, priorList){
 xi = particles$xi
 G = particles$G
 p = ncol(xi)
#
 log.detSigma = numeric(N)
 for(iN in 1:N){
  G.iN = matrix(G[iN,], p, p)
  log.detSigma[iN] = log(det(G.iN))
 }
 logPriorDens = # 0 +   					# xi
  (-(p+1)/2) * log.detSigma	# Sigma
 return(logPriorDens)
}
