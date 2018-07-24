logPriorDensN = function(N, p, particles, priorList){
# xi = particles$xi
# p = ncol(xi)
 G = particles$G
#
 log.detSigma = numeric(N)
 for(iN in 1:N){
  G.iN = matrix(G[iN,], p, p)
  log.detSigma[iN] = as.numeric(determinant(G.iN, logarithm=T)$modulus)  # log(det(G.iN))
 }
 logPriorDens = # 0 +   			# xi or B
  (-(p+1)/2) * log.detSigma	# Sigma
 return(logPriorDens)
}
