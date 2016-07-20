logPriorDensST = function(N, particles, priorList){
 p = ncol(particles$xi)
 #
 nu = particles$nu
 particlesStar = reparamSkew(N=N, particles=particles)
 #
 gdls.list = priorList$gdls
 nuIndices = numeric(N)
 for(iN in 1:N){
  nuIndices[iN] = which(gdls.list == nu[iN])
 }
 logPriorDens = # 0 +   						# xi
  - nSphereVolume(p, 1, LOG=T) - 0.5 * particlesStar$log.detOmega +	# delta
  (-(p+1)/2) * particlesStar$log.detSigma +		# Sigma
  #	(-2) * apply(particlesStar$h,1,sum) +		# |J1|
  (-1) * apply(log(particlesStar$h), 1, sum) +	# |J|
  priorList$nulogpriors[nuIndices]				# nu
 logPriorDens[which(particlesStar$pos == 0)] = rep(-Inf, sum(particlesStar$pos == 0))
 return(logPriorDens)
}
