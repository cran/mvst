logPriorDensSN = function(N, p, particles, priorList){
 particlesStar = reparamSkew(N=N, p=p, particles=particles)
 logPriorDens = # 0 +   							               # xi
  - nSphereVolume(p, 1, LOG=T) - 0.5 * particlesStar$log.detOmega +	# delta
  (-(p+1)/2) * particlesStar$log.detSigma +		# Sigma
  #	(-2) * apply(particlesStar$h,1,sum) +		# |J1|
  (-1) * apply(log(particlesStar$h), 1, sum)		# |J|
#  priorList$nulogpriors[nuIndices]			# nu
 logPriorDens[which(particlesStar$pos == 0)] = rep(-Inf, sum(particlesStar$pos == 0))
 return(logPriorDens)
}
