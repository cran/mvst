WUiterations = function(particles, N, Nnext, ti, Ti, y, modelType, priorList, parTypes, propFuncs, logPriorFunc, logPostFunc, saveParticles, outFolder, verbose){# parNames, otherNames,
# Given the arguments, this function performs a 'warm-up' MC step.

 n = nrow(y)
 p = ncol(y)

 # Proposal density matrix
 log.dQ = matrix(0, N, length(parTypes))
 colnames(log.dQ) = parTypes

 # Proposal step (with permuted sweeps)
 parPerm = sample(1:length(parTypes), size=length(parTypes), replace=F)
 for(ipar in 1:length(parTypes)){
  if(verbose) cat(parTypes[parPerm[ipar]],'\t')
  drawList = do.call(as.character(propFuncs[parTypes[parPerm[ipar]]]), list(y=y, N=N, particles=particles, priorList=priorList))
  particles[[parTypes[parPerm[ipar]]]] = drawList$values
  log.dQ[, parTypes[parPerm[ipar]]] = drawList$log.dq
 }

 # Importance weights
 log.pitilde = do.call(logPostFunc, list(y=y, N=N, particles=particles, priorList=priorList, logPriorFunc=logPriorFunc))
 log.q = apply(log.dQ, 1, sum)
 log.iw = log.pitilde - log.q
 K = max(log.iw)	# Constant needed for the computation of the marginal likelihood
 log.py.ti = K + log(sum(exp(log.iw-K))) - log(N)

 iw.bar = exp(log.iw)
 log.iw.b = log.iw - log(N) - log.py.ti
 iw.b = exp(log.iw.b) # unnormalized weights
 iw = iw.b/sum(iw.b)  # normalized weights

# perplexity.ti = exp(-sum(iw * log(iw))) / N
# iterEst = iterEstimates(particles, iw, parTypes)

 J = sample(1:N, size=Nnext, prob=iw, replace=T)
 nResampled.ti = length(unique(J))

 # Resampling step
 parClasses = lapply(particles, class)
 if(any((parClasses != 'numeric') & (parClasses != 'matrix'))) stop('Wrong class for one or more object(s) in particles.\n')
 for(iPar in which(as.character(parClasses) == 'numeric')){
  particles[[iPar]] = particles[[iPar]][J]
 }
 for(iPar in which(as.character(parClasses) == 'matrix')){
  particles[[iPar]] = particles[[iPar]][J,]
 }

 return(particles)
# return(list(particles=particles, iterEst=iterEst, log.py.ti=log.py.ti, nResampled.ti=nResampled.ti, perplexity.ti=perplexity.ti))
}
