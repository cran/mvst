WUiterations = function(particles, N, Nnext, ti, Ti, y, X, modelType, priorList, parInfo, propFuncs, logPriorFunc, logPostFunc, saveParticles, outFolder, verbose){
# Given the arguments, this function performs a 'warm-up' MC step.

 n = nrow(y)
 p = ncol(y)
 nPars = nrow(parInfo)
 parNames = parInfo$names
 parTypes = parInfo$type

 # Proposal density matrix
 log.dQ = matrix(0, N, nPars)
 colnames(log.dQ) = parNames

 # Proposal step (with permuted sweeps)
 parPerm = 1:nPars # sample(1:nPars, size=nPars, replace=F)
 for(ipar in 1:nPars){
  if(verbose) cat(parNames[parPerm[ipar]],'\t')
  drawList = do.call(as.character(propFuncs[parNames[parPerm[ipar]]]), list(y=y, X=X, N=N, particles=particles, priorList=priorList))
  particles[[parNames[parPerm[ipar]]]] = drawList$values
  log.dQ[, parNames[parPerm[ipar]]] = drawList$log.dq
 }

 # Importance weights
 if(verbose) cat('Resampling')
 log.pitilde = do.call(logPostFunc, list(y=y, X=X, N=N, particles=particles, priorList=priorList, logPriorFunc=logPriorFunc))
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
 for(iPar in 1:nPars){
  parType = parInfo$type[iPar]
  if(parType == 'u'){
   particles[[iPar]] = particles[[iPar]][J]
  }
  if(parType %in% c('m','SM')){
   particles[[iPar]] = particles[[iPar]][J,,drop=F]
  }
  if(parType == 'M'){
   particles[[iPar]] = particles[[iPar]][,,J,drop=F]
  }
 }

 return(particles)
}
