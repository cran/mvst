iterations = function(particles, N, ti, Ti, y, X, modelType, priorList, parInfo, propFuncs, logPriorFunc, logPostFunc, saveParticles, outFolder, verbose){
# Given the arguments, this function performs an MC step to update parameters, computes different quantities for the step and writes them in external files.

 n = nrow(y)
 p = ncol(y)
 nPars = nrow(parInfo)
 parNames = parInfo$names

 # Proposal density matrix
 log.dQ = matrix(0, N, nPars)
 colnames(log.dQ) = parNames

 # Proposal step (without permuted sweeps)
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
#	log.pitilde[which(is.infinite(log.dQ[,'z'] + log.dQ[,'v']))] = rep(0, sum((is.infinite(log.dQ[,'z'] + log.dQ[,'v']))))		# In realta' bisognerebbe modificare logPostDensT.R e logPostDensST.R
 log.q = apply(log.dQ, 1, sum)
 log.iw = log.pitilde - log.q
 K = max(log.iw)	# Constant needed for the computation of the marginal likelihood
 log.py.ti = K + log(sum(exp(log.iw-K))) - log(N)
 #
 # norm.cost[ti] = sum(exp(log.pitilde - log.q))/N

 # iw.bar = exp(log.iw)
 log.iw.b = log.iw - log(N) - log.py.ti
 iw.b = exp(log.iw.b) # unnormalized weights
 iw = iw.b/sum(iw.b)  # normalized weights

 perplexity.ti = exp(-sum((iw * log.iw.b)[which(log.iw.b > -Inf)])) / N
# perplexity.ti = exp(-sum(iw * log.iw.b)) / N
 iterEst = iterEstimates(particles, iw, parInfo)

 J = sample(1:N, size=N, prob=iw, replace=T)
 nResampled.ti = length(unique(J))

 # Save the proposed particles & the resampling vector
 if(saveParticles){
  saveRDS(particles, paste(outFolder, '/Proposed', ti, '.rds', sep=''))
  saveRDS(J, paste(outFolder, '/Resampled', ti, '.rds', sep=''))
#  saveRDS(iw, paste(outFolder, '/w', ti, '.rds', sep=''))
  saveRDS(log.iw.b, paste(outFolder, '/log.zeta', ti, '.rds', sep=''))
#  saveRDS(log.iw, paste(outFolder, '/log.w', ti, '.rds', sep=''))
 }

 # Resampling step
 # parClasses = lapply(particles, class)
 for(iPar in which(parInfo$type == 'u')){
  particles[[iPar]] = particles[[iPar]][J]
 }
 for(iPar in which(parInfo$type %in% c('m','SM'))){
  particles[[iPar]] = particles[[iPar]][J,,drop=F]
 }
 for(iPar in which(parInfo$type == 'M')){
  particles[[iPar]] = particles[[iPar]][,,J,drop=F]
 }

 return(list(particles=particles, iterEst=iterEst, log.py.ti=log.py.ti, nResampled.ti=nResampled.ti, perplexity.ti=perplexity.ti))
}
