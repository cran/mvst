finalEst = function(estimates, iterWeights, log.py){

 iter.w = iterWeights / sum(iterWeights)
  # Marginal likelihood
 K = max(log.py)
 log.margLike = K + log(sum(exp(log.py-K)*iter.w))

 #
 sampleEstimates = estimates
 for(iEst in 1:length(estimates)){
  estimate = estimates[[iEst]]
  for(iPar in 1:length(estimate)){
   parEstimate = estimate[[iPar]]
   if(is.null(dim(parEstimate))){
    sampleEstimates[[iEst]][[iPar]] = sum(estimates[[iEst]][[iPar]] * iter.w)
   } else {
    estVec = numeric(ncol(parEstimate))
    for(iCol in 1:ncol(parEstimate)){
    	estVec[iCol] = sum(estimates[[iEst]][[iPar]][,iCol] * iter.w)
    }
    sampleEstimates[[iEst]][[iPar]] = estVec
   }
  }
 }

 return(list(sampleEstimates=sampleEstimates, log.margLike=log.margLike))
}
