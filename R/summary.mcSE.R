summary.mcSE = function(object, ...){
 # Model's parameters info
 parInfo = object$parInfo

 # Iterations' weights
 iterWeights = object$perplexity
 iter.w = iterWeights / sum(iterWeights)

 #
 sampleEstimates = vector('list', length(object$estlist))
 names(sampleEstimates) = names(object$estlist)
 # sampleEstimates = object$estList

 parIndices = which(!(parInfo$names %in% c('z', 'v', 'parInfo')))
 parNames = parInfo$names[parIndices]
 parTypes = parInfo$type[parIndices]
 nPars = length(parNames)
 postSummaries = object$estlist # [1:5]
 for(iEst in 1:length(postSummaries)){
  postSummary = postSummaries[[iEst]]
  sampleEstimates[[iEst]] = vector('list', nPars)
  names(sampleEstimates[[iEst]]) = names(object$estlist[[iEst]])
  for(iPar in 1:nPars){
   parType = parTypes[iPar]
   values = postSummary[[iPar]]
   if(parType == 'u') sampleEstimates[[iEst]][[iPar]] = sum(object$estlist[[iEst]][[iPar]] * iter.w)
   if(parType %in% c('m', 'SM')){
    nCols = ncol(values)
    estVec = numeric(nCols)
    for(iCol in 1:nCols){
     estVec[iCol] = sum(values[,iCol] * iter.w)
    }
    sampleEstimates[[iEst]][[iPar]] = estVec   
   }
   if(parType == 'SM'){
    sampleEstimates[[iEst]][[iPar]] = matrix(sampleEstimates[[iEst]][[iPar]], object$p, object$p)
   }
   if(parType == 'M'){
    estVec = matrix(NA, dim(values)[1], dim(values)[2])
    for(iRow in 1:(dim(values)[1])){
     for(iCol in 1:(dim(values)[2])){
      estVec[iRow, iCol] = sum(values[iRow,iCol,] * iter.w)
     }
    }
    sampleEstimates[[iEst]][[iPar]] = estVec   
   }
   names(sampleEstimates[[iEst]][[iPar]]) = names(object$estlist[[iEst]][[iPar]])
  }
 }

 # Marginal likelihood
 log.py = object$log.py
 K = max(log.py)
 log.margLike = K + log(sum(exp(log.py-K)*iter.w))
 sampleEstimates$log.py = log.margLike
 
 # number of resampled particles
 sampleEstimates$nResampled = object$nResampled
 
 # perplexity
 sampleEstimates$perplexity = object$perplexity

 ans = list(sampleEstimates=sampleEstimates, log.margLike=log.margLike)
 class(ans) = 'mcSEsummary'
 print.summary.mcSE(ans)
 invisible(ans)
}
