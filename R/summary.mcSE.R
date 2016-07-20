summary.mcSE = function(estimates){

 # Iterations' weights
 iterWeights = estimates$perplexity
 iter.w = iterWeights / sum(iterWeights)

 # Marginal likelihood
 log.py = estimates$log.py
 K = max(log.py)
 log.margLike = K + log(sum(exp(log.py-K)*iter.w))

 #
 sampleEstimates = estimates$estList
 for(iEst in 1:length(estimates)){
  estimate = sampleEstimates[[iEst]]
  for(iPar in 1:length(estimate)){
   parEstimate = estimate[[iPar]]
   if(is.null(dim(parEstimate))){
    sampleEstimates[[iEst]][[iPar]] = sum(sampleEstimates[[iEst]][[iPar]] * iter.w)
   } else {
    estVec = numeric(ncol(parEstimate))
    for(iCol in 1:ncol(parEstimate)){
    	estVec[iCol] = sum(sampleEstimates[[iEst]][[iPar]][,iCol] * iter.w)
    }
    sampleEstimates[[iEst]][[iPar]] = estVec
   }
  }
 }

 # T marginal likelihood
 #	K = max(log.norm.cost)	# Costante usata per il calcolo della marginal likelihood
 #	log.py = K + log(sum(exp(log.norm.cost-K)*iter.w))

 # py = sum(norm.cost*iter.w)

 # py = sum(norm.cost*n.ricampionati) / sum(n.ricampionati)
 # cat(c(py,'\n'), file='Output/Estimates/pyST.txt', append=T)

 # # SN marginal likelihood
 # S = cov(y)*(n-1)/n
 # log.psi = 0
 # for(j in 1:p){
 #  log.psi = log.psi+lgamma((n-1)/2-(j-1)/2)
 # }
 # log.m0 = ((3-2*n)/2)*log(pi)+log.psi-n*log(n)-((n-1)/2)*log(det(S))
 # m0 = exp(log.m0)
 # cat(c(py,m0,'\n'), file='Output/Estimates/py.txt', append=T)
 # # BF = sum(n.ricampionati*exp(log.norm.cost-log.m0))/sum(n.ricampionati)	# Bayes factor
 # BF = sum(exp(log.norm.cost-log.m0)*iter.w)						# Bayes factor
 # cat(BF,'\n', file='Output/Estimates/BF.txt', append=T)

 # sampleEstimates = list(postMeans=post.means, otherMeans=other.means, postMode=estimates$post.mode[nrow(estimates$post.mode),], H=estimates$H, n.ricampionati=n.ricampionati, log.py=log.py)
 ans = list(sampleEstimates=sampleEstimates, log.margLike=log.margLike)
 class(ans) = 'summary.mcSE'
 return(ans)
}
