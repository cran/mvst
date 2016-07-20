iterEstimates = function(particles, iw, parTypes){# ti, particlesStar, iw.bar, otherNames, log.pitilde, log.post.mode.value, parNames
# iterEstimates = function(ti, particles, particlesStar, iw, iw.bar, log.pitilde, log.post.mode.value, parNames, otherNames){

# for(iItem in 1:length(particles)) assign(names(particles)[iItem], particles[[iItem]])

 N = length(iw)
 estNames = setdiff(parTypes, c('z', 'v'))
 estList = estListInit(parTypes)
 for(iPar in 1:length(estNames)){
  theta = particles[[iPar]]
  if(is.null(dim(theta))){
   thetaHat = sum(theta * iw)
   estList[['postMean']][[iPar]] = thetaHat
   estList[['postMeanSE']][[iPar]] = sqrt(sum(iw^2 * (theta - thetaHat)^2))
   estList[['postMedian']][[iPar]] = wQuantiles(theta, iw, 0.5)
   estList[['Q05']][[iPar]] = wQuantiles(theta, iw, 0.05)
   estList[['Q95']][[iPar]] = wQuantiles(theta, iw, 0.95)
  } else {
   nCols = ncol(theta)
   postMean = numeric(nCols)
   postMeanSE = numeric(nCols)
   postMedian = numeric(nCols)
   Q05 = numeric(nCols)
   Q95 = numeric(nCols)
   for(iCol in 1:nCols){
    thetaHat = sum(theta[,iCol] * iw)
    postMean[iCol] = thetaHat
    postMeanSE[iCol] = sum(iw^2 * (theta[,iCol] - thetaHat)^2)
    postMedian[iCol] = wQuantiles(theta[,iCol], iw, 0.5)
    Q05[iCol] = wQuantiles(theta[,iCol], iw, 0.05)
    Q95[iCol] = wQuantiles(theta[,iCol], iw, 0.95)
   }
   estList[['postMean']][[iPar]] = postMean
   estList[['postMeanSE']][[iPar]] = postMeanSE
   estList[['postMedian']][[iPar]] = postMedian
   estList[['Q05']][[iPar]] = Q05
   estList[['Q95']][[iPar]] = Q95
  }
 }

 # MODE
#  max.log.pitilde = which.max(log.pitilde)
#  if(log.pitilde[max.log.pitilde] > log.post.mode.value){
#   lista$post.mode = c(xi[max.log.pitilde,], psi[max.log.pitilde,],
#   matrix(G[max.log.pitilde,],p,p)[pmat.indices.w], nu[max.log.pitilde])
#   lista$log.post.mode.value = max(log.pitilde)
#   lista$updated.post.mode = TRUE
#  }
#  cat('\n',log.pitilde[max.log.pitilde], log.post.mode.value, lista$updated.post.mode,'\n')
 return(estList)
}
