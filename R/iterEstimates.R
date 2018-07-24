iterEstimates = function(particles, iw, parInfo){
# This function computes summary statistics on the
#  basis of a single iteration.

 N = length(iw)
 parIndices = which(!(parInfo$names %in% c('z', 'v')))
 parNames = parInfo$names[parIndices]
 parTypes = parInfo$type[parIndices]
 nPars = length(parNames)
 estList = vector('list', 5)
 names(estList) = c('postMean', 'postMeanSE', 'postMedian', 'Q05', 'Q95')
 for(iEst in 1:length(estList)){
  estList[[iEst]] = vector('list', nPars)
  names(estList[[iEst]]) = parNames
  for(iPar in 1:nPars){
   theta = particles[[iPar]]
   parType = parTypes[iPar]
   if(parType == 'u'){
    thetaHat = sum(theta * iw)
    estList[['postMean']][[iPar]] = thetaHat
#    estList[['postMeanSE']][[iPar]] = sqrt(sum(iw^2 * (theta - thetaHat)^2))
    estList[['postMeanSE']][[iPar]] = sqrt(sum(iw * (theta - thetaHat)^2))
    estList[['postMedian']][[iPar]] = wQuantiles(theta, iw, 0.5)
    estList[['Q05']][[iPar]] = wQuantiles(theta, iw, 0.05)
    estList[['Q95']][[iPar]] = wQuantiles(theta, iw, 0.95)
   }
   #
   if(parType %in% (c('m', 'SM'))){
    nCols = ncol(theta)
    postMean = numeric(nCols)
    postMeanSE = numeric(nCols)
    postMedian = numeric(nCols)
    Q05 = numeric(nCols)
    Q95 = numeric(nCols)
    for(iCol in 1:nCols){
     thetaHat = sum(theta[,iCol] * iw)
     postMean[iCol] = thetaHat
#     postMeanSE[iCol] = sum(iw^2 * (theta[,iCol] - thetaHat)^2)
     postMeanSE[iCol] = sqrt(sum(iw * (theta[,iCol] - thetaHat)^2))
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
   #
   if(parType == 'M'){
    dims = dim(theta[,,1,drop=F])
    postMean = matrix(0, nrow=dims[1], ncol=dims[2])
    postMeanSE = matrix(0, nrow=dims[1], ncol=dims[2])
    postMedian = matrix(0, nrow=dims[1], ncol=dims[2])
    Q05 = matrix(0, nrow=dims[1], ncol=dims[2])
    Q95 = matrix(0, nrow=dims[1], ncol=dims[2])
    for(iRow in 1:dims[1]){
     for(iCol in 1:dims[2]){
      thetaHat = sum(theta[iRow,iCol,] * iw)
      postMean[iRow, iCol] = thetaHat
#      postMeanSE[iRow, iCol] = sum(iw^2 * (theta[iRow,iCol,] - thetaHat)^2)
      postMeanSE[iRow, iCol] = sqrt(sum(iw * (theta[iRow,iCol,] - thetaHat)^2))
      postMedian[iRow, iCol] = wQuantiles(theta[iRow,iCol,], iw, 0.5)
      Q05[iRow, iCol] = wQuantiles(theta[iRow,iCol,], iw, 0.05)
      Q95[iRow, iCol] = wQuantiles(theta[iRow,iCol,], iw, 0.95)
     }
    }
    estList[['postMean']][[iPar]] = postMean
    estList[['postMeanSE']][[iPar]] = postMeanSE
    estList[['postMedian']][[iPar]] = postMedian
    estList[['Q05']][[iPar]] = Q05
    estList[['Q95']][[iPar]] = Q95
   }
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
