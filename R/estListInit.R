estListInit = function(parInfo, Ti){
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
   estList[[iEst]][[iPar]] = switch(parTypes[iPar],
    u = rep(NA, Ti),
    m = matrix(NA, Ti, parInfo$nCols[iPar]),
    M = array(NA, c(parInfo$nCols[iPar], parInfo$nRows[iPar], Ti)),
    SM = matrix(NA, Ti, parInfo$nCols[iPar]^2)
   )
  }
 }
 return(estList)
}