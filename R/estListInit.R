estListInit = function(parTypes){
 estNames = setdiff(parTypes, c('z', 'v'))
 estList = vector('list', 5)
 names(estList) = c('postMean', 'postMeanSE', 'postMedian', 'Q05', 'Q95')
 for(iEst in 1:length(estList)){
  estList[[iEst]] = vector('list', length(estNames))
  names(estList[[iEst]]) = estNames
 }
 return(estList)
}