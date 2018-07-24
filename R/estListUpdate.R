estListUpdate = function(oldEst, newEst, parInfo, ti){
 parIndices = which(!(parInfo$names %in% c('z', 'v')))
 parNames = parInfo$names[parIndices]
 parTypes = parInfo$type[parIndices]
 nPars = length(parNames)
 for(iEst in 1:length(oldEst)){
  for(iPar in 1:nPars){
   if(parTypes[iPar] == 'u'){
    oldEst[[iEst]][[iPar]][ti] = newEst[[iEst]][[iPar]]
   }
   if(parTypes[iPar] %in% c('m','SM')){
    oldEst[[iEst]][[iPar]][ti,] = newEst[[iEst]][[iPar]]
   }
   if(parTypes[iPar] == 'M'){
    oldEst[[iEst]][[iPar]][,,ti] = newEst[[iEst]][[iPar]]
   }
  }
 }
 return(oldEst)
}
