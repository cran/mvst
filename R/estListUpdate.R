estListUpdate = function(oldEst, newEst, ti){
 estimates = oldEst
 for(iEst in 1:length(estimates)){
  estimate = estimates[[iEst]]
  for(iPar in 1:length(estimate)){
   parEstimate = estimate[[iPar]]
   newEstimate = newEst[[iEst]][[iPar]]
   if(length(newEstimate) == 1){
    estimates[[iEst]][[iPar]] = c(parEstimate, as.numeric(newEstimate))
   } else {
    estimates[[iEst]][[iPar]] = rbind(parEstimate, newEstimate)
    rownames(estimates[[iEst]][[iPar]]) = NULL
   }
  }
 }
 return(estimates)
}
