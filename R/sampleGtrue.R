sampleGtrue = function(y, N, particles, priorList){
 p = ncol(y)
#
 theta = dget('Output/Samples/theta.txt')
 GTrue = theta$G
 G = matrix(rep(as.numeric(GTrue), each=N), ncol=p^2, byrow=F)
 log.dG = rep(0, N)
 return(list(values=G, log.dq=log.dG))
}
