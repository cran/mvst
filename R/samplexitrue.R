samplexitrue = function(y, N, particles, priorList){
 p = ncol(y)
#
 theta = dget('Output/Samples/theta.txt')
 xiTrue = theta$xi
 xi = matrix(rep(xiTrue, each=N), ncol=p, byrow=F)
 log.dxi = rep(0, N)
 return(list(values=xi, log.dq=log.dxi))
}
