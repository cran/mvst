samplepsitrue = function(y, X, N, particles, priorList){
 p = ncol(y)
 #
 theta = dget('Output/Samples/theta.txt')
 psiTrue = theta$psi
 psi = matrix(rep(psiTrue, each=N), ncol=p, byrow=F)
 log.dpsi = rep(0, N)

 return(list(values=psi, log.dq=log.dpsi))
}
