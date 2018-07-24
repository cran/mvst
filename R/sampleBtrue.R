sampleBtrue = function(y, X, N, particles, priorList){
 p = ncol(y)
 k = ncol(X)
 #
 theta = dget('Output/Samples/theta.txt')
 BTrue = theta$B
 B = array(NA, c(k, p, N))
 for(iN in 1:N) B[,,iN] = BTrue
 log.dB = rep(0, N)

 return(list(values=B, log.dq=log.dB))
}
