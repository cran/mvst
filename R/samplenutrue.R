samplenutrue = function(y, X, N, particles, priorList){
 theta = dget('Output/Samples/theta.txt')
 nuTrue = theta$nu
 nu = rep(nuTrue, N)
 log.dnu = rep(0, N)

 return(list(values=nu, log.dq=log.dnu))
}
