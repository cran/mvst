sampleztrue = function(y, X, N, particles, priorList){
 n = nrow(y)
 trueZ = dget('Output/Samples/z.txt')
 z = matrix(rep(trueZ,N), ncol=n, byrow=T)
 log.dz = rep(0, N)
 #
 return(list(values=z, log.dq=log.dz))
}
