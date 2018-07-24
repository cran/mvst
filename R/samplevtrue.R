samplevtrue = function(y, X, N, particles, priorList){
 n = nrow(y)
 trueV = dget('Output/Samples/v.txt')
 v = matrix(rep(trueV, N), ncol=n, byrow=T)
 log.dv = rep(0, N)

 return(list(values=v, log.dq=log.dv))
}
