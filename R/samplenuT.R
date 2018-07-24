samplenuT = function(y, X, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the latent variable Nu, in the p-variate T model.
# Notice: the full conditional is the same in the T and ST models
 nuList = samplenuST(y, X, N, particles, priorList)

 return(list(values=nuList$values, log.dq=nuList$log.dq))
}
