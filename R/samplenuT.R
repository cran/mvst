samplenuT = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the latent variable Nu, in the p-variate T model.
# Notice: the full conditional is the same in the T and ST models
 res = samplenuST(y, N, particles, priorList)
 return(list(values=res$values, log.dq=res$log.dq))
}
