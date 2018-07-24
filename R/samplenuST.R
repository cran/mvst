samplenuST = function(y, X, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the latent variable Nu, in the p-variate skew-t model.
 n = nrow(y)
#
 v = particles$v
 rm(particles)
 #
 gdls = priorList$gdls
 nulogpriors = priorList$nulogpriors
 sumlogv = apply(log(v), 1, sum)
 sumv = apply(v, 1, sum)
 logProbs = matrix(0, N, length(gdls))
 for(igdl in 1:length(gdls)){
  gdl = gdls[igdl]
  logProbs[,igdl] = nulogpriors[igdl] + n * gdl * log(gdl/2) / 2 - n * lgamma(gdl/2) +
   (gdl/2-1) * sumlogv - gdl * sumv / 2
 }
 values = propdens = nuIndices = numeric(N)
 for(inu in 1:N){
  logprobs = logProbs[inu,]
  nuprobs = exp(logprobs) / sum(exp(logprobs))
  nuIndices[inu] = sample(1:length(gdls), size=1, prob=nuprobs)
  values[inu] = gdls[nuIndices[inu]]
  propdens[inu] = log(nuprobs[nuIndices[inu]])
 }

 return(list(values=values, log.dq=propdens))
}
