initialPointsN = function(N, y, priorList){
# Given the arguments, initialPoints creates a population of particles from simple proposals.
 n = nrow(y)
 p = ncol(y)

 # theta
 xi = matrix(rep(apply(y, 2, mean),N), N, p, byrow=T)
 G = matrix(rep(as.numeric(var(y)),N), N, p^2, byrow=T)

 particles = list(xi=xi, G=G)
 return(particles)
}
