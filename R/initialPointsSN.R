initialPointsSN = function(N, y, priorList){
# Given the arguments, initialPoints creates a population of particles from simple proposals.
 n = nrow(y)
 p = ncol(y)

  # z
 z = matrix(rnorm(n*N), ncol=n)

 # zy-stats
 mean.y = apply(y, 2, mean)
 absz = abs(z)
 mean.absz = apply(absz, 1, mean)

 # theta
 zmat = absz - matrix(rep(mean.absz, n), N, n, byrow=F)
 ymat = y - matrix(rep(mean.y, n), n, p, byrow=T)
 devz = apply(zmat^2, 1, sum)
 psi = matrix(0, N, p)
 for(iN in 1:N){
  for(ip in 1:p){
   psi[iN,ip] = sum(zmat[iN,] * ymat[,ip]) / devz[iN]
  }
 }

 xi = matrix(rep(mean.y, N), N, p, byrow=T) - psi * matrix(rep(mean.absz, each=p), N, p, byrow=T)

 G = matrix(0, N, p^2)
 for(iN in 1:N){
  e = matrix(0, n, p)
  e = y - matrix(rep(xi[iN,], n), ncol=p, byrow=T) - matrix(rep(psi[iN,], n), ncol=p, byrow=T) * matrix(rep(abs(z[iN,]),each=p), ncol=p, byrow=T)
  Gsum = matrix(0, p, p)
  for(i in 1:n){
   Gsum = Gsum + e[i,] %*% t(e[i,])
  }
  G[iN,] = as.numeric(Gsum / n)
 }

 particles = list(xi=xi, G=G, psi=psi, z=z)
 return(particles)
}
