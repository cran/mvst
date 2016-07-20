initialPointsT = function(N, y, priorList){
# Given the arguments, initialPoints creates a population of particles from simple proposals.
 n = nrow(y)
 p = ncol(y)

  # nu
 nu = sample(priorList$gdls, size=N, prob=exp(priorList$nulogpriors), replace=T)
  # v
 v = matrix(0, N, n)
 for(iN in 1:N){
  v[iN,] = rgamma(n, nu[iN]/2, nu[iN]/2)
#  v[iN,] = rchisq(n, nu[iN]) / nu[iN]
 }
  # z
 z = matrix(rnorm(n*N), ncol=n)

#vy-stats
# sqrtv.absz = sqrt(v) * abs(z)
# sum.sqrtv.absz = apply(sqrtv.absz, 1, sum)
# sum.sqrtv.absz.y = sqrtv.absz %*% y
sum.v.y = v %*% y
# sum.z2 = apply(z^2,1,sum)
sum.v = apply(v,1,sum)

# psi = (sum.v * sum.sqrtv.absz.y - sum.sqrtv.absz * sum.v.y) / (sum.z2 * sum.v - sum.sqrtv.absz^2)

xi = sum.v.y / sum.v
G = matrix(0, N, p^2)
 for(iN in 1:N){
  e = matrix(0, n, p)
  e = y - matrix(rep(xi[iN,], n), ncol=p, byrow=T)
  sqrtv.e = matrix(rep(sqrt(v[iN,]),each=p), ncol=p, byrow=T) * e
  G.mat =  (n-1) * var(sqrtv.e) / n
  G[iN,] = as.vector(G.mat)
 }

 particles = list(xi=xi, G=G, nu=nu, v=v)
 return(particles)
}
