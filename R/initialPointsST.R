initialPointsST = function(N, y, X, priorList){
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
 sqrtv = sqrt(v)

  # z
 z = matrix(rnorm(n*N), ncol=n)
 abs.z = abs(z)

 # vzy-stats
 sqrtv.absz = sqrtv * abs.z
 sum.sqrtv.absz = apply(sqrtv.absz, 1, sum)
 sum.sqrtv.absz.y = sqrtv.absz %*% y
 sum.v.y = v %*% y
 sum.z2 = apply(z^2,1,sum)
 sum.v = apply(v,1,sum)

 # psi
 psi = (sum.v * sum.sqrtv.absz.y - sum.sqrtv.absz * sum.v.y) / (sum.z2 * sum.v - sum.sqrtv.absz^2)

 if(is.null(X)){
  # xi & G
  xi = (sum.v.y - psi * matrix(rep(sum.sqrtv.absz, p), N, p, byrow=F)) / sum.v
  G = matrix(0, N, p^2)
  for(iN in 1:N){
   e = matrix(0, n, p)
   e = y - matrix(rep(xi[iN,], n), ncol=p, byrow=T) - matrix(rep(psi[iN,], n), ncol=p, byrow=T) * matrix(rep(abs(z[iN,]),each=p), ncol=p, byrow=T) / matrix(rep(sqrt(v[iN,]),each=p), ncol=p, byrow=T)
   sqrtv.e = matrix(rep(sqrt(v[iN,]),each=p), ncol=p, byrow=T) * e
   G.mat = (n-1) * var(sqrtv.e) / n
   G[iN,] = as.vector(G.mat)
  }
  particles = list(xi=xi, G=G, psi=psi, nu=nu, z=z, v=v)
 } else {
  # B & G
  k = ncol(X)
  B = array(0, c(k, p, N))
  G = matrix(0, N, p^2)
  constFlag = identical(as.integer(X[,1]), as.integer(rep(1,n)))
  for(iN in 1:N){
   yStar = matrix(0, n, p)
   for(iCol in 1:p){
    yStar[,iCol] = y[,iCol] - rep(psi[iN,iCol], n) * abs.z[iN,] / sqrtv[iN,]
   }
   if(constFlag){
    lmFit = lm(yStar ~ X - 1)
   } else {
    lmFit = lm(yStar ~ X)
   }
   lmCoeff = lmFit$coefficients
   if(any(is.na(lmCoeff))) stop('initialPointsST.R: X is misspecified.\n')
   B.iN = matrix(lmCoeff, k, p)
   B[,,iN] = B.iN
   BX.iN = lmFit$fitted.values   # Conditional mean of y given X
   #
   e = yStar - BX.iN
   sqrtv.e = matrix(rep(sqrt(v[iN,]),each=p), ncol=p, byrow=T) * e
   G.mat = (n-1) * var(sqrtv.e) / n
   G[iN,] = as.vector(G.mat)
  }
  particles = list(B=B, G=G, psi=psi, nu=nu, z=z, v=v)
 }

 return(particles)
}
