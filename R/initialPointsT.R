initialPointsT = function(N, y, X, priorList){
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
 sqrt.v = sqrt(v)

 # vy-stats
 sum.v.y = v %*% y
 sum.v = apply(v,1,sum)

 if(is.null(X)){
  # xi & G
  xi = sum.v.y / sum.v
  G = matrix(0, N, p^2)
  for(iN in 1:N){
   e = y - matrix(rep(xi[iN,], n), ncol=p, byrow=T)
   sqrtv.e = matrix(rep(sqrt(v[iN,]),each=p), ncol=p, byrow=T) * e
   G.mat =  (n-1) * var(sqrtv.e) / n
   G[iN,] = as.vector(G.mat)
  }
  particles = list(xi=xi, G=G, nu=nu, v=v)
 } else {
  # B & G
  k = ncol(X)
  B = array(0, c(k, p, N))
  G = matrix(0, N, p^2)
  constFlag = identical(as.integer(X[,1]), as.integer(rep(1,n)))
  for(iN in 1:N){
   if(constFlag){
    lmFit = lm(y ~ X - 1)
   } else {
    lmFit = lm(y ~ X)
   }
   lmCoeff = lmFit$coefficients
   if(any(is.na(lmCoeff))) stop('initialPointsT.R: X is misspecified.\n')
   B.iN = matrix(lmCoeff, k, p)
   B[,,iN] = B.iN
   BX.iN = lmFit$fitted.values   # Conditional mean of y given X
   #
   e = y - BX.iN
   sqrtv.e = matrix(rep(sqrt(v[iN,]),each=p), ncol=p, byrow=T) * e
   G.mat =  (n-1) * var(sqrtv.e) / n
   G[iN,] = as.vector(G.mat)
  }
  particles = list(B=B, G=G, nu=nu, v=v)
 }

 return(particles)
}
