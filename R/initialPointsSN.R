initialPointsSN = function(N, y, X, priorList){
# Given the arguments, initialPoints creates a population of particles from simple proposals.
 n = nrow(y)
 p = ncol(y)

 # z
 z = matrix(rnorm(n*N), ncol=n)
 abs.z = abs(z)

 # zy-stats
 mean.y = apply(y, 2, mean)
 absz = abs(z)
 mean.absz = apply(absz, 1, mean)
 sum.absz = apply(absz, 1, sum)
 sum.absz.y = absz %*% y
 sum.y = matrix(rep(apply(y,2,sum), each=N), N, p)
 sum.z2 = apply(z^2,1,sum)

 # psi
 psi = (n * sum.absz.y - sum.absz * sum.y) / (sum.z2 * n - sum.absz^2)

 # psi
 # zmat = absz - matrix(rep(mean.absz, n), N, n, byrow=F)
 # ymat = y - matrix(rep(mean.y, n), n, p, byrow=T)
 # devz = apply(zmat^2, 1, sum)
 # psi = matrix(0, N, p)
 # for(iN in 1:N){
 #  for(ip in 1:p){
 #   psi[iN,ip] = sum(zmat[iN,] * ymat[,ip]) / devz[iN]
 #  }
 # }

 if(is.null(X)){
  # xi & G
  xi = matrix(rep(mean.y, N), N, p, byrow=T) - psi * matrix(rep(mean.absz, each=p), N, p, byrow=T)
  G = matrix(0, N, p^2)
  for(iN in 1:N){
   e = y - matrix(rep(xi[iN,], n), ncol=p, byrow=T) - matrix(rep(psi[iN,], n), ncol=p, byrow=T) * matrix(rep(abs(z[iN,]),each=p), ncol=p, byrow=T)
   Gsum = matrix(0, p, p)
   for(i in 1:n){
    Gsum = Gsum + e[i,] %*% t(e[i,])
   }
   G[iN,] = as.numeric(Gsum / n)
  }
  particles = list(xi=xi, G=G, psi=psi, z=z)
 } else {
  # B & G
  k = ncol(X)
  B = array(0, c(k, p, N))
  G = matrix(0, N, p^2)
  constFlag = identical(as.integer(X[,1]), as.integer(rep(1,n)))
  for(iN in 1:N){
   yStar = matrix(0, n, p)
   for(iCol in 1:p){
    yStar[,iCol] = y[,iCol] - rep(psi[iN,iCol], n) * abs.z[iN,]
   }
   if(constFlag){
    lmFit = lm(yStar ~ X - 1)
   } else {
    lmFit = lm(yStar ~ X)
   }
   lmCoeff = lmFit$coefficients
   if(any(is.na(lmCoeff))) stop('initialPointsSN.R: X is misspecified.\n')
   B.iN = matrix(lmCoeff, k, p)
   B[,,iN] = B.iN
   BX.iN = lmFit$fitted.values   # Conditional mean of y given X
   #
   e = yStar - BX.iN
   Gsum = matrix(0, p, p)
   for(i in 1:n){
    Gsum = Gsum + e[i,] %*% t(e[i,])
   }
   G[iN,] = as.numeric(Gsum / n)
  }

  particles = list(B=B, G=G, psi=psi, z=z)
 }

  return(particles)
}
