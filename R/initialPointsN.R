initialPointsN = function(N, y, X, priorList){
# Given the arguments, initialPoints creates a population of particles from simple proposals.
 n = nrow(y)
 p = ncol(y)
 # (xi | B) & G
 if(is.null(X)){
  xi = matrix(rep(apply(y, 2, mean),N), N, p, byrow=T)
  G = matrix(rep(as.numeric(var(y)),N), N, p^2, byrow=T)
  particles = list(xi=xi, G=G)
 } else {
  k = ncol(X)
  B = array(0, c(k, p, N))
  G = matrix(0, N, p^2)
  constFlag = identical(as.integer(X[,1]), as.integer(rep(1,n)))
  if(constFlag){
   lmFit = lm(y ~ X - 1)
  } else {
   lmFit = lm(y ~ X)
  }
  lmCoeff = lmFit$coefficients
  if(any(is.na(lmCoeff))) stop('initialPointsN.R: X is misspecified.\n')
  B.iN = matrix(lmCoeff, k, p)
  BX.iN = lmFit$fitted.values
  e = y - BX.iN
  for(iN in 1:N){
   B[,,iN] = B.iN
   G.mat = (n-1) * var(e) / n
   G[iN,] = as.vector(G.mat)
  }
  particles = list(B=B, G=G)
 }

 return(particles)
}
