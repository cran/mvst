logPostDensT = function(y, X, N, particles, priorList, logPriorFunc){
 n = nrow(y)
 p = ncol(y)
 #
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
 #
 # psi = particles$psi
 G = particles$G
 nu = particles$nu
 v = particles$v
 #
 RM = matrix(0, n*N, p) # xi or B'X
 Ry = y[rep(1:n, times=N),]
 # Rpsi = matrix(0, n*N, p)
 XFlag = !is.null(X)
 if(XFlag){
  k = ncol(X)
  B = particles$B
  RX = X[rep(1:n, times=N),]
 } else {
  xi = particles$xi
 }
 for(icol in 1:p){
  # Ry[,icol] = rep(y[,icol], N)
  #  Rpsi[,icol] = rep(psi[,icol], each=n)
  if(XFlag){
   Bcol = matrix(B[,icol,], ncol=k, byrow=T)
   RB = Bcol[rep(1:N, each=n),]
   if(k == 1){
    RM[,icol] = RB*RX
   } else {
    RM[,icol] = apply(RB*RX,1,sum)
   }
   # RM[,icol] = apply(RB*RX,1,sum)
  } else {
   RM[,icol] = rep(xi[,icol], each=n)
  }
 }
 # abs.z = abs(as.numeric(t(particles$z)))
 # Rabsz = matrix(rep(abs.z, p), ncol=p)
 # sums.z2 = apply((particles$z)^2, 1, sum)
 # sqrt.v = sqrt(as.numeric(t(v)))
 # Rsqrtv = matrix(rep(sqrt.v, p), ncol=p)
 #
 e = Ry - RM
 eList = vector('list', p)
 for(ilist in 1:p){
  eList[[ilist]] = matrix(e[,ilist], N, n, byrow=T)
 }
 LambdaMAT = matrix(0, N, n.pmat.indices.w)
 for(icol in 1:n.pmat.indices.w){
  LambdaMAT[,icol] = apply(v * eList[[pmat.indices.w[icol,1]]] * eList[[pmat.indices.w[icol,2]]], 1, sum)
 }
 detG = numeric(N)
 PsiVEC = numeric(N)
 for(iN in 1:N){
  G.iN = matrix(G[iN,], p, p)
  detG[iN] = det(G.iN)
  invG.iN = solve(G.iN)
  Lambda.iN = matrix(0, p, p)
  Lambda.iN[pmat.indices.w] = LambdaMAT[iN,]
  Lambda.iN[pmat.indices.w[,2:1]] = LambdaMAT[iN,]
  PsiVEC[iN] = sum(t(invG.iN) * Lambda.iN)
 }
 # log-prior density
 log.prior = do.call(logPriorFunc, list(N=N, p=p, particles=particles, priorList=priorList))
 # loglikelihood
 loglik.normalizing.constant = - (n*p)/2 * log(2*pi) # (log-)constants in the likelihood of non-Skew models
 loglikelihood = (- n/2) * log(detG) + (p/2) * apply(log(v),1,sum) - 0.5 * PsiVEC +
#  (- 0.5) * sums.z2 +					# z
  n * nu/2 * log(nu/2) - n * lgamma(nu/2) + (nu/2 - 1) * apply(log(v), 1, sum) - (nu/2) * apply(v, 1, sum) +		# v
  + loglik.normalizing.constant # - (n/2) * log(2*pi)
 # Numerator of the unnormalized importance weights
 log.pi = log.prior + loglikelihood
 # log.pi[which(particlesStar$pos == 0)] = rep(-Inf, sum(particlesStar$pos == 0))
 return(log.pi)
}
