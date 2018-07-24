logPostDensN = function(y, X, N, particles, priorList, logPriorFunc){
 n = nrow(y)
 p = ncol(y)
 #
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
 #
 G = particles$G
 #
 # Rxi = matrix(0, n*N, p)
 RM = matrix(0, n*N, p) # xi or B'X
 Ry = matrix(0, n*N, p)
 XFlag = !is.null(X)
 if(XFlag){
  k = ncol(X)
  B = particles$B
  RX = X[rep(1:n, times=N),]
 } else {
  xi = particles$xi
 }
 for(icol in 1:p){
  Ry[,icol] = rep(y[,icol], N)
  if(XFlag){
   Bcol = matrix(B[,icol,], ncol=k, byrow=T)
   RB = Bcol[rep(1:N, each=n),]
   if(k == 1){
    RM[,icol] = RB*RX
   } else {
    RM[,icol] = apply(RB*RX,1,sum)
   }
  } else {
   RM[,icol] = rep(xi[,icol], each=n)
  }
 }
 #
 e = Ry - RM
 eList = vector('list', p)
 for(ilist in 1:p){
  eList[[ilist]] = matrix(e[,ilist], N, n, byrow=T)
 }
 LambdaMAT = matrix(0, N, n.pmat.indices.w)
 for(icol in 1:n.pmat.indices.w){
  LambdaMAT[,icol] = apply(eList[[pmat.indices.w[icol,1]]] * eList[[pmat.indices.w[icol,2]]], 1, sum)
 }
 PsiVEC = numeric(N)
 detG = numeric(N)
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
 loglikelihood = (- n/2) * log(detG) - 0.5 * PsiVEC +
  loglik.normalizing.constant # - (n/2) * log(2*pi)
 # Numerator of the unnormalized importance weights
 log.num.iw.b = log.prior + loglikelihood
 return(log.num.iw.b)
}
