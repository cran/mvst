logPostDensN = function(y, N, particles, priorList, logPriorFunc){
 G = particles$G
 n = nrow(y)
 p = ncol(y)
 #
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
 Ry = matrix(0, n*N, p)
 Rxi = matrix(0, n*N, p)
# Rpsi = matrix(0, n*N, p)
 for(icol in 1:p){
  Ry[,icol] = rep(y[,icol], N)
  Rxi[,icol] = rep(particles$xi[,icol], each=n)
#  Rpsi[,icol] = rep(particles$psi[,icol], each=n)
 }
# abs.z = abs(as.numeric(t(particles$z)))
# Rabsz = matrix(rep(abs.z, p), ncol=p)
# sums.z2 = apply((particles$z)^2, 1, sum)
# sqrt.v = sqrt(as.numeric(t(particles$v)))
# Rsqrtv = matrix(rep(sqrt.v, p), ncol=p)
 #
 e = Ry - Rxi
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
  invG.iN = solve(matrix(particles$G[iN,], p, p))
  Lambda.iN = matrix(0, p, p)
  Lambda.iN[pmat.indices.w] = LambdaMAT[iN,]
  Lambda.iN[pmat.indices.w[,2:1]] = LambdaMAT[iN,]
  PsiVEC[iN] = sum(t(invG.iN) * Lambda.iN)
 }
 # log-prior density
 log.prior = do.call(logPriorFunc, list(N=N, particles=particles, priorList=priorList))
 # loglikelihood
 loglik.normalizing.constant = - (n*p)/2 * log(2*pi) # (log-)constants in the likelihood of non-Skew models
 loglikelihood = (- n/2) * log(detG) - 0.5 * PsiVEC +
  loglik.normalizing.constant # - (n/2) * log(2*pi)
 # Numerator of the unnormalized importance weights
 log.num.iw.b = log.prior + loglikelihood
 return(log.num.iw.b)
}
