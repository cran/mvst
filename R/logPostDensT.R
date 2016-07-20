logPostDensT = function(y, N, particles, priorList, logPriorFunc){
 n = nrow(y)
 p = ncol(y)
#
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
 #
 G = particles$G
 nu = particles$nu
 v = particles$v
 #
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
# Rabsz = matrix(rep(abs(particles$z), p), ncol=p)
# sqrt.v = sqrt(as.numeric(t(particles$v)))
# Rsqrtv = matrix(rep(sqrt.v, p), ncol=p)
# Rsqrtv = matrix(rep(sqrt(particles$v), p), ncol=p)
 #
 e = Ry - Rxi
# e = Ry - Rxi - Rabsz / Rsqrtv * Rpsi
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
 log.prior = do.call(logPriorFunc, list(N=N, particles=particles, priorList=priorList))
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
