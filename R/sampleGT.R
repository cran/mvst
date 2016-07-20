sampleGT = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable G, in the p-variate t model.
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
 v = particles$v
# sqrt.v = as.numeric(t(sqrt(v)))
# z = particles$z
# abs.z = as.numeric(t(abs(z)))

# e = Ry - Rxi - matrix(rep(abs.z, p), ncol=p) / matrix(rep(sqrt.v, p), ncol=p) * Rpsi
 e = Ry - Rxi
 eList = vector('list', p)
 for(ilist in 1:p){
  eList[[ilist]] = matrix(e[,ilist], N, n, byrow=T)
 }
 LambdaMAT = matrix(0, N, n.pmat.indices.w)
 for(icol in 1:n.pmat.indices.w){
  LambdaMAT[,icol] = apply(v * eList[[pmat.indices.w[icol,1]]] * eList[[pmat.indices.w[icol,2]]], 1, sum)
 }
 W.starMAT = matrix(rep(priorList$W[pmat.indices.w],N), N, n.pmat.indices.w, byrow=T) + LambdaMAT
 G = matrix(0, N, p^2)
 log.dG = numeric(N)
 for(iN in 1:N){
  W.star.iN = matrix(0, p, p)
  W.star.iN[pmat.indices.w] = W.starMAT[iN,]
  W.star.iN[pmat.indices.w[,2:1]] = W.starMAT[iN,]
  G.iN = riwish(v=(n+priorList$m), S=W.star.iN)
  G[iN,] = as.numeric(G.iN)
  log.dG[iN] = diwish.mia3(G.iN, n+priorList$m, W.star.iN, LOG=TRUE)
 }
 return(list(values=G, log.dq=log.dG))
}
