sampleGST = function(y, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable G, in the p-variate skew-t model.
 n = nrow(y)
 p = ncol(y)
#
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
#	pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
 n.pmat.indices.w = p * (p+1) / 2
#	n.pmat.indices.wo = p * (p-1) / 2
 Ry = matrix(0, n*N, p)
 Rxi = matrix(0, n*N, p)
 Rpsi = matrix(0, n*N, p)
 for(icol in 1:p){
  Ry[,icol] = rep(y[,icol], N)
  Rxi[,icol] = rep(particles$xi[,icol], each=n)
  Rpsi[,icol] = rep(particles$psi[,icol], each=n)
 }
 v = particles$v
 sqrt.v = as.numeric(t(sqrt(v)))
 z = particles$z
 abs.z = as.numeric(t(abs(z)))

 e = Ry - Rxi - matrix(rep(abs.z, p), ncol=p) / matrix(rep(sqrt.v, p), ncol=p) * Rpsi
 eList = vector('list', p)
 for(ilist in 1:p){
  eList[[ilist]] = matrix(e[,ilist], N, n, byrow=T)
 }
# e1 = matrix(e[,1], N, n, byrow=T)
# e2 = matrix(e[,2], N, n, byrow=T)
 LambdaMAT = matrix(0, N, n.pmat.indices.w)
 for(icol in 1:n.pmat.indices.w){
  LambdaMAT[,icol] = apply(v * eList[[pmat.indices.w[icol,1]]] * eList[[pmat.indices.w[icol,2]]], 1, sum)
 }
# LambdaMAT[,1] = apply(v*e1^2, 1, sum)
# LambdaMAT[,2] = apply(v*e1*e2, 1, sum)
# LambdaMAT[,3] = apply(v*e2^2, 1, sum)
 W.starMAT = matrix(rep(priorList$W[pmat.indices.w],N), N, n.pmat.indices.w, byrow=T) + LambdaMAT
# W.starMAT = matrix(rep(c(W[1,1],W[2,1],W[2,2]),N), N, n.pmat.indices.w, byrow=T) + LambdaMAT # W.star = W + Lambda
 G = matrix(0, N, p^2)
# detG = numeric(N)
 log.dG = numeric(N)
# PsiVEC = numeric(N)
# invG = matrix(0, N, p^2)
 for(iN in 1:N){
  W.star.iN = matrix(0, p, p)
  W.star.iN[pmat.indices.w] = W.starMAT[iN,]
  W.star.iN[pmat.indices.w[,2:1]] = W.starMAT[iN,]
  G.iN = riwish(v=(n+priorList$m), S=W.star.iN)	# ??? O LambdaMAT ???
  G[iN,] = as.numeric(G.iN)
#  detG[iN] = det(G.iN) # G[,1] * G[,4] - G[,2]^2
  log.dG[iN] = diwish.mia3(G.iN, n+priorList$m, W.star.iN, LOG=TRUE)
#  invG.iN = solve(G.iN)
#  invG[iN,] = as.numeric(invG.iN)
#  Lambda.iN = matrix(0, p, p)
#  Lambda.iN[pmat.indices.w] = LambdaMAT[iN,]
#  Lambda.iN[pmat.indices.w[,2:1]] = LambdaMAT[iN,]
#  PsiVEC[iN] = sum(t(invG.iN) * Lambda.iN)
 }
#  G.wo = G[,triangleIndices(p, side='u', dgn=T, dataframe=F)] # G without redundant columns
#  RG = matrix(0, n*N, n.pmat.indices.w)
#  for(icol in 1:n.pmat.indices.w){
#   RG[,icol] = rep(G.wo[,icol], each=n)
#  }
# invdetG = detG^(-1) # |G|^(-1)
# RinvdetG = rep(invdetG, each=n)
 return(list(values=G, log.dq=log.dG)) #, determinant=detG, invdet=invdetG, repinvdet=RinvdetG, inverse=invG, repG=RG, psi=PsiVEC)) #, Lambda=LambdaMAT))
}
