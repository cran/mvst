rmvSE = function(n, p, modelType, theta, seed=NULL){

 # Import parameter values & checks
 xi = theta$xi
 if(length(xi) != p) stop('(rmvse.R) length of xi != p.')
 G = theta$G
 if(any(dim(G) != c(p, p))) stop('(rmvse.R) non-conformable dimensions for G.')
 if((modelType == 'SN') | (modelType == 'ST')){
  psi = theta$psi
  if(length(theta$psi)!=p) stop('(rmvse.R) non-conformable dimension for psi.')
 } else psi = rep(0, p)
 if((modelType == 'T') | (modelType == 'ST')){
  nu = theta$nu
  if(length(theta$nu)!=1) stop('(rmvse.R) non-conformable dimension for nu.')
 }

 # Upper (with and without diagonal) indices of a pxp matrix
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
# n.pmat.indices.w = p * (p+1) / 2
# n.pmat.indices.wo = p * (p-1) / 2

 # Parameters in all the parameterizations
 Sigma = G + psi %*% t(psi)
 omega = diag(sqrt(diag(Sigma)))
 delta = as.numeric(solve(omega) %*% psi)
 Omega = cov2cor(Sigma)
 rho.ij = Omega[pmat.indices.wo]
 varcov = rbind(c(1, delta), cbind(delta, Omega))
 colnames(varcov) = NULL
 if(!is.null(seed)) set.seed(seed)
 ZX = rmvnorm(n, rep(0, p+1), varcov)
 z = ZX[,1]
 u = ZX[,-1] * sign(z)
# if(substr(modelType, nchar(modelType), nchar(modelType)) == 'T') v = rgamma(n, nu/2, nu/2)
 if((modelType=='T') | (modelType=='ST')) v = rgamma(n, nu/2, nu/2)
# alpha = as.numeric((1 - t(delta) %*% solve(Omega) %*% delta)^(-0.5)) * as.numeric(solve(Omega) %*% delta)
# omega = diag(sqrt(diag(Sigma)))
# psi = as.numeric(omega %*% delta)
# Sigma = omega %*% Omega %*% omega
# G = Sigma - psi %*% t(psi)

 #
# allNames = modelNames(p, modelType)
# parNames = allNames$parNames
# otherNames = allNames$otherNames

# theta = vector('list', 0)
# if(substr(modelType, 1, 1) == 'S') theta$z = z
# if(substr(modelType, nchar(modelType), nchar(modelType)) == 'T') theta$v = v
# theta$otherPars = c(delta, alpha, diag(omega), rho.ij) # , Sigma[pmat.indices.w]))
# names(theta[['otherPars']]) = otherNames

 y = matrix(0, n, p)
 if(substr(modelType, nchar(modelType), nchar(modelType)) == 'T'){
  for(i in 1:n) y[i,] = as.numeric(omega %*% u[i,]) * v[i]^(-0.5) + xi
 } else {
  for(i in 1:n) y[i,] = as.numeric(omega %*% u[i,]) + xi
 }

 completeData = list(y=y, z=NULL, v=NULL)
 if((modelType=='SN') | (modelType=='ST')) completeData$z = z
 if((modelType=='T') | (modelType=='ST')) completeData$v = v

 return(completeData)
}
