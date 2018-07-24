rmvSE = function(n, p, X=NULL, modelType, theta){

 # Import parameter values & checks
 XFlag = !is.null(X)
 if(XFlag){
  if(!('B' %in% names(theta))) stop('(rmvSE.R) cannot find the value of B in theta.\n')
  B = theta$B
  if(!is.matrix(B)){
   if(is.numeric(B) & length(B) == 1){
    B = matrix(B, 1, 1)
   } else {
    stop('(rmvSE.R) the object B should be a matrix.\n')
   }
  }
  k = nrow(B)
  if(ncol(B) != p) stop('(rmvSE.R) non-conformable dimensions for B.\n')
  if(!is.matrix(X)) stop('(rmvSE.R) the object X is not a matrix.\n')
  if(any(dim(X) != c(n, k))) stop('(rmvSE.R) non-conformable dimensions for X.\n')
 } else {
  if(!('xi' %in% names(theta))) stop('(rmvSE.R) cannot find the value of xi in theta.\n')
  xi = theta$xi
  if(length(xi) != p) stop('(rmvSE.R) non-conformable dimension for xi.\n')
 }
 if(!('G' %in% names(theta))) stop('(rmvSE.R) cannot find the value of G in theta.\n')
 G = theta$G
 if(any(dim(G) != c(p, p))) stop('(rmvSE.R) non-conformable dimensions for G.\n')
 if((modelType == 'SN') | (modelType == 'ST')){
  if(!('psi' %in% names(theta))) stop('(rmvSE.R) cannot find the value of psi in theta.\n')
  psi = theta$psi
  if(length(theta$psi)!=p) stop('(rmvSE.R) non-conformable dimension for psi.\n')
 } else psi = rep(0, p)
 if((modelType == 'T') | (modelType == 'ST')){
  if(!('nu' %in% names(theta))) stop('(rmvSE.R) cannot find the value of nu in theta.\n')
  nu = theta$nu
  if(length(theta$nu)!=1) stop('(rmvSE.R) non-conformable dimension for nu.\n')
 }

 # Upper (with and without diagonal) indices of a pxp matrix
 pmat.indices.w = triangleIndices(p, side='u', dgn=T, dataframe=T)
 pmat.indices.wo = triangleIndices(p, side='u', dgn=F, dataframe=T)
# n.pmat.indices.w = p * (p+1) / 2
# n.pmat.indices.wo = p * (p-1) / 2

 # Parameters in all the parameterizations
 Sigma = G + psi %*% t(psi)
 if(p > 1){
  omega = diag(sqrt(diag(Sigma)))
 } else {
  omega = matrix(sqrt(Sigma),1,1)
 }
 delta = as.numeric(solve(omega) %*% psi)
 Omega = cov2cor(Sigma)
 rho.ij = Omega[pmat.indices.wo]
 varcov = rbind(c(1, delta), cbind(delta, Omega))
 colnames(varcov) = NULL
 ZX = rmvnorm(n, rep(0, p+1), varcov)
 z = ZX[,1]
 # u = ZX[,-1] * sign(z)
 u = matrix(ZX[,-1], ncol=p) * sign(z)
 # if(substr(modelType, nchar(modelType), nchar(modelType)) == 'T') v = rgamma(n, nu/2, nu/2)
 if((modelType=='T') | (modelType=='ST')){
  v = rgamma(n, nu/2, nu/2)
 } else {
  v = rep(1, n)
 }
 # alpha = as.numeric((1 - t(delta) %*% solve(Omega) %*% delta)^(-0.5)) * as.numeric(solve(Omega) %*% delta)
 # omega = diag(sqrt(diag(Sigma)))
 # psi = as.numeric(omega %*% delta)
 # Sigma = omega %*% Omega %*% omega
 # G = Sigma - psi %*% t(psi)

 if(XFlag){
  M = X %*% B
 } else {
  M = matrix(rep(xi, each=n), n, p)
 }
 y = matrix(0, n, p)
 for(i in 1:n){
  y[i,] = as.numeric(omega %*% u[i,]) * v[i]^(-0.5) + M[i,]
 }
 if(p == 1) y = as.numeric(y)

 completeData = list(y=y, z=NULL, v=NULL, X=NULL)
 if((modelType=='SN') | (modelType=='ST')) completeData$z = z
 if((modelType=='T') | (modelType=='ST')) completeData$v = v
 if(XFlag) completeData$X = X

 return(completeData)
}
