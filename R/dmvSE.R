dmvSE = function(y, X=NULL, modelType, theta, LOG=FALSE){
 # Import data and parameter values
 if(is.vector(y)){
  p = nrow(theta$G)
  if(p == 1){
   y = matrix(y, nrow=length(y))
  } else {
   if(p == length(y)){
    y = matrix(y, ncol=length(y))
   } else stop('(dmvSE.R) non-conformable dimensions for y and theta.\n')
  }
 } else {
  if(!is.matrix(y)) y = as.matrix(y)
 }
 n = nrow(y)
 p = ncol(y)
 XFlag = !is.null(X)
 if(XFlag){
  if(!('B' %in% names(theta))) stop('(dmvSE.R) cannot find the value of B in theta.\n')
  B = theta$B
  k = nrow(B)
  if(ncol(B) != p) stop('(dmvSE.R) non-conformable dimensions for B.\n')
  if(any(dim(X) != c(n, k))) stop('(dmvSE.R) non-conformable dimensions for X.\n')
  M = X %*% B # Location parameter (one row per observation)
 } else {
  if(!('xi' %in% names(theta))) stop('(dmvSE.R) cannot find the value of xi in theta.\n')
  xi = theta$xi
  if(!is.vector(xi)) stop('(dmvSE.R) xi should be a numeric vector.\n')
  if(length(xi) != p) stop('(dmvSE.R) non-conformable dimension for xi.\n')
  M = matrix(rep(xi, each=n), n, p) # Location parameter (one row per observation)
 }

 if(!('G' %in% names(theta))) stop('(dmvSE.R) cannot find the value of G in theta.\n')
 G = theta$G
 if(any(dim(G) != rep(p,2))) stop('(dmvSE.R) non-conformable dimensions for the parameter G.\n')

 # Compute the log-density
 # N
 if(modelType == 'N') logDens = dmvnorm(y, M, G, log=T)

 # T
 if(modelType == 'T') logDens = dmvt(y, M, G, df=nu, log=T)

 # Skewed distributions
 if(modelType %in% c('SN', 'ST')){
  if(!('psi' %in% names(theta))) stop('(dmvSE.R) cannot find the value of psi in theta.\n')
  psi = theta$psi
  if(!is.vector(psi)) stop('(dmvSE.R) psi should be a numeric vector.\n')
  if(length(psi) != p) stop('(dmvSE.R) non-conformable dimensions for the parameter psi.\n')
  logDens = numeric(n)
  Sigma = G + psi %*% t(psi)
  invSigma = solve(Sigma)
  if(p > 1){
   omega = diag(sqrt(diag(Sigma)))
  } else {
   omega = matrix(sqrt(Sigma),1,1)
  }
  invomega = solve(omega)
  Omega = cov2cor(Sigma)
  invOmega = solve(Omega)
  delta = numeric(p)
  for(j in 1:p){
   delta[j] = (G[j,j] + psi[j]^2)^(-0.5) * psi[j]
  }
  alpha = (1-as.numeric(t(delta) %*% invOmega %*% delta))^(-0.5) * (invOmega %*% delta)
 }

 # SN
 if(modelType == 'SN'){
  for(i in 1:n){
   u = as.numeric(t(alpha) %*% invomega %*% (y[i,]-M[i,]))
   logDens[i] = log(2) + dmvnorm(y[i,], M[i,], Sigma, log=T) + pnorm(u, log.p=T)
  }
 }

 # ST
 if(modelType == 'ST'){
  if(!('nu' %in% names(theta))) stop('(dmvSE.R) cannot find the value of nu in theta.\n')
  nu = theta$nu
  if(length(nu) != 1) stop('(dmvSE.R) nu should be a scalar.\n')
  for(i in 1:n){
   Qy = as.numeric(t(y[i,] - M[i,]) %*% invSigma %*% (y[i,] - M[i,]))
   u = as.numeric(t(alpha) %*% invomega %*% (y[i,]-M[i,])) * ((nu+p) / (Qy+nu))^0.5
   logDens[i] = log(2) + dmvt(y[i,], M[i,], Sigma, df=nu, log=T) + pt(u, nu+p, log.p=T)
  }
 }

 # Return the (log-)density
 if(LOG == F) return(exp(logDens)) else return(logDens)
}
