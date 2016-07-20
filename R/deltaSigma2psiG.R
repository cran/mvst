deltaSigma2psiG = function(delta, Sigma){
 # Given a vector delta and a matrix Sigma, the function returns the corresponding values of psi and G.
 # Multiple input values can be entered as a pair of lists, in which a generic element of each list represents the value of the respective parameter. The output will be provided in the same form.
 dS2pG = function(d, S){
  # Reparameterization, when p > 1
  omega = diag(sqrt(diag(S)))
  psi = as.numeric(omega %*% d)
  G = S - psi %*% t(psi)
  return(list(psi, G))
 }

 ds2pg = function(d, s){
  # Reparameterization, when p = 1
  omega = sqrt(s)
  psi = omega * d
  G = s - psi^2
  return(list(psi, G))
 }

 if(is.list(delta) != is.list(Sigma)) stop('non-conformable arguments (errType=1).\n')
 if(!is.list(delta)){
  nValues = 1
  p = length(delta)
  if(any(dim(as.matrix(Sigma)) != c(p, p))) stop('non-conformable arguments (errType=2).\n')
  delta = list(delta=delta)
  Sigma = list(Sigma=Sigma)
 } else {
  nValues = length(delta)
  p = unlist(lapply(delta, length)) # Dimension of each element of the vector delta
  for(j in 1:nValues) if(any(dim(as.matrix(Sigma[[j]])) != c(p[j], p[j]))) stop('non-conformable arguments (errType=4).\n') # check: each element of delta is conformable to the respective element of Sigma.
 }

 psi = vector('list', nValues)
 G = vector('list', nValues)
 for(j in 1:nValues){
  d = delta[[j]]
  S = Sigma[[j]]
  if(p[j] == 1){
   values = ds2pg(d,S)
  } else {
   values = dS2pG(d,S)
  }
  if(nValues==1){
   psi = values[[1]]
   G = values[[2]]
  } else {
   psi[[j]] = values[[1]]
   G[[j]] = values[[2]]
  }
 }
 return(list(psi=psi, G=G))
}
