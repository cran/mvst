sampleBT = function(y, X, N, particles, priorList){
# Given the arguments, this function returns a population of MC draws for the values of the variable B, in the p-variate t regression model.
 n = nrow(y)
 p = ncol(y)
 k = ncol(X)
 B = array(0, c(k, p, N))
 log.dB = numeric(N)
 #
 G = particles$G
 # psi = particles$psi
 v = particles$v
 # z = particles$z
 for(iN in 1:N){
  viN = v[iN,]
  # ziN = z[iN,]
  GiN = matrix(G[iN,,drop=F], ncol=p)
  ViN = diag(viN)
  # absz.sqrtv = abs(ziN) / sqrt(viN)
  # HiN = y - matrix(rep(psi[iN,],n),ncol=p,byrow=T) / matrix(rep(sqrt(viN),p),ncol=p)
  S.iN = t(X) %*% ViN %*% X
  invS.iN = solve(S.iN)
  Cpsi.iN = t(X) %*% ViN %*% y
  M.iN = as.numeric(invS.iN %*% Cpsi.iN)
  # Sigma.temp = kronecker(invS.iN, GiN)
  Sigma.temp = kronecker(GiN, invS.iN)
  Sigma.iN = (Sigma.temp + t(Sigma.temp)) / 2 # force symmetry
  Bvec = as.numeric(rmnorm(1, M.iN, Sigma.iN))
#  Bvec = as.numeric(rmvnorm(1, M.iN, Sigma.iN))
  B[,,iN] = matrix(Bvec, ncol=p)
  log.dB[iN] = dmnorm(Bvec, M.iN, Sigma.iN, log=T)
#  log.dB[iN] = dmvnorm(Bvec, M.iN, Sigma.iN, log=T)
 }

 return(list(values=B, log.dq=log.dB))
}
