villaPrior = function(d, nuMax=30){
 # Prior distribution for the number of degrees of freedom for a d-variate t model.
 # See Villa, Rubio (2018) 'Objective priors for the number of degrees of freedom of a multivariate t distribution and the t-copula' arXiv

 # Checks
 if((d %% 1) != 0) stop('d should contain an integer value\n')
 if((nuMax %% 1) != 0) stop('nuMax should contain an integer value\n')

 # require(mvtnorm)
 nus = 1:nuMax
 DKLp = numeric(nuMax)

 ## Functions
 log.Kdnu = function(nu, dd=d){
  value = lgamma(0.5*(nu+dd)) - lgamma(0.5*nu) - 0.5*dd*log(pi*nu)
  return(value)
 }

 # Integrand function
 ft = function(tt, nu0, nuPrime, dd=d){
  value = (1+tt/nu0)^(-0.5*(nu0+dd)) * tt^(0.5*dd-1) * log(1 + tt/nuPrime)
  return(value)
 }

 for(iNu in 1:nuMax){
  nu = nus[iNu]
  # First expectation
  FE = digamma(0.5*(nu+d)) - digamma(0.5*nu)
  # Second expectation
  SE = exp(log.Kdnu(nu, d) + 0.5*d*log(pi) - lgamma(0.5*d)) * integrate(ft, 0, Inf, nu0=nu, nuPrime=nu+1)[[1]]
  # KL divergence among f_{nu} and f_{nu+1}
  DKLp[iNu] = log.Kdnu(nu, d) - log.Kdnu(nu+1, d) - 0.5*(nu+d) * FE + 0.5 * (nu+1+d) * SE
 }

 probsBar = exp(DKLp) - 1
 probs = probsBar / sum(probsBar)
 # return(DKLp)
 return(probs)
}
