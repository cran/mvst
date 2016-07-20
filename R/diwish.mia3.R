diwish.mia3 = function (W, v, S, LOG=F){
# Modified version of the function diwish() from the MCMCpack package (2016/07/18).
 if (!is.matrix(S)) S <- matrix(S)
 if (nrow(S) != ncol(S)) stop('W not square in diwish().\n')
 if (!is.matrix(W)) S <- matrix(W)
 if (nrow(W) != ncol(W)) stop('W not square in diwish().\n')
 if (nrow(S) != ncol(W)) stop('W and X of different dimensionality in diwish().\n')
 if (v < nrow(S))  stop('v is less than the dimension of S in  diwish().\n')
 k <- nrow(S)
 log.gammapart <- 0
 for (i in 1:k) {
  log.gammapart <- log.gammapart + lgamma((v + 1 - i)/2)
 }
 log.denom <- log.gammapart + (v * k/2) * log(2) + (k * (k - 1)/4) * log(pi)
 detS <- det(S)
 detW <- det(W)
 hold <- S %*% solve(W)
 tracehold <- sum(hold[row(hold) == col(hold)])
 log.num <- (v/2) * log(detS) + (-(v + k + 1)/2) * log(detW) - 0.5 * tracehold
 log.out = log.num - log.denom
 if(LOG==T) return(log.out)
 if(LOG==F) return(exp(log.out))
}
