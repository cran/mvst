nSphereVolume = function(n, r=1, LOG=TRUE){
 # Function that returns the volume of a n-sphere with radius r
 logV = (n/2) * log(pi) + n * log(r) - lgamma(n/2+1)
 if(LOG == TRUE) f = logV else f = exp(logV)
 return(f)
}
