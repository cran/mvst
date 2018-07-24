cmlSE = function(modelType, y, z=NULL, v=NULL, X=NULL){
 # Checks
 if(is.data.frame(y)) y = as.matrix(y)
 if(is.vector(y)) y = matrix(y, ncol=1)
 dimnames(y) = NULL
 #
 latentVars = NULL
 if(!is.null(v)){
  latentVars$v = v
 }
 if(!is.null(z)){
  latentVars$z = z
 }
 CML = do.call(paste('cml', modelType, sep=''), list(y=y, latentVars=latentVars, X=X))
 return(CML)
}
