cmlSE = function(modelType, y, z=NULL, v=NULL){
 n = nrow(y)
 p = ncol(y)
 #
 latentVars = NULL
 if(!is.null(v)){
#  modelType = 'T'
  latentVars$v = v
# } else {
#  modelType = 'N'
 }
 if(!is.null(z)){
#  modelType = paste('S', modelType, sep='')
  latentVars$z = z
 }
 CML = do.call(paste('cml', modelType, sep=''), list(y=y, latentVars=latentVars))
 return(CML)
}
