modelParInfo = function(y, X=NULL, modelType){
# This function returns a data.frame containing, the type ('u': scalar, 'm': multivariate, 'M': matrix, 'SM': symmetric matrix) and the dimension of each parameter of modelType.
 n = nrow(y)
 p = ncol(y)
 parInfo = switch(modelType,
  N = data.frame(names=c('xi', 'G'), type=c('m', 'SM'), nCols=c(p, p), nRows=c(1, p), stringsAsFactors=F),
  T = data.frame(names=c('xi', 'G', 'nu', 'v'), type=c('m', 'SM', 'u', 'm'), nCols=c(p, p, 1, n), nRows=c(1, p, 1, 1), stringsAsFactors=F),
  SN = data.frame(names=c('xi', 'G', 'psi', 'z'), type=c('m', 'SM', 'm', 'm'), nCols=c(p, p, p, n), nRows=c(1, p, 1, 1), stringsAsFactors=F),
  ST = data.frame(names=c('xi', 'G', 'psi', 'nu', 'z', 'v'), type=c('m', 'SM', 'm', 'u', 'm', 'm'), nCols=c(p, p, p, 1, n, n), nRows=c(1, p, 1, 1, 1, 1), stringsAsFactors=F)
 )

 if(!is.null(X)){
  k = ncol(X)
#  parInfo[1] = 'B'
  Bindex = which(parInfo['names'] == 'xi') 
  parInfo[Bindex,1:2] = c('B', 'M')
  parInfo[Bindex,3:4] = c(k, p)
 }

 return(parInfo=parInfo)
}