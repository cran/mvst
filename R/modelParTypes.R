modelParTypes = function(modelType){
 parTypes = switch(modelType,
  N = c('xi', 'G'),
  T = c('xi', 'G', 'nu', 'v'),
  SN = c('xi', 'G', 'psi', 'z'),
  ST = c('xi', 'G', 'psi', 'nu', 'z', 'v'))
 return(parTypes=parTypes)
}