triangleIndices = function(p, side='u', dgn = F, dataframe=T){
# Description: vector of indices of the upper triangle of a pxp matrix.
# Arguments:
#  p = dimension of the matrix
#  side = upper ('u') or lower ('l') triangle indices
#  dgn = logical flag indicating whether the indices of the diagonal should be returned
#  dataframe = should indices be returned as a vector (F) or as a data.frame (T)?
# Examples:
# Omega = matrix(9:1, 3) # 3x3 matrix
# p = nrow(Omega)
# # Extract the upper off-diagonal elements, as a data.frame
# indicesDF = triangleIndices(p, 'u', T)
# # or, as a vector
# indicesVec = triangleIndices(p, 'u', F)
# Omega[indicesDF]	 # 6 3 2
# Omega[indicesVec] # 6 3 2
# # To set both the upper AND the lower elements of the matrix, use dataframe=T
# Omega[indicesDF]	= -1
# Omega[indicesDF[,2:1]] = -2
 X = matrix(0, p, p)
 foo = function(x){
  switch(side,
   'u'=upper.tri(x, diag=dgn), 'l'=lower.tri(x, diag=dgn))
 }
 indices = which(foo(X) == T, arr.ind=dataframe)
 return(indices)
}
