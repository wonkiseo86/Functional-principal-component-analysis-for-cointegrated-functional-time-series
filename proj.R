function(basis_temp,f,grid){
  temp = 0 
  for (k in 1:ncol(basis_temp)){
    temp = temp + inner(basis_temp[,k],f,grid)*basis_temp[,k]
  }
  return(temp)
}