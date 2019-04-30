constructAdjMat <- function(x_tr){
  
  
  matched <- matchMatrices(x = x_tr, adjacency = adjMatrix, mapping = mapping)
  return(list(matched,ad.list,mapping))
  
}

nsvm_wrapper <- function(x_tr,y_tr){
  var <- constructAdjMat(x_tr)
  matched <- var[[1]]
  ad.list <- var[[2]]
  mapping <- var[[3]]
  nbSVM = fit.networkBasedSVM(matched$x, y_tr, DEBUG=F, adjacencyList = ad.list, lambdas = 10^(-1:2),sd.cutoff=sdCutoff)
  nbSVM
  return(nbSVM)
}