hhsvm_wrapper <- function(x_tr,y_tr){
  hhsvm = cv.gcdnet(as.matrix(x_tr),y_tr,method="hhsvm",pred.loss="misclass",nfolds=3)
  #browser()
  
  return(hhsvm)
}