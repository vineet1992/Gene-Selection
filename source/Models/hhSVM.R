hhsvm_wrapper <- function(x_tr,y_tr){
  sds = apply(x_tr,2,sd)
  x_tr = x_tr[,sds>sdCutoff]
  hhsvm = cv.gcdnet(as.matrix(x_tr),y_tr,method="hhsvm",pred.loss="misclass",nfolds=3)
  #browser()
  
  return(hhsvm)
}