rrfe_wrapper <- function(x_tr,y_tr){
  sds = apply(x_tr,2,sd)
  x_tr = x_tr[,sds>sdCutoff]
  var <- constructAdjMat(x_tr)
  matched <- var[[1]]
  ad.list <- var[[2]]
  mapping <- var[[3]]
  RRFE = fit.rrfe(x_tr, y_tr, DEBUG=F, mapping = mapping,Gsub = adjMatrix)
  RRFE
  return(RRFE)
}