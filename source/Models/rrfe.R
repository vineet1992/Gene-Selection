rrfe_wrapper <- function(x_tr,y_tr){
  #browser()
  var <- constructAdjMat(x_tr)
  matched <- var[[1]]
  ad.list <- var[[2]]
  mapping <- var[[3]]
  RRFE = fit.rrfe(x_tr, y_tr, DEBUG=T, mapping = mapping,Gsub = adjMatrix)
  RRFE
  return(RRFE)
}