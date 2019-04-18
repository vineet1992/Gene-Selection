##Pref-Div with no prior information
pdNP_wrapper = function(x_tr,y_tr)
{
  nc = ncol(x_tr)
  ###Subset by standard deviation
  sds = apply(x_tr,2,sd)
  
  x_tr = x_tr[,sds>sdCutoff]
  
  
  

  
  trainData = cbind(x_tr,y_tr)
  colnames(trainData)[ncol(trainData)] = "y"
  
  write.table(trainData,file="temp.txt",sep='\t',quote=F,row.names=F)
  
  runName = "PD"
  ###Submit jar file command
  cmd = paste("java -jar -Xmx4g PrefDiv.jar -t y -data temp.txt -outData data_summary.txt -outCluster clusters.txt -cv 3 1,3,5,10,15 -disc -name ",runName,sep="")
  system(cmd)
  
  ###Read in selected genes and train linear model (allow for summarization as well)
  newData = read.table(paste(runName,"data_summary.txt",sep='/'),sep='\t',header=T)
  col = colnames(newData)
  for(i in 1:length(col)) ## removes any leading X in the header names
  {
    if(substring(col[i], 1, 1) == 'X')
      col[i] <- substr(col[i], 2, nchar(col[i]))
  }
  colnames(newData) <- col
  
  newData$y = as.factor(newData$y)
  levels(newData$y) = c("1","2")
  mdl = glm(y~.,newData,family=binomial(link="logit"))
  
  ###return the linear model
  return(mdl)
}