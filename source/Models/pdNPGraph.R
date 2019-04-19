
###PrefDiv -> StEPS
pdNPGraph_wrapper = function(x_tr,y_tr)
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
  cmd = paste("java -jar -Xmx4g PrefDiv.jar -t y -data temp.txt -outData data_summary.txt -outCluster clusters.txt -outGraph graph.txt -useCausalGraph -numSelect ",graphSelect," -disc -name ",runName,sep="")
  system(cmd)
  
  ###Read in selected genes and train linear model (allow for summarization as well)
  newData = read.table(paste(runName,"graph.txt",sep='/'),sep=' ',header=F,skip=4)
  
  ###No variables were chosen to be connected to the target
  if(sum(newData$V2=="y"|newData$V4=="y")==0)
  {
    mdl = glm(y~1,data=trainData,family=binomial(link="logit"))
    return(mdl)
  }else
  {
    newData = newData[newData$V2=="y" | newData$V4 == "y",]
    predictors = unique(c(as.character(newData$V2),as.character(newData$V4)))
    
    newData = trainData[,colnames(trainData)%in%predictors]
    
    
    mdl = glm(y~.,data = newData,family=binomial(link="logit"))
    
    ###return the linear model
    return(mdl)
  }
  
  
}