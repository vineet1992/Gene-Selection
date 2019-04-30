
###PrefDiv -> StEPS
pdPCAGraph_wrapper = function(x_tr,y_tr)
{
  nc = ncol(x_tr)
  ###Subset by standard deviation
  sds = apply(x_tr,2,sd)
  
  x_tr = x_tr[,sds>sdCutoff]
  ###Load prior information and create temp directory to store
  dir = "./priors/"
  tempDir = "./tempPriors/"
  
  dir.create(tempDir)
  for(f in list.files(dir))
  {
    curr = matrix(scan(paste0(dir,f),sep="\t"),ncol=nc+1,nrow=nc+1,byrow=T)
    ###Subset prior information row-wise and column-wise based on sd cutoff
    
    ###Ensure that you keep the last row and column (for the target variable)
    curr = curr[c(sds>sdCutoff,TRUE),c(sds>sdCutoff,TRUE)]
    
    write.table(curr,paste0(tempDir,f),sep='\t',col.names=F,row.names=F,quote=F)
    
  }
  
  trainData = cbind(x_tr,y_tr)
  colnames(trainData)[ncol(trainData)] = "y"
  
  write.table(trainData,file="temp2.txt",sep='\t',quote=F,row.names=F)
  
  runName = "PD"
  ###Submit jar file command
  cmd = paste("java -jar -Xmx4g PrefDiv.jar -t y -data temp2.txt -ctype pca -priors tempPriors -outData data_summary.txt -outCluster clusters.txt -outGraph graph.txt -useCausalGraph -numSelect ",graphSelect," -disc -name ",runName,sep="")
  system(cmd)
  
  ###Read summarized dataset
  trainData = read.table(paste(runName,"data_summary.txt",sep='/'),sep='\t',header=T)
  
  
  ###Remove X's
  col = colnames(trainData)
  for(i in 1:length(col)) ## removes any leading X in the header names
  {
    if(substring(col[i], 1, 1) == 'X')
      col[i] <- substr(col[i], 2, nchar(col[i]))
  }
  colnames(trainData) <- col
  
  ###Read in selected genes and train linear model (allow for summarization as well)
  newData = read.table(paste(runName,"graph.txt",sep='/'),sep=' ',header=F,skip=4)
  
  ###No variables were chosen to be connected to the target
  if(sum(newData$V2=="y"|newData$V4=="y")==0)
  {
    trainData$y = factor(trainData$y)
    mdl = glm(y~1,data=trainData,family=binomial(link="logit"))
    return(mdl)
  }else
  {
    newData = newData[newData$V2=="y" | newData$V4 == "y",]
    predictors = unique(c(as.character(newData$V2),as.character(newData$V4)))
    
    predictors = gsub("|",".",predictors,fixed=T)
    
    newData = trainData[,colnames(trainData)%in%predictors]
    
    newData$y = factor(newData$y)
    
    
    mdl = glm(y~.,data = newData,family=binomial(link="logit"))
    
    ###return the linear model
    return(mdl)
  }
  
  
}