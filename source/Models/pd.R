

###PiPrefDiv 
pd_wrapper = function(x_tr,y_tr)
{
  nc = ncol(x_tr)
  ###Subset by standard deviation
  sds = apply(x_tr,2,sd)
  
  x_tr = x_tr[,sds>sdCutoff]
  
  ###Load prior information and create temp directory to store
  dir = "./priors/"
  tempDir = "./tempPriors/"
  
  dir.create(tempDir)
  
  print("Reading priors...")
  for(f in list.files(dir))
  {
    print(f)
    curr = matrix(scan(paste0(dir,f),sep="\t"),ncol=nc+1,nrow=nc+1,byrow=T)
    ###Subset prior information row-wise and column-wise based on sd cutoff
    
    ###Ensure that you keep the last row and column (for the target variable)
    curr = curr[c(sds>sdCutoff,TRUE),c(sds>sdCutoff,TRUE)]
    
    write.table(curr,paste0(tempDir,f),sep='\t',col.names=F,row.names=F,quote=F)
  
  }
  
  trainData = cbind(x_tr,y_tr)
  colnames(trainData)[ncol(trainData)] = "y"
  
  write.table(trainData,file="temp.txt",sep='\t',quote=F,row.names=F)
  
  print("Running pref-div")
  runName = "PD"
  ###Submit jar file command
  cmd = paste("java -jar -Xmx4g PrefDiv.jar -t y -data temp.txt -outData data_summary.txt -outCluster clusters.txt -cv 3 1,3,5,10 -priors tempPriors -disc -name ",runName,sep="")
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
  mdl = glm(y~.,data =newData,family=binomial(link="logit"))
  
  ###return the linear model
  return(mdl)
}