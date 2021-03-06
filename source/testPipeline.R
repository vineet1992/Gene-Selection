lib = "/net/maccutcheon/home/benos/vkr8/Rlib"

require(ggplot2,lib=lib)
library(tibble,lib=lib)
library(tidyr,lib=lib)
library(readr,lib=lib)
library(purrr,lib=lib)
library(dplyr,lib=lib)
library(tidyverse,lib=lib)


library(svmpath,lib=lib)
library(kernlab,lib=lib)
library(Biobase,lib=lib)
library(affy,lib=lib)
library(gplots,lib=lib)
library(igraph,lib=lib)
library(lpSolve,lib=lib)

if (!requireNamespace("clue", quietly = TRUE))
  install.packages(clue)

require(clue)

library(pathClass,lib=lib)
library(pROC,lib=lib)
require(gcdnet,lib=lib)
library(pipeliner,lib=lib)
library(modelr,lib=lib)
library(graphite,lib=lib)

library(org.Hs.eg.db,lib=lib)
library(hgu133a.db,lib=lib)
###Source helper files
source(paste0(getwd(),'/source/load_adjMat.R'))

###PARAMETERS
nFolds = 5
####Set standard deviation cutoff
sdCutoff = 0.5

graphSelect = 100

###NBSVM won't find shiz
#names = c("pd","pdNPPCA","pdNPPCAGraph","pdPCA","pdGraph","pdPCAGraph","rrfe","hhsvm","nbsvm")
#names = c("pd","pdPCA","pdGraph","pdNPPCAGraph","pdNPPCA","rrfe","hhsvm")
names = c("pdPCA")#,"pdNPPCAGraph","pdNPPCA","hhsvm","rrfe","pd","pdGraph","pdPCAGraph")
#names = c("pdNPGraph")


###Source Modles
modelDir = paste0(getwd(),'/source/Models/')
models = list.files(modelDir)
for(f in models)
{
  source(paste0(modelDir,f))
}


###Prepare dataset loading
dataDir = paste0(getwd(),'/data/')
datasets = list.files(dataDir)

datasets = datasets[startsWith(datasets,"GSE")]

#datasets = datasets[1:5]
num_rows <- (nfolds*(nfolds-1)/2)*length(datasets)
stability_scores <- matrix(ncol=2,nrow=num_rows)



####Load adjacency matrix
adj = load_adjMat()
adjMatrix <- adj[[1]]
ad.list <- adj[[2]]
mapping <- adj[[3]]


###Helper functions for the pipeline
pipeline_func <- function(df,FUN=lm) {
  
  ###This formula assumes that the first variable is the Target variable of interest 
  ### (you probably need to include target as a param)
  #form = as.formula(paste(target,"~ ."))
  
  ###Correctly runs the linear model (or any other model)
  x_tr = df[,-ncol(df)]
  y_tr = factor(df[,ncol(df)])
  
  mdl = FUN(x_tr,y_tr)
  
  ##Return the model
  return(mdl)
}


###Ensure all prediction types are the same
predict_use_model <- function(mdl,df){
  x_te <- df[,-ncol(df)]
  var <- constructAdjMat(x_te)
  matched <- var[[1]]
  ad.list <- var[[2]]
  mapping <- var[[3]]
  y_te <- factor(df[,ncol(df)])
  
  ###TODO only use matched$x if need be, check based upon the model
  
  
  if(class(mdl)=="cv.gcdnet")
  {
    x_te = x_te[,rownames(mdl$gcdnet.fit$beta)] 
    pred = predict(mdl,as.matrix(x_te),type="link",s="lambda.min")[,1]
  }else if(length(class(mdl))>1)
  {
    ###Check if this is PCA aggregated features
    if(sum(grepl(".",colnames(mdl$data),fixed=T)) > 0)
    {
      ###Extract the columns involved in each column of the data, apply PCA to them and include them in the final testing set
      testData = matrix(nrow=nrow(x_te),ncol=ncol(mdl$data))
      
    
      for(i in 1:ncol(mdl$data))
      {
        
        name = colnames(mdl$data)[i]
        
        if(name=="y")
        {
          testData[,i] = y_te
        }else
        {
          cols = unlist(strsplit(name,".",fixed=T))
          tempDat = x_te[,colnames(x_te)%in%cols]
          p = prcomp(tempDat,scale=T,center=T)
          testData[,i] = p$x[,1]
        }
      }
      colnames(testData) = colnames(mdl$data)
      pred = predict(mdl,as.data.frame(testData))
      pred = 1/(1+exp(-1*pred))

    }else
    {
      
      pred = predict(mdl,x_te)
      pred = 1/(1+exp(-1*pred))
    }
    
    
 
  }else if(class(mdl)!="networkBasedSVM")
  {
    pred = predict(mdl,x_te)[,1]
  }
  else
  {
    pred = predict(mdl,matched$x)[,1]
  }
  return (pred)
}


accuracy_wrapper <- function(x,pred){
  pred = pred$.y
  pred = exp(pred)/(1+exp(pred))
  #browser()
  return(as.numeric(auc(roc(x[,ncol(x)],pred))))
}

stability_wrapper = function(list1,list2){
  final_scores <- matrix(nrow=length(list1),ncol=length(list2))
  #print(ncol(final_scores))
  #print(nrow(final_scores))
  
  i<-1
  for(itemlist1 in list1){
    j<-1
    for(itemlist2 in list2){
      
      il1 = unlist(itemlist1)
      #print(il1[0:10])
      il2 = unlist(itemlist2)
      #print(il2[0:10])
      #print("intersect:")
      #print(length(intersect(il1,il2)))
      #print("union:")
      #print(length(union(il1,il2)))
      #print("score:")
      #print(length(intersect(il1,il2))/length(union(il1,il2)))
      #print("index:i")
      #print(i)
      #print("index:j")
      #print(j)
      final_scores[i,j] <- length(intersect(il1,il2))/length(union(il1,il2))
      #print(final_scores[i,j])
      #final_scores[j,i] = final_scores[i,j]
      j<-j+1
    }
    i<- i+1
  }
  if(nrow(final_scores) > ncol(final_scores))
    final_scores <- t(final_scores)
  #print(final_scores)
  rowscols <- solve_LSAP(final_scores,maximum=TRUE) 
  row<-1
  #print(c(rowscols))
  final_score<-0
  for(item in rowscols){
    #print(row)
    #print(as.numeric(item))
    final_score <- final_score + final_scores[row,item]
    row = row+1
  }
  final_score <- final_score/(row-1)
  
  return(final_score)
}


first = TRUE

for(d in datasets)
{
  
  FILE = paste0(getwd(),'/data/',d)
  
  print(paste0("Running dataset: ", d))
  
  ##Quickly determine the number of columns
  x <- max(count.fields(FILE, ","))
  
  ###Scan the file
  data <- scan(FILE,sep=",")
  
  ###Extract number of rows and columns
  nr = length(data)/x
  nc = x
  data = matrix(data,nrow=nr,ncol=nc,byrow=T)
  
  #data = data[,c(1:500,ncol(data))]
  
  col = scan(paste0(getwd(),"/data/headers.csv"),sep=",",what=character())[1:ncol(data)]
  ## remove the for loop ##
  
  
  ####Cannot use this for-loop, it should be abandoned
  ####A legal R variable cannot start with a number, so this causes problems when using a formula object
  for(i in 1:length(col)) ## removes any leading X in the header names
  {
    if(substring(col[i], 1, 1) == 'X')
      col[i] <- substr(col[i], 2, nchar(col[i]))
  }
  colnames(data) <- col
  
  
  
  ###Linear model doesn't work if  | features | > | samples | so take a small chunk of features
  target = colnames(as.data.frame(data[,ncol(data)]))[1]
  
  #change to precision recall
  #target variable sjo
  
  # 5-fold cross-validation using machine learning pipelines
  
  
  
  data[,ncol(data)] = factor(data[,ncol(data)])
  
  
  foldFile = paste0("data/Folds_",d,".Rdata")
  #.x = as.data.frame(data)
  ###Again here the target variable name should not be hardcoded
  folds = crossv_kfold(as.data.frame(data),nFolds)
  rm(data)
  
  
  if(file.exists(foldFile))
  {
    print("Loading folds...")
    load(foldFile)
  }else
  {
    print("Saving folds...")
    save(folds,file=foldFile)
    
  }
  
  
  ##x_tr = as.data.frame(folds$train[1])
  #y_tr= factor(x_tr[,ncol(x_tr)])
  #x_tr = x_tr[,-ncol(x_tr)]
  #x_te = as.data.frame(folds$test[1])
  #y_te = factor(x_te[,ncol(x_te)])
  #x_te = x_te[,-ncol(x_te)]
  
  for(i in 1:length(names))
  {
	fileName = paste0("Results/",d,"_",names[i],".Rdata")
    
    
    print(paste0("On Algorithm: ",names[i]))
	
    if(file.exists(fileName))
    {
      load(fileName)
      #cv_rmse$select_features=NULL
      cv_rmse$algorithm_name = rep(names[i],nFolds)
#      sf1 = strsplit(colnames(cv_rmse$select_features$`1`$data),".",fixed=TRUE)
#      sf2 = strsplit(colnames(cv_rmse$select_features$`2`$data),".",fixed=TRUE)
      #final_stability_scores <- stability_wrapper(sf1,sf2)
      
      
      for(l1 in (1:(nFolds-1))){
        #print(i)
        vari <- cv_rmse$select_features[[l1]]$data
        vari <- colnames(vari)
        vari <- vari[-length(vari)]
        sf1 = strsplit(vari,".",fixed=TRUE)
        for(l2 in ((l1+1):nFolds)){
          #print(j)
          varj <- cv_rmse$select_features[[l2]]$data
          varj <- colnames(varj)
          varj <- varj[-length(varj)]
          sf2 = strsplit(varj,".",fixed=TRUE)
          final_stability_scores <-  stability_wrapper(sf1,sf2)
          stability_scores[(l2+nFolds*(l1-1))-1,1] = rep(d)
          stability_scores[(l2+nFolds*(l1-1))-1,2] = final_stability_scores
          #print(final_stability_scores)
        }
      }
      save(stability_scores,file=paste0(names[i],".RData"))
      
    }else
    {
      cv_rmse <- folds %>% 
        mutate(select_features = map(train, ~ pipeline_func(as.data.frame(.x),eval(parse(text=paste0(names[i],"_wrapper"))))),
               predictions = map2(select_features, test, ~ predict_use_model(.x,as.data.frame(.y))),
               algorithm_name = rep(names[i],nFolds),
               accuracy = map2(test, predictions , ~ accuracy_wrapper(as.data.frame(.x),as.data.frame(.y))),
               dataset = rep(d,nFolds))
      
      cv_rmse$train=NULL
      cv_rmse$test=NULL
      
     
      
      
      save(cv_rmse,file=fileName)
      
    }
  
    if(first)
    {
      
      cv_rmse$train=NULL
      cv_rmse$test=NULL
      results = cv_rmse
      first = FALSE
    }else
    {
      
      cv_rmse$train=NULL
      cv_rmse$test=NULL
      results = rbind(results,cv_rmse)
    }
gc()    
  }
}

#save(results,file="Results/Full_Results.Rdata")

df = as.data.frame(cbind(results$dataset,results$algorithm_name,results$accuracy))
rownames(df) = seq(1,nrow(df))
colnames(df) <- c('dataset','algorithm_name','accuracy')

df$algorithm_name = as.factor(as.character(df$algorithm_name))
df$accuracy = as.numeric(df$accuracy)
df$dataset = as.factor(as.character(df$dataset))

row = df$dataset%in%datasets[1:3]
row = ifelse(row,1,2)
col = ifelse(df$dataset%in%datasets[c(1,4)],1,ifelse(df$dataset%in%datasets[c(2,5)],2,3))

df$row = factor(row)
df$col = factor(col)

g = ggplot(df, aes(x=algorithm_name, y=accuracy,fill=algorithm_name)) + 
  geom_boxplot() + theme_bw()+ ylim(0.4,1) + facet_grid(row ~ col)


