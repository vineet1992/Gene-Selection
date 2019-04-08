if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(pathClass)) install.packages('../dependencies/pathClass_0.9.4.tar.gz')
library(pathClass)
if(!require(pROC)) install.packages('pROC')
library(pROC)

if(!require(pipeliner)) devtools::install_github('tidyverse/tidyverse')
if(!require(pipeliner)) devtools::install_github('tidyverse/modelr')

library(pipeliner)
library(modelr)
library(tidyverse)
library(graphite)

if(!require(hgu133a.db)) BiocManager::install(hgu133a.db)
library(hgu133a.db)
  
FILE = paste0(getwd(),'/data/GSE11121_2.csv')
source(paste0(getwd(),'/source/load_adjMat.R'))


####Load adjacency matrix
adj = load_adjMat()
adjMatrix <- adj[[1]]
ad.list <- adj[[2]]
mapping <- adj[[3]]

##Quickly determine the number of columns
x <- max(count.fields(FILE, ","))

###Scan the file
data <- scan(FILE,sep=",")

###Extract number of rows and columns
nr = length(data)/x
nc = x
data = matrix(data,nrow=nr,ncol=nc,byrow=T)


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

pipeline_func <- function(df,FUN=lm) {
  
  ###This formula assumes that the first variable is the Target variable of interest 
  ### (you probably need to include target as a param)
  #form = as.formula(paste(target,"~ ."))
  
  ###Correctly runs the linear model (or any other model)
  x_tr = df[,-ncol(df)]
  y_tr = factor(df[,ncol(df)])
  print(table(y_tr))
  
  mdl = FUN(x_tr,y_tr)
  
  ##Return the model
  return(mdl)
}

predict_use_model <- function(mdl,df){
  x_te <- df[,-ncol(df)]
  var <- constructAdjMat(x_te)
  matched <- var[[1]]
  ad.list <- var[[2]]
  mapping <- var[[3]]
  y_te <- factor(df[,ncol(df)])
  pred = predict(mdl, matched$x)
  return (pred)
}

###Linear model doesn't work if  | features | > | samples | so take a small chunk of features
temp = data[,-ncol(data)]
target = colnames(as.data.frame(data[,ncol(data)]))[1]

#change to precision recall
#target variable sjo

# 5-fold cross-validation using machine learning pipelines

constructAdjMat <- function(x_tr){

  
matched <- matchMatrices(x = x_tr, adjacency = adjMatrix, mapping = mapping)
return(list(matched,ad.list,mapping))

}

nsvm_wrapper <- function(x_tr,y_tr){
  #browser()
  var <- constructAdjMat(x_tr)
  matched <- var[[1]]
  ad.list <- var[[2]]
  mapping <- var[[3]]
  nbSVM = fit.networkBasedSVM(matched$x, y_tr, DEBUG=T, adjacencyList = ad.list, lambdas = 10^(-1:2),sd.cutoff=0.5)
  nbSVM
  mapping[mapping[,2] %in% nbSVM$features,]
  return(nbSVM)
}

data[,ncol(data)] = factor(data[,ncol(data)])

#.x = as.data.frame(data)
###Again here the target variable name should not be hardcoded

cv_rmse <- crossv_kfold(as.data.frame(data[,c(1:3000,ncol(data))]), 5) %>% 
  mutate(select_features = map(train, ~ pipeline_func(as.data.frame(.x),nsvm_wrapper)),
         predictions = map2(select_features, test, ~ predict_use_model(.x,as.data.frame(.y))))
         #residuals = map2(predictions, test, ~ .x - as.data.frame(.y)[,target]),
#[,c(1:100,ncol(as.data.frame(.x)))]
         #rmse = map_dbl(residuals, ~ sqrt(mean(.x ^ 2)))) %>% summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse))
         #pred = predict(nbSVM, mdev$x)
         #print(pred)	
         #A = roc(y_te,pred)
         #AUC.nbsvm[i+1] = auc(A)
cv_rmse