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
if(!require(gcdnet)) install.packages("gcdnet")


require(gcdnet)
library(pipeliner)
library(modelr)
library(tidyverse)
library(graphite)

if(!require(hgu133a.db)) BiocManager::install(hgu133a.db)
library(hgu133a.db)


install.packages('ggplot2')
require(ggplot2)  

FILE = paste0(getwd(),'/data/GSE11121_2.csv')

###Source helper files
source(paste0(getwd(),'/source/load_adjMat.R'))



###Source Modles
modelDir = paste0(getwd(),'/source/Models/')
models = list.files(modelDir)
for(f in models)
{
  source(paste0(modelDir,f))
}




####Set standard deviation cutoff
sdCutoff = 0.25

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

###Check for low standard deivation
sds = apply(data,2,sd)
data = data[,sds>sdCutoff | col=="y"]
col = col[sds>sdCutoff | col=="y"]


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
  
  ###TODO only use matched$x if need be, check based upon the model
  if(class(mdl)=="cv.gcdnet")
  {
    #browser()
    pred = predict(mdl,as.matrix(x_te),type="link",s="lambda.min")
  }else if(class(mdl)!="nbSVM")
  {
    pred = predict(mdl,x_te)
  }else
  {
    pred = predict(mdl,matched$x)
  }
  return (pred)
}

###Linear model doesn't work if  | features | > | samples | so take a small chunk of features
temp = data[,-ncol(data)]
target = colnames(as.data.frame(data[,ncol(data)]))[1]

#change to precision recall
#target variable sjo

# 5-fold cross-validation using machine learning pipelines



data[,ncol(data)] = factor(data[,ncol(data)])


#.x = as.data.frame(data)
###Again here the target variable name should not be hardcoded

accuracy_wrapper <- function(x,pred){
  pred = pred$"1"
  pred = exp(pred)/(1+exp(pred))
  #browser()
  return(auc(roc(x[,ncol(x)],pred)))
}

cv_rmse <- crossv_kfold(as.data.frame(data[,c(1:3000,ncol(data))]), 5) %>% 
  mutate(algorithm_name = 'nbsvm',
          select_features = map(train, ~ pipeline_func(as.data.frame(.x),hhsvm_wrapper)),
         predictions = map2(select_features, test, ~ predict_use_model(.x,as.data.frame(.y))),
         #residuals = map2(predictions, test, ~ .x - as.data.frame(.y)[,target]),
#[,c(1:100,ncol(as.data.frame(.x)))]
         #rmse = map_dbl(residuals, ~ sqrt(mean(.x ^ 2)))) %>% summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse))
         #pred = predict(nbSVM, mdev$x),
        #print(pred)	
        accuracy = map2(test, predictions , ~ accuracy_wrapper(as.data.frame(.x),as.data.frame(.y))))
cv_rmse


g = ggplot(cv_rmse, aes(x=algorithm_name, y=accuracy, fill=grouping_variable)) + 
  geom_boxplot() + theme_bw() + ylim(0,1)
