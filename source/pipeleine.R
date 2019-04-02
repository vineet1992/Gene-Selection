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

if(!require(hgu133a.db)) BiocManager::install(hgu133a.db)
library(hgu133a.db)
  

data <- read.csv(paste0(getwd(),'/data/GSE11121_2.csv'),header=FALSE)

col = read.csv(paste0(getwd(),"/data/headers.csv"),header=TRUE)[1:ncol(data)]
col <- colnames(col)
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
  y_tr = df[,ncol(df)]
  
  mdl = FUN(x_tr,y_tr)
  
  ##Return the model
  return(mdl)
}

predict <- function(mdl,df){
  x_te <- df[,-ncol(df)]
  y_te <- df[,ncol(df)]
  pred = predict(mdl, x_te)
  return (pred)
}

###Linear model doesn't work if  | features | > | samples | so take a small chunk of features
temp = data[,-ncol(data)]
target = colnames(as.data.frame(data[,ncol(data)]))[1]

#change to precision recall
#target variable sjo

# 5-fold cross-validation using machine learning pipelines

custom_func <- function(x_tr){
  
  
  #browser()
  mapped.probes <- mappedkeys(hgu133aENTREZID)
  refseq <- as.list(hgu133aENTREZID[mapped.probes])
  times <- sapply(refseq, length)
  mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq), row.names=NULL, stringsAsFactors=FALSE)
  mapping <- unique(mapping)
  #put if conditional
  
  paths <- pathways("hsapiens", "kegg")
  alledges <- NULL
  
  ## creates an adjacency matrix
  ## remove for loops
  
  for(i in 1:length(paths)) ## loop through each of the pathways
  {
    curr2 <- paths[[i]]
    edges <-graphite::edges(curr2)
    edges <- edges[,c("src","dest")]
    alledges <- rbind(alledges, edges)
    alledges = unique(alledges) ## Get the unique edges in the pathway just as a list of gene-gene interactions
  }
  
  src <- as.matrix(sort(as.numeric(unique(alledges[,c("src")]))))
  dest <- as.matrix(sort(as.numeric(unique(alledges[,c("dest")]))))
  vals <- unique(rbind(src, dest))
  
  adjMatrix <- matrix(0, length(vals), length(vals))
  rownames(adjMatrix) <- vals
  colnames(adjMatrix) <- vals
  ### remove for loops
  for(i in 1:nrow(alledges))
  {
    col1 <- alledges[i,"src"]
    col2 <- alledges[i,"dest"]
    adjMatrix[toString(col1),toString(col2)] <- 1
    adjMatrix[toString(col2),toString(col1)] <- 1
  }
  
  
  
  ad.list <- as.adjacencyList(adjMatrix)
  
  
matched <- matchMatrices(x = x_tr, adjacency = adjMatrix, mapping = mapping)
return(list(matched,ad.list,mapping))

}

custom_nsvm <- function(x_tr,y_tr){
  #browser()
  var <- custom_func(x_tr)
  matched <- var[[1]]
  ad.list <- var[[2]]
  mapping <- var[[3]]
  x = x_tr
  nbSVM = fit.networkBasedSVM(matched$x, y_tr, DEBUG=FALSE, adjacencyList = ad.list, lambdas = 10^(-1:2),sd.cutoff=0.5)
  mapping[mapping[,2] %in% nbSVM$features,]
  return(nbSVM)
}


#.x = as.data.frame(data)
###Again here the target variable name should not be hardcoded

cv_rmse <- crossv_kfold(as.data.frame(data), 5) %>% 
  mutate(select_features = map(train, ~ pipeline_func(as.data.frame(.x),custom_nsvm)),
         predictions = map2(select_features, test, ~ predict(select_features,as.data.frame(.x))))
         #residuals = map2(predictions, test, ~ .x - as.data.frame(.y)[,target]),
#[,c(1:100,ncol(as.data.frame(.x)))]
         #rmse = map_dbl(residuals, ~ sqrt(mean(.x ^ 2)))) %>% summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse))
         #pred = predict(nbSVM, mdev$x)
         #print(pred)	
         #A = roc(y_te,pred)
         #AUC.nbsvm[i+1] = auc(A)
cv_rmse