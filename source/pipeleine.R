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

pipeline_func <- function(df,FUN = lm,target) {
  
  ###This formula assumes that the first variable is the Target variable of interest 
  ### (you probably need to include target as a param)
  form = as.formula(paste(target,"~ ."))
  
  ###Correctly runs the linear model (or any other model)
  
  mdl = FUN(form,df)
  
  ##Return the model
  return(mdl)
}

###Linear model doesn't work if  | features | > | samples | so take a small chunk of features
temp = data[,1:100]
target = colnames(as.data.frame(data))[1]

#change to precision recall
#target variable sjo

# 5-fold cross-validation using machine learning pipelines

custom_fun <- function(x_tr){
  
  
  
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
return(matched)

}

custom_nsvm <- function(y_tr,x_tr){
  m_tr <- custom_func(x_tr)
  model <- fit.networkBasedSVM(m_tr$x, y_tr, DEBUG=FALSE, adjacencyList = ad.list, lambdas = 10^(-1:2),sd.cutoff=1)
  
  return(model)
}

###Again here the target variable name should not be hardcoded
cv_rmse <- crossv_kfold(temp, 5) %>% 
  mutate(model = map(train, ~ pipeline_func(as.data.frame(.x),custom_nsvm,colnames(as.data.frame(.x))[1])),
         predictions = map2(model, test, ~ predict(.x, as.data.frame(.y))),
         #residuals = map2(predictions, test, ~ .x - as.data.frame(.y)[,target]),
         
         rmse = map_dbl(residuals, ~ sqrt(mean(.x ^ 2)))) %>% 
  summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

cv_rmse