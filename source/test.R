#this will load each dataset, and set up the prior knowledge adjacency matrix
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(hgu133a.db))
  BiocManager::install('hgu133a.db')
library(hgu133a.db)
if(!require(graphite))
  BiocManager::install(graphite)
library(graphite)
setwd("./")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(pathClass)) install.packages('../dependencies/pathClass_0.9.4.tar.gz')
library(pathClass)
if(!require(pROC)) install.packages('pROC')
library(pROC)
setwd("./")

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


x_tr <- data[,-ncol(data)]
y_tr <- data[,ncol(data)]



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

matched <- matchMatrices(x = x_tr, adjacency = adjMatrix, mapping = mapping)
#dimad.List <- as.adjacencyList(matched$adjacency)
nbSVM = fit.networkBasedSVM(matched$x, y_tr, DEBUG=FALSE, adjacencyList = ad.list, lambdas = 10^(-1:2),sd.cutoff=0.3)
paste(mapping[mapping[,2] %in% nbSVM$features,])
m <- matchMatrices(x=x_te, adjacency = adjMatrix, mapping=mapping)
pred = predict(nbSVM, mdev$x)
print(pred)	
A = roc(y_te,pred)
AUC.nbsvm[i+1] = auc(A)