## gets .csv file with the header manually added in
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(hgu133a.db)) install.packages('hgu133a.b')
library(hgu133a.db)
if(!require(graphite)) install.packages('graphite')
library(graphite)
if(!require(pathClass)) install.packages('../dependencies/pathClass_0.9.4.tar.gz')
library(pathClass)
setwd("./")

data <- scan("./data/GSE1456_2.csv", sep = ',') ## read in data
col <- colnames(data)
for(i in 1:length(col)) ## removes any leading X in the header names
{
  if(substring(col[i], 1, 1) == 'X')
    col[i] <- substr(col[i], 2, nchar(col[i]))
}
colnames(data) <- col
x <- (data)
y <- factor(data$y)

## creates a mapping using ENTREZID

mapped.probes <- mappedkeys(hgu133aENTREZID)
refseq <- as.list(hgu133aENTREZID[mapped.probes])
times <- sapply(refseq, length)
mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq), row.names=NULL, stringsAsFactors=FALSE)
mapping <- unique(mapping)

paths <- pathways("hsapiens", "kegg")
alledges <- NULL

## creates an adjacency matrix
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

for(i in 1:nrow(alledges))
{
  col1 <- alledges[i,"src"]
  col2 <- alledges[i,"dest"]
  adjMatrix[toString(col1),toString(col2)] <- 1
  adjMatrix[toString(col2),toString(col1)] <- 1
}

matched <- matchMatrices(x = x, adjacency = adjMatrix, mapping = mapping)
ad.list <- as.adjacencyList(adjMatrix)


## networkBasedSVM
res.nBSVM <- crossval(matched$x, y, theta.fit=fit.networkBasedSVM, folds=3, repeats=1, DEBUG=TRUE,
                      parallel=FALSE, adjacencyList=ad.list, lambdas=10^(-1:2), sd.cutoff=1)


## RRFE
res.rrfe <- crossval(x, y, DEBUG=TRUE, theta.fit=fit.rrfe, folds=3, repeats=1, parallel=TRUE,
                     Cs=10^(-3:3), mapping=mapping, Gsub=adjMatrix, d=1/2)

