if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require(hgu133a.db)) BiocManager::install(hgu133a.db)
library(hgu133a.db)
if(!require(pathClass)) install.packages('../dependencies/pathClass_0.9.4.tar.gz')
library(pathClass)
setwd("./")
data <- scan("./data/GSE11121_2.csv", sep='\t',skip=1)
x <- (data)
y <- factor(data$Response)

mapped.probes <- mappedkeys(hgu133aREFSEQ)
refseq <- as.list(hgu133aREFSEQ[mapped.probes])
times <- sapply(refseq, length)
mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq),
row.names=NULL, stringsAsFactors=FALSE)

data(adjacency.matrix)
matched <- matchMatrices(x=x, adjacency=adjacency.matrix, mapping=mapping)
ad.list <- as.adjacencyList(matched$adjacency)
res.nBSVM <- crossval(matched$x, y, theta.fit=fit.networkBasedSVM, folds=3, repeats=1, DEBUG=TRUE,
                      parallel=FALSE, adjacencyList=ad.list, lambdas=10^(-1:2), sd.cutoff=50)
