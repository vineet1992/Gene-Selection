data <- read.table("C:/Users/Marcus/Desktop/eset_highest_variance.txt", header=TRUE, sep='\t')
x <- (data)
y <- factor(data$Response)

library(hgu133a.db)
mapped.probes <- mappedkeys(hgu133aREFSEQ)
refseq <- as.list(hgu133aREFSEQ[mapped.probes])
times <- sapply(refseq, length)
mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq),
row.names=NULL, stringsAsFactors=FALSE)

library(pathClass)
data(adjacency.matrix)
matched <- matchMatrices(x=x, adjacency=adjacency.matrix, mapping=mapping)
ad.list <- as.adjacencyList(matched$adjacency)
res.nBSVM <- crossval(matched$x, y, theta.fit=fit.networkBasedSVM, folds=3, repeats=1, DEBUG=TRUE,
                      parallel=FALSE, adjacencyList=ad.list, lambdas=10^(-1:2), sd.cutoff=50)
