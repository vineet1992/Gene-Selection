data <- read.table("C:/Users/Marcus/Desktop/eset_highest_variance.txt", header=TRUE, sep='\t')
x <- (data)
y <- factor(data$Response)

library(hgu133a.db)
mapped.probes <- mappedkeys(hgu133aREFSEQ)
refseq <- as.list(hgu133aREFSEQ[mapped.probes])
times <- sapply(refseq, length)
mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq), row.names=NULL, stringsAsFactors=FALSE)
mapping <- unique(mapping)

library(graphite)
paths <- pathways("hsapiens", "kegg")
alledges <- NULL

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

library(pathClass)
matched <- matchMatrices(x=x, adjacency=adjMatrix, mapping=mapping)
ad.list <- as.adjacencyList(adjMatrix)
res.nBSVM <- crossval(matched$x, y, theta.fit=fit.networkBasedSVM, folds=3, repeats=1, DEBUG=TRUE,
                      parallel=FALSE, adjacencyList=ad.list, lambdas=10^(-1:2), sd.cutoff=50)
