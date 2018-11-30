#this will load each dataset, and set up the prior knowledge adjacency matrix

GSE1456 <- read.csv('GSE1456_2.csv', header = FALSE, sep=',') ##loads the data
GSE2034 <- read.csv('GSE2034_2.csv', header = FALSE, sep=',')
GSE2990 <- read.csv('GSE2990_2.csv', header = FALSE, sep=',')
GSE4922 <- read.csv('GSE4922_2.csv', header = FALSE, sep=',')
GSE7390 <- read.csv('GSE7390_2.csv', header = FALSE, sep=',')
GSE11121 <- read.csv('GSE11121_2.csv', header = FALSE, sep=',')

col = read.csv("headers.csv",header=TRUE,sep=",")[1:22284]
col <- colnames(col)
for(i in 1:length(col)) ## removes any leading X in the header names
{
  if(substring(col[i], 1, 1) == 'X')
    col[i] <- substr(col[i], 2, nchar(col[i]))
}

## creates a mapping using ENTREZID
library(hgu133a.db)
mapped.probes <- mappedkeys(hgu133aENTREZID)
refseq <- as.list(hgu133aENTREZID[mapped.probes])
times <- sapply(refseq, length)
mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq), row.names=NULL, stringsAsFactors=FALSE)
mapping <- unique(mapping)

library(graphite)
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