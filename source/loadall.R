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

## use scan with flag skip=1 instead of read.csv() but tell read.csv -> double
GSE1456 <- read.csv(paste0(getwd(),'/data/GSE1456_2.csv'), header = FALSE) ##loads the data
GSE2034 <- read.csv(paste0(getwd(),'/data/GSE2034_2.csv'), header = FALSE)
GSE2990 <- read.csv(paste0(getwd(),'/data/GSE2990_2.csv'), header = FALSE)
GSE4922 <- read.csv(paste0(getwd(),'/data/GSE4922_2.csv'), header = FALSE)
GSE7390 <- read.csv(paste0(getwd(),'/data/GSE7390_2.csv'), header = FALSE)
GSE11121 <- read.csv(paste0(getwd(),'/data/GSE11121_2.csv'), header = FALSE)

col = read.csv(paste0(getwd(),"/data/headers.csv"),header=TRUE)[1:22284]
col <- colnames(col)
## remove the for loop ##
for(i in 1:length(col)) ## removes any leading X in the header names
{
  if(substring(col[i], 1, 1) == 'X')
    col[i] <- substr(col[i], 2, nchar(col[i]))
}

## creates a mapping using ENTREZID

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