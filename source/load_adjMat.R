
load_adjMat = function()
{
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
    
    ###Remove Metabolic Pathways
    #browser()
    if(as.numeric(unlist(strsplit(curr2@id,":"))[[2]]) >1230)
    {
      edges <-graphite::edges(curr2)
      edges <- edges[,c("src","dest")]
      alledges <- rbind(alledges, edges)
      alledges = unique(alledges) ## Get the unique edges in the pathway just as a list of gene-gene interactions
    }
    
    
   
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
  
  return (list(adjMatrix,ad.list,mapping))

}