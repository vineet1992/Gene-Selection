####Conert graph file with selected genes to gene symbols or directly to top pathways?###

source("http://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
library(clusterProfiler)

graph = "C:/Users/vinee_000/Desktop/Gene-Selection/Interpretability_Results/GSE1456_WP/graph.txt"



data = read.table(graph,sep=" ", skip=4)

connections = data[data$V2=="y" | data$V4=="y",]


variables = connections$V2

variables = as.character(variables)

require(hgu133a.db)

mapped.probes <- mappedkeys(hgu133aSYMBOL)
names <- as.list(hgu133aSYMBOL[mapped.probes])
names = unlist(names)

for(i in 1: length(variables))
{
  curr = unlist(strsplit(variables[i],"|",fixed=T))
  symbols = unlist(names)[curr]
  symbols = unique(symbols[!is.na(symbols)])
  eg = bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ggo <- enrichGO(gene     = eg$ENTREZID,
                  keyType="ENTREZID",
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 pvalueCutoff=0.1,
                 readable = TRUE)
  head(ggo)
  
  kk <- enrichKEGG(gene         = eg$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  
  head(kk)
  mkk <- enrichMKEGG(gene = eg$ENTREZID,
                     organism = 'hsa')
  
  head(mkk)
  do <- enrichDO(gene     = eg$ENTREZID,
                  pvalueCutoff=0.1,
                  readable = TRUE)
  
  head(do)

}

