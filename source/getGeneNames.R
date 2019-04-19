require(hgu133a.db)

mapped.probes <- mappedkeys(hgu133aSYMBOL)
names <- as.list(hgu133aSYMBOL[mapped.probes])
col = scan(paste0(getwd(),"/data/headers.csv"),sep=",",what=character())[1:ncol(data)]
for(i in 1:length(col)) ## removes any leading X in the header names
{
  if(substring(col[i], 1, 1) == 'X')
    col[i] <- substr(col[i], 2, nchar(col[i]))
}
names = unlist(names)
symbols = unlist(names)[col]

names(symbols)=col

symbols[length(symbols)] = "y"

write.table(symbols,file="Affy_Gene_Symbols.txt",quote=F,sep='\t',row.names=F,col.names=T)

symbols = symbols[c(1:500,length(symbols))]
write.table(symbols,file="Affy_Symbols_500.txt",quote=F,sep='\t',row.names=F,col.names=T)
