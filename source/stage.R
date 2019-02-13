## working through the links Vineet provided to download, unzip, convert to expression data for the 6 datasets, then save as a .csv 
## sometimes when I it it tells me that one study or another doesnt have supplemental files, but if I run it again I'll usually get them?

##I havent yet included the clinical data from each study
##one issue I am still having/ having to work through is an error I get for GSE1456 ---
##		Error in read.celfile.probeintensity.matrices(filenames = filenames, cdfInfo = cdfInfo,  : 
##  		Cel file C:/Users/Dan/Documents/GSE1456_data/GSM107232.CEL does not seem to be of HG-U133A type

##I guess the .cel files are of different types and the function I found only runs on one of those types? Looking at them on NCBI, indeed half of them are HG-U133B

library(GEOquery)
library(affy)
datasets <- c('./data/GSE2034','./data/GSE1456','./data/GSE2990','./data/GSE4922','./data/GSE7390','./data/GSE11121')
for (i in datasets){
	getGEOSuppFiles(i)
	
	FP <- paste(i,"/",i,"_RAW.tar",sep="")
	untar(FP, exdir=paste(i,"_data/",sep=""))
	
	FP <- paste(i,"_data/",sep="")
	cels <- list.files(FP, pattern="[gz]")
	
	FP <- paste(i, "_data/", cels, sep="")
	sapply(FP, gunzip)	
}

for (i in datasets){
	setwd(paste("./",i,"_data",sep=""))
	eset <- justRMA()
	write.csv(eset, paste(i,".csv",sep=""))
	setwd("..")
}