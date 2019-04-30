lib = "/net/maccutcheon/home/benos/vkr8/Rlib"

require(ggplot2,lib=lib)
library(tibble,lib=lib)
library(tidyr,lib=lib)
library(readr,lib=lib)
library(purrr,lib=lib)
library(dplyr,lib=lib)
library(tidyverse,lib=lib)


library(svmpath,lib=lib)
library(kernlab,lib=lib)
library(Biobase,lib=lib)
library(affy,lib=lib)
library(gplots,lib=lib)
library(igraph,lib=lib)
library(lpSolve,lib=lib)




library(pathClass,lib=lib)
library(pROC,lib=lib)
require(gcdnet,lib=lib)
library(pipeliner,lib=lib)
library(modelr,lib=lib)
library(graphite,lib=lib)

library(org.Hs.eg.db,lib=lib)
library(hgu133a.db,lib=lib)



###PARAMETERS
nFolds = 5
####Set standard deviation cutoff
sdCutoff = 0.5

graphSelect = 150


model = "pdNPPCAGraph"
#names = c("pdNPGraph")


###Source Modles
modelDir = paste0(getwd(),'/source/Models/')
source(paste0(modelDir,model,".R"))

###Prepare dataset loading
dataDir = paste0(getwd(),'/data/')
datasets = list.files(dataDir)

datasets = datasets[startsWith(datasets,"GSE")]

d = datasets[2]


FILE = paste0(getwd(),'/data/',d)

print(paste0("Running dataset: ", d))

##Quickly determine the number of columns
x <- max(count.fields(FILE, ","))

###Scan the file
data <- scan(FILE,sep=",")

###Extract number of rows and columns
nr = length(data)/x
nc = x
data = matrix(data,nrow=nr,ncol=nc,byrow=T)

#data = data[,c(1:500,ncol(data))]

col = scan(paste0(getwd(),"/data/headers.csv"),sep=",",what=character())[1:ncol(data)]
## remove the for loop ##


####Cannot use this for-loop, it should be abandoned
####A legal R variable cannot start with a number, so this causes problems when using a formula object
for(i in 1:length(col)) ## removes any leading X in the header names
{
  if(substring(col[i], 1, 1) == 'X')
    col[i] <- substr(col[i], 2, nchar(col[i]))
}
colnames(data) <- col



###Linear model doesn't work if  | features | > | samples | so take a small chunk of features
target = colnames(as.data.frame(data[,ncol(data)]))[1]

#change to precision recall
#target variable sjo

# 5-fold cross-validation using machine learning pipelines



data[,ncol(data)] = factor(data[,ncol(data)])





FUN = eval(parse(text=paste0(model,"_wrapper")))

y_tr = factor(data[,ncol(data)])
x_tr = as.data.frame(data[,-ncol(data)])


mdl = FUN(x_tr,y_tr)
