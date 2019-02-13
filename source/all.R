if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(pathClass)) install.packages('../dependencies/pathClass_0.9.4.tar.gz')
library(pathClass)
if(!require(pROC)) install.packages('pROC')
library(pROC)
setwd("./")
##before running this, run loadall.R to set up the datasets and the prior knowledge adjacency matrix

#sampleData = GSE1456
#sampleData = GSE2034
#sampleData = GSE2990
#sampleData = GSE4922
#sampleData = GSE7390
#sampleData = GSE11121
sampleData = rbind(GSE1456,GSE2034,GSE2990,GSE4922,GSE7390,GSE11121)

n = dim(sampleData)[1]
f = floor(n/10)
sampleData = sampleData[sample(1:n,n),]

AUC.nbsvm = c(1:10)
AUC.rrfe = c(1:10)
AUC.diffex = c(1:10)

ad.list <- as.adjacencyList(adjMatrix)

for (i in 0:9){
print(i)
	OUT = sampleData[c(1:f+f*i),] #leave one out for testing
	IN = sampleData[-c(1:f+f*i),] #use the rest as a training set
	
	colnames(IN) <- col
	colnames(OUT) <- col

	x_tr <- (IN)
	y_tr <- factor(IN$y)
	x_te <- (OUT)
	y_te <- factor(OUT$y)
	
	
##-------------------NBSVM-------------------------------
	matched <- matchMatrices(x = x_tr, adjacency = adjMatrix, mapping = mapping)
	nbSVM = fit.networkBasedSVM(matched$x, y_tr, DEBUG=FALSE, adjacencyList = ad.list, lambdas = 10^(-1:2),sd.cutoff=1)
	
	m <- matchMatrices(x=x_te, adjacency = adjMatrix, mapping=mapping)
	pred = predict(nbSVM, m$x)
print(pred)	
	A = roc(y_te,pred)
	AUC.nbsvm[i+1] = auc(A)
##-------------------RRFE--------------------------------
	rrfe = fit.rrfe(x_tr, y_tr, Cs = 10^(-3:3), mapping=mapping, Gsub = adjMatrix, d=1/2)
	pred = predict(rrfe, x_te, type="response")
print(pred)
	A = roc(y_te,pred)
	AUC.rrfe[i+1] = auc(A)
##-------------------t test------------------------------
	ind = matrix(,nrow=5,ncol=50) #this will store the top 50 predictors with greatest difference in means between positive and negative samples
	for(j in 0:4){
		TEST = IN[c(1:s+j*s),] #similarly, leave one out for testing. this time we treat only the training set as if it were the whole data set
		TRAIN = IN[-c(1:s+s*j),] #and use rest for training
		POS = TRAIN[which(TRAIN[,22284]==1),]
		NEG = TRAIN[which(TRAIN[,22284]==0),]

		m = c(1:22283)
		for (d in 1:22283){	#calculate the difference in means between pos and neg for each predictor
			m[d] = abs(mean(POS[,d])-mean(NEG[,d]))		
		}
		ind[j+1,] = order(m, 1:22283, decreasing=TRUE)[1:50] #order their indices highest to lowest, and select only the top 50
	}

	best_k=0
	best_acc=0
	## remove this loop,vectorize ##
	for(k in 12:20){
		acc = c(0,0,0,0,0) #init a vector to store the AUC for each fold
		for (j in 0:4){	
			m = ind[j+1,1:k]
			TRAIN = IN[-c(1:s+s*j),] #these will be the same sets we computed difference of means on, we just did it above so that we didnt have to do it inside the loop. this shiz takes long enough already
			TEST = IN[c(1:s+j*s),]
			
			xx = data.frame(TRAIN[,c(22284,m)])	#glm wants a data frame, and predicts the first column on the remaining columns, so put y first and the best predictors (indexed by m) of x
			model = glm(xx, family=binomial(link="logit")) #train! #i'm not sure if these are appropriate parameters
			
			xx = data.frame(TEST[,m]) #now we need to exclude the labels; just use the indexed best predictors from the cv test fold
			yy = TEST[,22284]
			pred = predict(model, xx, type = "response")
			
			A = roc(yy,pred)	#create an ROC object
			acc[j+1]=auc(A)	#get area under curve from it, and record in the array we init'd
		}
			
		if(mean(acc)>best_acc){	#once we have an AUC for each cv test fold, compute their average and update the best as needed
			best_k=k
			best_acc=mean(acc)}
	}
			
	#equipped with out best_k, we basically do the same model building and testing as above, but now using the 9 folds from the whole data set to predict the 1 left out fold
	POS = IN[which(IN[,22284]==1),] 
	NEG = IN[which(IN[,22284]==0),]
	m = c(1:22283)
	for (d in 1:22283)	m[d] = abs(mean(POS[,d])-mean(NEG[,d]))
	m = order(m, 1:22283, decreasing=TRUE)[1:best_k]
	
	xx = data.frame(IN[,c(22284,m)])
	model = glm(xx, family=binomial(link="logit"))
	
	xx = data.frame(OUT[,m])
	pred = predict(model, xx, type="response")
print(pred)
	yy = OUT[,22284]
	A = roc(yy,pred)	
	AUC.diffex[i+1] = auc(A)
}

AUC = cbind(AUC.nbsvm, AUC.rrfe, AUC.diffex)

png(paste0(getwd() + "/images/plot.png"))
boxplot(AUC, ylab = "AUC", ylim = c(0,1), range=1, names = c("nbSVM","RRFE","t-test"))
dev.off()