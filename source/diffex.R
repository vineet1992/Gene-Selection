if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(gcdnet)) install.packages('gcdnet')
library(gcdnet)
if(!require(pROC)) install.packages('pROC')
library(pROC)
setwd("./")

#sampleData <- read.csv('./data/GSE1456_2.csv', header = FALSE, sep=',') ##loads the data
sampleData <- scan('./data/GSE11121_2.csv',  sep=',')
#sampleData <- read.csv('./data/GSE2990_2.csv', header = FALSE, sep=',')
#sampleData <- read.csv('./data/GSE2034_2.csv', header = FALSE, sep=',')
#sampleData <- read.csv('./data/GSE4922_2.csv', header = FALSE, sep=',')
#sampleData <- read.csv('./data/GSE7390_2.csv', header = FALSE, sep=',')


n = dim(sampleData)[1]
rand = sample(1:n,n)
sampleData = sampleData[rand,] #randomly shuffle data
f = floor(n/10)

AUC = c(1:10) #a vector to store the AUC for each test fold
selected = matrix(,nrow=10,ncol=50) #a matrix to store indices of selected predictors for each fold; used for comparison in computing stability

for (i in 0:9){ #for each of 10 folds
	OUT = sampleData[c(1:f+f*i),] #leave one out for testing
	IN = sampleData[-c(1:f+f*i),] #use the rest as a training set
	
	s = floor(dim(IN)[1]/5)	 #to pick best k, will do 5-fold cross validation
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
	for(k in 5*1:10){
		acc = c(0,0,0,0,0) #init a vector to store the AUC for each fold
		for (j in 0:4){	
			m = ind[j+1,1:k]
			TRAIN = IN[-c(1:s+s*j),] #these will be the same sets we computed difference of means on, we just did it above so that we didnt have to do it inside the loop. this shiz takes long enough already
			TEST = IN[c(1:s+j*s),]
			
			x = data.frame(TRAIN[,c(22284,m)])	#glm wants a data frame, and predicts the first column on the remaining columns, so put y first and the best predictors (indexed by m) of x
			model = glm(x, family=binomial(link="logit")) #train! #i'm not sure if these are appropriate parameters
			
			x = data.frame(TEST[,m]) #now we need to exclude the labels; just use the indexed best predictors from the cv test fold
			y = TEST[,22284]
			p = predict(model, x, type = "response")
			
			A = roc(y,p)	#create an ROC object
			acc[j+1]=auc(A)	#get area under curve from it, and record in the array we init'd
		}
		
		if(mean(acc)>best_acc){		#once we have an AUC for each cv test fold, compute their average and update the best as needed
			best_k=k
			best_acc=mean(acc)}
	}

	
	#equipped with out best_k, we basically do the same model building and testing as above, but now using the 9 folds from the whole data set to predict the 1 left out fold
	POS = IN[which(IN[,22284]==1),] 
	NEG = IN[which(IN[,22284]==0),]
	m = c(1:22283)
	for (d in 1:22283){
		m[d] = abs(mean(POS[,d])-mean(NEG[,d]))		
	}
	m = order(m, 1:22283, decreasing=TRUE)[1:best_k]
	
	x = data.frame(IN[,c(22284,m)])
	model = glm(x, family=binomial(link="logit"))
	
	x = data.frame(OUT[,m])
	p = predict(model, x, type="response")
	y = OUT[,22284]
	A = roc(y,p)
	AUC[i+1] = auc(A)
	for(z in 1:length(m)){ #each fold may not have the same best_k, so we would get length mismatches if we tried something like AUC[i+1] = m
		selected[i+1,z] = m[z]
	}
}


print(mean(AUC)) #this gives the avergae AUC. we could probably do more statistics on it, unless whatever plot we use can just accept this as is

I = selected[1,] #we'll intersect all the selected predictors. we start with the intersection of a single set
U = selected[1,] #similarly, we union 
for(i in 2:10){		#go through each set and intersect and union
	I = intersect(I, selected[i,])
	U = union(U, selected[i,])
}
I = I[which(I!='NA')] #remove any NA
U = U[which(U!='NA')]

print(length(I)/length(U))