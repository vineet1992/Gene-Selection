#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

data <- read.csv('1456_2.csv', header = FALSE, sep=',') ##loads the data

#separates into positive and negative groups
pos <- data[which(data[,22284]==1),]
neg <- data[which(data[,22284]==0),]

#calculates difference of means between each group, for each feature
val= c(1:22283)
for (i in val){
	t = t.test(pos[,i], neg[,i])
	val[i] = abs(t$estimate[1] - t$estimate[2])
}
i = c(1:22283)

#orders the features indices high to low, according the difference of means
m = order(val,1:22283,decreasing=TRUE)

n = dim(data)[1]
s = floor(n/5)
x = data[,1:22284] #matrix of predictors and classes
y = data[,22284] #vector of labels
rand = sample(1:n,n) #creates a random permutation (for indexing samples randomly)
x = x[rand,]





	best_k=1 #assume using 5 features has the best accuracy
	best_acc = 0 #assume the accuracy of the best value of i is 0
	best_accSet = c(0,0,0,0,0)
	
	#performs 5-fold cross validation using i features with greatest difference in mean
	for (i in 5:50){ #set to 7 for testing script; change to desired number when actually running
		acc = c(0,0,0,0,0)
		for (j in 0:4){
			xtr = x[-c(1:s+j*s),] #excludes (j+1)th fold for cross validation
									#do note that we may leave out some samples entirely (up to 4)
			xtr = data.frame(xtr[,c(22284, m[1:i])]) #glm() wants a dataframe; glm() predicts the first column on the remaining columns
			xte = x[c(1:s+j*s),] #select the values from the left out fold
			xte = data.frame(xte[,m[1:i]]) #create a dataframe of the left out fold, for testing
			yte = x[c(1:s+j*s),22284] #vector of actual classes for the leftout fold, used for comparison to predicted classes
			model = glm(xtr, family=binomial(link="logit")) #i'm not sure if these are appropriate parameters
			p= predict(model, xte, type = "response")

			p[which(p>=.5)]=1
			p[which(p<.5)]=0
			acc[j+1] = length(which(p==yte))/s
		}
		if (mean(acc)>best_acc){
			best_k = i
			best_acc = mean(acc)
			best_accSet = acc
		}
	}

	print(best_k)
	print(best_acc)
	


s = floor(n/3) #now split 2/3 train, 1/3 test
rand = sample(1:n,n) #get a new random permutation
x = x[rand,]
xi = x[,c(22284,m[1:best_k])] #use k features according to how many previously determined best 
xtr = xi[(s+1):n,] #set up dataframe to train on
xte = xi[1:s,2:(best_k+1)] #set up dataframe to test on; do not include class labels
yte = xi[1:s,1] #vector of class labels for testing samples

model = glm(xtr, family=binomial(link="logit")) #build
p = predict(model, xte, type = "response") #predict
p[which(p>=.5)] = 1 #convert probabilities to labels
p[which(p<.5)] = 0 #convert probabilities to labels
print(length(which(p==yte))) #compare predicted and actual labels
print(length(which(p==yte))/s)




#1456:
#11121: 86.4%
#2990: 80.9%
#2034:
#4922:
#7390: 68.2%