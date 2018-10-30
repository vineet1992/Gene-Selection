#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

data <- read.csv('GSE7390_2.csv', header = FALSE, sep=',') ##loads the data

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

best_i=1 #assume using 5 features has the best accuracy
overall_acc = 0 #assume the accuracy of the best value of i is 0
#performs 5-fold cross validation using i features with greatest difference in mean
for (i in 5:7){ #set to 7 for testing script; change to desired number when actually running
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
		print(length(which(p>=.5)))
		#1
		#2
	}
	#3
	print("----")
}

#4



##work still to do:

#1will need to convert probabilities to classes, and compare predicted and actual classes for accuracy
#2record accuracy, and average at the end. 
#3compare this accuracy to the as of yet best accuracy from using i features, and update best_i and overall_acc as needed
#4split data into training and test set
##train glm on train set
##use to get predicted classes
##compare to actual classes