if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(graphite)) install.packages('graphite')
library(gcdnet)
setwd("./")


data <- scan('./data/GSE2990_2.csv', sep=',') ##loads the data
n <- dim(data)[1] ##this is a count of how many samples
lim <- floor(n/3) ##this gives a value 1/3 of the count, to be used to create training and testing sets
x <- data[,1:22283] ##separates out the predictors
y <- data[,22284] ##creates a vector of labels; already is of type integer?

rand <- sample(1:n, n) ##creates a random permutation on the number of samples, to be used as indices to get random training/testing data
x_tr <- x[rand[(lim+1):n],] ##split into train and test sets
x_te <- x[rand[1:lim],]

#this uses data exactly as read in from file
#y_tr <- y[rand[(lim+1):n]] ##split into corresponding train and test sets
#y_te <- y[rand[1:lim]]
#m <- cv.gcdnet(x_tr, y_tr, nfolds=10, pred.loss = "misclass") ##trains a model
#pclass <- predict(m$gcdnet.fit, x_te, s=m$lambda.1se) ##uses the trained model to make predictions on unseen data

#this converts y to two-level factor before training model
yf <- as.factor(y)
levels(yf) = c("A","B")
yf_tr <- yf[rand[(lim+1):n]] ##split into corresponding train and test sets
yf_te <- yf[rand[1:lim]]

x_tr = as.matrix(x_tr)

m2 =cv.gcdnet(x_tr, yf_tr, method ="hhsvm",lambda2=0.01, pred.loss="misclass", nfolds=5)

x_te = as.matrix(x_te)
plot(m2)
preds = predict(m2$gcdnet.fit,x_te,s=m2$lambda.min,type="link")
preds
probs = 1/(1+exp(-preds)) ##These contain the probability of predicting the positive class