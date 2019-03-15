matched <- matchMatrices(x = x_tr, adjacency = adjMatrix, mapping = mapping)
#dimad.List <- as.adjacencyList(matched$adjacency)
nbSVM = fit.networkBasedSVM(matched$x, y_tr, DEBUG=FALSE, adjacencyList = ad.list, lambdas = 10^(-1:2),sd.cutoff=0.3)
mapping[mapping[,2] %in% nbSVM$features,]
m <- matchMatrices(x=x_te, adjacency = adjMatrix, mapping=mapping)
pred = predict(nbSVM, mdev$x)
print(pred)	
A = roc(y_te,pred)
AUC.nbsvm[i+1] = auc(A)