stability_wrapper = function(list1,list2){
  final_scores <- matrix(nrow=length(list1),ncol=length(list2))
  #print(ncol(final_scores))
  #print(nrow(final_scores))
  
  i<-1
  for(itemlist1 in list1){
    j<-1
    for(itemlist2 in list2){
      
      il1 = unlist(itemlist1)
      #print(il1[0:10])
      il2 = unlist(itemlist2)
      #print(il2[0:10])
      #print("intersect:")
      #print(length(intersect(il1,il2)))
      #print("union:")
      #print(length(union(il1,il2)))
      #print("score:")
      #print(length(intersect(il1,il2))/length(union(il1,il2)))
      #print("index:i")
      #print(i)
      #print("index:j")
      #print(j)
      final_scores[i,j] <- length(intersect(il1,il2))/length(union(il1,il2))
      #print(final_scores[i,j])
      #final_scores[j,i] = final_scores[i,j]
      j<-j+1
    }
    i<- i+1
  }
  if(nrow(final_scores) > ncol(final_scores))
    final_scores <- t(final_scores)
  #print(final_scores)
  rowscols <- solve_LSAP(final_scores,maximum=TRUE) 
  row<-1
  #print(c(rowscols))
  final_score<-0
  for(item in rowscols){
    #print(row)
    #print(as.numeric(item))
    final_score <- final_score + final_scores[row,item]
    row = row+1
  }
  final_score <- final_score/(row-1)
  
  return(final_score)
}

nfolds<-5
for(i in (1:(nfolds-1))){
  #print(i)
  vari <- cv_rmse$select_features[[i]]$data
  vari <- colnames(vari)
  vari <- vari[-length(vari)]
  sf1 = strsplit(vari,".",fixed=TRUE)
  for(j in ((i+1):nfolds)){
    #print(j)
    varj <- cv_rmse$select_features[[j]]$data
    varj <- colnames(varj)
    varj <- varj[-length(varj)]
    sf2 = strsplit(varj,".",fixed=TRUE)
    final_stability_scores <-  stability_wrapper(sf1,sf2)
    #stability[d,1] = rep(d)
    #stability[d,2] = final_stability_scores
    #print(final_stability_scores)
  }
}