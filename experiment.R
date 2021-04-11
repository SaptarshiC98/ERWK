library(igraph)
library(pracma)
source('functions.R')

#### Example Run on Iris ####

data(iris)
X=iris
X=data.matrix(X)
toss=X[,5]
X=X[,-5]

p=dim(X)[2]
for(i in 1:p){
  if(sd(X[,i])==0){
    X[,i]=0
  }else{
    X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i])
  }
}
n=dim(X)[1]
k=3
S=matrix(0,20,k)
for(i in 1:20){
  S[i,]=sample(n,k)
}
indexing=matrix(rep(0,20*2),ncol=2)
for(i in 1 : 20){
  sa=S[i,]
  M=X[sa,]
  #l=kmeans(X,M,30)
  #M=l$centers
  l=erwkmeans(X,M,40,30)
  indexing[i,1]=compare(toss,l[[1]],method='nmi')
  indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
  cat(i)
  cat('\n')
}
colMeans(indexing)
