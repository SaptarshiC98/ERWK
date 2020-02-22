msq.euc.dist= function(x1, x2,a) (sum((abs(x1 - x2)) ^ a))^(1/a)


mwt.euc.dist.sq=function(x1,x2,a,w){
  p=(abs(x1-x2))^a
  p=w*p
  return(sum(p))
}

mvec.wt.euc.dist.sq=function(x1,x2,a,w){
  p=(abs(x1-x2))^a
  p=w*p
  return(p)
}

mE=function(mu,Y,b){
  n=dim(Y)[1]
  s=0
  for(i in 1:n){
    s=s+msq.euc.dist(Y[i,],mu,b)^b
  }
  return(s)
}


mwkmeans=function(X,M,a,beta,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  weight=rep(1/d,d)
  label=numeric(n)
  dist=numeric(c)
  t=0
  D=numeric(d)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=mwt.euc.dist.sq(X[i,],M[j,],a,weight^beta)
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=colMeans(X[I,])#optim(colMeans(X[I,]),mE,Y=X[I,],b=a)$par
    }
    
    #update weights
    for(j in 1:d){
      D[j]=0
    }
    for(i in 1:c){
      I=which(label==i)
      for(k in I){
        D=D+mvec.wt.euc.dist.sq(X[k,],M[i,],a,rep(1,d))
      }
    }
    
    for(i in 1:d){
      if(D[i]!=0){
        D[i]=1/D[i]
        D[i]=D[i]^(1/(beta-1))
      }
    }
    sum=sum(D)
    weight=D/sum
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}

wkmeans=function(X,M,beta,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  weight=rep(1/d,d)
  label=numeric(n)
  dist=numeric(c)
  t=0
  D=numeric(d)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=wt.euc.dist.sq(X[i,],M[j,],weight^beta)
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=colMeans(X[I,])
    }
    
    #update weights
    for(j in 1:d){
      D[j]=0
    }
    for(i in 1:c){
      I=which(label==i)
      for(k in I){
        D=D+vec.wt.euc.dist.sq(X[k,],M[i,],rep(1,d))
      }
    }
    
    for(i in 1:d){
      if(D[i]!=0){
        D[i]=1/D[i]
        D[i]=D[i]^(1/(beta-1))
      }
    }
    sum=sum(D)
    weight=D/sum
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}
euc.dist.sq=function(x1,x2){
  p=(x1-x2)^2
  return(sum(p))
}

cpdkmeans=function(X,M,lambda,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  weight=rep(1/d,d)
  label=numeric(n)
  dist=numeric(c)
  t=0
  D=numeric(d)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=wt.euc.dist.sq(X[i,],M[j,],weight)
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=colMeans(X[I,])
    }
    
    #update weights
    for(j in 1:d){
      D[j]=0
    }
    for(i in 1:c){
      I=which(label==i)
      for(k in I){
        D=D+vec.wt.euc.dist.sq(X[k,],M[i,],rep(1,d))
      }
    }
    
    for(i in 1:d){
      
      D[i]=exp(-D[i]/lambda)
      
      
    }
    sum=sum(D)
    weight=D/sum
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}



HWkmeans=function(X,M,lambda,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  weight=matrix(1/d,nrow=c,ncol=d)
  label=numeric(n)
  dist=numeric(c)
  t=0
  D=matrix(0,nrow=c,ncol=d)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=wt.euc.dist.sq(X[i,],M[j,],weight)
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=colMeans(X[I,])
    }
    
    #update weights
    for(l in 1:c){
      for(j in 1:d){
        D[l,j]=0
      }
    }
    for(i in 1:c){
      I=which(label==i)
      for(k in I){
        D[i,]=D[i,]+vec.wt.euc.dist.sq(X[k,],M[i,],rep(1,d))
      }
    }
    
    for(i in 1:d){
      for(l in 1:c){
        D[l,i]=exp(-D[l,i]/lambda)
      }
    }
    sum=numeric(c)
    for(l in 1:c){
      sum[l]=0
      for(j in 1:d){
        sum[l]=sum[l]+D[l,j]
        j=j+1
      }
    }
    for(l in 1:c){
      for(j in 1:d){
        weight[l,j]=D[l,j]/sum[l]
      }
    }
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}

euc.dist.sq=function(x1,x2){
  p=(x1-x2)^2
  return(sum(p))
}


k.means=function(X,M,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  label=numeric(n)
  dist=numeric(c)
  t=0
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=euc.dist.sq(X[i,],M[j,])
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=colMeans(X[I,])
    }
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M))
  
}
sq.euc.dist= function(x1, x2) sum((x1 - x2) ^ 2)

wt.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(sum(p))
}

vec.wt.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(p)
}

k=8
p=dim(X)[2]
for(i in 1:p){
  if(sd(X[,i])==0){
    X[,i]=0
  }else{
    X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i])
  }
}
n=dim(X)[1]
S=matrix(0,20,k)
for(i in 1:20){
  S[i,]=sample(n,k)
}

toss=c(rep(1,100),rep(2,100))#,rep(3,100),rep(4,100))

# k-means
#toss=c(rep(1,64),rep(2,64),rep(3,64),rep(4,64))#rep(5,64),rep(6,64),rep(7,64),rep(8,64))
indexing=matrix(rep(0,10*2),ncol=2)
for(i in 1 : 10){
  sa=S[i,]
  M=X[sa,]
  l=kmeans(X,M,30)
  indexing[i,1]=compare(toss,l[[1]],method='nmi')
  indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
  cat(i)
  cat('\n')
}
colMeans(indexing)

# CPD
indexing=matrix(rep(0,10*2),ncol=2)
for(i in 1 : 10){
  sa=S[i,]
  M=X[sa,]
  #l=kmeans(X,M,30)
  #M=l$centers
  l=cpdkmeans(X,M,40,30)
  indexing[i,1]=compare(toss,l[[1]],method='nmi')
  indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
  cat(i)
  cat('\n')
}
colMeans(indexing)
plot(X,col=l[[1]],pch=l[[1]]+14,cex=1.1,xlab = 'Feature 1',ylab='Feature 1')
plot(X,col=l$cluster,pch=l$cluster+14,cex=1.1)
hist(l[[3]])
l# Huang
indexing=matrix(rep(0,10*2),ncol=2)
for(i in 1 : 10){
  sa=S[i,]
  M=X[sa,]
  l=HWkmeans(X,M,1,30)
  indexing[i,1]=compare(toss,l[[1]],method='nmi')
  indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
  cat(i)
  cat('\n')
}
colMeans(indexing)
#   WK-means
indexing=matrix(rep(0,10*2),ncol=2)
for(i in 1 : 10){
  sa=S[i,]
  M=X[sa,]
  #l=kmeans(X,M,30)
  #M=l$centers
  l=wkmeans(X,M,1.5,20)
  indexing[i,1]=compare(toss,l[[1]],method='nmi')
  indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
  cat(i)
  cat('\n')
}
colMeans(indexing)

indexing=matrix(rep(0,10*2),ncol=2)
for(i in 1 : 10){
  sa=S[i,]
  M=X[sa,]
  #l=kmeans(X,M,30)
  #M=l$centers
  l=mwkmeans(X,M,4,1.5,30)
  indexing[i,1]=compare(toss,l[[1]],method='nmi')
  indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
  cat(i)
  cat('\n')
}
colMeans(indexing)
