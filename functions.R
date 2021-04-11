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
euc.dist.sq=function(x1,x2){
  p=(x1-x2)^2
  return(sum(p))
}

erwkmeans=function(X,M,lambda,tmax){
  
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

