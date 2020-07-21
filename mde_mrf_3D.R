
f.test=function(x){
  fun.prepare=function(t){ return(ker(x,t,sig)) }
  comp.vec=apply(mu,1,fun.prepare)
  return(sum(theta.hat*comp.vec))
}

ker=function(x,mu,sig){
  yseq <- sapply((x-mu)/sig,dnorm)
  return((sig^{-length(x)})*prod(yseq))
}

rao = function(int){
  temp = table(int)
  freq = temp/sum(temp)
  int.unique = sort(unique(int))
  dist.mat = as.matrix(dist(int.unique, method = "euclidean", 
                            diag = TRUE, p = 2))
  return(t(freq)%*%dist.mat%*%freq)
}

# Euclidean norm
norm.euclidean=function(x){ return(norm(matrix(x,ncol=1),type = "F")) }
# Distance function 
dist.euclidean=function(x,y){ return(norm.euclidean(x-y)) }


### This function computes the weights in the GMM 
weight.seq.mrf=function(x.obs, coord, beta = 1, coord.neigh, mu, sig, epsilon=0.001, max.iter=1000){ 
  
  n.D=dim(x.obs)                  
  n=n.D[1]
  D=n.D[2]
  N=dim(mu)[1]  
  d.c = dim(coord)[2]
  
  A=matrix(NA,ncol = N,nrow = n)
  for(j in 1:N){
    A.prepare=function(x){ return(ker(x,mu[j,],sig))  }
    A[,j]=apply(x.obs,1,A.prepare) 
  }
  
  theta.old=matrix(1/N, n, N)            
  abs.diff=10*epsilon                
  count=0                      
  lambda.hat.old=rep(1,n)      
  
  while((abs.diff>epsilon)&(count<=max.iter)){   
    
    W=A*theta.old            
    W=W/apply(W,1,sum)       
    Gjv = matrix(0, n, N)
    
    ### Calculate the G_j,v
    for(j in 1:N){
      
      smwp = 0
      mask.w = array(0, apply(coord,2,max)+10)
      mask.w[coord] = W[,j]
      mask.pi = array(0, apply(coord,2,max)+10)
      mask.pi[coord] = theta.old[,j]
      
      for(vpr in 2:(dim(coord.neigh)[2]/d.c)){
        temp = coord.neigh[,((vpr-1)*d.c+1):(d.c*vpr)]
        smwp = smwp + (mask.w[temp] + mask.pi[temp])
      }
      Gjv[,j] = exp(beta/(2*vpr)*smwp)
    }
    
    #### Update lambda.v
    lambda.hat.old = 1 + apply(Gjv,1,sum)
    
    theta.new = (W + Gjv)/lambda.hat.old
    abs.diff=dist.euclidean(theta.new,theta.old)                        
    if(is.na(abs.diff)){ abs.diff=0; theta.new=theta.old }                 
    
    theta.old=pmax(theta.new,0)     
    theta.old=pmin(theta.old,1)    
    count=count+1
  }
  
  theta.hat=pmax(theta.new,0)      
  theta.hat=pmin(theta.hat,1)      
  return(theta.hat)                
}



### A function for choosing the N and estimated GMM weights 
est.func=function(x.obs, coord, beta = 1, coord.neigh, N0=4, sig, max.knots=100,alpha=0.05,mu=NULL){
  
  zalpha = qnorm(1-alpha/2)
  n=length(x.obs)
  ind.full = 1:max.knots
  
  ########## First set of the means and standard deviation vector
  N=N0
  if(is.na(mu[1])) {mu=seq(min(x.obs)+0.5,max(x.obs)-0.5,l=N)}
  mu.full = matrix(mu[order(mu)],length(mu),1)
  
  #### Start with N0 means out of the full mean sequence 
  ind.mu = round(seq(max.knots/(N0+2), max.knots-max.knots/(N0+2),l=N0))
  ### Current set of means
  mu = matrix(mu.full[ind.mu,],N0,1)
  ### indices of the remaining means
  ind.temp = ind.full[-ind.mu]
  ### remaining means
  mu.temp = mu.full[ind.temp,]
  #### Find the weight sequence by using the constrained EM algorithm
  theta.hat=weight.seq.mrf(x.obs, coord, beta = beta, coord.neigh, mu=mu, sig, max.iter = max.iter)
  f.test=function(x){
    fun.prepare=function(t){return(ker(x,t,sig)) }
    return(apply(mu,1,fun.prepare))
  }
  temp = t(sapply(x.obs,f.test))
  delta.old=log(apply(temp*theta.hat, 1, sum))
  
  st=1
  while((st==1)&(length(ind.temp)>0)){
    # Add means to the sequence at a time
    if(length(ind.temp)>2){
      ind.mu = ind.temp[round(seq(1, length(ind.temp),l=min(4,length(ind.temp))))[c(2,3)]]
    }
    if(length(ind.temp)<=2){
      ind.mu = ind.temp
    }
    mu = c(mu.full[ind.mu], mu)
    mu=matrix(mu[order(mu)],length(mu),1)
    # Remove these from the list of remaining indices
    ind.temp = ind.temp[-match(ind.mu, ind.temp)]
    theta.hat=weight.seq.mrf(x.obs, coord, beta = beta, coord.neigh, mu=mu, sig, max.iter = max.iter)
    f.test=function(x){
      fun.prepare=function(t){ return(ker(x,t,sig)) }
      return(apply(mu,1,fun.prepare))
    }
    temp = sapply(x.obs,f.test)
    ## Calculate the delta temp
    delta.temp = log(apply(t(temp)*theta.hat, 1, sum))
    delta.new=delta.temp-delta.old
    z=sqrt(n)*mean(delta.new)/sd(delta.new)
    if(abs(z)<zalpha){st=0}
    delta.old=delta.temp
  }
  
  resp=list(theta.hat,mu)
  return(resp)
}
