
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

# Norm function in an Euclidean space of any dimension
norm.euclidean=function(x){ return(norm(matrix(x,ncol=1),type = "F")) }
# Distance function in an Euclidean space of any dimension
dist.euclidean=function(x,y){ return(norm.euclidean(x-y)) }
## Subsection 1.2, Kernels for minimization in a semi-normed space of Sobolev type

weight.seq.mrf=function(x.obs, coord, beta = 1, coord.neigh, mu, sig, epsilon=0.001, max.iter=1000){ 
  
  # "x.obs" is the data set of interest. 
  #         There are n observations, and each observation is a D-dimensional point.
  #         x.obs is a n-by-D matrix.
  # "mu" is a vector of the knots in a mixture density estimation.
  # "sigma" is the bandwidth of this density estimation.
  # "epsilon" is a predetermined tolerance of the Euclidean distance between thetas in two consecutive steps.
  # "max.iter" is a predetermined upper bound of the number of steps in this iteration.
  
  n.D=dim(x.obs)                  
  n=n.D[1]
  D=n.D[2]
  N=dim(mu)[1]  
  d.c = dim(coord)[2]
  
  A=matrix(NA,ncol = N,nrow = n)
  for(j in 1:N){
    A.prepare=function(x){ return(ker(x,mu[j,],sig))  }
    A[,j]=apply(x.obs,1,A.prepare)   # A[i,j] is \f_{j,K} (x_i).
  }
  
  theta.old=matrix(1/N, n, N)               # The initial guess for weights, say \theta_j's.
  abs.diff=10*epsilon                # Absolute value of the difference between "theta.new" and "theta.old".
  count=0                            # Counting the number of steps in this iteration.
  lambda.hat.old=rep(1,n)      # The initial guess of the Lagrangian multipliers
  
  while((abs.diff>epsilon)&(count<=max.iter)){   # The iteration for computing desired theta's
    
    W=A*theta.old              # \theta_j^{(k)} \times \psi_\sigma(x_i-mu_j)
    W=W/apply(W,1,sum)               # W[i,j] is the posterior probability of Z=j|X=x_i, say w_{i,j}(\theta.old).
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
    
    # We set the Lagrangian multipliers in the previous step 
    # as the initial guess in this step.
    theta.new = (W + Gjv)/lambda.hat.old
    abs.diff=dist.euclidean(theta.new,theta.old)                              # The Euclidean distance between the old and new theta vectors.
    if(is.na(abs.diff)){ abs.diff=0; theta.new=theta.old }                    # It helps us avoid "NA trouble".
    
    theta.old=pmax(theta.new,0)      # Set the new theta as the old theta for the next iteration step.   
    theta.old=pmin(theta.old,1)      # pmax() and pmin() guarantee that theta_j's are in [0,1]. 
    count=count+1
  }
  
  theta.hat=pmax(theta.new,0)        # The iteration has stopped.  
  theta.hat=pmin(theta.hat,1)        # Again, we guarantee that theta_j's are in [0,1].
  return(theta.hat)                  # It returns the estimation of weights \theta_j's.
}




# A function that returns the mean vector, weight vector and sigmaN for a given sample from the true density.
# The starting value of the sample size is given by the user as well as the corresponding standard deviation.
# The estimated density can be computed by using the results of this function in a mixture density.

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
    # Add two means to the sequence at a time
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
