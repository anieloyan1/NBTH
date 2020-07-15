
### Compare clustering methods
library(spatstat)
library(moments)
library(ClusterR)
library(mclustcomp)
library(png)
library(furrr)
library(dplyr) # for the pipe
plan(multiprocess)

source('~/Dropbox/2020Radiomics/RcodeMay/FeatureExtractionCode/MRF/mde_mrf_3D.R')
temp = readPNG("~/Dropbox/2020Radiomics/RcodeMay/Simulations/nkar2.png")
out.dir = "~/Dropbox/2020Radiomics/RcodeMay/Simulations/ClusterComparison"

### Generate the data
Y1 = round(temp[3:190,190:3,2]*10)
Y1[Y1==9] = 10
Y1[Y1==8] = 10
Y1[Y1==7] = 6
Y1[Y1==5] = 6
Y1[Y1==4] = 3
Y1[Y1==2] = 3
Y1[Y1==1] = 0
ddim = dim(Y1)[1]

### Extract the coordinates and neighborhood coordinates
coord = which(!is.na(Y1), arr.ind = TRUE)
xx = coord
### Add column corresponding to the immediate neighbors on all sides
for(k in c(-1,1)){
  xx = cbind(xx, xx[,1]+k, xx[,2])
  xx = cbind(xx, xx[,1], xx[,2]+k)
}
### Add column corresponding to the second neighbors on all sides
for(k in c(-1, 1)){
  for(j in c(-1, 1)){
    xx = cbind(xx, xx[,1]+k, xx[,2]+j)
  }
}
xx[xx==0] <- 1
xx[xx==ddim+1] <- ddim
coord.neigh = xx

ave.dice = function(y1, y2){
  els = unique(y1)
  res = 0
  for(i in 1:length(els)){
    y1.temp = y1
    y1.temp[y1.temp != els[i]] = 0
    res = res + mclustcomp(y1.temp, y2, type = "overlap")[2]
  }
  return(res/length(els))
}

epsilon=0.001
max.iter=100
beta = 4
N0 = 4
alpha=0.05
max.N = 10

run.simulation = function(x){
  
  #### Add the random noise to the data
  set.seed(x)
  rand = Y1+matrix(rnorm(ddim*ddim, 0, 3),ddim,ddim)  
  bl.Y1 = t(as.matrix(blur(as.im(t(rand)), sigma = 1)))
  
  ### Create a matrix of intensities
  ints = matrix(bl.Y1, dim(coord)[1], 1)

  ### Calculate the size of clusters for k-means
  wss <- function(k) {
    kmeans(ints, k, nstart = 10 )$tot.withinss
  }
  # Compute wss for k = 1 to k = max.cluster.size
  k.values <- 1:max.N
  # extract wss for all clusters
  wss_values <- sapply(k.values, wss)
  # smallest cluster size that explains at least 80% of within sum of squares
  cl.size = max(which(cumsum(wss_values)/sum(wss_values)<0.80),1)
  
  max.N = max(N0, cl.size)
  ### Use k-means to estimate the means and variance for mixture model
  fit.clust = kmeans(ints, max.N)  
  ## Retrieve the clusters
  temp.clusters1 = fit.clust$cluster
  ### Extract the means of the clusters to set as mu compoments
  mu = fit.clust$centers
  
  ### Estimate the variance
  sigma.vec=rep(NA, max.N)               
  for(j in 1:max.N){                      
    index.temp=which(fit.clust$cluster==j)   
    xi.j = matrix(ints[index.temp,],length(ints[index.temp,]),1)           
    sig.prepare=function(x){ return((dist.euclidean(x,mu[j,]))^2) }
    s=apply(xi.j,1,sig.prepare)
    sigma.vec[j]=mean(s)
  }
  sig=sqrt(mean(sigma.vec)/dim(ints)[2])
  
  ##### Estimate the N and weights of the MRF model
  res.mrf = est.func(ints, coord, beta = beta, coord.neigh, N0=N0,sig = sig,max.knots=max.N,alpha=0.05,mu=mu)
  theta.mrf = res.mrf[[1]]
  
  clust.mrf1 = rep(NA, length(ints))
  for(i in 1:length(clust.mrf1)){
    clust.mrf1[i] = which(theta.mrf[i,]==max(theta.mrf[i,]))[1]
  }
  res.mrf1 = matrix(clust.mrf1, ddim, ddim)

  cl.size.gmm = Optimal_Clusters_GMM(ints, max.N, plot_data = F)
  j = 1
  d = 1
  while((d>0.01)&(j<max.N)){
    d = abs(cl.size.gmm[j+1]/sum(cl.size.gmm)-cl.size.gmm[j]/sum(cl.size.gmm))
    j=j+1
  }
  gmm = GMM(ints, j, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10, 
            em_iter = 10, verbose = F) 
  gmm.res = predict_GMM(ints, gmm$centroids, gmm$covariance_matrices, gmm$weights)
  res.gmm1 = matrix(gmm.res$cluster_labels,ddim,ddim)
  
  fit.clust = kmeans(ints, max(cl.size,2))  
  ## Retrieve the clusters
  temp.clusters1 = fit.clust$cluster
  
  m1.gmm1 = mclustcomp(as.vector(Y1),as.vector(res.gmm1),type="all")[,2]
  m1.gmm1[14] = ave.dice(as.vector(Y1), as.vector(res.gmm1))
  m1.kmeans1 = mclustcomp(temp.clusters1,as.vector(Y1),type="all")[,2]
  m1.kmeans1[14] = ave.dice(as.vector(Y1), temp.clusters1)
  m1.mrf1 = mclustcomp(as.vector(Y1),as.vector(res.mrf1),type="all")[,2]
  m1.mrf1[14] = ave.dice(as.vector(Y1), as.vector(res.mrf1))
  
  return(cbind(m1.gmm1, m1.kmeans1, m1.mrf1))

}

temp = 1:10
results = temp%>%future_map(run.simulation, .progress = T)
res.dat = as.data.frame(do.call(cbind, results))
sim.res = res.dat

for(runs in 2:10){
  temp = ((runs-1)*10+1):(runs*10)
  results = temp%>%future_map(run.simulation, .progress = T)
  res.dat = as.data.frame(do.call(cbind, results))
  sim.res = cbind(sim.res, res.dat)
  print(runs)
}
setwd(out.dir)
save(sim.res, file = "sim3out.rda")
