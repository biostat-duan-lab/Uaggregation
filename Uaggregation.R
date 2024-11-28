library(tidyverse)
library(mvtnorm)
library(RSpectra)
library(irlba)

####AMP function#########################################################
amp.svd <- function(data, v0=v0, n.ite=10, percent_cut=0.9){
  
  n=dim(data)[1]
  d=dim(data)[2]
  
  Pj <- function(x, percent_cut){
    temp = if_else(abs(x) >= quantile(abs(x),percent_cut), abs(x), 0)
    temp = sign(x)*temp
    return(temp)
  }
  
  u.amp=matrix(ncol=n.ite, nrow = n)
  u.th.amp=matrix(ncol=n.ite, nrow = n)
  v.amp=matrix(ncol=n.ite, nrow = d)
  
  u.amp[,1] = rep(0,n)
  v.amp[,1] = v0*sqrt(d)
  v.amp[,2] = v.amp[,1]
  
  
  for(tt in 2:(n.ite-1)){
    
    u.amp[,tt] = data %*% v.amp[,tt]-(d/n)*Pj(u.amp[,tt-1],percent_cut)
    v.amp[,tt+1] = t(data) %*% Pj(u.amp[,tt],percent_cut) - sum(Pj(u.amp[,tt],percent_cut)!=0)*v.amp[,tt]/n
    u.th.amp[,tt] = Pj(u.amp[,tt],percent_cut)
  }
  
  return(list(v.amp, u.th.amp))
}


#cross validation for choosing sparsity 
cross_validate_percent_cut = function(Y,n.folds=5,r=1){
  K=n.folds
  d = dim(Y)[2]
  folds = sample(1:n.folds,d,replace = T)
  sparseall = (1:9)/10
  loss = rep(NA,length(sparseall))
  for (s in 1:length(sparseall)) {
    loss.k = 0
    for (k in 1:K) {
      Y_train = Y[,folds!=k]
      Y_test = Y[,folds==k]
      percent_cut = 1-sparseall[s]
      Ybar = diag(1/sqrt(rowSums(Y_train^2)))%*%Y_train
      Y_std = Ybar
      Y_svd <- svd(Y_std)
      theta <- Y_svd$d
      eta <- median(theta)
      U <- Y_svd$u
      V <- Y_svd$v
      g1 <- map_dbl(1:nrow(U), function(x){sum(eta/(theta^2 + eta^2)*U[x,]^2)})
      g2 <- map_dbl(1:nrow(V), function(x){1/eta+sum((eta/(theta^2 + eta^2)-1/eta)*V[x,]^2)})
      h <- (1/g1-eta)/sqrt(n-eta*sum(abs(g1)))
      f <- (1/g2-eta)/sqrt(d-eta*sum(abs(g2)))
      Sigmainverse = diag(d^(-1/4)*(h)^(-1/2))
      Gammainverse = diag(d^(-1/4)*(f)^(-1/2))
      Ytilde = Sigmainverse%*%Ybar%*%Gammainverse
      svd.Y = svds(Ytilde,k=r)
      Yr = (svd.Y$d[1]) * svd.Y$u[,1] %*%  t(svd.Y$v[,1])
      hatY = diag((h)^(1/2))%*% Yr %*% diag((f)^(1/2))
      hatv = svd(hatY)$v[,1]
      hatu = hatY%*%hatv
      hatu = if_else(abs(hatu) >= quantile(abs(hatu),percent_cut), abs(hatu), 0)
      hatu = hatu/sqrt(sum(hatu^2))
      Y_test = diag(1/sqrt(rowSums(Y_test^2)))%*%Y_test
      #evaluate testing reconstruction error
      v_test = t(Y_test)%*%(hatu)
      loss.k = sum(((hatu)%*%t(v_test)-Y_test)^2)+loss.k
    }
    loss[s] = loss.k
  }
  sparseall[which.min(loss)]
}



##############Generate Data############################
d = 1000
r=1
set.seed(1234)
perset = 1
heterset = 1
noiseset = 0
if(noiseset==1){set.seed(1000)}
v0 = runif(d,-1,1)
v0 = v0/sqrt(sum(v0^2)) #true V
n = 50 #n=50,100,200
sparsity = 0.3 #0.7,0.5,0.3,0.1
f0 = c(runif(d/2,0.2,1),runif(d/2,1,2))
u0 = c(rep(1,sparsity*n),rep(0,n-sparsity*n))
u0 = u0/sqrt(sum(u0)*d/n)
sigma0 = rep(1,n)
m0 = sum(u0==0)
sigma0[u0==0] = c(runif(m0,1,5))
m1 = sum(u0!=0)
sigma0[u0!=0] = c(runif(2,0.1,0.9),runif(m1-2,0.9,2))
if(heterset == 0){sigma0 =rep(1,n)}
if(noiseset == 0){w0 <- matrix(rnorm(d*n,0,1/sqrt(d)),nrow = n)}
if(noiseset == 1){w0 <- matrix(runif(d*n,-0.11/2,0.11/2),nrow = n)}
Y = u0%*%t(v0)*3+diag(sigma0)%*%w0%*%diag(f0)
Ybar = diag(1/sqrt(rowSums(Y^2)))%*%Y

#simple average
vhatnaive = colSums(Ybar)


##PCA based ensemble 
w = (t(Ybar))
G = t(w)%*%w
u = eigs(G,k=1)$vectors[,1]
u = abs(u)
v_pca = w%*%u

# perform HeteroSVD
## initial with diagonal-deletion rank-one SVD
pgs.similarity <- G
diag(pgs.similarity) <- 0
rst_SVD_woDiag <- irlba(pgs.similarity, 
                        # right_only = TRUE,
                        nv = 1)
SVD_niters <- 0
while (1|SVD_niters<=100) {
  SVD_niters <- SVD_niters + 1
  v_pre <- rst_SVD_woDiag$v
  # replace the diagonal entries with the ones of rank-one SVD
  diag(pgs.similarity) <- diag(
    rst_SVD_woDiag$u %*% t(rst_SVD_woDiag$v) * rst_SVD_woDiag$d)
  rst_SVD_woDiag <- irlba(pgs.similarity, 
                          # right_only = TRUE,
                          nv = 1)
  v_pos <- rst_SVD_woDiag$v
  if(norm(v_pos - v_pre, type = '2')<1e-5)break
}
v_heteropca <- w %*% v_pos

#estimate sparsity 
estimated_s = cross_validate_percent_cut(Y)
#######uaggregation-o############################################
percent_cut = 1-sparsity
####Biwhitening######
Ybar = diag(1/sqrt(rowSums(Y^2)))%*%Y
Y_std = Ybar
Y_svd <- svd(Y_std)
theta <- Y_svd$d
eta <- median(theta)

U <- Y_svd$u
V <- Y_svd$v
g1 <- map_dbl(1:nrow(U), function(x){sum(eta/(theta^2 + eta^2)*U[x,]^2)})
g2 <- map_dbl(1:nrow(V), function(x){1/eta+sum((eta/(theta^2 + eta^2)-1/eta)*V[x,]^2)})
h <- (1/g1-eta)/sqrt(n-eta*sum(abs(g1)))
f <- (1/g2-eta)/sqrt(d-eta*sum(abs(g2)))
Sigmainverse = diag(d^(-1/4)*(h)^(-1/2))
Gammainverse = diag(d^(-1/4)*(f)^(-1/2))
Ytilde = Sigmainverse%*%Ybar%*%Gammainverse
svd.Y = svds(Ytilde,k=r)
Yr = (svd.Y$d[1]) * svd.Y$u[,1] %*%  t(svd.Y$v[,1])
hatY = diag((h)^(1/2))%*% Yr %*% diag((f)^(1/2))
hatv = svd(hatY)$v[,1]
hatu = hatY%*%hatv
hatu = if_else(abs(hatu) >= quantile(abs(hatu),percent_cut), abs(hatu), 0)
hatu = hatu/sqrt(sum(hatu^2))
amp.out = amp.svd(Ytilde, v0=hatv*d^(-1/4)*(f)^(-1/2), n.ite=5, percent_cut=percent_cut)
hatu_uagg = ((amp.out[[2]][,4])*(h)^(1/2))
hatv_uagg = t(Ybar) %*% hatu

#######uaggregation-o############################################
#cross validation
percent_cut=1-estimated_s
amp.out = amp.svd(Ytilde, v0=hatv*d^(-1/4)*(f)^(-1/2), n.ite=5, percent_cut=percent_cut)
hatu_uagg.cv = ((amp.out[[2]][,4])*(h)^(1/2))
hatv_uagg.cv = t(Ybar) %*% hatu

#performance
abs(c(cor(hatv_uagg,v0),cor(hatv_uagg.cv,v0),cor(v_pca,v0),cor(v_heteropca,v0), cor(vhatnaive,v0),max(apply(Ybar,1, function(x)cor(x,v0))))


