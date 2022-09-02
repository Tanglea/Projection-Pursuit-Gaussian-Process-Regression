#This is the code of PPGPR in the paper "Projection Pursuit Gaussian Process Regression".
#This code is for the experiment for Table 4.
.
# learning_rate = e-9, nu=2.5, train_epoch= 150, number of nodes =28, and Matern kernel with nu=2.5.
library(fields)
library(MASS)
library(SDraw)
library(RandomFieldsUtils)
library(DiceKriging)


set.seed(3)


kappa=2.5


ts0=14 #testing sample size



materncorr<-function(d,nu)
  #matern(x,nu)
  2/gamma(nu)*(sqrt(nu)*d)^nu*besselK(2 * sqrt(nu)*d,nu)

diff<-function(x,nu)
{
  if(x==0)res=0
  else 
    res=(2*x)*nu/(nu-1)*matern(abs(x)*sqrt(nu/(nu-1)),nu-1)
  res
}

Cov<-function(X1)
{
  n1=ncol(X1)
  n2=nrow(X1)
  Cov=matrix(0,n1,n1)
  for(i in 1:n1)
    for(j in 1:n1)
      for(k in 1:n2)
        Cov[i,j]=Cov[i,j]+matern(abs(X1[k,i]-X1[k,j]),nu=kappa)
  
  Cov/n2
}



GD<-function(X,W,Y,sl)
{
  X1=W%*%X
  
  C=Cov(X1)
  
  
  u=solve(C+diag(nrow(C))*1E-8,Y)
  
  n1=ncol(X1)
  n2=nrow(X1)
  
  d=ncol(W)
  
  grad=matrix(0,n2,d)
  
  div=rep(0,n2)
  
  
  
  for(k in 1:n2)
  {
    for(i in 1:n1)
      for(j in 1:n1)
        
        grad[k,]=grad[k,]+u[i]*u[j]*diff(X1[k,i]-X1[k,j],nu=kappa)*t(X[,i]-X[,j])/n2
      
      
      
      grad[k,] = grad[k,] + sum(diag(solve(C)))
      
      #	grad[k,]=grad[k,]-(t(grad[k,])%*%W[k,])[1]*W[k,]
      
      div[k]=sqrt((t(grad[k,])%*%grad[k,])[1])
      
      
      
      
      
      
      
  }
  
  
  
  dmax=max(div)
  
  for(k in 1:n2)
  {
    W[k,]=W[k,]-0.000000001 * grad[k,]
    #		W[k,]=W[k,]/sqrt(t(W[k,])%*%W[k,])[1]		
  }
  
  res<-list(W=W,norm=(t(Y)%*%u)[1],grad=mean(div))
  res
}

##data




## data

ss=64

dimX=4

data_gene<-function()
{
  
  
  #X=matrix(0,dimX,ss)
  #X = halton(ss,dimX)
  train_data = as.matrix(HeatExchangerTrain[,2:6])
  X = train_data[,1:4]
  X=t(X)
  #X=t(X)
  #X=2*X-1
  
  
  
  
  Y = train_data[,5]
  res<-list(X=X,Y=Y)
  
  res$X=X
  res$Y=Y
  
  res
}

##Architecture

initW<-function(X,d)
{
  n=nrow(X)
  W=matrix(0,d,n)
  for(i in 1:d)
  {
    W1=rnorm(n)
    W1=W1/sqrt(t(W1)%*%W1)[1]
    W[i,]=t(W1)
  }
  W
}

training<-function(dimFt,X,Y)
{
  
  ##dimFt=14
  
  W=initW(X,dimFt)
  
  ##training
  
  j=1
  
  prenorm=0
  
  for(i in 1:100)
  {
    a=GD(X,W,Y,1/j)
    if(a$norm>prenorm)
    {
      j=j*2
    }
    prenorm=a$norm
    W=a$W
    
    if(j>1000)break
  }
  
  res=W
  
  res
  
}

##prediction

prediction<-function(X,Y,W,Xtest)
{
  
  X1=W%*%X
  
  res<-list(E1=0,E2=0)
  
  
  kriging<-function(X1,Y)
  {
    C=Cov(X1)
    
    u0=solve(C+diag(nrow(C))*1E-6,Y)
    
    n1=ncol(X1)
    n2=nrow(X1)
    
    #=matrix(0,n2,n1)
    #for(i in 1:n1)
    #{
    #  for(j in 1:n2)
    #    for(k in 1:n1)
    #      resY[j,i]=resY[j,i]+Matern(d=abs(X1[j,k]-X1[j,i]),nu=kappa)*u0[k]
    #}  
    #res<-list(u=u0,Y=resY,norm=t(Y)%*%u0)
    res<-list(u=u0,norm=t(Y)%*%u0)
    res
  }
  
  
  kriging_res=kriging(X1,Y)
  
  
  
  f_pred<-function(xnew, X1, u, W)
  {
    res=0
    xnew=W%*%xnew
    n1=ncol(X1)
    n2=nrow(X1)
    for(i in 1:n1)
    {
      for(j in 1:n2)
        res=res+matern(abs(X1[j,i]-xnew[j]),nu=kappa)*u[i]/n2
    }
    res
  }
  
  
  
  
  Z=rep(0,ts0)
  ZT=Z
  
  test_data = as.matrix(HeatExchangerTest[,2:6])
  ZT = test_data[,5]
  
  
  for (i in 1:ts0)
  {
    Z[i]=f_pred(X_test[,i],X1,kriging_res$u,W)
  }
  
  
  #print(sum(abs(ZT-Z)/abs(ZT))/14)
  print(sqrt(mean(((ZT-Z)/ZT)^2)))
  #res$E1=sum((ZT-Z)*(ZT-Z))
  
  
  
}

Er=matrix(0,7,10)
Erold=rep(0,10)

for(p in 1:1)
{
  
  data=data_gene()
  
  X=data$X
  Y=data$Y
  test_data = as.matrix(HeatExchangerTest[,2:6])
  #X_test=matrix(runif(ts0*dimX),dimX,ts0)
  #X_test = 2*X_test -1
  X_test = test_data[,1:4]
  X_test = t(X_test)
  
  
  for(q in 1:1)
  {
    
    t1 = Sys.time()
    W=training(28,X,Y)
    
    
    
    E=prediction(X,Y,W,X_test)
    t2 = Sys.time()
    
    print(t2-t1)
    #Er[j,i]=E$E1
    
  }
  t3 = Sys.time()
  X = t(X)
  X_test = t(X_test)
  dice_model = km(design = data.frame(X),response = data.frame(Y))
  Y_pred_dice = predict(dice_model,newdata=data.frame(X_test),type="UK")
  Y_true = test_data[,5]
  #print(sum(abs(Y_pred_dice$mean-Y_true)/abs(Y_true))/14)
  print(sqrt(mean(((Y_pred_dice$mean-Y_true)/Y_true)^2)))
  t4 = Sys.time()
  print(t4-t3)
  #Erold[i]=E$E2
}
