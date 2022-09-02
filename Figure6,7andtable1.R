#This is the code of PPGPR in the paper "Projection Pursuit Gaussian Process Regression".
#This code is for the experiment for Figure 6,7 and Table 1.
.
# learning_rate = e-9, nu=2.5, train_epoch= 150, number of nodes =35, Matern kernel with nu=1.5, 2.5, 3.5, 4 or Gaussian kernel with phi=0.4, 0.5, 0.6, 0.7.
library(fields)
library(MASS)
library(SDraw)
library(RandomFieldsUtils)
library(DiceKriging)

rm(list=ls())
set.seed(3)


kappa=2.5


ts0=500 #testing sample size



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


detpep10curv <- function(xx)
{

  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  rw=rw*0.05+0.1
  r=r*20000+30000
  Tu=100000+Tu*15000
  Hu=1000+10*Hu
  Tl=90+20*Tl
  Hl=750+50*Hl
  L=1400+200*L
  Kw=1000+5000*Kw
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

## data

ss=40

dimX=8

data_gene<-function()
{
  
  
  X=matrix(0,dimX,ss)
  X = halton(ss,dimX)
  X=t(X)
  #X=t(X)
  X=2*X-1
  
  
  
  for (i in 1:ss)
  {
    #X[,i]=2*runif(dimX)-1
  }
  
  Y=rep(0,ss)
  for(i in 1:ss)
  {
    ##Y[i]=f(X[1,i],X[2,i])
    Y[i]=detpep10curv(X[,i])
  }
  
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
  
  for(i in 1:150)
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
  

  for (i in 1:ts0)
  {
    Z[i]=f_pred(X_test[,i],X1,kriging_res$u,W)
    ZT[i]=detpep10curv(X_test[,i])
  }


  print(sum(abs(ZT-Z)/abs(ZT))/500)
  #print(sqrt(mean(((ZT-Z)/ZT)^2)))
  #res$E1=sum((ZT-Z)*(ZT-Z))
  
  
  
}

Er=matrix(0,7,10)
Erold=rep(0,10)

for(p in 1:1)
{
  
  data=data_gene()
  
  X=data$X
  Y=data$Y
  
  X_test=matrix(runif(ts0*dimX),dimX,ts0)
  #X_test = 2*X_test -1
  
  for(q in 1:1)
  {
    
    t1 = Sys.time()
    W=training(35,X,Y)

    

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
  Y_true = rep(0,500)
  for (i in 1:500)
    {
          
         Y_true[i]=detpep10curv(X_test[i,])
  }
  print(sum(abs(Y_pred_dice$mean-Y_true)/abs(Y_true))/500)
  t4 = Sys.time()
  print(t4-t3)
  #Erold[i]=E$E2
}
