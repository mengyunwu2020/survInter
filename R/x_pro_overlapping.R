#' @export



x_pro_overlapping<-function(n,n_var=1000,M=100,Sigma){
  #Generating main genetic effects
  #@param n sample size
  #
  #@param n_var dimension
  #@param M number of pathways
  #@param Sigma covariance matrix
  p0=n_var/M
  mean <-rep(0,p0)
  temp<- mvrnorm((n*M),mean,Sigma)
  X1=NULL
  x=list()
  for(i in 1:M){
    x[[i]]=temp[(1+(i-1)*n):(n*i),]
    X1=cbind(X1,x[[i]])
  }
  X=X1

  x1=x

  temp=x1[[1]][,8:10]
  x[[2]]=cbind(x1[[2]],temp)
  x[[3]]=cbind(x1[[3]],temp)
  x[[4]]=cbind(x1[[4]],temp)
  x[[5]]=cbind(x1[[5]],temp)
  temp=x1[[6]][,8:10]
  for(m in 7:11) x[[m]]=cbind(x1[[m]],temp)
  temp=x1[[12]][,8:10]
  for(m in 13:17) x[[m]]=cbind(x1[[m]],temp)
  temp=x1[[18]][,8:10]
  for(m in 19:23) x[[m]]=cbind(x1[[m]],temp)
  temp=x1[[24]][,8:10]
  for(m in 25:29) x[[m]]=cbind(x1[[m]],temp)
  temp=x1[[30]][,8:10]
  for(m in 31:35) x[[m]]=cbind(x1[[m]],temp)
  temp=x1[[36]][,8:10]
  for(m in 37:40) x[[m]]=cbind(x1[[m]],temp)
  temp=x1[[41]][,10]
  for(m in 42) x[[m]]=cbind(x1[[m]],temp)


  X1=matrix(NA,n,1)
  for(i in 1:M){
    X1=cbind(X1,x[[i]])
  }
  X1=X1[,-1]
  return(list(x=x,X=X1,X1=X,x1=x1))
}