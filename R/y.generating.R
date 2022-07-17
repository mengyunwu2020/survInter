#' @export

#Generating survival response
y.generating<-function(n,x,X,rate,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,M,p,c){
  p0=p/M
  Mm=M+M*(M-1)/2
  Phix=list()
  for(m in 1:M){
    r=0
    p_hat=p0+p0*(p0-1)/2
    temp1=matrix(NA,n,p_hat)
    temp1[,1:p0]=x[[m]]
    for(i in 1:(p0-1)){
      for(j in (1+i):p0){
        r=r+1
        temp1[,p0+r]=(x[[m]][,i])*(x[[m]][,j])
      }
    }
    Phix[[m]]=temp1
  }


  cp=0
  R1=matrix(0,M,M)
  for(m in 2:4){
    cp=cp+1
    temp1=matrix(NA,n,(p0*p0))
    r=1
    for(i in 1:p0){
      temp1[,r:(r+p0-1)]=(x[[m]][,1:p0])*(x[[1]][,i])
      r=r+p0
    }
    Phix[[M+cp]]=temp1
  }


  b_s_int=p*(p-1)/2+p
  pp_int=p_hat*M
  p_sum_int=p
  p_int=rep(p0,M)
  # M_int=M
  # r=1
  # temp2=0
  # s=0
  # Phi2_int=matrix(0,n,(b_s_int-pp_int))
  # for(m in 1:(M_int-1)){
  #   temp=p_sum_int-temp2-p_int[m]
  #   for(i in 1:p_int[m]){
  #     Phi2_int[,((i-1)*temp+1+s):(i*temp+s)]=X[,(temp2+i)]*X[,(p_int[m]+temp2+1):p_sum_int]
  #   }
  #   temp2=temp2+p_int[m]
  #   s=temp*p_int[m]+s
  # }

  inter.group=matrix(0,p_sum_int,p_sum_int)
  r=s=ii=1
  for(j in 1:(M-1)){
    s=sum(p_int[1:j])+1
    for(i in (j+1):M){
      inter.group[s:(s+p_int[i]-1),r:(r+p_int[j]-1)]=ii
      s=s+p_int[i]
      ii=ii+1
    }
    r=r+p_int[j]
  }
  inter.group.vec=inter.group[which(inter.group!=0)]+M





  real_w=list()
  for(i in 1:M)
    real_w[[i]]=rep(0,times=p0+(p0-1)*p0/2)
  for(i in (1+M):Mm){
    real_w[[i]]=rep(0,times=p0*p0)
  }
  real_w[[1]][1:length(a1)]=a1
  real_w[[1]][(p0+1):(p0+length(b1))]=b1
  real_w[[2]][1:length(a2)]=a2
  real_w[[2]][(p0+1):(p0+length(b2))]=b2
  real_w[[3]][1:length(a3)]=a3
  real_w[[3]][(p0+1):(p0+length(b3))]=b3
  real_w[[4]][1:length(a4)]=a4
  real_w[[4]][(p0+1):(p0+length(b4))]=b4
  real_w[[5]][1:length(a5)]=a5
  real_w[[5]][(p0+1):(p0+length(b5))]=b5
  real_w[[M+1]][1:length(c[,1])]=c[,1]
  real_w[[M+2]][1:length(c[,2])]=c[,2]
  real_w[[M+3]][1:length(c[,3])]=c[,3]



  real_y=0
  for(m in c(1:5,(M+1),(M+2),(M+3))){
    real_y=Phix[[m]]%*%real_w[[m]]+real_y
  }



  t = real_y + rnorm(n)

  T=exp(t)
  delta=rep(0,n)

  niter=1
  rrrr=NULL
  if(rate==0.8){
    while((abs(sum(delta)/n-rate)>0.02)&niter<=50){
      niter=niter+1
      censoring.time<-rgamma(round(n*0.6), shape =10, scale =10)  # censoring time
      cc<-sample(c(censoring.time,rep(exp(9999),n-round(n*0.6))),n)
      time<- pmin(T, cc)  # observed time is min of censored and true
      delta = time == T   # set to 1 if event is observed

    }
    if(niter==51){

      ss=c(0.005,seq(0.05,0.1,0.01),seq(0.1,1,0.1),seq(1,1000,1))
      for(i in ss){
        censoring.time<-rgamma(round(n*0.6), shape =i, scale =i)  # censoring time
        cc<-sample(c(censoring.time,rep(exp(9999),n-round(n*0.6))),n)
        time<- pmin(T, cc)  # observed time is min of censored and true
        delta = time == T   # set to 1 if event is observed
        #print(sum(delta)/n)
        rrrr=append(rrrr,sum(delta)/n)
        if(abs(sum(delta)/n-rate)<=0.02){
          break;
        }
      }
    }
    if(i==1000){
      temp=which.min(abs(rate-rrrr))
      censoring.time<-rgamma(n, shape =ss[temp], scale =ss[temp])  # censoring time
      time<- pmin(T, censoring.time)  # observed time is min of censored and true
      delta = time == T   # set to 1 if event is observed
    }

  }

  ct=cbind(time,delta)

  return(list(t=t,ct=ct,real_w=real_w,observed.rate=sum(delta)/n))

}