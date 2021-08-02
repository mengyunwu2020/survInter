#' Two-level Bayesian interaction analysis for survival data incorporating pathway information
#'
#' One of the main functions in the survInter package. Fits a path of survInter models over different values of the tunning parameters.
#'

#' @importFrom stats var
#' @importFrom stats rnorm
#' @importFrom stats rgamma
#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' @importFrom MASS mvrnorm
#' @export
#'
#' @param X a matrix of predictor variables
#' @param ct a two-column matrix with the first column being the survival time and the second column being the censoring indicator.
#' The indicator is a binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param pathway a vector of pathway information, for instance,
#' value 1 represent the corresponding node location belongs to pathway 1
#' @param s_init a vector of tunning parameter related to the variance
#' @param a_init a vector of the initial selecting probability
#' @param r_init a vector of tunning parameter related to the variance
#' @param tol convergence tolerance
#'
#' @return a list, each element corresponding to its tunning parameters contains values  as follows:
#' \item{mean_w_red}{The coefficients vector}
#' \item{beta_interaction}{The coefficients of interaction terms}
#' \item{beta_main}{The coefficients of main effects}
#' \item{BIC}{Bayesian Information Criterion}
#' \item{df}{The number of nonzero coefficients}
#' \item{index_path}{Selected networks' number}
#' @examples
#' set.seed(11)
#' u=runif(15,0.8,1.2)
#' u=round(u,2)
#'
#' rho=0.6
#' rate=0.8
#' p=n_var=1000
#' n=400
#' a1=rep(u[1],5)
#' a2=rep(u[2],5)
#' a3=rep(u[3],5)
#' a4=rep(u[4],5)
#' a5=rep(0,5)
#' #interaction within networks
#' b1=rep(u[5],4)
#' b2=rep(u[6],4)
#' b3=rep(u[7],4)
#' b4=rep(u[8],4)
#' b5=rep(0,4)
#' #interaction among networks
#' c=matrix(NA,nrow =4 ,ncol=3)
#' c[,1]=rep(u[9],4)
#' c[,2]=rep(u[10],4)
#' c[,3]=rep(0,4)
#' M=100
#' p0=p/M
#' Sigma=matrix(NA,p0,p0)
#' for(i in 1:p0){
#' for(j in 1:p0)  Sigma[i,j]=rho^(abs(i-j))
#' }
#'
#' temp=x_pro_overlapping(n,n_var=1000,M=M,Sigma)
#' x=temp$x
#' X=temp$X
#' x1=temp$x1
#' X1=temp$X1
#' output<-y.generating(n,x1,X1,rate,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,M,p,c)
#' ct<-as.matrix(output$ct)
#' real_w=output$real_w
#' cat('rate=',output$observed.rate)
#' rm(output)
#'
#' p_temp=as.numeric(lapply(x,function(x) dim(x)[2]))
#' pathway=rep(1:M,p_temp)
#' a_init=seq(0.96,0.97,by=0.01)
#' s_init=seq(2e-3,3e-3,1e-3)
#' temp_1<-survInter(X,ct,pathway,s_init=s_init,a_init=a_init)


survInter<-function(X,ct,pathway,s_init=1e-4,a_init,r_init=NULL,tol=1e-5){

  t=ct[,1]
  delta=ct[,2]
  n=dim(X)[1]
  n_var=dim(X)[2]
  M0=max(pathway)
  p_total=matrix(0,M0,1)
  for(i in 1:M0){
    p_total[i]=sum(pathway==i)
  }
  t_int=log(t)
  censoring.time=log(t)

  p_int=p_total
  q_int=p_int*(p_int-1)/2
  M_int=M0
  p_sum_int=sum(p_int)              #node
  p_hat_int=p_int*(p_int-1)/2+p_int
  pp_int=sum(p_hat_int[1:M0])      #pathway
  b_s_int=p_sum_int+p_sum_int*(p_sum_int-1)/2


  R_int=lapply(p_int,location)

  node_location_int=matrix(1,p_sum_int,1)
  node_location_int[1:p_int[1]]=1:p_int[1]
  for (m in 2:M_int){
    node_location_int[(sum(p_int[1:(m-1)])+1):(sum(p_int[1:(m-1)])+p_int[m])]<-(sum(p_hat_int[1:(m-1)])+1):(sum(p_hat_int[1:(m-1)])+p_int[m])
  }
  within_interaction_location_int=setdiff(1:pp_int,node_location_int)
  within_interaction_location2_int=matrix(FALSE,p_sum_int,p_sum_int)
  within_interaction_location2_int[1:p_int[1],1:p_int[1]]=TRUE
  for (m in 2:M_int){
    within_interaction_location2_int[(sum(p_int[1:(m-1)])+1):(sum(p_int[1:(m)])),(sum(p_int[1:(m-1)])+1):(sum(p_int[1:(m)]))]=TRUE
  }


  out_interaction_location2_int=!(within_interaction_location2_int)
  RR_int=matrix(0,p_sum_int,p_sum_int)
  temp=(lower.tri(out_interaction_location2_int))*out_interaction_location2_int
  RR_int[which(temp==1)]=c(1:(b_s_int-pp_int))

  Phi_int=matrix(0,n,pp_int)
  Phi_int[,node_location_int]=X

  s=0
  r1=1
  for(m in 1:M_int){
    if(q_int[m]>0){
      temp=matrix(0,n,q_int[m])
      r=0
      for(i in 1:(p_int[m]-1)){
        temp1=p_int[m]-i
        temp[,(r+1):(r+temp1)]=(X[,(s+i)])*(X[,(s+i+1):(p_int[m]+s)])
        r=r+temp1
      }
      Phi_int[,(r1+p_int[m]):(r1+p_hat_int[m]-1)]=temp
      s=s+p_int[m]
      r1=r1+p_hat_int[m]
    }
  }


  r=1
  temp2=0
  s=0
  Phi2_int=matrix(0,n,(b_s_int-pp_int))
  for(m in 1:(M_int-1)){
    temp=p_sum_int-temp2-p_int[m]
    for(i in 1:p_int[m]){
      Phi2_int[,((i-1)*temp+1+s):(i*temp+s)]=X[,(temp2+i)]*X[,(p_int[m]+temp2+1):p_sum_int]
    }
    temp2=temp2+p_int[m]
    s=temp*p_int[m]+s
  }

  inter.group=matrix(0,p_sum_int,p_sum_int)
  r=s=ii=1
  for(j in 1:(M0-1)){
    s=sum(p_total[1:j])+1
    for(i in (j+1):M0){
      inter.group[s:(s+p_total[i]-1),r:(r+p_total[j]-1)]=ii
      s=s+p_total[i]
      ii=ii+1
    }
    r=r+p_total[j]
  }
  inter.group.vec=inter.group[which(inter.group!=0)]+M0


  cols_PP1_int=colSums(Phi_int*Phi_int)
  cols_PP2_int=colSums((Phi2_int*Phi2_int))
  pathway_new_int=matrix(0,pp_int,1)
  pathway_new_int[node_location_int]=pathway


  ls=length(s_init)
  la=length(a_init)
  result_temp=list()
  lll=1
  temp=rep(0,times=(M0-1)*M0/2)
  s=1
  for(i in 1:(M0-1)){
    temp[s:(s+M0-i-1)]=p_int[i]*p_int[(i+1):M0]
    s=s+M0-i
  }

  p_hat_int=append(p_hat_int,temp)
  censored.id=which(delta==0)
  n.observed=n-length(censored.id)


  for(iii in 1:ls){
    ss=s_init[iii]
    if(is.null(r_init)){
      s4=ss*0.8
    }else s4=r_init

    for(jjj in 1:la){
      t=t_int
      a=a_init[jjj]
      p=p_int
      M=M_int
      p_sum=p_sum_int
      p_hat=p_hat_int
      pp=pp_int
      b_s=b_s_int
      R=R_int
      node_location=node_location_int
      within_interaction_location=within_interaction_location_int
      within_interaction_location2=within_interaction_location2_int
      out_interaction_location2=out_interaction_location2_int
      RR=RR_int
      Phi=Phi_int
      q=q_int
      Phi2=Phi2_int
      cols_PP1=cols_PP1_int
      cols_PP2= cols_PP2_int

      Mm=M+(M-1)*M/2
      gamma=a*rep(1,times=Mm)
      eta_ori=matrix(a,b_s,1)

      mean_w_ori=rep(0,times=b_s)
      tau =as.numeric(sqrt(1/var(t)))
      r1=s1=r3=1
      r2=r4=s2=ss
      zeta1=zeta2=zeta3=zeta4=a

      niter = 0

      diff_sum=1
      index_path_ori=c(1:M)

      sigma_vec=sigma_vec_ori=matrix(0,b_s,1)
      mean_w=mean_w_ori
      pathway_inter_includ=matrix(0,pp,1)
      s=1
      for (i in 1:M){
        pathway_inter_includ[s:(s+p_hat[i]-1),1]=i
        s=s+p_hat[i]
      }
      locat_within=c(1:pp)
      eta=eta_ori
      uz=t
      w_mul_Phi=rep(0,n)
      while (diff_sum>tol&niter<1000){
        niter = niter+1

        for(i in censored.id){
          temp=w_mul_Phi[i]
          uz[i]=temp+dnorm(censoring.time[i],mean=temp,sd=sqrt(1/tau))/(tau*(1-pnorm(censoring.time[i],mean=temp,sd=sqrt(1/tau))))
        }
        t[censored.id]=uz[censored.id]


        gamma_new=rep(0,pp)
        gamma_new[1:p_hat[1]]=gamma[1]
        if(M>1){
          for (m in 2:M){
            gamma_new[(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:m]))]=gamma[m]
          }
        }
        eta_temp=matrix(eta[node_location],p_sum,1)
        temp=eta_temp%*%t(eta_temp)
        temp2=as.vector(t(upper.tri(temp)*within_interaction_location2))
        eta_prod=(as.vector(temp))[temp2==1]
        A=1/s1*gamma_new+1/s2*(1-gamma_new)+eta[1:pp]/s1+(1-eta[1:pp])/s2
        temp=rep(0,pp)
        temp[within_interaction_location]=eta_prod/s1+(1-eta_prod)/s2
        A=A+temp
        rm(gamma_new)
        sigma_vec[1:pp]=1/(tau*cols_PP1[1:pp]+A)
        #----------------------------------------------------------------------
        mean_w_old=mean_w_ori

        if(M>1){
          y1=t-Phi%*%mean_w[1:pp]-Phi2%*%mean_w[(pp+1):b_s]
        }else{
          y1=t-Phi%*%mean_w[1:pp]
        }
        for(i in 1:pp){
          y1=y1+Phi[,i]*mean_w[i]
          mean_w[i]=sigma_vec[i]*(sum(Phi[,i]*y1)*tau)
          y1=y1-Phi[,i]*mean_w[i]
        }



        rm(A)
        mean_w_ori[1:pp_int]=sigma_vec_ori[1:pp_int]=rep(0,pp_int)
        sigma_vec_ori[locat_within]=sigma_vec[1:pp]
        mean_w_ori[locat_within]=mean_w[1:pp]


        #update Q(gamma)
        s=1
        for(m in 1:M){
          temp=log((1-zeta1)/zeta1)+0.5*(p_hat[m])*log(s1/s2)+0.5*(sum(mean_w[s:(s+p_hat[m]-1)]^2+sigma_vec[s:(s+p_hat[m]-1)]))*(1/s1-1/s2) ##----different from the first one
          gamma[m]= 1/(1+exp(temp))
          s=p_hat[m]+s
        }

        #seleted pathway
        gamma_trunc=ifelse(gamma[1:M]>0.5,1,0)
        index_path1=which(gamma[1:M]>0.5)
        if(length(index_path1)==0) break;
        index_path_ori=gamma_trunc*(index_path_ori[index_path_ori!=0])

        gamma_path=gamma[index_path1]
        if(M>1){
          temp=location(M)
          temp=temp[index_path1,index_path1]
          bi_index=temp[temp!=0]
          de=gamma[(bi_index+M)]
          gamma=c(gamma_path,de)
        }else gamma=gamma_path



        p=p_total[index_path_ori]
        M=length(p)
        p_sum=sum(p)
        q=p*(p-1)/2
        p_hat=p*(p-1)/2+p


        Mm=M+(M-1)*M/2
        if(M>1){
          temp=rep(0,times=(M-1)*M/2)
          s=1
          for(i in 1:(M-1)){
            temp[s:(s+M-i-1)]=p[i]*p[(i+1):M]
            s=s+M-i
          }
          p_hat=append(p_hat,temp)
        }


        pp=sum(p_hat[1:M])      #pathway
        b_s=p_sum+p_sum*(p_sum-1)/2  #total
        pathway_second=matrix(0,p_sum,1)
        s=1
        for (i in 1:M){
          pathway_second[s:(s+p[i]-1),1]=i
          s=s+p[i]
        }

        locat_within=which(is.element(pathway_inter_includ,index_path_ori))#######this only find the location of node and inter within pathway
        temp=which(is.element(pathway,index_path_ori))###the node rank only in 1000 node for the location of inter among pathway
        temp1=RR_int[temp,temp]





        locat_inter_out=as.vector(temp1[temp1!=0])+pp_int
        mean_w=mean_w_ori[c(locat_within,locat_inter_out)]
        sigma_vec=sigma_vec_ori[c(locat_within,locat_inter_out)]
        eta=eta_ori[c(locat_within,locat_inter_out)]
        R=lapply(p,location)#here p is a vector


        if(M>1){
          inter.group=matrix(0,p_sum,p_sum)
          r=s=ii=1
          for(j in 1:(M-1)){
            s=sum(p[1:j])+1
            for(i in (j+1):M){
              inter.group[s:(s+p[i]-1),r:(r+p[j]-1)]=ii
              s=s+p[i]
              ii=ii+1
            }
            r=r+p[j]
          }
          inter.group.vec=inter.group[which(inter.group!=0)]+M
        }else inter.group.vec=NULL

        node_location=matrix(1,p_sum,1)
        node_location[1:p[1]]=1:p[1]
        if(M>1){
          for (m in 2:M){
            node_location[(sum(p[1:(m-1)])+1):(sum(p[1:(m-1)])+p[m])]<-(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:(m-1)])+p[m])
          }
        }
        within_interaction_location=setdiff(1:pp,node_location)
        within_interaction_location2=matrix(FALSE,p_sum,p_sum)
        within_interaction_location2[1:p[1],1:p[1]]=TRUE
        if(M>1){
          for (m in 2:M){
            within_interaction_location2[(sum(p[1:(m-1)])+1):(sum(p[1:(m)])),(sum(p[1:(m-1)])+1):(sum(p[1:(m)]))]=TRUE
          }
        }
        out_interaction_location2=!(within_interaction_location2)
        RR=matrix(0,p_sum,p_sum)
        temp=(lower.tri(out_interaction_location2))*out_interaction_location2
        RR[which(temp==1)]=c(1:(b_s-pp))


        Phi=matrix(Phi_int[,locat_within],n,pp)

        Phi2=matrix(0,n,(b_s-pp))
        if(M>1){
          Phi2=as.matrix(Phi2_int[,locat_inter_out-pp_int])
        }


        cols_PP1=colSums(Phi*Phi)
        cols_PP2=colSums((Phi2*Phi2))

        pathway_new=matrix(0,pp,1)
        pathway_new[node_location]=pathway_second



        if(M>1){

          B=gamma[inter.group.vec]/s1+(1-gamma[inter.group.vec])/s4
        }


        if(M>1){
          eta_temp=matrix(eta[node_location],p_sum,1)
          temp=eta_temp%*%t(eta_temp)
          temp2=as.vector(t(upper.tri(temp)*within_interaction_location2))
          eta_prod=(as.vector(temp))[temp2==1]
          temp2=as.vector(t(upper.tri(temp)*out_interaction_location2))
          eta_prod_out=(as.vector(temp))[temp2==1]
          summ=(eta[(pp+1):b_s]/r1+eta_prod_out/r3)+(1-eta[(pp+1):b_s])/r2+(1-eta_prod_out)/r4+B
          sigma_vec[(pp+1):b_s]=as.vector(1/(tau*cols_PP2+summ))
          y1=t-Phi%*%mean_w[1:pp]-Phi2%*%mean_w[(pp+1):b_s]
          for(i in (1+pp):b_s){
            y1=y1+Phi2[,(i-pp)]*mean_w[i]
            mean_w[i]=(sum(Phi2[,(i-pp)]*y1))*tau*sigma_vec[i]
            y1=y1-Phi2[,(i-pp)]*mean_w[i]
          }
        }else{
          y1=t-Phi%*%mean_w[1:pp]
        }


        mean_w_ori[(pp_int+1):b_s_int]=sigma_vec_ori[(pp_int+1):b_s_int]=rep(0,(b_s_int-pp_int))
        if(M>1){
          mean_w_ori[locat_inter_out]=mean_w[(pp+1):b_s]
          sigma_vec_ori[locat_inter_out]=sigma_vec[(pp+1):b_s]
          wwmtp_path_out=sigma_vec[(pp+1):b_s]+mean_w[(pp+1):b_s]^2
        }else{
          wwmtp_path_out=0
        }



        if(M>1){
          for(m in (1+M):(Mm))
          {
            index=which(inter.group.vec==m)
            ttemp=sum(wwmtp_path_out[index])
            temp=log((1-zeta4)/zeta4)+0.5*p_hat[(m)]*log(s1/s4)+0.5*ttemp*(1/s1-1/s4)  ############
            gamma[m]=1/(1+exp(temp))
          }
        }

        wwmtp_path_within=sigma_vec[1:pp]+mean_w[1:pp]^2
        for(i in node_location){
          m=pathway_new[i]
          id_node=which(node_location==i)
          if((p[m]!=1)){
            id_temp=as.vector(setdiff(which(pathway_new==pathway_new[i]),i))#different node in its pathway,location in the whole pathway including inter
            id_nodep=ifelse(m>1,i-sum(p_hat[1:(m-1)]),i)  #ith node location in its pathway valued 1:p[m] only, for the sum_kj calculation

            id_temp1=ifelse(rep(m,(p[m]-1))>1,id_temp-sum(p_hat[1:(m-1)]),id_temp)#id_temp corresponding nodes location in its pathway
            temp=pmax((R[[m]][id_nodep,id_temp1]),(R[[m]][id_temp1,id_nodep]))#### interaction within its pathway correlated with node i
            temp=ifelse(rep(m,(p[m]-1))>1,temp+sum(p_hat[1:(m-1)]),temp)+p[m]## interaction location in wwmtp_path_within
            sum_kj=sum(eta[id_temp]*(0.5*log(r3/r4)+0.5*(wwmtp_path_within[temp])*(1/r3-1/r4)))#0.5*diag(wwmtp[[m]][(temp+p[m]),(temp+p[m])])*(1/s1-1/s2)))
          }else{
            sum_kj=0
          }
          id_temp=as.vector(setdiff(1:p_sum,which(pathway_second==m)))#different node all pathway
          id_temp1=as.vector(which((pathway_new!=m)&(pathway_new!=0)))
          temp=pmax(RR[id_node,id_temp],RR[id_temp,id_node])  ##interation location all pathway
          sum_mkj=sum(eta[id_temp1]*(0.5*log(r3/r4)+0.5*wwmtp_path_out[temp]*(1/r3-1/r4)))
          eta[i]=1/(1 +exp(0.5*log(r1/r2)+ 0.5*wwmtp_path_within[i]*(1/r1-1/r2)+log((1-zeta2)/zeta2)+sum_kj+sum_mkj))
        }
        eta[within_interaction_location]=1/(1 +exp(0.5*log(r1/r2)+ 0.5*(wwmtp_path_within[within_interaction_location])*(1/r1-1/r2)+ log((1-zeta3)/zeta3)))

         if(M>1) eta[(pp+1):b_s]=1/(1 +exp(0.5*log(r1/r2)+ 0.5*wwmtp_path_out[1:(b_s-pp)]*(1/r1-1/r2)+log((1-zeta3)/zeta3))) ################
        rm(wwmtp_path_within)


        eta_ori=rep(0,b_s_int)
        temp=which(is.element(pathway,index_path_ori))
        eta_ori[locat_within]=eta[1:pp]
        if(M>1){
          eta_ori[locat_inter_out]=eta[(pp+1):b_s]
        }

        #M-step
        w_mul_Phi=t-y1
        trace0=sum(w_mul_Phi^2)+sum(cols_PP1_int*sigma_vec_ori[1:pp_int])#-

         if(M>1) trace0=trace0+sum(cols_PP2_int*sigma_vec_ori[(pp_int+1):b_s_int])#trace0=trace0+sum(cols_PP1*sigma_vec_ori[locat_within])+sum(cols_PP2*sigma_vec_ori[locat_inter_out])

        tau=n.observed/(sum((t^2-2*t*w_mul_Phi)[which(delta==1)])+trace0-sum((w_mul_Phi^2)[censored.id]))
        temp=t(location(M0))
        temp=temp[index_path_ori,index_path_ori]
        bi_index=temp[temp!=0]

        gamma_ori=rep(0,M0+M0*(M0-1)/2)
        gamma_ori[index_path_ori]=gamma[1:M]
        gamma_ori[bi_index+M0]=gamma[(M+1):Mm]

        zeta1=zeta4=mean(gamma_ori)
        zeta2=mean(eta_ori[node_location_int])
        inter_location=setdiff(1:b_s_int,node_location_int)
        zeta3=mean(eta_ori[inter_location])

        temp=which((mean_w_old-mean_w_ori)!=0)
        diff_sum = mean(c(abs((mean_w_old-mean_w_ori)[temp]/mean_w_old[temp])))
      }


      if(length(index_path1)==0) {
        eta_ori=mean_w_ori=rep(0,b_s_int)
        eta=mean_w=rep(0,length(eta))
        gamma_ori=rep(0,(M0+M0*(M0-1)/2))
      }
      index_pred1=which(eta>0.5)
      temp=rep(0,b_s)
      temp[index_pred1]=mean_w[index_pred1]

      if(M>1){
         temp=((t_int-Phi%*%temp[1:pp]-Phi2%*%temp[(pp+1):b_s])^2)[delta==1]
        error=mean(temp)
      }else{
        temp=((t_int-Phi%*%temp[1:pp])^2)[delta==1]
        error=mean(temp)
      }



      index_pred=which(eta_ori>0.5)
      mean_w_red=rep(0,b_s_int)
      mean_w_red[index_pred]=mean_w_ori[index_pred]

      df=length(which(eta>0.5))
      BIC=log(error)+df*log(n.observed)/n.observed
      beta_interaction=c(mean_w_red[c(within_interaction_location_int,c((pp_int+1):b_s_int))])#
      beta_main=mean_w_red[node_location_int]
      temp_2=list(mean_w_red=mean_w_red,beta_interaction=beta_interaction,beta_main=beta_main,BIC=BIC,df=df,index_path=index_path_ori)
      result_temp[[lll]]<-temp_2
      lll=lll+1




    }

  }


  return(result_temp)

}
