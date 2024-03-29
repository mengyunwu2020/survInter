% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survInter.R
\name{survInter}
\alias{survInter}
\title{Two-level Bayesian interaction analysis for survival data incorporating pathway information}
\usage{
survInter(X, ct, pathway, s_init = 1e-04, a_init, r_init = NULL, tol = 0.001)
}
\arguments{
\item{X}{a matrix of predictor variables}

\item{ct}{a two-column matrix with the first column being the survival time and the second column being the censoring indicator.
The indicator is a binary variable, with "1" indicating dead, and "0" indicating right censored}

\item{pathway}{a vector of pathway information, for instance,
value 1 represent the corresponding node location belongs to pathway 1}

\item{s_init}{a vector of tuning parameter related to the variance}

\item{a_init}{a vector of the initial selecting probability}

\item{r_init}{a vector of tuning parameter related to the variance}

\item{tol}{convergence tolerance}
}
\value{
a list, each element corresponding to its tuning parameters contains the following values:
\item{mean_w_red}{The coefficients vector}
\item{beta_interaction}{The coefficients of interaction terms}
\item{beta_main}{The coefficients of main effects}
\item{df}{The number of nonzero coefficients}
\item{index_path}{Selected pathways' number}
}
\description{
One of the main functions in the survInter package. Fits a path of survInter models over different values of the tuning parameters.
}
\examples{
set.seed(11)
u=runif(15,0.8,1.2)
u=round(u,2)

rho=0.6
rate=0.8
p=n_var=1000
n=400
a1=rep(u[1],5)
a2=rep(u[2],5)
a3=rep(u[3],5)
a4=rep(u[4],5)
a5=rep(0,5)
#interaction within networks
b1=rep(u[5],4)
b2=rep(u[6],4)
b3=rep(u[7],4)
b4=rep(u[8],4)
b5=rep(0,4)
#interaction among networks
c=matrix(NA,nrow =4 ,ncol=3)
c[,1]=rep(u[9],4)
c[,2]=rep(u[10],4)
c[,3]=rep(0,4)
M=100
p0=p/M
Sigma=matrix(NA,p0,p0)
for(i in 1:p0){
for(j in 1:p0)  Sigma[i,j]=rho^(abs(i-j))
}

temp=x_pro_overlapping(n,n_var=1000,M=M,Sigma)
x=temp$x
X=temp$X
x1=temp$x1
X1=temp$X1
output<-y.generating(n,x1,X1,rate,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,M,p,c)
ct<-as.matrix(output$ct)
real_w=output$real_w
cat('rate=',output$observed.rate)
rm(output)

p_temp=as.numeric(lapply(x,function(x) dim(x)[2]))
pathway=rep(1:M,p_temp)
a_init=seq(0.96,0.97,by=0.01)
s_init=seq(2e-3,3e-3,1e-3)
temp_1<-survInter(X,ct,pathway,s_init=s_init,a_init=a_init,tol=1e-5)
}
