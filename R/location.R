location<-function(p){
  if(p>1){
    temp=matrix(0,p,p)
    r=1
    for(i in 1:(p-1)){
      temp[i,(i+1):p]=c(r:(r+p-i-1))
      r=r+p-i
    }
  }else{
    temp=0
  }
  return(temp)
}