
std<-function(x){
	x<-x[!is.na(x)]
	(x-mean(x))/sd(x)
}




#fold change function
#positive: class1/class0
#negative: class0/class1
FC <- function(x,cl,logTransformed=TRUE,log.base=2,average.method=c('geometric','arithmetic'),pseudoCount=0){
	x.class0 <- x[(cl == 0)]+pseudoCount
	x.class1 <- x[(cl == 1)]+pseudoCount
	if(missing(average.method))
		average.method<-'geometric'
	if(logTransformed){
		if(is.na(log.base)|log.base<0)
			stop('You must specify log.bsae !\n')
		logFC<-mean(x.class1)-mean(x.class0)
		FC.val<-sign(logFC)*log.base^abs(logFC)
	}else{
		logFC<-ifelse(average.method=='arithmetic',log(mean(x.class1))-log(mean(x.class0)),mean(log(x.class1)-mean(log(x.class0))))
		FC.val<-sign(logFC)*exp(abs(logFC))
	}
	FC.val[FC.val==0 | is.na(FC.val)]<-1
	FC.val
}



es <- function(z,es.method="mean"){
  if(es.method=="maxmean"){
    n<-length(z)
    m1<-ifelse(sum(z>0)>0,sum(z[z>0])/n,0)
    m2<-ifelse(sum(z<0)>0,sum(z[z<0])/n,0)
    if(m1>-m2) es<-m1
    else es<-m2
  }
  else if(es.method=='absmean'){
    es<-mean(abs(z))
  }
  else if(es.method=='mean'){
    es<-mean(z)
  }
  return(es)
}

