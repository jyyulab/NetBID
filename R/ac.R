

#######get activity for expression in table format
###gsc: gene set collections
###eset: eset of expression profiles
###es.method: mean (default), mean, basmean
#' @export
getActivity<-function(gsc,dat,size.min=2,es.method='mean'){

  if(!is.list(gsc))
    stop('gene sets must in list format!\n')

  gsc<-gsc[sapply(gsc,length)>=max(2,size.min)]

  if(class(dat)=='ExpressionSet'){
    ac<-getActivity4eset(gsc,dat,es.method=es.method)
  }

  else if(is.data.frame(dat)|is.matrix(dat)){

    ac<-getActivity4tbl(gsc,dat,es.method=es.method)

  }else{
    stop('Woops, unsupported data format! We only supports data.frame or matrix or ExpressionSet format so far!\n')
  }

  ac

}



########scale activity
#' @export
scaleActivity<-function(dat){

  if(class(dat)=='ExpressionSet'){
    ac<-dat
    exprs(ac)<-as.matrix(scaleActivity4tbl(exprs(dat)))
  }

  else if(is.data.frame(dat)|is.matrix(dat)){

    ac<-scaleActivity4tbl(dat)

  }else{
    stop('Woops, unsupported data format! We only supports data.frame or matrix or ExpressionSet format so far!\n')
  }

  ac

}


getActivity4tbl<-function(gsc,exp,es.method='mean'){

  ac<-matrix(NA,nrow=ncol(exp),ncol=length(gsc),dimnames = list(colnames(exp),names(gsc)))

  exp<-apply(exp,2,std)

  for(i in 1:length(gsc)){
    gsc[[i]]<-unique(intersect(gsc[[i]],row.names(exp)))
    n<-length(gsc[[i]]);n
    colnames(ac)[i]<-paste(colnames(ac)[i],n,sep='_')

    if(n>0){
      ac[,i]<-apply(exp[gsc[[i]],],2,es,es.method)
    }

  }

  ac<-data.frame(t(ac))
  ac
}



getActivity4eset<-function(gsc,eset,es.method='mean'){

  require(Biobase)

  ac<-getActivity4tbl(gsc,exprs(eset),es.method=es.method)


  fd<-data.frame(id=names(gsc),netSize=as.numeric(gsub('(.*)_([0-9]+$)','\\2',row.names(ac))))
  row.names(fd)<-row.names(ac)<-fd$id

  pd<-pData(eset[,names(ac)])
  if(!all(row.names(pd)==names(ac)))
    stop("sample names don\'t match!\n")

  aset<-new("ExpressionSet",phenoData = new("AnnotatedDataFrame",pd),featureData=new("AnnotatedDataFrame",fd),annotation='',exprs=as.matrix(ac))

  aset
}


scaleActivity4tbl<-function(dat){

  #	ds<-data.frame(t(apply(dat,1,std)))
  ds<-data.frame(t(apply(dat,1,scale)))
  colnames(ds)<-colnames(dat)
  ds
}
