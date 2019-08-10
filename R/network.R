#' Transform network ids
#'
#' Transform network ids given a reference, e.g. from gene symbols to probe ids
#'
#' @param net.list network in list format
#' @param ref data.frame with two columns: id in net.list, new id to tranform
#'
#' @return network in list format
#' @examples
#'
#' dci.kinase.wprotein<-transformNetworkListIds(dci.kinase.symbols,ref=fData(wprotein)[,c('geneSymbol','proteinId')])
#'
#' @export
transformNetworkListIds<-function(net.list,ref){
  names(ref)<-c('idInList','newId')
  for(i in 1:length(net.list)){
    net.list[[i]]<-unique(subset(ref,idInList%in%net.list[[i]])$newId)
  }
  net.list
}



#' Create network list format from table
#'
#' Transform network from table format to list format
#'
#' @param net network in tbl format, the first two columns: source, target
#' @param annotation if not NULL, it has two columns: the first column is ids that are mapped to source and target columns; the second column is the annotation such as geneSymbol, etc
#' @param annotateSize if TRUE, add size of the network
#'
#' @return network in list format
#'
#' @examples
#'
#' dci.kinase.probes<-network.tbl2list(dci.kinase.tbl[,c('source.symbol','target')])
#'
#' @export
network.tbl2list<-function(net,annotation=NULL,annotateSize=FALSE){

  names(net)[1:2]<-c('source','target')

  gs<-unique(net$source)

  gsc<-list()
  length(gsc)<-length(gs)

  if(is.null(annotation)){
    names(gsc)<-gs
  }else{
    annotation<-unique(annotation)
    names(gsc)<-paste(annotation[match(gs,annotation[,1]),2],gs,sep='_')
  }

  i<-1

  for(i in 1:length(gsc)){
    gsc[[i]]<-unique(subset(net,source==gs[i])$target)
    if(annotateSize){
      names(gsc)[i]<-paste(names(gsc)[i],length(gsc[[i]]),sep='_')
    }
  }

  gsc

}
