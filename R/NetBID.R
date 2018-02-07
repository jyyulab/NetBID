#' Network-based Bayesian Inference of Drivers
#'
#' Driver inference of of case phenotype (vs. ctrl) from one source data and a pre-defined network
#'
#' @param gsc gsc (gene set collection) is a input network in list format with drivers as names and network targets as values that can be mapped to featureNames(eset)
#' @param eset ExpressionSet object of the data used to make inference: pData(eset)$group must exist to indicate the groups which you will compare to make the inference; featureNames(eset) are ids that can be mapped by gsc values.
#' @param group.case case name in pData(eset)$group
#' @param group.ctrl control name in pData(eset)$group
#' @param signed If TURE, consider directions of case vs. ctrl, If FALSE, ignore the directions
#' @return data.frame of case vs. ctrl inference with statistics including netSize (network size), z, pval, FDR.BH (FDR by BH method), log2FC (if signed)
#' @examples
#'
#' ###transform network into list of probes
#' dci.kinase.probes<-network.tbl2list(dci.kinase.tbl[,c('source.symbol','target')])
#' dmrna<-netbid(dci.kinase.probes,mrna,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')
#'
#'
#'
#' @export
netbid<-function(gsc,eset,group.case,group.ctrl,signed=FALSE){

	aset<-getActivity(gsc,eset)
	rs<-bid(aset,group.case,group.ctrl,signed)
	rs

}


#' NetBID Integration
#'
#' Driver inference by intergrating multiple NetBID results
#'
#' @param netbid.list list of netbid outputs, e.g. list(mRNA=dmrna,wProtein=dwprotein,pProtein=dpprotein)
#' @param signed If TURE, consider directions when combining NetBID results, If FALSE, ignore the directions
#' @return data.frame of case vs. ctrl inference with statistics including cobmbined and individual netSize (network size), z, pval, FDR.BH (FDR by BH method), log2FC (if signed)
#' @examples
#'
#' ###transform network into list of probes
#' dci.kinase.probes<-network.tbl2list(dci.kinase.tbl[,c('source.symbol','target')])
#' dmrna<-netbid(dci.kinase.probes,mrna,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')
#' ##whole proteomcis
#' dci.kinase.wprotein<-transformNetworkListIds(dci.kinase.symbols,ref=fData(wprotein)[,c('geneSymbol','proteinId')])
#' dwprotein<-netbid(dci.kinase.wprotein,wprotein,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')

#' dcomb<-netbidi(list(mRNA=dmrna,wProtein=dwprotein))
#'
#'
#'
#' @export
netbidi<-function(netbid.list,signed=FALSE){


	n<-length(netbid.list)
	if(n<2) stop('At least two netbid outputs are required for netbidi analysis!\n')

	rs<-NULL
	for(i in 1:n){

		if(names(netbid.list)[i]==''|is.null(names(netbid.list)[i]))
			names(netbid.list)[i]<-i

		names(netbid.list[[i]])[-1]<-paste(names(netbid.list[[i]])[-1],names(netbid.list)[i],sep='_')

		netbid.list[[i]]$id<-as.character(netbid.list[[i]]$id)
		if(i==1){
			rs<-netbid.list[[i]]
		}else{
			rs<-dplyr::full_join(rs,netbid.list[[i]],by='id')
		}
	}


	conds<-names(netbid.list)
	if(sum(grepl('log2FC',names(rs)))>0){
		FC.cols<-paste('log2FC',conds,sep='_')
	}else{
		FC.cols<-NULL
	}
	comb<-combinePvalues(rs,pvalue.cols=paste('pval',conds,sep='_'),sign.cols=paste('z',conds,sep='_'),FC.cols =FC.cols,logTransformed=TRUE,log.base=2,method=c('Stouffer'),twosided=TRUE,signed=signed,byRow=TRUE)
	comb<-comb[,-1]
	names(comb)<-paste(names(comb),'comb',sep='_')

	rs<-cbind(comb,rs)

	col.sel<-c('id',
			grep('^netSize_',names(rs),value=TRUE)
	,
			grep('^z_',names(rs),value=TRUE),
			grep('^log2FC_',names(rs),value=TRUE),
			grep('^pval_',names(rs),value=TRUE),
			grep('^FDR.BH_',names(rs),value=TRUE)
	)

	rs<-rs[,c(col.sel,setdiff(names(rs),col.sel))]
	rs<-dplyr::arrange(rs,pval_comb)

	rs

}
