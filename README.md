# NetBID
Network-based Bayesian Inference of Drivers


## Installation
```r
devtools::install_github("jyyulab/NetBID")
```


## Examples

> Data

* dci.kinase.tbl: DC kinase network in tbl format
* dci.kinase.symbols: DC kinase network in list format with symbol as values
* mrna: ExpressionSet of mRNA epxression data
* wprotein: ExpressionSet of whole proteomics data


> R Codes

```r

library(NetBID)

#####infer drivers from mRNA data
###transform network into list of probes
dci.kinase.probes<-network.tbl2list(dci.kinase.tbl[,c('source.symbol','target')])
dmrna<-netbid(dci.kinase.probes,mrna,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')

#####infer drivers from whole proteomcis data
dci.kinase.wprotein<-transformNetworkListIds(dci.kinase.symbols,ref=fData(wprotein)[,c('geneSymbol','proteinId')])
dwprotein<-netbid(dci.kinase.wprotein,wprotein,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')

#####integrate mRNA and whole proteomcis results
dcomb<-netbidi(list(mRNA=dmrna,wProtein=dwprotein))
head(dcomb)


#####diff-exp
de.mrna<-bid(mrna[fData(mrna)$geneSymbol%in%dcomb$id,],group.case='CD8pos.DC',group.ctrl='CD8neg.DC',signed=TRUE)
de.mrna<-filter(de.mrna,!duplicated(geneSymbol))%>%mutate(id.old=id,id=geneSymbol)

de.wprotein<-bid(wprotein[fData(wprotein)$geneSymbol%in%dcomb$id,],group.case='CD8pos.DC',group.ctrl='CD8neg.DC',signed=TRUE)
de.wprotein<-filter(de.wprotein,!duplicated(geneSymbol))%>%mutate(id.old=id,id=geneSymbol)

####comb diff-exp
decomb<-netbidi(list(mRNA=de.mrna,wProtein=de.wprotein,pProtein=de.pprotein),signed=TRUE)
head(decomb)

####
dfinal<-full_join(dcomb,decomb,by='id)
head(dfinal)

```

