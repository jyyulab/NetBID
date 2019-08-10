# NetBID
data-driven Network-based Bayesian Inference of Drivers

## Reference: 
Du X, Wen J, Wang Y, Karmaus PWF, Khatamian A, Tan H, Li Y, Guy C, Nguyen TM, Dhungana Y, Neale G, Peng J, **Yu J***, Chi H*. Hippo/Mst signalling couples metabolic state and immune function of CD8&alpha;<sup>+</sup> dendritic cells. Nature. 2018;558(7708):141-5. doi: 10.1038/s41586-018-0177-0. PubMed PMID: [29849151.](https://www.nature.com/articles/s41586-018-0177-0) * Corresponding authors

## Installation
```r
devtools::install_github("jyyulab/NetBID")
```


## Tutorial with Examples

> Data

* dci.kinase.tbl: DC kinase network in tbl format
* dci.kinase.symbols: DC kinase network in list format with symbol as values
* mrna: ExpressionSet of mRNA epxression data
* wprotein: ExpressionSet of whole proteomics data


> Step 1

**Network Reconstruction:** 
Please refer to [SJARACNe](https://github.com/jyyulab/SJARACNe).

> Step 2 

**Driver inference from single data source:**

**mRNA data**

```r
library(NetBID)
###transform mRNA network into list of probes
dci.kinase.probes<-network.tbl2list(dci.kinase.tbl[,c('source.symbol','target')])
###get driver activty and differential activity analysis
dmrna<-netbid(dci.kinase.probes,mrna,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')
```
**Proteomic data**

```r
###transform proteomic network into list of genes
dci.kinase.wprotein<-transformNetworkListIds(dci.kinase.symbols,ref=fData(wprotein)[,c('geneSymbol','proteinId')])
###get driver activty and differential activity analysis
dwprotein<-netbid(dci.kinase.wprotein,wprotein,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')
```
>Step 3

**Integrate driver inference from multiple data sources:**

```r
dcomb<-netbidi(list(mRNA=dmrna,wProtein=dwprotein))
head(dcomb)
```
>Additional analysis

**Combine driver inference with differential expression analysis:**

```r
#####differential expression
de.mrna<-bid(mrna[fData(mrna)$geneSymbol%in%dcomb$id,],group.case='CD8pos.DC',group.ctrl='CD8neg.DC',signed=TRUE)
de.mrna<-filter(de.mrna,!duplicated(geneSymbol))%>%mutate(id.old=id,id=geneSymbol)

de.wprotein<-bid(wprotein[fData(wprotein)$geneSymbol%in%dcomb$id,],group.case='CD8pos.DC',group.ctrl='CD8neg.DC',signed=TRUE)
de.wprotein<-filter(de.wprotein,!duplicated(geneSymbol))%>%mutate(id.old=id,id=geneSymbol)

####comb diff-exp
decomb<-netbidi(list(mRNA=de.mrna,wProtein=de.wprotein),signed=TRUE)
head(decomb)

####
dfinal<-full_join(dcomb,decomb,by='id')
head(dfinal)

```

