---
layout: default
title: Network construction
nav_order: 2
has_children: false
permalink: /docs/network_construction
---

# Network Construction

The purpose for this part: 

**generate a gene regulatory network based on a transcriptomic datasets**.

The full demo script for this part could be found in [pipeline_network_demo1.R](https://github.com/jyyulab/NetBID-dev/blob/master/demo_scripts/pipeline_network_demo1.R).

----------
## Quick Navigation for this page

- [Step0: preparations](#step0-preparations)
- [Step1: load in gene expression datasets for network construction (exp-load)](#step1-load-in-gene-expression-datasets-for-network-construction-exp-load)
   - [Q&A: The choice of expression dataset for network construction](#the-choice-of-expression-dataset-for-network-construction)
   - [Q&A: Input RNASeq dataset](#input-rnaseq-dataset)
   - [Q&A: Input expression matrix](#input-expression-matrix)
- [Step2: normalization for the expression dataset (exp-QC)](#step2-normalization-for-the-expression-dataset-exp-qc)
   - [Q&A: QC for RNASeq dataset](#qc-for-rnaseq-dataset)
   - [Q&A: Combine two datasets](#combine-two-datasets)
- [Step3: check sample cluster information, optional (exp-cluster)](#step3-check-sample-cluster-information-optional-exp-cluster)
- [Step4: prepare SJARACNE (sjaracne-prep)](#step4-prepare-sjaracne-sjaracne-prep)
   - [Q&A: ID conversion](#id-conversion)
   
----------


## Step0: preparations

Before start, we need to library the installed NetBID2. 

```R
# library the package
library(NetBID2)
```

This tutorial is based on the suggested pipeline design in NetBID2. 
So, we need to set up a main directory for the project in `project_main_dir`. 
The main directory could have multiple projects under the path, separeted by the `project_name` as the name for the sub-directory. 
So user could use the same main directory for another project. 
But project related files with the same `project_main_dir` and `project_name` will be covered. 
So, in the script below, user could add a time tag in the `project_name` to avoid this by accident. 

```R
#set up paramters
project_main_dir <- 'test/' ### user defined main directory for the project, one main directory could have multiple projects, separeted by project name
current_date <- format(Sys.time(), "%Y-%m-%d") ## current date for project running, suggested to add the time tag
project_name <- sprintf('project_%s',current_date) ## project name for the project
```

Once decided the `project_main_dir` and `project_name`, user could run `NetBID.network.dir.create()` to generate the sub-directories for the working directory, including QC/ to save QC related files, DATA/ to save RData and SJAR/ to save data for running [SJARACNe](https://github.com/jyyulab/SJARACNe). Besides, a global list variable `network.par` will be returned by the function. 
Attention, if the current environment already has this variable, the function will do nothing, report a warning message and return the original `network.par`.  

```R
## network.par is very essential in the analysis, if the environment already has this parameter, strongly suggest to delete it first
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,prject_name=project_name) ## create working directory
```

## Step1: load in gene expression datasets for network construction (exp-load)

Here, we use same demo dataset for network construction and following analysis (Check the ***The choice of expression dataset for network construction*** section below). 
This dataset could be directly downloaded from GEO database by input the GSE ID and GPL ID.
If set `getGPL=TRUE`, will download the gene annotation file. 
The output of this function will be the [eSet](https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet) class object and will save the RData into `out.dir`.
Next time by running this function, it will try to load the RData in the `out.dir/{GSE}_{GPL}.RData` first (if `update=FALSE`).

```R
# from GEO, need to provide GSE and GPL
net_eset <- load.exp.GEO(out.dir=network.par$out.dir.DATA,GSE='GSE116028',GPL='GPL6480',getGPL=TRUE,update=FALSE)
```

*Optional:*
User could choose to update the feature data frame in the eSet object by `update_eset.feature()`. 
This function allows user to change the main ID from `from_feature` to `to_feature` by inputting the transfer table `use_feature_info` and choosing the `merge_method`.

```R
# ID conversion, or merge transcript level to expression level, use_feature_info can be other dataframe info; optional;
net_eset <- update_eset.feature(use_eset=net_eset,use_feature_info=fData(net_eset),from_feature='ID',to_feature='GENE_SYMBOL',merge_method='median') ### !!! need modify
```

*Optional:*
User could choose to update the phenotype data frame in the eSet object by `update_eset.phenotype()`. 
This function allows user to get phenotype information from the data frame `use_phenotype_info` by indicating the column that matched the sample name.
The `use_col` could be used to tell the function which column in `use_phenotype_info` will be kept. If set to 'auto', wil extract columns with unique sample feature ranges from 2 to sample size-1. 
If set to 'GEO-auto', will extract columns: 'geo_accession','title','source_name_ch1',and columns end with ':ch1'.

```R
# select phenotype columns or user add phenotype info; optional
net_eset <- update_eset.phenotype(use_eset=net_eset,use_phenotype_info=pData(net_eset),use_sample_col='geo_accession',use_col='GEO-auto')
```

Now, we need to check the data quality of the eSet by `draw.eset.QC()`. 
The QC report html mainly contains four parts, *heatmap*, *pca*, *density* and *meansd*.
The `intgroup` could be used to indicate which column from the `fData(eset)` will be used in the plot (*heatmap*, *pca*, *density*).
If set to NULL, will automatcially extract all possible groups by `get_int_group()`.

For the usage of `draw.eset.QC`, pandoc is required for `generate_html=TRUE`. 
If `pandoc_available()=FALSE`, please install pandoc and set environment of pandoc by `Sys.setenv(RSTUDIO_PANDOC=---installed path---)`.

```R
# QC for the raw expdatasets
# if intgroup==NULL, will auto get intgroup
# prefix: prefix for the pdf files
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='beforeQC_')
```

- What to check for the QC report html [before_QC.html](beforeQC_QC.html) ? 
The number of samples and genes (probes/transcripts/...);
All samples will be clustered by the expression pattern from all genes, check possible mis-labeled samples; 
The density plots show the range and distribution for the expression values, may judge whether the original dataset has been log transformed;
The meansd plots show the relationship between the mean and standard deviation of the genes, used to verify whether there is a dependence of the standard deviation (or variance) on the mean. 

Now, basic processing steps are finished for the input expression dataset. Save it to `network.par`.

```R
# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
```

Save `network.par` into the RData file and give the step name to it, e.g `exp-load`. 
The RData could be found in `network.par$out.dir.DATA/network.par.Step.{exp-load}.RData`.

```R
# save to RData
NetBID.saveRData(network.par = network.par,step='exp-load')
```

----------
### *The choice of expression dataset for network construction*

- For a NetBID2 project, the analysis expression dataset is selected first to assist the investigation of a biological story. 
The network construction expression dataset could be the same as the analysis expression dataset but need to consider some other factors. 

   - The theory of using expression dataset to infer gene regulatory networks is based on [SJARACNe](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty907/5156064). 
   It uses an information theoretic approach to eliminate the majority of indirect interactions inferred by co-expression methods. More samples, higher sensitivity and precision will be obtained by experiment. 
   Typically, more than 100 samples is a better choice. 
   - Large size public datasets from the same tissue, cell line or biological background as the analysis dataset are recommended. User could search the databases such as GEO, TCGA. 
   - Computational inferred networks will surely have false positive edges, especially for those with relative small mutual information (MI) score. Functions related with network processing will be described in the [Driver estimation](../driver_estimation) part. 
   - Once a high quality network is generated, user could put them into a common shared place that multiple projects with similar biological background could rely on that. 

- The demo used in the tutorial actually is not a good network construction expression dataset in real practice. Just used to assist user get familar with the procedure of NetBID2. 

----------

### *Input RNASeq dataset*   

- Another two functions can be applied to load expression dataset from RNASeq, `load.exp.RNASeq.demo()` and `load.exp.RNASeq.demoSalmon()`. 
**BUT** this two function are just demo functions, which do not support complicated options in `tximport()` and `DESeq()`. 
Besides, the output format may be different by using different dataset as reference that these two load functions may not work well.
Suggest to use the original functions if have some experience of coding.

- If try to use `load.exp.RNASeq.demo()` and `load.exp.RNASeq.demoSalmon()`, be **ATTENTION** to the `return_type` option in the functions. 
    - 'txi' is the output of `tximport()`, a simple list containing matrices: abundance, counts, length.
    - 'counts' is the output of raw count matrix.
    - 'tpm' is the output of raw tpm.
    - 'dds' is the DESeqDataSet class object, the data has been processed by `DESeq()`.
    - 'eset' is the ExpressionSet class object, the expression data matrix has been processed by `DESeq(), vst()`.
    
    Default is 'tpm'. If user do not choose 'eset', the output will not be directly used in the following scripts in the tutorial. Check the ***Input expression matrix*** section below.

----------
### *Input expression matrix*   
- If the user decide to prepare the expression matrix by themselves. The eSet object could be directly obtained by using `generate.eset()`.
For example for RNASeq dataset with 'tpm' as the output:
```R
#tpm <- load.exp.RNASeq.demo(XXX)
tmp_mat  <- log2(tpm)
tmp_eset <- generate.eset(exp_mat = tmp_mat, phenotype_info = NULL,feature_info = NULL, annotation_info = "")
```
For the option of `generate.eset()`, 
if `phenotype_info = NULL`, a one column named with 'group' will be automatically generated. 
if `feature_info = NULL`, a one column named with 'gene' will be automatically generated. 

----------

## Step2: normalization for the expression dataset (exp-QC)

Before start, load the RData from the previous step if user want to re-run the following steps or re-open the R session. 
Remember to set the temporary `network.par` if re-open the R session and `network.par$out.dir.DATA` will be used to find the correct RData file. 
No need to run this if continue working from the previous step. 

```R
# load from RData
#network.par <- list()
#network.par$out.dir.DATA <- 'test//project_2019-05-02//DATA/'
NetBID.loadRData(network.par = network.par,step='exp-load')
```

Following basic QC steps are suggested procedure for **microarray** dataset, all of the steps are ***optional***.

Firstly, check the `NA` values in the expression dataset by counting the number of `NA` values for each sample and for each gene (or probes/transcripts/...). 
If one sample or gene with too many `NA` values, user could choose to remove that gene or sample or do imputation by `impute.knn()`. 

```R
## following QC steps are optional !!!!
mat <- exprs(network.par$net.eset)
# remove possible NA? or imputation ? ## need to user-decide
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data
```

Secondly, the log2 transformation. Sometimes it is hard to know whether the original dataset has been log2 transformed or not. 
Here, we provide an experienced judging threshold for the median value. This is may not be suitable for all conditions. 

```R
# log2 transformation
med_val <- median(apply(mat,2,median));print(med_val)
if(med_val>16){mat <- log2(mat)}
```

Thirdly, the quantile normalization between samples. This is suggested for dealing with microarray dataset and not for RNASeq, even the log2tpm etc.  

```R
# quantile normalization
mat <- normalizeQuantiles(mat) ## limma quantile normalization
```

Fourthly, remove low expressed genes across nearly all samples. 
The suggested threshold is shown below, try to remove genes that in more than 90% samples, the expression value is lower than 5%. 

```R
# remove low expressed genes, such as whether in more than 90% samples, the expression value is lower than 5%
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]
```

Now, the expression matrix has been updated, need to save into the eSet class object by using `generate.eset()`. 
Update the `network.par$net.eset`, generate the QC report html by `draw.eset.QC()` and save to RData by `NetBID.saveRData()`. 

```R
# update eset
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),],
                          annotation_info=annotation(network.par$net.eset))
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='afterQC_')
# save to RData
NetBID.saveRData(network.par = network.par,step='exp-QC')
```

- What to check for the QC report html after QC steps [after_QC.html](afterQC_QC.html)? 
The number of samples and genes (probes/transcripts/...), whether large amount of genes/samples are removed;
All samples will be clustered by the expression pattern from all genes, check possible mis-labeled samples; 
The density plots show the range and distribution for the expression values, whether the low expressed genes have been removed;
The meansd plots show the relationship between the mean and standard deviation of the genes, used to verify whether there is a dependence of the standard deviation (or variance) on the mean. 


----------
### *QC for RNASeq dataset*  
- No matter what strategy is used to input the RNASeq dataset, only the fourth step 'remove low expressed genes' is suggested. 
For example, if use `load.exp.RNASeq.demo()` or `load.exp.RNASeq.demoSalmon()` and output `dds`, no previous normlization is required.
- If user use the raw count as the expression matrix, the `RNASeqCount.normalize.scale()` could be used to normalize the count data, followed by 'log2 transformation'.
- For the fpkm(Fragments per kilobase of exon per million reads mapped ), tpm(Transcripts Per Million), cpm(Counts Per Million), 
the second step 'log2 transformation' is suggested. 
- User is strongly suggested to judge which QC step to use for their own dataset or follow the pipeline suggested by the different calling software. 

----------
### *Combine two datasets*  
- If user want to combine two datasets, `merge_eset()` could be used. 
- If the original two datasets are generated from the same platform with the same expression gene list, no Z-transformation will be performed; 
otherwise Z-transformation will be performed before merging the dataset. 
- The merged eSet will automatically generate one phenotype column named by `group_col_name`. By default, the function will not remove batches between the two datasets. 
Strongly suggest to remove the batch for microarray dataset but not for RNASeq dataset. 
For the RNASeq dataset, user could follow the tutorial in the Step3 for detailed sample clustering checking and re-run the code in this step. 

----------


## Step3: check sample cluster information, optional (exp-cluster)

Before start, load the RData from the previous step if user want to re-run the following steps or re-open the R session. 
Remember to set the temporary `network.par` if re-open the R session and `network.par$out.dir.DATA` will be used to find the correct RData file. 
No need to run this if continue working from the previous step. 


```R
# load from RData
#network.par <- list()
#network.par$out.dir.DATA <- 'test//project_2019-05-02//DATA/'
NetBID.loadRData(network.par = network.par,step='exp-QC')
```

Select the most variable genes for the sample clustering analysis by `IQR.filter()`. 
In the script below, the most 50% variable genes will be used. 
For the `IQR.filter()` function, it allows user to input a list of genes (`loose_gene`), which can be applied to `loose_thre`. 
This is applicable when user need to keep more interested genes (e.g transcription factors) by using a looser threshold for filteration. 

```R
# use most variable genes for cluster
mat <- exprs(network.par$net.eset)
choose1 <- IQR.filter(exp_mat=mat,use_genes=rownames(mat),thre = 0.5,loose_gene=NULL,loose_thre=0.1)
print(table(choose1))
mat <- mat[choose1,]
```

Generate a temporary eSet and get the html QC report [Cluster_QC.html](Cluster_QC.html). 
Here give a first galance of sample clustering results Vs. pre-defined sample groups. 

```R
# generate tmp eset
tmp_net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),], annotation_info=annotation(network.par$net.eset))
# QC plots
draw.eset.QC(tmp_net_eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='Cluster_')
```

In the following scripts, user could get lots of plots. 
Firstly, get the phenotype information data frame `pData(network.par$net.eset)` and all possible phenotype classes `intgroup`.

For each `intgroup`, use could choose to use `draw.pca.kmeans()` or `draw.umap.kmeans()` to display the sample clustering results between the observed label and predicted label. The prediction is performed by kmeans based on pca or umap dimension reduction results. 
If user use [MICA](https://github.com/jyyulab/scMINER/tree/master/MICA) for clustering, `draw.MICA()` also could be used for result display.

The output for those functions will be the predicted label for the best `k` if setting `return_type='optimal'` or the results for `all_k` if `return_type='all'`.

```R
# more cluster functions (will not directly save to file, but actively layout)
phe <- pData(network.par$net.eset)
intgroup <- get_int_group(network.par$net.eset)
# pca+kmeans in 2D
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]))
}
```

Take `subgroup` as an example, the following scripts will generate:

```R
use_int <- 'subgroup'
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D')
```

![sample_cluster_1](sample_cluster_1.png)

This is the basic scatter plot to display the samples with color coded by observed and predicted label. 
The statistics in the right figure is the score between predicted label and observed label by `get_clustComp()`. 
ARI stands for 'adjusted rand index', which ranges from 0 to 1 with higher value indicates higher similarity. 

```R
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.ellipse')
```

![sample_cluster_2](sample_cluster_2.png)

This is the scatter plot with ellipse to cover the points belong to one class. 

```R
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.text')
```

![sample_cluster_3](sample_cluster_3.png)

This is the scatter plot with sample name directly labelled on the plot, which is useful for outlier checking. 

```R
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='3D')
```

![sample_cluster_4](sample_cluster_4.png)

This is the 3D scatter plot.

```R
print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe, use_int))))
draw.clustComp(pred_label,obs_label=get_obs_label(phe,use_int),outlier_cex=1,low_K=10) ## display the comparison in detail
```

![sample_cluster_4](sample_cluster_5.png)

This is the table to display the detailed difference between predicted label and observed label. We can see from the table here, 4 WNTs are further separated into two groups.

In this demo dataset, no clear outlier samples are observed. If user find some somes need to remove, please remove the samples and re-run the exp-QC steps. 

## Step4: prepare SJARACNE (sjaracne-prep)

Before start, load the RData from the 'exp-QC' step if user want to re-run the following steps or re-open the R session. 
Remember to set the temporary `network.par` if re-open the R session and `network.par$out.dir.DATA` will be used to find the correct RData file. 
No need to run this if continue working from the previous step. 

```R
# load from RData
#network.par <- list()
#network.par$out.dir.DATA <- 'test//project_2019-05-02//DATA/'
NetBID.loadRData(network.par = network.par,step='exp-QC')
```

Load the transcription facotrs (**TF**) and signaling factors (**SIG**) list from database by `db.preload()`.NetBID2 has prepared both 'gene' and 'transcript' level RData files for human.

If user's input is not human, could run `db.preload()` to prepare the database files. 
If leave `main.dir=NULL`, the RData will be saved to `system.file(package = "NetBID2")/db/`. 
If NetBID2 is installed in a public place with no permission to user, just set `main.dir` to another place and remember to use the same path next time using it.

For the TF and SIG list, NetBID2 has provided the files in 'external_gene_name' and 'ensembl_gene_id' ID type for human and mouse. (e.g `MOUSE_SIG_ensembl_gene_id.txt` in `system.file(package = "NetBID2")/db/`). 
The function will automatically use those files if set `TF_list=NULL` or `SIG_list=NULL`. 
User could also input their own list and input by setting `TF_list` or `SIG_list`.

```R
# load database
db.preload(use_level='gene',use_spe='human',update=FALSE)
```

After loading the database, user need to set the ID attribute type `use_gene_type` for the input expression matrix. Check ***ID conversion*** section below for detailed description of ID conversion issue. 

```R
# ID convertion, get TF/SIG list !!!!
use_gene_type <- 'external_gene_name' ## this should user-defined !!!
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
```

The above two steps are *not required* if user could get the `TF_list` and `SIG_list` with the same ID type as the expression matrix, just input them in the `SJAracne.prepare()`. 

The final step is to prepare the input for running SJAracne. 

User could choose to use part of the samples or use all. 
And in one project, multiple networks could be generated by setting different `prj.name`. For example, if want to generate Group4 specific network, user could choose to use samples in Group4 by setting the `use.samples` and give it an easily identified project name such as `prj.name='Group4_net'`. This `prj.name` is important in the [Driver estimation](../driver_estimation) part. 

The `IQR.thre` and `IQR.loose_thre` will be passed to `IQR.filter()`. The `loose_gene` in this function will be the genes in `TF_list` and `SIG_list` as we want to keep more possible drivers in the network construction.
In the demo network of NetBID2, in order to control the file size, the `IQR.thre=0.9` and `IQR.loose_thre=0.7`. In real practice, `IQR.thre=0.5` and `IQR.loose_thre=0.1` is recommended.

```R
# select sample for analysis
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples
prj.name <- network.par$project.name # can use other names, if need to run different use samples
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                    TF_list=use_list$tf,SIG_list=use_list$sig,
                    IQR.thre = 0.5,IQR.loose_thre = 0.1,
                    SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)
```

Next is to follow the message to run [SJAracne](https://github.com/jyyulab/SJARACNe). 
As SJAracne will consume lots of memory in running, which may be not suitable in R session, user need to follow the instructions to run SJAracne.

----------
### *ID conversion*  

We will use the ID name from [biomaRt](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html).
Some common attribute names are 'ensembl_transcript_id', 'ensembl_gene_id', 'external_transcript_name', 'external_gene_name', 'hgnc_symbol', 'entrezgene', 'refseq_mrna'. 
If the original input is 'ensembl_gene_id_version' or 'ensembl_transcript_id_version', user could set `ignore_version=TRUE` to neglect the version number. 

**ATTENTION!** 
- biomaRt will use the newest version number of [GENCODE](https://www.gencodegenes.org) and all the ID conversion related functions `db.preload(), get.TF_SIG.list(), get_IDtransfer(), get_IDtransfer2symbol2type(), get_IDtransfer_betweenSpecies()` will remotely call the database from biomaRt through the web link.
So, the version number may be different when running the same code at different time. 
- `get_IDtransfer(), get_IDtransfer2symbol2type(), get_IDtransfer_betweenSpecies()` will output a transfer table used for `get_name_transfertab()`. User could use their curated one. 

----------

Finish network construction part !!! Cheers !!!




