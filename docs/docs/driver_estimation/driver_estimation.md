---
layout: default
title: Driver estimation
nav_order: 3
has_children: false
permalink: /docs/driver_estimation
---

# Driver Estimation

The purpose for this part: 

**estimate potential drivers in a biological process and generate a master table**.

The full demo script for this part could be found in [pipeline_analysis_demo1.R](https://github.com/jyyulab/NetBID-dev/blob/master/demo_scripts/pipeline_analysis_demo1.R).

----------
## Quick Navigation for this page

- [Step0: preparations](#step0-preparations)
- [Step1: load in gene expression dataset for analysis (exp-load,exp-cluster,exp-QC)](#step1-load-in-gene-expression-datasets-for-analysis-exp-load-exp-cluste-exp-qc)
   - [Q&A: What to do if the ID type is different between the network construction dataset and analysis dataset ?](#what-to-do-if-the-id-type-is-different-between-the-network-construction-dataset-and-analysis-dataset-)
- [Step2: read in network files and activity calculation (act-get)](#step2-read-in-network-files-and-activity-calculation-act-get)
    - [Q&A: Why study driver’s activity ?](#why-study-drivers-activity-) 
- [Step3: get differentiated expression/differentiated activity for all possible drivers (act-DA)](#step3-get-differentiated-expressiondifferentiated-activity-for-all-possible-drivers-act-da)
- [Step4: generate master table (ms-tab)](#step4-generate-master-table-ms-tab)
    - [Q&A: How to read and use the master table ?](#how-to-read-and-use-the-master-table-)
   
---------

## Step0: preparations

Library the installed NetBID2. 

```R
# library the package
library(NetBID2)
```

This tutorial is based on the suggested pipeline design in NetBID2, and before start, some parameters need to be set.

First, user need to know which network to use. If followed the tutorial in the [Network construction](../docs/network_construction) part, user could specify the `network.dir` and `network.project.name`.

```R
network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo network in the package
network.project.name <- 'project_2019-02-14' #
```

The `network.dir` is the path for the network project, i.e the `$project_main_dir/$project_name` when calling the `NetBID.network.dir.create()` function. 
The `network.project.name` must be the project name whening calling the `SJAracne.prepare()` function. 
Under the `network.dir` multiple `network.project.name` is allowed, and by specifying these two paramters simultaneously, unique network files (TF, SIG networks will be generated separately) could be obtained. 

In the above scripts, to ensure the testing of demo scripts in this part independently, user could use the demo network in the NetBID2 package. 

Next, we need to set up a main directory for saving the analysis results for this part, `project_main_dir`. User could choose to use the same directory as in the network construction. 
And the project name `project_name` should also be settled. Similarly, the main directory could have multiple projects under the path, separeted by the `project_name` as the name for the sub-directory. 
So user could use the same main directory for another project. But project related files with the same `project_main_dir` and `project_name` will be covered. 
So, in the script below, user could add a time tag in the `project_name` to avoid this by accident. 

```R
#set up paramters
project_main_dir <- 'test/' ### user defined main directory for the project, one main directory could have multiple projects, separeted by project name
current_date <- format(Sys.time(), "%Y-%m-%d") ## current date for project running, suggested to add the time tag
project_name <- sprintf('driver_%s',current_date) ## project name for the project
```

Once decided the `network.dir`, `network.project.name`, `project_main_dir` and `project_name`, user could run `NetBID.analysis.dir.create()` to generate the sub-directories for the working directory, including QC/ to save QC related files, DATA/ to save RData and PLOT/ to save the visulization plots in the [Advanced analysis](../docs/advanced_analysis) part. 
Besides, a global list variable `analysis.par` will be returned by the function. 
Attention, if the current environment already has this variable, the function will do nothing, report a warning message and return the original `analysis.par`.  

```R
## analysis.par is very essential in the analysis, if the environment already has this parameter, strongly suggest to delete it first
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, prject_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)
```

## Step1: load in gene expression dataset for analysis (exp-load,exp-cluster,exp-QC)

Here, the analysis expression dataset need to be loaded for analysis. For demo, we used the same dataset for network construction and analysis. 
If so, the expression dataset has already been processed by QC steps, just load from the RData in the `network.par$net.eset` and save it to `analysis.par$cal.eset`. 

```R
# if use same expression datasets as in the network construction, directly do the followings:
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) ## RData from network construction
analysis.par$cal.eset <- network.par$net.eset
```

However, in most of the cases, the datasets are not the same. User need to follow the similar steps in [Network construction](../docs/network_construction) for expression dataset loading, QC and sample cluster checking.

Save `analysis.par` into the RData file and give the step name to it, e.g `exp-QC`. 
The RData could be found in `analysis.par$out.dir.DATA/analysis.par.Step.{exp-QC}.RData`.

```R
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')
```

----------
### *What to do if the ID type is different between the network construction dataset and analysis dataset ?*

- It is strongly recommended that **the main ID should be the ID type in the network construction dataset !**

- The purpose of NetBID2 is to infer potential drivers in a biological process. 
The driver is defined as the source node in the network and its activity is calculated based on the network structure (its target nodes). 
If one gene does not exist in the pre-generated network, it will not be used in the following analysis. 
On the contrary, if one driver does not exist in the analyis expression dataset, it can still has its activity and could be used in the following analysis.

- So, we strongly recommend to align the ID type from the expression dataset to the network construction dataset by using `update_eset.feature()` and ID transfer table could be obtained from GPL files or use `get_IDtransfer()` or `get_IDtransfer_betweenSpecies()`. 

- The most complicated condition is that the level of network construction ID type and the analysis expression dataset is different (e.g transcript level Vs. gene level). `update_eset.feature()` could deal with such conditions by setting the `distribute_method` and `merge_method`. But be careful by doing so and try to avoid such condition !

----------

## Step2: read in network files and activity calculation (act-get)

Before start, load the RData from the previous step if user want to re-run the following steps or re-open the R session. 
Remember to set the temporary `analysis.par` if re-open the R session and `analysis.par$out.dir.DATA` will be used to find the correct RData file. 
No need to run this if continue working from the previous step. 

```R
## load from RData
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')
```

Firstly, get the network generated by [SJAracne](https://github.com/jyyulab/SJARACNe). 
If followed the pipeline, the file path will be `analysis.par$tf.network.file` and `analysis.par$sig.network.file`.
Use `get.SJAracne.network()` to read in the network information. If not followed the pipeline, it is also working by directly input the file path here.

```R
# get network info ! three list(network_dat, target_list, igraph_obj)
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
```

This function aims to read in network files generated by SJAracne and save the network information into three lists, `network_dat` is a data.frame to save all the information in the network file; `target_list` is a list for containing the target genes' information for the drivers. The names for the list is the driver name and each object in the list is a data.frame to save the target genes. `igraph_obj` is an igraph object to save the network, it is a directed, weighted network. The function will set two edge attributes to the `igraph_obj`, `weight` is the mutual information (MI) values and sign is the `sign` for the spearman value to indicate positive regulation (1) or negative regulation (-1).

The network dataset could be updated by using `update_SJAracne.network()`. 
This function allows user to filter drivers or to add drivers without edges by setting `all_possible_drivers`. Similar for targets by `all_possible_targets`.
Statistical filteration is also allowed by setting `min_MI`, `max_p.value`, `min_spearman_value`, `min_pearson_value` and the attribute of the network by setting `directed` and `weighted`.

Generate an html QC report file for the network file with the `igraph_obj`.

```R
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_')
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_')
```

The QC report files are generated for [TF network](TF_net_netQC.html) and [SIG network](SIG_net_netQC.html). 
Basic statistics including `size`, `diameter` and several kinds of centrality will be calculated. 
Scale free distrbution will be checked as well. 

Secondly, merge these two networks by `merge_TF_SIG.network()`. User could choose to merge the network first by `merge_TF_SIG.network()` and calculate driver activity based on the merged network or calculate the activity separately for TF and SIG network first and merge the activity results later by `merge_TF_SIG.AC()`. Both strategy will get the same final results. Drivers in the `analysis.par$merge.network` will have suffix of '_TF' and '_SIG' to indicate the driver type. One driver may have both '_TF' and '_SIG', but sometimes may have slightly different following results due to the network source. 

```R
# merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
```

Thirdly, get the activity matrix for all possible drivers by `cal.Activity()`. Activity will be used to estimate the driver's performance in a biological process, beside the expression level. 
Basic idea to get driver's activity is to estimate the cumulative effect on its targets. 
The summarization strategy could be 'mean', 'weighted mean', 'maxmean' and 'absmean', in which the weighted in the weighted mean is the MI (mutual information) value * sign of correlation (use the spearman correlation sign). If set `std=TRUE`, the expression matrix will be Z-transformed before calculating the activity. Higher expression of its positively regulated genes and lower expression of its negatively regulated genes indicate higher activity for this driver in the sample (if set 'weighted mean' as the calculation strategy). 

```R
# get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
```

Consider that the activity matrix is similar as the expression matrix that each row represents a driver with each column a sample. The value in the matrix is the activity level for the driver in the corresponding sample. Thus, we will also use the eSet class to save the data information. The phenotype information for the activity matrix is the same as for the analysis expression dataset. 

```R
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')
```

Finally, use `draw.eset.QC()` to get the QC report for the `analysis.par$merge.ac.eset`. Check [QC for AC](AC_QC.html). 
All functions related with eSet class manipulation are also applicable for `analysis.par$merge.ac.eset`. 

```R
# QC plots
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_')
```

Save `analysis.par` to the RData file by `NetBID.saveRData()`. 

```R
# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')
```

----------
### *Why study driver's activity ?*

- Drivers, such as transcription factors (TF) bind to the enhancer or promoter regions of the target genes and regulate their expression level. Once the TF has been synthesized, there are still many steps between mRNA translation of a TF and the actual transcriptional regulation of target genes. And the activity will be controlled by lots of following processes, such as: nuclear localization (need to be directed into nucleus), activation through signal-sensing domain (e.g ligand binding or post-translational modification: methylated, ubiquinated or phosphorylated,  phosphorylations are often necessary for dimerization and binding to the target gene’s promoter), access to DNA-binding site (epigenetic feature of the genome), and interaction with other cofactors or TFs (form a complex). Thus, sometimes the expression trend and the activity pattern of the driver may be violated. Due to some 'hidden' effect, the activity of a driver may be essential in regulating the process while the expression level may not show significant difference in the process.

----------

## Step3: get differentiated expression/differentiated activity for all possible drivers (act-DA)

Load the RData from the previous step if user want to re-run the following steps or re-open the R session. 
Remember to set the temporary `analysis.par` if re-open the R session and `analysis.par$out.dir.DATA` will be used to find the correct RData file. 
No need to run this if continue working from the previous step. 

```R
## load from RData
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')
```

The main purpose of this step is to get significantly differentiated expression (DE) and activity (DA) for all possible drivers between a specific process. 
For the RNASeq dataset, user could choose to use the DE output from DESeq2 but need to be attention about the column names (e.g pvalue) and prepare one column of `Z-statistics` in the following functions. 
Here, NetBID2 provides two functions `getDE.BID.2G()` and `getDE.limma.2G()` to assist the analysis between two groups (`G1` Vs. `G0`). 
User only need to get the sample list for `G1` and `G0` and specificy the name of `G1_name` and `G0_name` to use these two functions. 
`getDE.BID.2G()` will use the bayesian inference strategy to get the DE and DA, if user choose `method='Bayesian'`, the calculation will take more time compared with `method='MLE'`.

Ordinal phenotype condition is complicated and user could directly call `bid()` or functions in `limma()` to analyze. 

In the demo script below, user could get DE/DA for drivers between `G4.Vs.WNT` and `G4.Vs.SHH`. Each comparison must specify a clear name, and results will be saved into `analysis.par$DE` and `analysis.par$DA`. 
The comparison names will be displayed in the final master table. 

```R
# generate a list for multiple comparisons
analysis.par$DE <- list()
analysis.par$DA <- list()
# get compared sample names
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # get sample list for G1
# G4 Vs. WNT
comp_name <- 'G4.Vs.WNT'
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='WNT')] # get sample list for G0
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
# save to analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# G4 Vs. SHH
comp_name <- 'G4.Vs.SHH'
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='SHH')] # get sample list for G0
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
# save to analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid
```

If user want to combine the results from multiple comparisons, `combineDE()` could be applied just by creating a `DE_list` containing all the comparisons ready for merge. 
The output will be a list of DE/DA results, with one more component named "combine" that include the combined results. Give a `comp_name` and save into `analysis.par$DE` and `analysis.par$DA`. 

```R
## G4 Vs. others
# combine the above two comparisons
comp_name <- 'G4.Vs.otherTwo'
DE_gene_comb <- combineDE(DE_list=list(WNT=analysis.par$DE$`G4.Vs.WNT`,SHH=analysis.par$DE$`G4.Vs.SHH`))
DA_driver_comb <- combineDE(DE_list=list(WNT=analysis.par$DA$`G4.Vs.WNT`,SHH=analysis.par$DA$`G4.Vs.SHH`))
analysis.par$DE[[comp_name]] <- DE_gene_comb$combine
analysis.par$DA[[comp_name]] <- DA_driver_comb$combine
```

User could use `draw.combineDE()` to visualize the top significant DE/DA combining results compared with previous ones.

```R
draw.combineDE(DE_gene_comb)
draw.combineDE(DE_gene_comb,pdf_file=sprintf('%s/combineDE.pdf',analysis.par$out.dir.PLOT))
```

![`combineDE`](combineDE.png)

```R
draw.combineDE(DA_driver_comb)
draw.combineDE(DA_driver_comb,pdf_file=sprintf('%s/combineDA.pdf',analysis.par$out.dir.PLOT))
```

![`combineDA`](combineDA.png)


Also, user could choose to combine the sample list first and call DE/DA between the two sample lists. The results may be different from the above `combineDE` strategy as they are different in the statistical hypothesis stating. 

```R
# combine the sample list
comp_name <- 'G4.Vs.others'
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # get sample list for G0
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
# save to analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid
```

When finished, save to RData.

```R
# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')
```
Now, for the calculation part, we have obtained the top DA/DE list, we could directly draw the statistics for them (by default the top 30 drivers will be displayed):

```R
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others')
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others',pdf_file=sprintf('%s/NetBID_TOP.pdf',analysis.par$out.dir.PLOT),text_cex=0.8)
```

![`NetBID_TOP`](NetBID_TOP.png)

User could choose which column to display in DA/DE results, check `?draw.NetBID` for detailed instruction. 
**ATTENTION** for real practice, user may choose to display one comparison result, and the input of `DA_list` and `DE_list` must be the list class with names on it.

## Step4: generate master table (ms-tab)

Load the RData from the previous step if user want to re-run the following steps or re-open the R session. 
Remember to set the temporary `analysis.par` if re-open the R session and `analysis.par$out.dir.DATA` will be used to find the correct RData file. 
No need to run this if continue working from the previous step. 

```R
## load from RData
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')
```

Here, the main purpose is to generate the final master table for all possible drivers. 
User need to prepare:
- Use `db.preload()` to load in pre-saved TF/SIG gene's information by setting the `use_level`.
- Get the list of comparison names `use_comp` for output in the master table, if all, could use `all_comp <- names(analysis.par$DE)` to get them.
- The `DE` and `DA` results, if followed the above pipeline, just the `analysis.par$DE` and `analysis.par$DA`.
- The `network` target list file, i.e `analysis.par$merge.network$target_list`.
- The TF/SIG list `tf_sigs`, will be in the global environment if run `db.preload()`.
- The column name for Z statistics (`z_col`), for the output of `getDE.limma.2G()` and `getDE.BID.2G()`, it will be `Z-statistics`.
- Other column names wish to display in the final master table, e.g `logFC`, `P.Value`. 
- `main_id_type`, this may be the *only* important parameter need to set here if followed the above pipeline. Check [ID conversion section](../network_construction#id-conversion) for the detailed description of the id types.
- `transfer_tab`, the transfer table for ID conversion, could be obtained by `get_IDtransfer2symbol2type()`. If NULL, will automatically get the transfer table within the function if the `main_id_type` is not in the column names of 'tf_sigs'. User could also generate their own transfer table if do not want to use BioMart.
- `column_order_stratey`, an option to order the columns in the mater table, if set to `type`, the columns will be ordered according to the column type; if set to `comp`, the columns will be ordered according to the comparisons. 

```R
# load db
db.preload(use_level='gene',use_spe='human',update=FALSE)
# generate master table
all_comp <- names(analysis.par$DE) ## get all comparison name for output
# prepare transfer table (optinal)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes) ## get transfer table !!!
analysis.par$transfer_tab <- transfer_tab
# generate master table
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                               network=analysis.par$merge.network$target_list,
                                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                               main_id_type='external_gene_name')
```

Note: There may exist some drivers only have activity but no expression level, this is due to the difference between the network-construction dataset and analysis dataset mentioned above. 

Finally, output the master table into an excel file. Basic usage is just input the master table data frame into `out2excel()`. Marker genes `mark_gene` could also be highlighted. `mark_strategy` could be chosen from 'color' and 'add_column'. 'Color' means the mark_gene will be displayed by its background color; 'add_column' means the mark_gene will be displayed in separate columns with content TRUE/FALSE indicating whether the genes belong to each mark group. Check `?out2excel` for detailed usage description. Here, we create a list of MB subtype specific marker genes. The color could be random generated by `get.class.color()` or to user specified ones. Here, we could set up the `mark_col` as MB subtypes have pre-definied color codes.

```R
# output into excel files
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# can add marker gene
mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),
                  SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
                  G4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
#mark_col <- get.class.color(names(mark_gene)) # this will randomly generate color code
mark_col <- list(G4='green','WNT'='blue','SHH'='red')
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col)
```

Download the master table excel file [ms_tab.xlsx](driver_ms_tab.xlsx) to study the organization of the master table.

Save the `analysis.par` to the RData file. This file contains all information for this project, which could be used for the [Advanced analysis](../docs/advanced_analysis) part and [NetBID2 shiny server](https://github.com/jyyulab/NetBID_shiny). 
The `analysis.par` list need to include 13 components (`main.dir`, `project.name`, `out.dir`, `out.dir.QC`, `out.dir.DATA`, `out.dir.PLOT`, `merge.network`, `cal.eset`, `merge.ac.eset`, `DE`, `DA`, `final_ms_tab`, `transfer_tab`) in order to run in the shiny server. 

**Strongly suggest to save the RData file in this step !**

```R
# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')
```

-------
### *How to read and use the master table ?*

- `out2excel()` could output multiple master tables in different excel sheet. 
- In each master table, it is divided into three parts:
  - The first six columns are `gene_label`, `geneSymbol`, `originalID`, `originalID_label`, `funcType` and `Size`
    - `gene_label` is the driver's gene symbol or transcript symbol with suffix '_TF' or '_SIG' to indicate the driver's type. This column is often used for display as gene/transcript symbol is the most common accpeted gene ID for communication.
    - `geneSymbol` is the driver's gene symbol or transcript symbol.
    - `originalID` is the original ID type used in network construction, which should match the ID type in `analysis.par$cal.eset`, `analysis.par$DE`.
    - `originalID_label` is the original ID type with suffix '_TF' or '_SIG', which should match the the ID type in `analysis.par$merge.network`, `analysis.par$merge.ac.eset`,`analysis.par$DA`.
    - **originalID_label** is the only ensured unique ID for each row !!!
    - `funcType` is the 'TF' or 'SIG' used to indicate the driver's type.
    - `Size` is the target size for the driver. 
  - The main columns are named by $prefix.$comp_name_{DA or DE}, in which the prefix could be `Z`, `P.Value`, `logFC`, `AveExpr` to indicate the column's data type. `comp_name` is the comparison name in DA/DE. Columns of Z statistics will be automatically marked by the background color to indicate their significance.
  - The next 13 columns (start from `ensembl_gene_id` to `refseq_mrna`) are the information for the genes.
  - The last columns (optional) will be the marker information if set `mark_strategy='add_column'`.
- User could filter by the target size, sort the columns with Z-statistics to get top significant drivers. OR, follow the tutorial in [Advanced analysis](../docs/advanced_analysis) part for analyze and visualization.

-------

Finish driver estimation part !!! Cheers !!!

