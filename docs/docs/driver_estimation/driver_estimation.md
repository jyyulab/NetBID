---
layout: default
title: "- Driver estimation"
nav_order: 5
permalink:  /docs/driver_estimation
---

## Driver Estimation

The purpose of this part: 

**retrieve potential drivers for interested phenotype and generate a master table for drivers**.

One "lazy mode" function without flexiable options is available for this part of analysis: `NetBID.lazyMode.DriverEstimation()`. 
User could check the manual for this function and try the demo code for usage.

The complete step-by-step demo script for network construction can be found here, 
[pipeline_analysis_demo1.R](https://github.com/jyyulab/NetBID-dev/blob/master/demo_scripts/pipeline_analysis_demo1.R).

----------
## Quick Navigation for this page

- [Step 0: Preparations](#step-0-preparations)
- [Step 1: Load in the expression dataset for analysis (exp-load, exp-cluster, exp-QC)](#step-1-load-in-the-expression-dataset-for-analysis-exp-load-exp-cluster-exp-qc)
   - [Q&A: What to do if the ID types from network-construction dataset and analysis dataset are different?](#what-to-do-if-the-id-types-from-network-construction-dataset-and-analysis-dataset-are-different)
- [Step 2: Read in network files and calcualte driver activity (act-get)](#step-2-read-in-network-files-and-calcualte-driver-activity-act-get)
    - [Q&A: Why study driver’s activity ?](#why-study-drivers-activity-) 
- [Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA)](#step-3-get-differential-expression-de--differential-activity-da-for-drivers-act-da)
- [Step 4: Generate a master table for drivers (ms-tab)](#step-4-generate-a-master-table-for-drivers-ms-tab)
    - [Q&A: How to interpret and use the master table ?](#how-to-interpret-and-use-the-master-table-)
   
---------

## Step 0: Preparations
**Purpose: create an organized working directory for the driver estimation step in NetBID2 analysis.**

Make sure you have NetBID2 package. 

```R
library(NetBID2)
```

**First, retrieve a constructed network for driver analysis.**
If user followed the [Network construction](../docs/network_construction) tutorial, 
the path of the network project `network.dir` should have been set. For details, please check [Network construction: Step 0](../docs/network_construction#step-0-preparations).
The `network.project.name` is used to distinguish different network construction jobs when using `SJAracne.prepare()`, under the main path of `network.dir`.
By specifying `network.dir` and `network.project.name`, user should be able to retrieve the target network constructed by network construction part in NetBID2 and SJARACNe.

For the online tutorial, we have already prepared the demo dataset's network constructed by SJARACNe. So users don't need to run SJARACNe for the demo, and can direclty proceed to the following pipeline.

```R
network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo network in the package
network.project.name <- 'project_2019-02-14' 
```
**Next, create directories and folders to save and organize your analysis results.**

We have designed a function `NetBID.analysis.dir.create()` to handle the working directories, so users can have a better data organization. Similar to the one in network construction.
This function needs users to define the main working directory `project_main_dir` and the project’s name `project_name`. 
To prevent previous project with the same `project_main_dir` and `project_name` from being rewrite, it is highly suggested to add a time tag to your `project_name`.

```R
# Define main working directory and project name
project_main_dir <- 'test/' # user defined main directory for the project, one main directory could
current_date <- format(Sys.time(), "%Y-%m-%d") # optional, if user like to add current date to name the project folder
project_name <- sprintf('driver_%s',current_date) # project name for the project folders under main directory
```

`NetBID.analysis.dir.create()` creates a main working directory with a subdirectory of the project. It also automatically creates three subfolders (QC, DATA and PLOT) within the project folder. QC/, storing Quality Control related plots; DATA/, saving data in RData format; PLOT/, storing output plots. 
It also returns a list object, here named `analysis.par` with directory information wrapped inside. 
This list is an ESSENTIAL variable for driver estimation step, all the important intermediate data generated later will be wrapped inside.

```R
# Create a hierarchcial working directory and return a list contains the hierarchcial working directory information
# This list object (analysis.par) is an ESSENTIAL variable in driver estimation pipeline
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)
```

## Step 1: Load in the expression dataset for analysis (exp-load, exp-cluster, exp-QC)
**Purpose: load in the expression dataset for driver estimation analysis step.**

For driver estimation, we will need an expression dataset containing the interested phenotype with proper control samples. Compared with the dataset used for network construction, there is no strong requirement for large sample size. User could choose to use the same dataset in network construction for this part but not necessary. For demo, we choose to use the same dataset here. 

If user choose to use the same dataset, no pre-processing of the dataset is required, since we saved the expression dataset after quality control as RData from network construction part in NetBID2, we can load it back directly. Users need to assign `network.par$net.eset` to `analysis.par$cal.eset`, cause for driver estimation step, the ESSENTIAL variable is `analysis.par`.

```R
# If use the same expression dataset as in the network construction, just reload it directly
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) # RData saved after QC in the network construction step
analysis.par$cal.eset <- network.par$net.eset
```

If user choose to use a different dataset, please complete the first three steps in [Network construction](../docs/network_construction) first to process the dataset before further analysis.

For persistent storage of data, and prevent the re-run of all previous steps. Users can checkout and save `analysis.par` as RData for this part.
The function `NetBID.saveRData()` provide easier pipeline step checkout and reference. 
Giving pipeline step name to `step`, the `analysis.par` will be saved with step name as `analysis.par$out.dir.DATA/analysis.par.Step.{exp-QC}.RData`.

```R
# Save Step 1 network.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')
```

----------
### *What to do if the ID types from network-construction dataset and analysis dataset are different?*

- **It is highly suggested to use the ID type from network-construction dataset as the main ID type.**

- The purpose of NetBID2 is to find potential drivers in a biological process of interest based on the data-driven gene regulatory network.  Each drivers' activity is evaluated based on the network structure (its directly targeted genes). If a driver doesn't exist in the pre-generated network, it will not be used in activity calculation. On the contrary, if a driver doesn't exist in the analysis dataset, it can still have its activity been calculated (if its target genes' expression values are available).

- We strongly recommend to align the ID type from the analysis dataset to the network-construction dataset using `update_eset.feature()`. The ID conversion table can be obtained from GPL files or `get_IDtransfer()` and `get_IDtransfer_betweenSpecies()`. 

- The most complicated situation is, the levels of network-construction ID type and the analysis expression dataset are different (e.g. transcript level vs. gene level).
NetBID2 provides `update_eset.feature()` to solve this problem, by assigning `distribute_method` and `merge_method` parameters.

----------

## Step 2: Read in network files and calcualte driver activity (act-get)

Please skip the following line if you didn't close R session after completed Step 1.

Don't skip, if you have checked out and closed R session after completed the Step 1. Before start Step 2, please reload `analysis.par` RData from Step 1.
`NetBID.loadRData()` reloads RData saved by `NetBID.saveRData()`. It prevents user from repeating former pipeline steps.
If the re-opened R session doesn't have `analysis.par` in the environment, please comment off the first two command lines. It will create a temporary `analysis.par`
with path of the saved Step 1 RData, `analysis.par$out.dir.DATA`. The path `test//driver_2019-05-06//DATA/` here is just an example, 
users need to give their own path used to save `analysis.par` RData from Step 1.

```R
# Reload network.par RData from Step 1
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')
```

**Firstly, get the network constructed by [SJAracne](https://github.com/jyyulab/SJARACNe).** If one followed the pipeline, the file path should be `analysis.par$tf.network.file` and `analysis.par$sig.network.file`.
Use `get.SJAracne.network()` to read in the network information. 
If one didn't follow the pipeline, just pass the directory of network to `get.SJAracne.network()`.

```R
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
```

**More about `get.SJAracne.network()` and related function.** It reads SJARACNe network construction result and returns a list object contains three elements, `network_data`, `target_list` and `igraph_obj`. `network_dat` is a data.frame, contains all the information of the network SJARACNe constructed. `target_list` is a driver-to-target list object. Please check details in `get_net2target_list()`. `igraph_obj` is an igraph object used to save this directed and weighted network. Each edge of the network has two attributes, weight and sign. Weight is the "MI (mutual information)" value and sign is the sign of the spearman correlation coefficient (1, positive regulation; -1, negative regulation). To updata the network dataset, user can call `update_SJAracne.network()`. It updates the network object created by `get.SJAracne.network`, using constraints like statistical thresholds and interested gene list.

**Generate an HTML QC report for the constructed network, using `igraph_obj`.**

```R
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_',html_info_limit=FALSE)
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_',html_info_limit=TRUE)
```

Two QC reports have been created. One for the transcription factor [TF network](TF_net_netQC.html), the other one is for signaling factor [SIG network](SIG_net_netQC.html). 

**- What information can you get from the HTML QC report of network?**
  - A table for network. It shows some basic statistics to characterize the target network, including size and centrality etc.
  - A table for drivers. It shows detailed statistics for all drivers in the network.
  - A density over histogram to show the distribution of the degree of nodes and the target size of all drivers. The average target size around several hundreds may be good by experience. 
  - A scatter plot to check if the network is scale-free. It is good but not necessary for a gene regulatory network to have the scale-free topology feature. 

**Second, merge TF-network and SIG-network.**
There are two ways to merge the networks. (1) Merge the networks first using `merge_TF_SIF.network()`, then calculate driver's activity value. (2) Calculate the driver's activity value in both TF-network and SIG-network, then merge them together using `merge_TF_SIG.AC()`. Both ways give the same result. 
Drivers in the `analysis.par$merge.network` will have suffix of '_TF' and '_SIG' to distinguish the type of the driver. It is possible to see a driver with both '_TF' and '_SIG'.

```R
# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
```

**Third, get the activity matrix for all possible drivers using `cal.Activity()`.** 
Driver's activity value quantifies its influence in a biological process. It is evaluated from the driver's accumulative effects to its targets.
The evaluation strategies are "mean"", "weighted mean", "maxmean" and "absmean". "Weighted mean" is the MI (mutual information) value with the sign of the spearman correlation. 
For example, if users choose "weighted mean" to calcualte the activity of driver. The higher expression value of its positively-regulated genes and the lower expression value of its negatively-regulated genes, the higher activity value of that driver will be.
If user would like to perform Z-transformation to the expression matrix before calculating the activity values, he can set `std=TRUE`.

```R
# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
```
Now, we have an activity matrix `ac_mat` for drivers. Rows are drivers, columns are samples. 
Due to the similar way of data display, and same phenotype information sharing, we can wrap the activity matrix into the ExpressionSet class object, just like expression matrix.

```R
# Create eset using activity matrix
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')
```

**Fourth, create QC report for the activity matrix.**
Use `draw.eset.QC()` to create the HTML QC report [QC for AC](AC_QC.html) of `analysis.par$merge.ac.eset`. 

```R
# QC plot for activity eset
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_',
             pre_define=c('WNT'='blue','SHH'='red','G4'='green'),pca_plot_type='2D.interactive')
```

**Last, save `analysis.par`.**
For persistent storage of data, and prevent the re-run of all previous steps. Users can checkout now and save `analysis.par` as RData for this part.

```R
# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')
```

----------
### *Why study driver's activity ?*

- Drivers, such as transcription factors (TF) bind to the enhancer or promoter regions of the target genes and regulate their expression level. Once the TF has been synthesized, there are still many steps between mRNA translation of a TF and the actual transcriptional regulation of target genes. The activity is controlled by most of the following processes, nuclear localization (import into the cell nucleaus), activation through signal-sensing domain (e.g. ligand binding or post-translational modifications: methylation, ubiquitination and phosphorylation, which is essential for dimerization and promoter binding), access to the DNA-binding site (epigenetic features of the genome), and interaction with other cofactors or TFs (form a complex). Thus, sometimes the expression trend and the activity pattern of a driver may be contradicted. Due to these "hidden effects, the analysis of the activity of a driver maybe more fruitful than the analysis its expression level.

----------

## Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA)
**Purpose: get significantly differential expression (DE) and activity (DA) for all possible drivers between two phenotype groups.**

Please skip the following line if you didn't close R session after completed Step 1 and Step 2.

Don't skip, if you have checked out and closed R session after completed the Step 1 and Step 2. Before start Step 3, please reload `analysis.par` RData from Step 2.
`NetBID.loadRData()` reloads RData saved by `NetBID.saveRData()`. It prevents user from repeating former pipeline steps.
If the re-opened R session doesn't have `analysis.par` in the environment, please comment off the first two command lines. It will create a temporary `analysis.par`
with path of the saved Step 2 RData, `analysis.par$out.dir.DATA`. The path `test//driver_2019-05-06//DATA/` here is just an example, 
users need to give their own path used to save `analysis.par` RData from Step 2.

```R
# Reload network.par RData from Step 2
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')
```

**To compare the drivers' activity between two phenotype groups.**
In the demo dataset, one column of the phenotype data frame is "subgroup". It contains 3 phenotype groups, `WNT`, `SHH` and `G4`. To compare the driver's DE/DA between any two groups (`G1` vs. `G0` in the following functions), NetBID2 provides two major functions `getDE.BID.2G()` and `getDE.limma.2G()`. Users need to assign sample names from each group to `G1_name` and `G0_name`.
In the script below, we compared between `G4.Vs.WNT` and `G4.Vs.SHH`. Each comparison has a name (highly suggested, will be displayed in the final master table), and the results are saved in `analysis.par$DE` and `analysis.par$DA`. 

**More detail.**
`getDE.BID.2G()` uses Bayesian Inference method to calculate DE and DA values, by setting `method='Bayesian'`. If set `method='MLE'`, it will take shorter time to calculate.
If the input is RNA-Seq dataset, user can use `DE` output from `DESeq2`. Just pay attention to the column names, one column should be `Z-statistics`.
If the phenotype is ordinal or more complicated, user can use `bid()` or `limma()` to analyze. 

```R
# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# First comparison: G4 vs. WNT
comp_name <- 'G4.Vs.WNT' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='WNT')] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# Second comparison: G4 vs. SHH
comp_name <- 'G4.Vs.SHH' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='SHH')] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

```

**To combine multiple comparison results.** Wrap all the comparison results into the `DE_list`, and pass it to `combineDE()`.
It returns a list contains the combined DE/DA analysis. A data frame named "combine" inside the list is the combined analysis. Rows are genes/drivers, columns are combined statistics (e.g. "logFC", "AveExpr", "t", "P.Value" etc.).

```R
## Third comparison: G4 vs. others
# Combine the comparison results from `G4.Vs.WNT` and `G4.Vs.SHH`
comp_name <- 'G4.Vs.otherTwo' # Each comparison must has a name
DE_gene_comb <- combineDE(DE_list=list('G4.Vs.WNT'=analysis.par$DE$`G4.Vs.WNT`,'G4.Vs.SHH'=analysis.par$DE$`G4.Vs.SHH`))
DA_driver_comb <- combineDE(DE_list=list('G4.Vs.WNT'=analysis.par$DA$`G4.Vs.WNT`,'G4.Vs.SHH'=analysis.par$DA$`G4.Vs.SHH`))
analysis.par$DE[[comp_name]] <- DE_gene_comb$combine
analysis.par$DA[[comp_name]] <- DA_driver_comb$combine
```
NetBID2 also provides `draw.combineDE()` to visualize the top drivers with significant DE/DA from the combined comparison, with DE/DA values from seperate comparisons listed as well.
Drivers are sorted based on the P-values from the third column, which is the combined comparison (`G4 vs. WNT+SHH`). 

- Driver table of top DE:

```R
# Driver table of top DE
draw.combineDE(DE_gene_comb)
draw.combineDE(DE_gene_comb,pdf_file=sprintf('%s/combineDE.pdf',analysis.par$out.dir.PLOT)) # Save it as PDF
```

![`combineDE`](combineDE.png)

- Driver table of top DA:

```R
# Driver table of top DA
draw.combineDE(DA_driver_comb)
draw.combineDE(DA_driver_comb,pdf_file=sprintf('%s/combineDA.pdf',analysis.par$out.dir.PLOT)) # Save it as PDF
```

![`combineDA`](combineDA.png)

Another way to perform the comparison between one group versus multiple groups, is to choose the sample names from multiple groups as `G0_name` parameter in the `getDE.BID.2G()`. 
In this way, the user doesn't need to call `combineDE`. But the final result may vary, due to the different statistical hypothesis.

```R
# Another way to do the third comparison: G4 vs. others
comp_name <- 'G4.Vs.others'
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # Combine other groups as the Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid
```

Save the comparison results as RData.

```R
# Save Step 3 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')
```

Now, we have get the top differential expression (DE) and activity (DA) drivers from different comparisons. We can use `draw.NetBID()` to visualize the top drivers (top 30 by default). Here, the displayed statistics for NetBID is set at `P.Value` and `logFC` for differential expression (set by `DA_display_col` and `DE_display_col`). 

```R
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others')
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others',pdf_file=sprintf('%s/NetBID_TOP.pdf',analysis.par$out.dir.PLOT),text_cex=0.8) # Save as PDF
```

![`NetBID_TOP`](NetBID_TOP.png)

Users can customize the table above, by choosing which column to display DE or DA, or by choosing which comparison to display. 

## Step 4: Generate a master table for drivers (ms-tab)
**Purpose: gather previous calculated data to create a final master table for all possible drivers.**

Please skip the following line if you didn't close R session after completed Step 1-3.

Don't skip, if you have checked out and closed R session after completed the Step 1-3. Before start Step 4, please reload `analysis.par` RData from Step 3.
`NetBID.loadRData()` reloads RData saved by `NetBID.saveRData()`. It prevents user from repeating former pipeline steps.
If the re-opened R session doesn't have `analysis.par` in the environment, please comment off the first two command lines. It will create a temporary `analysis.par`
with path of the saved Step 3 RData, `analysis.par$out.dir.DATA`. The path `test//driver_2019-05-06//DATA/` here is just an example, 
users need to give their own path used to save `analysis.par` RData from Step 3.

```R
# Reload analysis.par RData from Step 3
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')
```

**To create the final master table, users need to gather certian previous data and pass them to the parameters in `generate.masterTable()` function.**
Here are the details,
- `use_comp`, the vector of multiple comparison names. It will be used as the columns of master table. If use all the comparisons, just call `names(analysis.par$DE)`.
- `DE` and `DA`, if following the previous pipeline, just use `analysis.par$DE` and `analysis.par$DA`.
- `network`, the driver-to-target list. The names of the list elements are drivers. Each element is a data frame, usually contains three columns. "target", target gene names; "MI", mutual information; "spearman", spearman correlation coefficient. If following the previous pipeline, just use `analysis.par$merge.network$target_list`.
- `tf_sigs`, contains all the detailed information of TF and Sig. Users need to call `db.preload()` for access.
  - `db.preload()` reloads the TF/SIG gene lists into R workspace and saves it locally under db/ directory with specified species name and analysis level.
- `z_col`, name of the column in `getDE.limma.2G()` and `getDE.BID.2G()` output data frame contains the Z-statistics. By default, it is "Z-statistics".
- `display_col`, other driver's statistical values for display. This must columns from the `getDE.limma.2G()` and `getDE.BID.2G()` output data frame. For example, "logFC" and "P.Value".
- `main_id_type`, the type of driver’s ID, **IMPORTANT**. It comes from the attribute name in `biomaRt` package. Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna". For details, please heck [ID conversion section](../network_construction#id-conversion).
- `transfer_tab`, the data frame for ID conversion. If NULL and `main_id_type` is not in the column names of `tf_sigs`, it will use the conversion table within the function. Users can also create their own ID conversion table. If user want to save it into `analysis.par`, we suggest to call `get_IDtransfer2symbol2type()` instead of `get_IDtransfer()` as the output of the former function could be used in more visualization functions (e.g `draw.bubblePlot()`).
- `column_order_stratey`, an option to order to the columns in the mater table. If set as `type`, the columns will be ordered by column type; If set as `comp`, the columns will be ordered by comparison. 

```R
# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones
# Prepare the conversion table (OPTIONAL)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes)
analysis.par$transfer_tab <- transfer_tab
# Creat the final master table
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                               target_list=analysis.par$merge.network$target_list,
                                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                               main_id_type='external_gene_name')
```

**Please note**, there may exist some drivers with only activity values but no expression values. This is due to the fact that the network-construction dataset is different from the analysis dataset.

**Save the final master table as EXCEL file.** 
NetBID2 provides `out2excel()` to save the master table as EXCEL, with multiple options to highlight interested data (e.g. marker genes, drivers with significant Z-values). 
For more options, please check `?out2excel()`. **Download the master table EXCEL file [ms_tab.xlsx](driver_ms_tab.xlsx) here to see.**

```R
# Path and file name of the output EXCEL file
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# Highlight marker genes
mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),
                  SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
                  G4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
# Customize highlight color codes
#mark_col <- get.class.color(names(mark_gene)) # this randomly assign color codes
mark_col <- list(G4='green','WNT'='blue','SHH'='red')
# Save the final master table as EXCEL file
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col)
```

We've finished the driver estimation part of NetBID2. We need to save the `analysis.par` as RData, cause it contains all the important results of driver estimation.
This is **ESSENTIAL** to run [Advanced analysis](../docs/advanced_analysis) part of NetBID2 and [NetBID2 shiny server](https://github.com/jyyulab/NetBID_shiny). 
The `analysis.par` list includes 13 elements (`main.dir`, `project.name`, `out.dir`, `out.dir.QC`, `out.dir.DATA`, `out.dir.PLOT`, `merge.network`, `cal.eset`, `merge.ac.eset`, `DE`, `DA`, `final_ms_tab`, `transfer_tab`). They will be used to run NetBID shiny server. 

```R
# Save Step 4 analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')
```

-------
### *How to interpret and use the master table ?*

- `out2excel()` can save multiple master tables as excel sheets in one EXCEL file.
- For one master table, it consists of three parts:
- For one master table, it consists of three parts:
  - The first six columns are `gene_label`, `geneSymbol`, `originalID`, `originalID_label`, `funcType` and `Size`.
    - `gene_label` is the driver's gene symbol or transcript symbol, with suffix "_TF" or "_SIG" to show driver's type. 
    - `geneSymbol` is the driver's gene symbol or transcript symbol, without suffix.
    - `originalID` is the original ID type used in network construction, which should match the ID type in `analysis.par$cal.eset`, `analysis.par$DE`.
    - `originalID_label` is the original ID type with suffix "_TF" or "_SIG", which should match the the ID type in `analysis.par$merge.network`, `analysis.par$merge.ac.eset`,`analysis.par$DA`.
    - **`originalID_label`** is the only column to ensure unique ID for row record.
    - `funcType` is either "TF" or "SIG" to mark driver's type. 
    - `Size` is number of target genes for the driver. 
  - The statistical columns are named as `prefix.comp_name_{DA or DE}`. The `prefix` can be `Z`, `P.Value`, `logFC` or `AveExpr` to indicate which statistical value is stored. The `comp_name` is the comparison name. For example, `Z.G4.Vs.WNT_DA` means the Z-statistics of the differential activity (DA) calculated from comparison between phenotype G4 and phenotype WNT. The color shade of the background indicated the significance of Z-statistics.
  - The next 13 columns (from `ensembl_gene_id` to `refseq_mrna`) are detailed information of genes.
  - The last columns (optional) are the detailed information of marker genes, users use `mark_strategy='add_column'` to set.
- Users can fileter drivers by target size, or sort the Z-statistics to get top significant drivers. NetBID2 also provides [Advanced analysis](../docs/advanced_analysis) for advanced
analysis and visualization.

-------


