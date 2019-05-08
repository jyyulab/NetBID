---
layout: default
title: Advanced analysis
nav_order: 4
has_children: true
permalink: /docs/advanced_analysis
---

# Advanced analysis

The purpose for this part: 

**further analyze and visualize the results for drivers**.

The full demo script for this part could be found in [analysis_and_plot_demo1.R](https://github.com/jyyulab/NetBID-dev/blob/master/demo_scripts/analysis_and_plot_demo1.R).

In this part, we will use the question-induced strategy to assist user for the choice and usage of those visualization functions. 
The main purpose is that we want to **find potential hidden drivers in Group4 compared with the other subtypes by using NetBID2**, for which may be related with the specific clinical feature for Group4 MB.

## Preparations

Library the installed NetBID2, set the directory to get the pre-saved RData in the `ms-tab` step. Here, we could use the demo in NetBID2 package.
Print the `analysis.par$out.dir.PLOT` to get know where the plot figures will be saved. User could also choose to change to another local writable path.

```R
# library the package
library(NetBID2)
# set the directory to get the pre-saved RData
analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
print(analysis.par$out.dir.PLOT)
#analysis.par$out.dir.PLOT <- 'test/driver/PLOT' ## directory to save the plot figures
```

## Part I: for the top list of significant drivers 

### QI.1: How to get the top list of drivers with significantly different activity (DA) in G4 Vs. other subtypes ?

For simplicity, save the final master table to another variable `ms_tab` and set the comparison name `comp_name`. 
Here we choose the `G4.Vs.others` as the demo.
Before detailed analysis, the target size for the driver could be used as a filteration, drivers with too small target size (e.g <30) or too large target size (e.g >1000) will not be used.

```R
ms_tab <- analysis.par$final_ms_tab ## get the master table data frame
ms_tab <- ms_tab[which(ms_tab$Size>=30 & ms_tab$Size <=1000),] ## target size filteration
comp_name <- 'G4.Vs.others' ## get the comparison name
```

We will use `draw.volcanoPlot()` to get the list of top DA drivers and plot the results. 
The input is the data frame containing the columns of `label_col`, `logFC_col` and `Pv_col`. 
User need to input the threshold for logFC (`logFC_thre`) and P.Value (`Pv_thre`).
This function could choose not to display figures by setting `show_plot=FALSE` and only output the list of significant top drivers.
If choose `show_label=TRUE`, the number of the drivers passed the threshold may not be very high as too many labels displayed on the figure will become a mess.

For the activity value, the threshold for logFC may not be very high especially setting `std=TRUE` in `cal.Activity()`. 

- Top drivers with DA passed absolute logFC&ge;0.4 and P-value&le;1e-8, with detailed label information are shown:

```R
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                               Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.4,Pv_thre=1e-8,
                               main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                               pdf_file=sprintf('%s/vocalno_label_DA.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
```                               

![vocalno_label_DA](vocalno_label_DA.pdf)

- Volcano plot is also applicable for get the top DE genes. Top genes with DE passed absolute logFC&ge;1.5 and P-value&le;1e-4 are shown:

```R
sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                             Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-4,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                             pdf_file=sprintf('%s/vocalno_label_DE.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
```

![vocalno_label_DE](vocalno_label_DE.pdf)

`sig_driver` and `sig_gene` are the data frame containing drivers and genes passed the threshold with detailed logFC and P.value statistics.

### QI.2: How to understand the significance of those top DA drivers ?

Those drivers are significant due to their activity difference in Group4 Vs. other subtypes. 
The activity is calculated based on the expression pattern of their target genes.
If the driver's activity is significantly up-regulated in Group4, which means that its positively regulated target genes are significantly up-regulated in Group4 and negative target genes are down-regulated.
Based on this, we could use `draw.GSEA.NetBID()` to visualize this trend for top drivers. This function has more than 10 options, check the manual for detail by `?draw.GSEA.NetBID`.

First, the differentiated expression (DE) profile in Group4 Vs. others should be used to estimate the performance of the target genes. 
And `driver_list` need to be the `originalID_label` in the master table as it is the only unique column. In the above volcano plot, the rownames for the output data frame will be the rownames from the original master table, i.e the `originalID_label`. 

```R
DE <- analysis.par$DE[[comp_name]]
driver_list <- rownames(sig_driver) ## the rownames is the originalID_label
```

The following figures differ in the option of `profile_trend`, `profile_col`, `target_nrow`, `target_col`, `target_col_type`. 

- In the first demo plot, the profile (above part of the figure) is from negative values to positive values (`profile_trend='neg2pos'`), indicating the `logFC` (set in `profile_col`) value for genes.
The position for the genes are consistent in the profile figure and the GSEA-pattern figure (below part of the figure). 
Each vertical line in each row stands for one target gene for the corresponding driver. 
The target genes for one driver will be separated into positive-regulated sets and negative-regulated sets (`target_nrow=2`).
The color for those target genes are the significance in the profile (`target_col_type='DE'`), with value passed the `profile_sig_thre` color in blue (negative value) and red (positive value).
We could study from the figure that, the positively regulated target genes for the *down* DA drivers tend to lay in the left part of the profile and the negatively regulated in the right part. Contrary trend for *up* DA drivers. 

```R
draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',name_col='ID',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in others',right_annotation = 'high in G4',
                 main=comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo1.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo1](NetBID_GSEA_demo1.pdf)

- In the second demo plot,the profile is from positive values to negative values (`profile_trend='pos2neg'`), indicating the `t` (set in `profile_col`) value for genes.
The color for those target genes are the same for positively-regulaion or negatively-regulation (`target_col_type='PN'`)

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main=comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo2.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo2](NetBID_GSEA_demo2.pdf)

- In the third demo plot, the color for those target genes are just black (`target_col='black'`).

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='black',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main=comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo3.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo3](NetBID_GSEA_demo3.pdf)

- In the forth demo plot, the positively-regulaion or negatively-regulation is not separately into two rows (`target_nrow = 1`), `target_col='RdBu'` and `target_col_type='DE'`.

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='RdBu',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main=comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo4.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo4](NetBID_GSEA_demo4.pdf)

- In the fifth demo plot, the positively-regulaion or negatively-regulation is not separately into two rows (`target_nrow = 1`) and `target_col='black'`.

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='black',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main=comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo5.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo5](NetBID_GSEA_demo5.pdf)

### QI.3: What is the expression/activity pattern of those top DA drivers in samples with different subtype ?

Next, we'd like to know the expression/activity pattern for those top DA drivers in different group of samples, the most straightforward strategy is to use heatmap.
In NetBID2, `draw.heatmap()` is designed to assist simple usage of heatmap. It is based on `Heatmap()` function in [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html). 

Get the expression matrix `exp_mat` and activity matrix `ac_mat` from `analysis.par`, and phenotype data frame `phe_info`.

```R
exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames must be the originalID
ac_mat <- exprs(analysis.par$merge.ac.eset) ## ac,the rownames must be the originalID_label
phe_info <- pData(analysis.par$cal.eset)
# could use addtional paramters in Heatmap()
```

- Draw the expression value (scaled by samples `scale='row'`)  for top DA drivers across all samples and display the sample class label above. 
Here, the `use_genes` must be the `originalID`.
User could directly input `phenotype_info=phe_info` and choose which columns to display in `use_phe`.

```R
draw.heatmap(mat=exp_mat,use_genes=ms_tab[driver_list,'originalID'],use_gene_label=ms_tab[driver_list,'geneSymbol'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo1.pdf',analysis.par$out.dir.PLOT))
```

![heatmap_demo1](heatmap_demo1.pdf)

- Draw the activity value for top DA drivers, here, the `use_genes` must be the `originalID_label` but user could to display label in the `gene_label` column. 
In this demo, the original label is also gene symbol, but in some cases the original label maybe the ensemble_gene_id, the originalID will be different and user could choose to display the gene symbol in `gene_label` (set by `use_gene_label`).

```R
draw.heatmap(mat=ac_mat,use_genes=driver_list,use_gene_label=ms_tab[driver_list,'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo2.pdf',analysis.par$out.dir.PLOT))
```

![heatmap_demo2](heatmap_demo2.pdf)

Compared with the two figures above, the activity pattern is much more clean than the expression pattern for top DA drivers, and also for top DE drivers (figures not shown here, but demo script is below), suggesting that the activity of drivers may be more robust in separating specific sample groups. 

```R
# try for the top DE genes             
gene_list <- rownames(sig_gene) 
draw.heatmap(mat=exp_mat,use_genes=ms_tab[gene_list,'originalID'],use_gene_label=ms_tab[gene_list,'geneSymbol'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo3.pdf',analysis.par$out.dir.PLOT))

draw.heatmap(mat=ac_mat,use_genes= gene_list,use_gene_label=ms_tab[gene_list,'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo4.pdf',analysis.par$out.dir.PLOT))
```

### QI.4: What is the biological function of those top DA drivers ?

Function enrichment plot for top drivers

```R

res_up <- funcEnrich.Fisher(input_list=ms_tab[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP'),
                           Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
                           
res_down <- funcEnrich.Fisher(input_list=ms_tab[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP'),
                           Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       
                           
# draw barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers',pdf_file=sprintf('%s/funcEnrich_bar_nogene.pdf',analysis.par$out.dir.PLOT))
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_withgene.pdf',analysis.par$out.dir.PLOT))

draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.75)
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGS.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = FALSE)
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGENE.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = TRUE)
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 1.5,gene_cex=1.4,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterNO.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = FALSE)
```


gs.preload()
merge_gs()
funcEnrich.Fisher()
draw.funcEnrich.bar()
draw.funcEnrich.cluster()


### QI.5: What is the biological function of the target genes of those top DA drivers ?
Bubble plot for top drivers

## Part II: for a selected interested driver

### QII.1: How to understand the significance of the selected driver ?
GSEA plot for each driver

### QII.2: How to visualize the network structure of the selected driver ?
target network structure for each driver, or two drivers (overlap testing)

### QII.3: What is the expression/activity for this selected driver in samples with different subtypes ?
category plot

## Part III: further analysis

### QIII.1: What are the activities of the curated gene sets in all samples and what are the top significantly differed gene sets ?
gene set-based activity analysis, including vocalno, heatmap, category and GSEA plot

### QIII.2: How to find drivers with significantly overlapped target genes ?
bubble plot for target gene list


-------
### *How to modify the figures by adjusting the paramters in the draw.* functions ?*


-------
