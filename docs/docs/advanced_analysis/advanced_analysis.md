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

Top drivers with DA passed absolute logFC&ge;0.4 and P-value&le;1e-8, with detailed label information is shown:

```R
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                               Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.4,Pv_thre=1e-8,
                               main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                               pdf_file=sprintf('%s/vocalno_label_DA.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
```                               

![vocalno_label_DA](vocalno_label_DA.pdf)

Top drivers with DE passed absolute logFC&ge;1.5 and P-value&le;1e-4 is shown:

```R
sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                             Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-4,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                             pdf_file=sprintf('%s/vocalno_label_DE.pdf',analysis.par$out.dir.PLOT))
```

![vocalno_label_DE](vocalno_label_DE.pdf)

`sig_driver` and `sig_gene` are the data frame containing drivers and genes passed the threshold with detailed logFC and P.value statistics.

### QI.2: How to understand the significance of those top DA drivers ?

Those drivers are significant due to their activity difference in Group4 Vs. other subtypes. 
The activity is calculated based on the expression pattern of their target genes.
If the driver's activity is significantly up-regulated in Group4, which means that its positively regulated target genes are significantly up-regulated in Group4 and negative target genes are down-regulated.
Based on this, we could use `draw.GSEA.NetBID()` to visualize this trend for top drivers.

First, the differentiated expression (DE) profile in Group4 Vs. others should be used to estimate the performance of the target genes. And `driver_list` is used to extract out the driver's name from above.

```R
DE <- analysis.par$DE[[comp_name]]
driver_list <- rownames(sig_driver)
```

```R
draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',name_col='ID',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in others',right_annotation = 'high in G4',
                 main='test',target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo1.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo1](NetBID_GSEA_demo1.pdf)

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main='test',target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo2.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo2](NetBID_GSEA_demo2.pdf)

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='black',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main='test',target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo3.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo3](NetBID_GSEA_demo3.pdf)

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='RdBu',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main='test',target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo4.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo4](NetBID_GSEA_demo4.pdf)

```R
draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='black',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main='test',target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo5.pdf',analysis.par$out.dir.PLOT))
```

![NetBID_GSEA_demo5](NetBID_GSEA_demo5.pdf)


### QI.3: What is the expression/activity pattern of those top DA drivers ?
Heatmap for top drivers


### QI.4: What is the biological function of those top DA drivers ?
Function enrichment plot for top drivers


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


