# NetBID2
Network-based Bayesian Inference of Drivers, Version II

# Version 0.1.2 Update Notes:

1. Create two "lazy mode" functions, NetBID.lazyMode.DriverVisualization(),NetBID.lazyMode.DriverEstimation()

2. Add warning message for all functions

3. Let color code set by user-defined, with two options: use_color,pre_define passed to get.class.color(). 
The modified functions are: 
draw.2D(), draw.2D.interactive(), draw.2D.text(), draw.3D(), draw.2D.ellipse(), draw.eset.QC(), draw.pca.keans(), draw.umap.kmeans(), draw.heatmap(), draw.categoryValue(), 

4. Add package information in calling functions to avoid possible conflict of function names

5. Re-write the cal.Activity(), add in Matrix Cross Products, which will accelerate calculation time but it is memory consuming. memory_constrain option could be set.

6. Modify the draw.eset.QC(), add correlation plot. Modify draw.network.QC(), add html_info_limit. 

7. Add draw.2D.interactive(), and add option "2D.interactive" to draw.pca.kmeans(), draw.umap.kmeans()

8. Add functions to judge abnormal values.

9. modify option name from network --> target_list in generate.masterTable()

10. add option geneSymbol_column for SJAracne.prepare()

11. rebuild on R 3.6.0

# Install

## remote install (not available yet)

```R
library(devtools)
library(BiocManager)
# set repos, for R version 3.6.0, Bioconductor version 3.9
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.rstudio.com/"
  r["BioCsoft"] <- "https://bioconductor.org/packages/3.9/bioc"
  r["BioCann"] <- "https://bioconductor.org/packages/3.9/data/annotation"
  r["BioCexp"] <- "https://bioconductor.org/packages/3.9/data/experiment"
  options(repos = r)
})

devtools::install_github("jyyulab/NetBID-dev",ref='master',dependencies='Depends') 
```

or download the release version from https://github.com/jyyulab/NetBID-dev/releases/download/0.1.2/NetBID2_0.1.2.tar.gz

## local install

pull the repos from github and install locally:

```R
devtools::install(pkg='.',dependencies=TRUE) ## Install the package with dependencies.
devtools::install_deps(pkg = ".", dependencies = TRUE) ## Install package dependencies if needed.
```

download the directory to your workspace and then run:

```R
devtools::install_local('NetBID2_0.1.2.tar.gz') ## 
```

# Manual & Tutorial

manual: NetBID2_0.1.2.pdf

tutorial: https://jyyulab.github.io/NetBID-dev/

# Demo
in demo_scripts/ directory

## demo scripts for network generation 
pipeline_network_demo1.R
* Step1: load in gene expression datasets for network construction (exp-load)
* Step2: normalization for the exp dataset (exp-QC)
* Step3: check sample cluster info, optional (exp-cluster)
* Step4: prepare SJARACNE (sjaracne-prep)

## demo scripts for network-based analysis
pipeline_analysis_demo1.R
* Step1: load in gene expression datasets for analysis (exp-load,exp-cluster,exp-QC)
* Step2: activity calculation (act-prep,act-get)
* Step3: get DE/DA (act-DA)
* Step4: generate master table (ms-tab)

## demo scripts for following analysis, mainly focus on visualization
analysis_and_plot_demo1.R

Part I: More details about the top drivers

QI.1: How to get the top drivers with significant differential activity (DA) in the comparison between G4 vs. other subtypes ?

QI.2: How to interpret the significance of top DA drivers ?

QI.3: What is the expression/activity pattern of these top DA drivers across sample subtypes?

QI.4: What are the biological functions of these top DA drivers ?

QI.5: What are the biological functions of the target genes of these top DA drivers ?

Part II: More details about the selected driver

QII.1: How to interpret the significance of the selected driver ?

QII.2: How to visualize the network structure of the selected driver ?

QII.3: What is the expression/activity of this selected driver across subtypes of sample ?

QII.4: What are the functions of the target genes of this selected driver ?

Part III: Other analyses NetBID2 can do

QIII.1: What are the activities of the curated gene sets across all samples ?

QIII.2: How to find drivers share significantly overlapped target genes ?

Q&A: How to modify the figure size created by draw. functions ?

