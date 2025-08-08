---
layout: default
title:  Installation
nav_order: 2
permalink: /docs/installation
---

## Installation

### 1. Install from github main branch

```R
# install.packages("devtools")
devtools::install_github("jyyulab/NetBID") 
```

### 2. Install a released R package

Download a released version from [https://github.com/jyyulab/NetBID/releases](https://github.com/jyyulab/NetBID/releases) and run:

```R
# install.packages("devtools")
devtools::install_local('NetBID2_2.2.0.tar.gz')
```
NOTE:
1) For MAC users, please make sure the XCode and Xquartz have been installed properly before you start installing NetBID2.
2) Some dependencies might not be available in your situation. If so, please run the commands below to install them:

```R
# install.packages("BiocManager")
BiocManager::install(c('SummarizedExperiment', 'Biobase', 'GEOquery', 'limma', 'impute', 'tximport', 'DESeq2', 'ConsensusClusterPlus', 'ComplexHeatmap', 'biomaRt', 'GSVA', 'rhdf5')) # Modify the list per your requirements
```

---
