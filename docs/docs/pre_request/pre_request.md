---
layout: default
title: "- Dependency"
nav_order: 3
permalink: /docs/pre_request
---

## Pre-requested R packages

**R version >= 3.6.0**

The following is a description of the dependent R packages, mainly to help developers better understand the code in NetBID2.

The users does not need to digest the information below. 


|R package|Functions[^1]|Category|Purpose|
|---------|---------|--------|-------|
|Biobase|-|Data processing|ExpressionSet class|
|GEOquery|getGEO|Data processing|Get expression dataset from GEO database|
|limma|-|Data processing|Expression data normalization|
|impute|impute.knn|Data processing|Data imputation|
|tximport|Data processing|Data import|
|DESeq2|DESeqDataSetFromTximport,DESeq|Data processing|Data import from txi and normalization for RNASeq data|
|ConsensusClusterPlus|-|Clustering|Get consensus clustering results|
|aricode|clustComp|Clustering|For cluster comparison statistic calculation|
|igraph|-|Visualization|Igraph class and basic network-based calculation|
|RColorBrewer|brewer.pal|Visualization|Get color bar|
|plot3D|scatter3D|Visualization|3D plot|
|plotrix|draw.ellipse|Visualization|Drawing ellipse|
|umap|umap|Visualization|Data dimension reduction and visualization|
|ComplexHeatmap|Heatmap|Visualization|Heatmap drawing|
|plotly|-|Visualization|For interactive plot|
|ordinal|clm,clmm|bid|Cumulative link mixed models|
|MCMCglmm|MCMCglmm|bid|Multivariate generalized linear mixed models|
|arm|bayesglm|bid|Bayesian generalized linear models|
|reshape|melt|bid|Melt an object into a form suitable for easy casting|
|biomaRt|-|ID conversion|ID conversion|
|GSVA|gsva|Gene Set|Gene set activity calculation|
|msigdbr|-|Gene Set|MSigDB database|
|openxlsx|-|IO|Output into EXCEL file (master table)|
|rhdf5|H5Fopen,H5Fclose|IO|For hd5 formation data processing (MICA)|
|rmarkdown|render|Report|For generating HTML report file|
|kableExtra|-|Report|For table layout in the report file|


[^1]:'-' in the Function column represents multiple functions in the package were used in NetBID2.

