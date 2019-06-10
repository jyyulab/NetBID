## Complete Pipeline for Network Construction from Transcriptome Dataset ##

############### Step 0: Preparation ###############

# Load NetBID2 package
library(NetBID2)

# Define main working directory and project name
project_main_dir <- './test' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
current_date <- format(Sys.time(), "%Y-%m-%d") # optional, if user like to add current date to name the project folder.
project_name <- sprintf('project_%s',current_date) # project name for the project folders under main directory.

# Create a hierarchcial working directory and return a list contains the hierarchcial working directory information
# This list object (network.par) is an ESSENTIAL variable in network construction pipeline
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

############### Step 1: Load in gene expression datasets for network construction (exp-load) ###############

# Download expression dataset from GEO, need to provide GSE ID and GPL ID
net_eset <- load.exp.GEO(out.dir=network.par$out.dir.DATA,GSE='GSE116028',GPL='GPL6480',getGPL=TRUE,update=FALSE)

# ID conversion, or merge transcript level to expression level, use_feature_info can be other dataframe instead of fData; optional;
net_eset <- update_eset.feature(use_eset=net_eset,use_feature_info=fData(net_eset),from_feature='ID',to_feature='GENE_SYMBOL',merge_method='median')

# Select phenotype columns or user added phenotype info; optional
net_eset <- update_eset.phenotype(use_eset=net_eset,use_phenotype_info=pData(net_eset),use_sample_col='geo_accession',use_col='GEO-auto')

# Add the variable into network.par. ESSENTIAL STEP.
network.par$net.eset <- net_eset

# QC for the raw eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='beforeQC_')

# Save Step 1 network.par as RData
NetBID.saveRData(network.par = network.par,step='exp-load')

############### Step 2: Normalization for the expression dataset (exp-QC) ###############

# Reload network.par RData from Step 1
NetBID.loadRData(network.par = network.par,step='exp-load')

###  The following QC steps are optional
## Firstly, the handling of missing data.
# Get the expression matrix from ExpressionSet object
mat <- exprs(network.par$net.eset)
# Count and show number of NAs across samples and genes
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
# Perform imputation
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data

## Secondly, the log2 transformation.
med_val <- median(apply(mat,2,median)); print(med_val)
if(med_val>16){mat <- log2(mat)}

## Thirdly, the quantile normalization across samples.
# Perform limma quantile normalization
mat <- normalizeQuantiles(mat)

## Fourthly, filter out genes with very low expression values (bottom 5%) in most samples (more than 90%).
# Filter out low-expression genes
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]

# Update eset with normalized expression matrix
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                            feature_info=fData(network.par$net.eset)[rownames(mat),],
                            annotation_info=annotation(network.par$net.eset))
# Updata network.par with new eset
network.par$net.eset <- net_eset

# QC for the normalized eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='afterQC_')

# Save Step 2 network.par as RData
NetBID.saveRData(network.par = network.par,step='exp-QC')


############### Step 3: Check sample cluster analysis, optional (exp-cluster) ###############

# Reload network.par RData from Step 2
NetBID.loadRData(network.par = network.par,step='exp-QC')

# Select the most variable genes across samples
mat <- exprs(network.par$net.eset)
choose1 <- IQR.filter(exp_mat=mat,use_genes=rownames(mat),thre = 0.5)
print(table(choose1))
mat <- mat[choose1,]

# Generate temporary eset
tmp_net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                              feature_info=fData(network.par$net.eset)[rownames(mat),], annotation_info=annotation(network.par$net.eset))
# QC plot for IQR filtered eset
draw.eset.QC(tmp_net_eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='Cluster_')

# The following scripts provide various ways to visualize and check if the IQR filter selected genes
# can be used to perform good sample cluster analysis (predicted labels vs. real labels). Figures will be displayed instead of saving as files.

# Extract phenotype information data frame from eset
phe <- pData(network.par$net.eset)
# Extract all "cluster-meaningful" phenotype columns
intgroup <- get_int_group(network.par$net.eset)
# Cluster analysis using Kmean and plot result using PCA biplot (pca+kmeans in 2D)
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]))
}
# Cluster analysis using Kmeans and plot result using PCA (pca+kmeans in 3D)
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]),plot_type='3D')
  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
}
# Pick one phenotype column "subgroup" from the demo eset and show various plots NetBID2 can create
use_int <- 'subgroup'
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D')
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.ellipse')
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.text')
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='3D')
print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe, use_int))))
draw.clustComp(pred_label,obs_label=get_obs_label(phe,use_int),outlier_cex=1,low_K=10)

############### Step 4: Prepare files to run SJARACNe (sjaracne-prep) ###############

# Reload network.par RData from Step 2
NetBID.loadRData(network.par = network.par, step='exp-QC') ## do not load file from exp-cluster

# Load database
db.preload(use_level='gene',use_spe='human',update=FALSE)

# Converts gene ID into the corresponding TF/SIG list
use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)

# Select samples for analysis
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) # here is using all samples, users can modify
prj.name <- network.par$project.name # if use different samples, need to change the project name
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.1,
                 SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)
