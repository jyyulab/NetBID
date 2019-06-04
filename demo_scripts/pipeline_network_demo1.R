## This script aims to generate network ##

############# Step0: preparations ###############
######### library the package
library(NetBID2)

######### set up paramters
project_main_dir <- 'test/' ### user defined main directory for the project, one main directory could have multiple projects, separeted by project name
current_date <- format(Sys.time(), "%Y-%m-%d") ## current date for project running, suggested to add the time tag
project_name <- sprintf('project_%s',current_date) ## project name for the project

## network.par is very essential in the analysis, if the environment already has this parameter, strongly suggest to delete it first
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name) ## create working directory

########### !!!!!!!!!! pipelines !!!!!!!!!! ############

######################################################### Step1: load in gene expression datasets for network construction (exp-load)
## from GEO, need to provide GSE and GPL
net_eset <- load.exp.GEO(out.dir=network.par$out.dir.DATA,GSE='GSE116028',GPL='GPL6480',getGPL=TRUE,update=FALSE)

# ID conversion, or merge transcript level to expression level, use_feature_info can be other dataframe info; optional;
net_eset <- update_eset.feature(use_eset=net_eset,use_feature_info=fData(net_eset),from_feature='ID',to_feature='GENE_SYMBOL',merge_method='median') ### !!! need modify
# select phenotype columns or user add phenotype info; optional
net_eset <- update_eset.phenotype(use_eset=net_eset,use_phenotype_info=pData(net_eset),use_sample_col='geo_accession',use_col='GEO-auto')

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset

# QC for the raw expdatasets
# if intgroup==NULL, will auto get intgroup
# prefix: prefix for the pdf files
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='beforeQC_')

# save to RData
NetBID.saveRData(network.par = network.par,step='exp-load')

######################################################### Step2: normalization for the expression dataset (exp-QC)
# load from RData
NetBID.loadRData(network.par = network.par,step='exp-load')

## following QC steps are optional !!!!
mat <- exprs(network.par$net.eset)
# remove possible NA? or imputation ? ## need to user-decide
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data
# log2 transformation
med_val <- median(apply(mat,2,median));print(med_val)
if(med_val>16){mat <- log2(mat)}
# quantile normalization
mat <- normalizeQuantiles(mat) ## limma quantile normalization
# remove low expressed genes, such as whether in more than 90% samples, the expression value is lower than 5%
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]
# update eset
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),],
                          annotation_info=annotation(network.par$net.eset))
network.par$net.eset <- net_eset
# QC for the QC expdatasets
# if intgroup==NULL, will auto get intgroup
# prefix: prefix for the pdf files
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='afterQC_')

# save to RData
NetBID.saveRData(network.par = network.par,step='exp-QC')

######################################################### Step3: check sample cluster information, optional (exp-cluster)
# load from RData
NetBID.loadRData(network.par = network.par,step='exp-QC')
# use most variable genes for cluster
mat <- exprs(network.par$net.eset)
choose1 <- IQR.filter(exp_mat=mat,use_genes=rownames(mat),thre = 0.5,loose_gene=NULL,loose_thre=0.1)
print(table(choose1))
mat <- mat[choose1,]
# generate tmp eset
tmp_net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),], annotation_info=annotation(network.par$net.eset))
# QC plots
draw.eset.QC(tmp_net_eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='Cluster_')

# more cluster functions (will not directly save to file, but actively layout)
phe <- pData(network.par$net.eset)
intgroup <- get_int_group(network.par$net.eset)
# pca+kmeans in 2D
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]))
  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
}
# pca+keans in 3D
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]),plot_type='3D')
  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
}
use_int <- 'subgroup'
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D')
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.ellipse')
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.text')
pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='3D')
print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe, use_int))))
draw.clustComp(pred_label,obs_label=get_obs_label(phe,use_int),outlier_cex=1,low_K=10) ## display the comparison in detail


######################################################### Step4: prepare SJARACNE (sjaracne-prep)
# load from RData
NetBID.loadRData(network.par = network.par,step='exp-QC') ## do not load file from exp-cluster

# load database
db.preload(use_level='gene',use_spe='human',update=FALSE)

# ID convertion, get TF/SIG list !!!!
use_gene_type <- 'external_gene_name' ## this should user-defined !!!
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
##
# select sample for analysis
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples
prj.name <- network.par$project.name # can use other names, if need to run different use samples

#### New version
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                    TF_list=use_list$tf,SIG_list=use_list$sig,
                    IQR.thre = 0.5,IQR.loose_thre = 0.1,
                    SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)

###################################### finish network construction part !!! Cheers !!! #########################################

