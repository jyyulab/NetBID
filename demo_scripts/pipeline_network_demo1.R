## This script aims to generate network ##
############# preparations ###############

library(NetBID2)

######## global parameters
# currently some steps will be run on server, the fist two parameters is used for changing the path (local to server)
# RP(run parameters)
RP.main_dir      <- '/Volumes/project_space' ## local path
RP.bash.main_dir <- '/research/projects/yu3grp' ## server path
RP.python3.path  <- '/hpcf/apps/python/install/3.6.1/bin/python3.6' ## required
RP.load          <- 'module load python/3.6.1' ## required
RP.MICA.main     <- '/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/scMINER-master/' ## required for MICA
RP.SJAR.main     <- '/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/SJARACNe-master/' ## required for SJAracne
RP.SJAR.pre_cmd  <- sprintf("export SJARACNE_PATH=%s \nexport PYTHON_PATH=%s \n%s",RP.SJAR.main,RP.python3.path,RP.load) ## required

######### preload all knowledge info, choose species, use_level (gene/transcript)
db.preload(use_level='gene',use_spe='human',update=FALSE)

######### write in paramters
project_main_dir <- 'test/' ### user definied !!!!
current_date <- format(Sys.time(), "%Y-%m-%d") ## current date for project running
project_name <- sprintf('project_%s',current_date) ## recommend to add time tag

## network.par is very essential in the analysis!!
if(exists('network.par')==TRUE) rm(network.par)
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,prject_name=project_name)

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

######################################################### Step2: normalization for the exp dataset (exp-QC)
# load from RData
NetBID.loadRData(network.par = network.par,step='exp-load')

## following QC steps are optional !!!!
mat <- exprs(network.par$net.eset)
# remove possible NA? or imputation ? ## need to user-decide
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
mat <- impute.knn(mat)$data
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

######################################################### Step3: check sample cluster info, optional (exp-cluster)
# load from RData
NetBID.loadRData(network.par = network.par,step='exp-QC')
# use most variable genes for cluster
mat <- exprs(network.par$net.eset)
choose1 <- IQR.filter(exp_mat=mat,use_genes=rownames(mat),thre = 0.5,loose_gene=NULL,loose_thre=0.1)
print(table(choose1))
mat <- mat[choose1,]
# update eset
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),], annotation_info=annotation(network.par$net.eset))
network.par$net.eset <- net_eset
# QC plots
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='Cluster_')
# save to RData
NetBID.saveRData(network.par = network.par,step='exp-cluster')

# more cluster functions (will not directly save to file, but avtively layout)
# load from RData
NetBID.loadRData(network.par = network.par,step='exp-cluster')
#
mat <- exprs(network.par$net.eset)
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
draw.clustComp(pred_label,obs_label=get_obs_label(phe,intgroup[i])) ## display the comparison in detail
# perform clustering based on MICA (current run on server !!)
mat <- exprs(network.par$net.eset)
outdir <- file.path(network.par$out.dir,'MICA') ## set the output directory, this should be run on server
prj.name <- network.par$project.name
SJ.MICA.prepare(mat,outdir=outdir,prjname=prj.name,all_k=2:6,retransformation="False",perplexity=5)
# then run the two bash files on SJ server by the output instructions
# when finished, plot functions:
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.MICA(outdir=outdir,prjname=prj.name,all_k = 2:6,obs_label=get_obs_label(phe,intgroup[i]))
  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
}

######################################################### Step4: prepare SJARACNE (sjaracne-prep)
# load from RData
NetBID.loadRData(network.par = network.par,step='exp-QC') ## do not load file from exp-cluster
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

#### Old version for internal use in SJ
result_info <- SJ.SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                                   TF_list=use_list$tf,SIG_list=use_list$sig,
                                   IQR.thre = 0.5,IQR.loose_thre = 0.1,
                                   SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR,mem=10240) ## memory shoud be definied !!!
# check result_info$bash.tf result_info$bash.sig; and run on cluster !!!
# or auto generate bash files for Step1, Step2, Step3, Step4 for all bash files under one directory
SJ.SJAracne.step(network.par$out.dir.SJAR)
# !!! then run four steps one by one

###################################### finish network generation part !!! Cheers !!! #########################################

