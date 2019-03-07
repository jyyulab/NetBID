## This script aims to get activity and ms table ##
############# preparations ########################

library(NetBID2)

network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#network.dir <- 'inst/demo1/network'
network.project.name <- 'project_2019-02-14' #

######### write in paramters
project_main_dir <- 'test'
project_name <-  'driver'

## analysis.par is very essential in the analysis!!
if(exists('analysis.par')==TRUE) rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, prject_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

########### !!!!!!!!!! pipelines !!!!!!!!!! ############

######################################################### Step1: load in gene expression datasets for analysis (exp-load,exp-cluster,exp-QC)

################## if use same expression datasets as in the network generation, directly do the followings:
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) ## RData from network generation
# can make modifications to the network.par$net.eset and save to analysis.par$cal.eset
analysis.par$cal.eset <- network.par$net.eset
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

################################################################################################ OR:
################## otherwise can follow the similar steps in pipeline_network_Demo1.R for expression dataset QC, cluster
# Step1: exp-load
#cal_eset <- load.exp.GEO(out.dir=network.par$out.dir.DATA,GSE='GSE116028',GPL='GPL6480',getGPL=TRUE,update=FALSE)
#cal_eset <- update_eset.feature(use_eset=cal_eset,use_feature_info=fData(net_eset),from_feature='ID',to_feature='GENE_SYMBOL',merge_method='median') ### !!! need modify
#cal_eset <- update_eset.phenotype(use_eset=cal_eset,use_phenotype_info=pData(net_eset),use_sample_col='geo_accession',use_col='GEO-auto')
#analysis.par$cal.eset <- cal_eset
#draw.eset.QC(analysis.par$cal.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='beforeQC_')
#NetBID.saveRData(analysis.par=analysis.par,step='exp-load')
# Step2: exp-QC
# load from RData
#NetBID.loadRData(analysis.par=analysis.par,step='exp-load')
#mat <- exprs(analysis.par$cal.eset)
#mat <- impute.knn(mat)$data
#med_val <- median(apply(mat,2,median));print(med_val)
#if(med_val>16){mat <- log2(mat)}
#mat <- normalizeQuantiles(mat) ## limma quantile normalization
#choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
#mat <- mat[choose1,]
# update eset
#cal_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(analysis.par$cal.eset)[colnames(mat),],
#                          feature_info=fData(analysis.par$cal.eset)[rownames(mat),], annotation_info=annotation(analysis.par$cal.eset))
#analysis.par$cal.eset <- cal_eset
#draw.eset.QC(analysis.par$cal.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='afterQC_')
#NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')
# Step3: exp-cluster
# load from RData
#NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')
# use most variable genes for cluster
#mat <- exprs(analysis.par$cal.eset)
#choose1 <- IQR.filter(exp_mat=mat,use_genes=rownames(mat),thre = 0.5,loose_gene=NULL,loose_thre=0.1)
#print(table(choose1))
#mat <- mat[choose1,]
## update eset
#cal_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(analysis.par$cal.eset)[colnames(mat),],
#                          feature_info=fData(analysis.par$cal.eset)[rownames(mat),], annotation_info=annotation(analysis.par$cal.eset))
#analysis.par$cal.eset <- cal_eset
## QC plots
#draw.eset.QC(analysis.par$cal.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='Cluster_')
## save to RData
#NetBID.saveRData(analysis.par=analysis.par,step='exp-cluster')
# more cluster functions (will not directly save to file, but actively layout) , will not show demo here, please refer pipeline_network_Demo1.R
################################################################################################

######################################################### Step2: Activity Calculation (act-prep,act-get)
## load the RData for calculation
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')
#### Step2.1
# get network info ! three list(network_dat, target_list, igraph_obj)
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_')
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_')

#### Step2.2
# get activity from network + expression data ### !!! ID issue again, prefer to use ID in networks;
#   ensure more genes in network
# if ID in network and expression data is same, just move to the next step
# otherwise (following are demo scripts):
  # can get transfer tab from biomart
  #use_network_genes  <- union(names(V(analysis.par$tf.network$igraph_obj)),names(V(analysis.par$sig.network$igraph_obj)))
  #db.preload(use_level='transcript',use_spe='human',update=FALSE)
  #network_gene_type  <- 'ensembl_transcript_id' ## this should user-defined !!!
  #analysis_gene_type <- 'refseq_mrna' ## this should user-defined !!!
  #transfer_tab <- get_IDtransfer(from_type=network_gene_type,to_type=analysis_gene_type,use_genes=use_network_genes)
  ## need to check whether the gene_type is correct, or input user-definied gene-id transfer table
  #cal_eset <- analysis.par$cal.eset
  #cal_eset <- update_eset.feature(use_eset=cal_eset,use_feature_info=transfer_tab,from_feature=analysis_gene_type,to_feature=network_gene_type,merge_method='median') ### !!! need modify
  # OR get from GPL, such as the GPL from GPL6480
  #net_eset <- load.exp.GEO(out.dir=sprintf('%s/DATA',network.dir),GSE='GSE28245',GPL='GPL6480',getGPL=TRUE,update=FALSE)
  #transfer_tab <- fData(net_eset)[,c('GB_ACC','ENSEMBL_ID')]
  #cal_eset <- analysis.par$cal.eset
  #cal_eset <- update_eset.feature(use_eset=cal_eset,use_feature_info=transfer_tab,from_feature='GB_ACC',to_feature='ENSEMBL_ID',merge_method='median') ### !!! need modify
  # update eset
  #analysis.par$cal.eset <- cal_eset

# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='act-prep')

#### Step2.3
NetBID.loadRData(analysis.par=analysis.par,step='act-prep')

## You can get activity matrix separately for TF/SIG and merge
#ac_mat <- cal.Activity(target_list=analysis.par$tf.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
#analysis.par$ac.tf.eset  <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
#                                           feature_info=NULL,annotation_info='activity in net-dataset')
#ac_mat <- cal.Activity(target_list=analysis.par$sig.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
#analysis.par$ac.sig.eset  <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
#                                           feature_info=NULL,annotation_info='activity in net-dataset')
# merge tf/sig info (network + ac eset)
#analysis.par$merge.ac.eset <- merge_TF_SIG.AC(TF_AC=analysis.par$ac.tf.eset,SIG_AC=analysis.par$ac.sig.eset)
# QC plots
#draw.eset.QC(analysis.par$ac.tf.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_TF_')
#draw.eset.QC(analysis.par$ac.sig.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_SIG_')

### OR merge network first and get activity matrix (results will be same~)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')
# QC plots
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_')

# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

######################################################### Step3: get DE/DA (act-DA)
## load RData
NetBID.loadRData(analysis.par=analysis.par,step='act-get')
phe_info <- pData(analysis.par$cal.eset)
# generate a list for multiple comparisons
analysis.par$DE <- list()
analysis.par$DA <- list()
# get compared sample names
all_subgroup <- unique(phe_info$subgroup) ## for example need to get subgroup-specific drivers
for(each_subtype in all_subgroup){
  comp_name <- sprintf('%s.Vs.others',each_subtype) ## each comparison must give a name !!!
  G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
  G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
  # get DE, between samples in G1 to samples in G0, random_effect is optional (two strategies to get DE/DA, limma and bid, just choose one)
  #DE_gene_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name=each_subtype,G0_name='other')
  DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name=each_subtype,G0_name='other')
  # get DA
  #DA_driver_limma <- getDE.limma.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name=each_subtype,G0_name='other')
  DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name=each_subtype,G0_name='other')
  # save to analysis.par
  analysis.par$DE[[comp_name]] <- DE_gene_bid
  analysis.par$DA[[comp_name]] <- DA_driver_bid
}
# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

######################################################### Step4: generate master table (ms-tab)
## load RData
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')
# load db
db.preload(use_level='gene',use_spe='human',update=FALSE)
# generate master table
all_comp <- names(analysis.par$DE) ## get all comparison name for output
### attention !!! there may exist some drivers only use AC but no expression level (this is due to platform difference between the network-generation dataset and analysis dataset)
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                               network=analysis.par$merge.network$target_list,
                                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                               main_id_type='external_gene_name')

# output into excel files
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# can add marker gene
mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
                  Group3=c('IMPG2','GABRA5','EGFL11','NRL','MAB21L2','NPR3','MYC'),
                  Group4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
mark_col <- get.class.color(names(mark_gene))
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col)

# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')

###################################### finish basic analysis part !!! Cheers !!! #########################################


