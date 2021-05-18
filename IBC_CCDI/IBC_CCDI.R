## A test run of NetBID on HDAC6 data, for CCDI.
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org) on Apr. 20, 2021
## R-3.6

library(NetBID2) # version: 0.1.2

#setwd("/Users/qpan/Desktop/IBC_CCDI")
## 0 configuration
exp.eset <- "GEP.195.quantileNormalized.geneLevel_log2.eset" # path to expression eset
network.tf <- "TCGABRCA_log2fpkm.TF.txt" # path to TF network by SJARACNe
network.sig <- "TCGABRCA_log2fpkm.SIG.txt" # path to SIG network by SJARACNe
outdir <- "./" # path of output directory
project_name <- sprintf('CCIDTestRun_HDAC6ofIBC.%s',format(Sys.time(), "%Y_%m_%d")) # project name should be readable

analysis.par  <- NetBID.analysis.dir.create(project_main_dir=outdir, project_name=project_name, tf.network.file = network.tf, sig.network.file = network.sig)
load(exp.eset)
analysis.par$exp.eset <- GEP.195.quantileNormalized.geneLevel_log2.eset # add expression eset to analysis.par
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

## 1 calculate the activity
#NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_',html_info_limit=FALSE) # QC of network
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
#draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_',html_info_limit=TRUE) # QC of network

# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$exp.eset),es.method='weightedmean')
# calculate activity
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$exp.eset)[colnames(ac_mat),], feature_info=NULL,annotation_info='activity in net-dataset')
# QC plot for activity eset
#draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_',emb_plot_type = '3D')
# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

## 2 differential analysis
# NetBID.loadRData(analysis.par=analysis.par,step='act-get')
analysis.par$DE <- list(); analysis.par$DA <- list() # Create empty list to store comparison result
# compare IBC and non-IBC
comp_name <- 'IBC.Vs.NIBC'
phe_info <- pData(analysis.par$exp.eset)
G1  <- rownames(phe_info)[which(phe_info$`IBC`=='IBC')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`IBC`=='NIBC')] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$exp.eset,G1=G1,G0=G0,G1_name='IBC',G0_name='non-IBC')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='IBC',G0_name='non-IBC')
analysis.par$DE[[comp_name]] <- DE_gene_bid; analysis.par$DA[[comp_name]] <- DA_driver_bid
# save differential analysis results
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

## 3 highlight the top drivers
draw.NetBID(DA_list=analysis.par$DA, DE_list=analysis.par$DE, main_id='IBC.Vs.NIBC',pdf_file=sprintf('%s/NetBID_TOP.pdf',analysis.par$out.dir.PLOT),text_cex=0.8) # Save as PDF

## 4 generate the masster tables
#NetBID.loadRData(analysis.par=analysis.par,step='act-DA')

# Get all comparison names
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones

# Creat the final master table
db.preload(use_level='gene',use_spe='human',update=FALSE)
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp, DE=analysis.par$DE, DA=analysis.par$DA, target_list=analysis.par$merge.network$target_list,
                                                  tf_sigs=tf_sigs, z_col='Z-statistics', display_col=c('logFC','P.Value'), main_id_type='external_gene_name')

# Path and file name of the output EXCEL file
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
mark_gene <- list(IBC=c('WIF1','TNC','GAD1','DKK2','EMX2'), NIBC=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1')) # highlight marker genes
mark_col <- list(IBC='green',NIBC='blue') # customize highlight color codes
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col) # Save the final master table as EXCEL file

# Save Step 4 analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')



