## This script aims to get activity and master table ##
############# preparations ############################

library(NetBID2)

######### write in paramters
network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo network in the package
network.project.name <- 'project_2019-02-14' #

project_main_dir <- 'test/' ### user defined main directory for the project, one main directory could have multiple projects, separeted by project name
current_date <- format(Sys.time(), "%Y-%m-%d") ## current date for project running, suggested to add the time tag
project_name <- sprintf('driver_%s',current_date) ## project name for the project

## analysis.par is very essential in the analysis!!
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, prject_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

########### !!!!!!!!!! pipelines !!!!!!!!!! ############

######################################################### Step1: load in gene expression dataset for analysis (exp-load,exp-cluster,exp-QC)

################## if use same expression datasets as in the network construction, directly do the followings:
# if use same expression datasets as in the network construction, directly do the followings:
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) ## RData from network construction
analysis.par$cal.eset <- network.par$net.eset
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

######## otherwise can follow the similar steps in pipeline_network_Demo1.R for expression dataset QC, cluster

######################################################### Step2: Step2: read in network files and activity calculation (act-get)
## load the RData for calculation
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')

# get network info ! three list(network_dat, target_list, igraph_obj)
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_')
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_')

### merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

# get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# QC plots
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_')

# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

######################################################### Step3: get differentiated expression/differentiated activity for all possible drivers (act-DA)
## load RData
NetBID.loadRData(analysis.par=analysis.par,step='act-get')

# generate a list for multiple comparisons
analysis.par$DE <- list()
analysis.par$DA <- list()
# get compared sample names
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # get sample list for G1
# G4 Vs. WNT
comp_name <- 'G4.Vs.WNT'
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='WNT')] # get sample list for G0
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
# save to analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# G4 Vs. SHH
comp_name <- 'G4.Vs.SHH'
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='SHH')] # get sample list for G0
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
# save to analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

## G4 Vs. others
# combine the above two comparisons
comp_name <- 'G4.Vs.otherTwo'
DE_gene_comb <- combineDE(DE_list=list(WNT=analysis.par$DE$`G4.Vs.WNT`,SHH=analysis.par$DE$`G4.Vs.SHH`))
DA_driver_comb <- combineDE(DE_list=list(WNT=analysis.par$DA$`G4.Vs.WNT`,SHH=analysis.par$DA$`G4.Vs.SHH`))
analysis.par$DE[[comp_name]] <- DE_gene_comb$combine
analysis.par$DA[[comp_name]] <- DA_driver_comb$combine
# draw combineDE/DA
draw.combineDE(DE_gene_comb)
draw.combineDE(DE_gene_comb,pdf_file=sprintf('%s/combineDE.pdf',analysis.par$out.dir.PLOT))

draw.combineDE(DA_driver_comb)
draw.combineDE(DA_driver_comb,pdf_file=sprintf('%s/combineDA.pdf',analysis.par$out.dir.PLOT))

# combine the sample list
comp_name <- 'G4.Vs.others'
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # get sample list for G0
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
# save to analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# draw results for NetBID
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others')
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others',pdf_file=sprintf('%s/NetBID_TOP.pdf',analysis.par$out.dir.PLOT),text_cex=0.8)

# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

######################################################### Step4: generate master table (ms-tab)
## load RData
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')
# load db
db.preload(use_level='gene',use_spe='human',update=FALSE)
# generate master table
all_comp <- names(analysis.par$DE) ## get all comparison name for output
# prepare transfer table (optinal)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes) ## get transfer table !!!
analysis.par$transfer_tab <- transfer_tab
# generate master table
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                               network=analysis.par$merge.network$target_list,
                                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                               main_id_type='external_gene_name')

# output into excel files
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# can add marker gene
mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
                  G4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
#mark_col <- get.class.color(names(mark_gene))
mark_col <- list(G4='green','WNT'='blue','SHH'='red')
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col)

# save to RData
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')

###################################### finish driver estimation part !!! Cheers !!! #########################################


