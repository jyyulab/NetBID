## Complete Pipeline for Driver Estimation and Master Table Creation ##

############### Step 0: Preparation ###############

# Load NetBID2 package
library(NetBID2)

# Get the demo's constructed network data
network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo network in the package
network.project.name <- 'project_2019-02-14'
# Define main working directory and project name
project_main_dir <- 'test/' # user defined main directory for the project, one main directory could
current_date <- format(Sys.time(), "%Y-%m-%d") # optional, if user like to add current date to name the project folder
project_name <- sprintf('driver_%s',current_date) # project name for the project folders under main directory

# Create a hierarchcial working directory and return a list contains the hierarchcial working directory information
# This list object (analysis.par) is an ESSENTIAL variable in driver estimation pipeline
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

############### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############

# If use the same expression dataset as in the network construction, just reload it directly
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) # RData saved after QC in the network construction step
analysis.par$cal.eset <- network.par$net.eset

# Save Step 1 network.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

# Reload network.par RData from Step 1
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')

# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

# Creat QC report for the network
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_')
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_')

# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

# Create eset using activity matrix
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# QC plot for activity eset
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_')

# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############

# Reload network.par RData from Step 2
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# First comparison: G4 vs. WNT
comp_name <- 'G4.Vs.WNT' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='WNT')] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# Second comparison: G4 vs. SHH
comp_name <- 'G4.Vs.SHH' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='SHH')] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

## Third comparison: G4 vs. others
# Combine the comparison results from `G4.Vs.WNT` and `G4.Vs.SHH`
comp_name <- 'G4.Vs.otherTwo' # Each comparison must has a name
DE_gene_comb <- combineDE(DE_list=list(WNT=analysis.par$DE$`G4.Vs.WNT`,SHH=analysis.par$DE$`G4.Vs.SHH`))
DA_driver_comb <- combineDE(DE_list=list(WNT=analysis.par$DA$`G4.Vs.WNT`,SHH=analysis.par$DA$`G4.Vs.SHH`))
analysis.par$DE[[comp_name]] <- DE_gene_comb$combine
analysis.par$DA[[comp_name]] <- DA_driver_comb$combine

# Driver table of top DE
draw.combineDE(DE_gene_comb)
draw.combineDE(DE_gene_comb,pdf_file=sprintf('%s/combineDE.pdf',analysis.par$out.dir.PLOT)) # Save it as PDF

# Driver table of top DA
draw.combineDE(DA_driver_comb)
draw.combineDE(DA_driver_comb,pdf_file=sprintf('%s/combineDA.pdf',analysis.par$out.dir.PLOT)) # Save it as PDF

# Another way to do the third comparison: G4 vs. others
comp_name <- 'G4.Vs.others'
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # Combine other groups as the Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# Save Step 3 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

# Visualize top drivers
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others')
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others',pdf_file=sprintf('%s/NetBID_TOP.pdf',analysis.par$out.dir.PLOT),text_cex=0.8) # Save as PDF

############### Step 4: Generate a master table for drivers (ms-tab) ###############

# Reload analysis.par RData from Step 3
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones
# Prepare the conversion table (OPTIONAL)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes)
analysis.par$transfer_tab <- transfer_tab
# Creat the final master table
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                                  network=analysis.par$merge.network$target_list,
                                                  tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                                  main_id_type='external_gene_name')

# Path and file name of the output EXCEL file
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# Highlight marker genes
mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),
                  SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
                  G4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
# Customize highlight color codes
#mark_col <- get.class.color(names(mark_gene)) # this randomly assign color codes
mark_col <- list(G4='green','WNT'='blue','SHH'='red')
# Save the final master table as EXCEL file
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col)

# Save Step 4 analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')


