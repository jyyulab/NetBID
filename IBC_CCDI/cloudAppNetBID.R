#################
##Jingjing.liu@stjude.org
##2021-07-26
#################
library(NetBID2)
library(optparse)
####input data required####

option_list = list(
  make_option(
   c("-e", "--expression-set"), 
    type="character", 
    default=NULL, 
    help="File containing gene expression data.", 
    metavar="character"
  ),
  make_option(
    c("-t", "--tf-set"), 
    type="character", 
    default=NULL,
    help="Filename of the transcription factor network from SJARACNe.", 
    metavar="character"
  ),
  make_option(
    c("-s", "--sig-set"),
    type="character",
    default=NULL,
    help="Filename of SIG network from SJARACNe",
    metavar="character"
  ),
  make_option(
    c("-m", "--metadata"),
    type="character",
    default=NULL,
    help="Filename of metadata describing samples",
    metavar="character"
  ),
  make_option(
    c("-p", "--project"),
    type="character",
    default="project",
    help="Output project name",
    metavar="character"
  )
); 
 
opt_parser = OptionParser(prog = "cloudAppNetBID.R", 
                          description = "Analyze an expression set with NetBID.",
                          option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt) || is.null(opt$`expression-set`) || is.null(opt$`tf-set`) || is.null(opt$`sig-set`) || is.null(opt$`metadata`)) {
  print_help(opt_parser)
  q(status=1)
}

exp_mat_path <- opt$`expression-set` # path to expression matrix first column with unique gene name/probeID
pd_path<- opt$`metadata` # path to metadata file
network.tf_path <- opt$`tf-set` # path to TF network by SJARACNe
network.sig_path <- opt$`sig-set`  # path to SIG network by SJARACNe

####input data  optional####
outdir <- "./" # path of output directory 
project_name <- opt$`project` # user define or default 

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE) #default use gene levle and human species 

#####step0 load data####
exp_mat<-read.csv(exp_mat_path,row.names = 1)
pd<-read.csv(pd_path)
rownames(pd)<-pd$sampleID
cal.eset<-generate.eset(exp_mat = exp_mat,phenotype_info = pd)
  
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=outdir, project_name=paste("NetBID_", project_name, sep=""), tf.network.file = network.tf_path, sig.network.file = network.sig_path)

analysis.par$cal.eset <- cal.eset # add expression eset to analysis.par
NetBID.saveRData(analysis.par=analysis.par,step='exp-load')

####step1. build network####
# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

####step2. calculate activity####
# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

# Create eset using activity matrix
analysis.par$ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],feature_info=NULL)

###step3. Get differential expression (DE) / differential activity (DA) for drivers####

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# the comparison group
pd<-pData(analysis.par$cal.eset)
levels<-as.character(unique(pd$comparison))

g1_name<-levels[1];g0_name<-levels[2];comp_name<-sprintf("%s.Vs.%s",g1_name,g0_name)

G1<-rownames(pd)[which(pd$comparison==g1_name)];G0<-rownames(pd)[which(pd$comparison==g0_name)]

DE_gene_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name=g1_name,G0_name=g0_name)
DA_driver_limma   <- getDE.limma.2G(eset=analysis.par$ac.eset,G1=G1,G0=G0,G1_name=g1_name,G0_name=g0_name)

# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_limma
analysis.par$DA[[comp_name]] <- DA_driver_limma

####step4. generate master table####
# Get all comparison names
all_comp <- names(analysis.par$DE)

analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                                  target_list=analysis.par$merge.network$target_list,
                                                  tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                                  main_id_type='external_gene_name')

out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file)

# Save analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')

####plot_TOP30 NetBID drivers####
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,pdf_file =sprintf("%s/NetBID_top30.pdf",analysis.par$out.dir.PLOT),text_cex = 0.8,col_srt = 0)
