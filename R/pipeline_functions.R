############################# NetBID2 Functions
######## pre-request packages
# R >= 3.4.0
library(Biobase) ## basic functions for bioconductor
library(GEOquery) ## for samples from GEO
library(limma) ## for data normalization for micro-array
library(DESeq2) ## for data normalization for RNASeq
library(tximport) ## for data import from Salmon/sailfish/kallisto/rsem/stringtie output
library(arrayQualityMetrics) ## for QC

library(RColorBrewer) ## for color scale
library(colorspace) ## for color scale
library(plot3D) ## for 3D plot
library(vsn) ## for QC plot, require hexbin
library(gplots) ## for some plot like venn
library(igraph) ## for network related functions
library(plotrix) ## for draw.ellipse

library(biomaRt) ## for gene id conversion
library(openxlsx) ## for output into excel
library(impute) ## for impute

library(msigdbr) ## for msigDB gene sets
library(ComplexHeatmap) ## for complex heatmap
library(umap) ## for umap visualization
library(rhdf5) ## for read in MICA results

###### preload database into R
# gene info from GENCODE, support human and mouse, other species may try
# required: use_level, use_spe
# use_level: transcript, gene
# use_spe: suggested human / mouse, others may work
# TF_list, SIG_list: can define byself or use in db/TF.txt db/SIG.txt (external_gene_symbol/hgnc_symbol)
# input_attr_type: applied when TF_list/Sig_list not equal to NULL
# update: 0, judge existence, 1, always update
db.preload <- function(use_level='transcript',use_spe='human',update = 0,
                       TF_list=NULL,SIG_list=NULL,input_attr_type='external_gene_name'){
  ## load annotation info, including: TF/Sig list, gene info
  if(is.null(main.dir)){message('main.dir required, please check and re-try !');return(FALSE)}
  db.dir <- sprintf("%s/db/",main.dir)
  use_spe <- toupper(use_spe)
  output.db.dir <- sprintf('%s/%s',db.dir,use_spe)
  if(!file.exists(output.db.dir)){
    dir.create(output.db.dir)
  }
  RData.file <- sprintf('%s/%s_%s.RData', output.db.dir,use_spe,use_level)
  if (update == 1 | !file.exists(RData.file)) {
    ## get info from use_spe
    ensembl <- useMart("ensembl")
    all_ds  <- listDatasets(ensembl)
    w1 <- grep(sprintf("^%s GENES",use_spe),toupper(all_ds$description))
    if(length(w1)==1){
      w2 <- all_ds[w1,1]
      mart <- useMart(biomart="ensembl", dataset=w2) ## get id for input spe
      message(sprintf('Read in ensembl annotation file for %s and output all db files in %s/%s !',use_spe,db.dir,use_spe))
    }
    if(length(w1)==0){
      message(sprintf('Check input use_spe parameter: %s, not included in the ensembl database',use_spe))
      return(FALSE)
    }
    if(length(w1)>1){
      w2 <- paste(all_ds[w1,2],collapse=';')
      message(sprintf('Check input use_spe parameter: %s, more than one species match in ensembl database : %s,
                      please check and re-try',use_spe,w2))
      return(FALSE)
    }
    ## get attributes for mart
    filters <- listFilters(mart)
    attributes <- listAttributes(mart)
    ensembl.attr.transcript <- c('ensembl_transcript_id_version','ensembl_gene_id_version',
                                 'external_transcript_name','external_gene_name',
                      'transcript_biotype','gene_biotype',
                      'chromosome_name','strand','start_position','end_position','band','transcript_start','transcript_end',
                      'description','phenotype_description','hgnc_symbol','entrezgene','refseq_mrna')
    ensembl.attr.gene <- c('ensembl_gene_id_version','external_gene_name',
                          'gene_biotype',
                          'chromosome_name','strand','start_position','end_position','band',
                          'description','phenotype_description','hgnc_symbol','entrezgene','refseq_mrna')
    ## do not output hgnc in non-human species
    if(use_spe != 'HUMAN')
      ensembl.attr.transcript <- setdiff(ensembl.attr.transcript,'hgnc_symbol')
    if(use_spe != 'HUMAN')
      ensembl.attr.gene <- setdiff(ensembl.attr.gene,'hgnc_symbol')
    ## if too much: Query ERROR: caught BioMart::Exception::Usage: Too many attributes selected for External References
    # judge input type if TF_list, Sig_list not equal to NULL
    if(is.null(TF_list)==FALSE | is.null(SIG_list)==FALSE){
      if(!input_attr_type %in% filters$name){
        message(sprintf('%s not in the filter name, please retry !',input_attr_type));return(FALSE)
      }
      if(!input_attr_type %in% ensembl.attr.transcript)
        ensembl.attr.transcript <- c(ensembl.attr.transcript,input_attr_type)
      if(!input_attr_type %in% ensembl.attr.gene)
        ensembl.attr.gene <- c(ensembl.attr.gene,input_attr_type)
    }
    ## get TF/SIG list and output to output.db.dir, if not defined by user, will use in db/ (human)
    # for TF list
    filter_attr <- input_attr_type
    if(is.null(TF_list)){ ## use TF.txt in db/
      filter_attr <- 'hgnc_symbol' ## for human
      TF_list <- read.delim(sprintf('%s/TF.txt',db.dir),stringsAsFactors=FALSE,header=F)$V1
      if(use_spe != 'HUMAN'){ ## if spe not human
        filter_attr <- 'external_gene_name'
        tmp1 <- getBM(attributes=c('hsapiens_homolog_associated_gene_name','external_gene_name'),values=TRUE,mart=mart,filters='with_hsapiens_homolog')
        TF_list <- unique(tmp1[which(tmp1[,1] %in% TF_list),2])
      }
    }
    if(use_level=='transcript'){
      message(sprintf('Begin read TF list information from ensembl for %s !',use_spe))
      TF_info  <- getBM(attributes = ensembl.attr.transcript,values=TF_list, mart=mart, filters=filter_attr)
    }
    if(use_level=='gene'){
      message(sprintf('Begin read TF list information from ensembl for %s !',use_spe))
      TF_info  <- getBM(attributes = ensembl.attr.gene,values=TF_list, mart=mart, filters=filter_attr)
    }
    # for SIG list
    if(is.null(SIG_list)){
      filter_attr <- 'hgnc_symbol' ## for human
      SIG_list <- read.delim(sprintf('%s/SIG.txt',db.dir),stringsAsFactors=FALSE,header=F)$V1
      if(use_spe != 'HUMAN'){ ## if spe not human
        filter_attr <- 'external_gene_name'
        tmp1 <- getBM(attributes=c('hsapiens_homolog_associated_gene_name','external_gene_name'),values=TRUE,mart=mart,filters='with_hsapiens_homolog')
        SIG_list <- unique(tmp1[which(tmp1[,1] %in% SIG_list),2])
      }
    }
    if(use_level=='transcript'){
      message(sprintf('Begin read SIG list information from ensembl for %s !',use_spe))
      SIG_info  <- getBM(attributes = ensembl.attr.transcript,values=SIG_list, mart=mart, filters=filter_attr)
    }
    if(use_level=='gene'){
      message(sprintf('Begin read SIG list information from ensembl for %s !',use_spe))
      SIG_info  <- getBM(attributes = ensembl.attr.gene,values=SIG_list, mart=mart, filters=filter_attr)
    }
    # check input not in the list
    miss_TF  <- unique(setdiff(TF_list,TF_info[[filter_attr]]))
    miss_SIG <- unique(setdiff(SIG_list,SIG_info[[filter_attr]]))
    if(length(miss_TF)>0){message(sprintf("%d TFs could not match,please check and choose to re-try : %s",
                                          length(miss_TF),paste(sort(miss_TF),collapse=';')))}
    if(length(miss_SIG)>0){message(sprintf("%d SIGs could not match,please check and choose to re-try : %s",
                                           length(miss_SIG),paste(sort(miss_SIG),collapse=';')))}
    ####### output full info
    tf_sigs <- list();tf_sigs$tf <- list();tf_sigs$sig <- list();
    tf_sigs$tf$info  <- TF_info; tf_sigs$sig$info <- SIG_info;
    for(each_id_type in intersect(c('ensembl_transcript_id_version','ensembl_gene_id_version',
                            'external_transcript_name','external_gene_name','hgnc_symbol',
                            'entrezgene','refseq_mrna'),colnames(TF_info))){
      tf_sigs$tf[[each_id_type]] <- setdiff(unique(TF_info[[each_id_type]]),'')
      tf_sigs$sig[[each_id_type]] <- setdiff(unique(SIG_info[[each_id_type]]),'')
    }
    db_info <- all_ds[w1,]
    save(tf_sigs,db_info=db_info,file = RData.file)
  }
  load(RData.file,.GlobalEnv)
  return(TRUE)
}
## get TF SIG list for the input genes
get.TF_SIG.list <- function(use_genes,use_gene_type='ensembl_gene_id',dataset=db_info[1]){
  n1 <- names(tf_sigs$tf)[-1]
  if(use_gene_type %in% n1){
    TF_list <- unique(intersect(use_genes,tf_sigs$tf[[use_gene_type]]))
    SIG_list <- unique(intersect(use_genes,tf_sigs$sig[[use_gene_type]]))
  }else{
    mart <- useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
    #filters <- listFilters(mart)
    attributes <- listAttributes(mart)
    if(!use_gene_type %in% attributes$name){
      message(sprintf('%s not in the attributes for %s, please check and re-try !',use_gene_type,db_info[1]));return(FALSE)
    }
    tmp1 <- getBM(attributes=c(n1[1],use_gene_type),values=tf_sigs$tf[n1[1]][[1]],mart=mart,filters=n1[1])
    TF_list <- tmp1[,2]
    tmp1 <- getBM(attributes=c(n1[1],use_gene_type),values=tf_sigs$sig[n1[1]][[1]],mart=mart,filters=n1[1])
    SIG_list <- tmp1[,2]
    TF_list <- unique(intersect(use_genes,TF_list))
    SIG_list <- unique(intersect(use_genes,SIG_list))
  }
  message(sprintf('%d TFs and %s SIGs are included in the expression matrix !',length(TF_list),length(SIG_list)))
  return(list(tf=TF_list,sig=SIG_list))
}
## get Id transfer table, could be input info for update.eset.feature
get_IDtransfer <- function(from_type=NULL,to_type=NULL,use_genes=NULL,dataset=db_info[1]){
  mart <- useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
  attributes <- listAttributes(mart)
  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,db_info[1]));return(FALSE)
  }
  if(!to_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',to_type,db_info[1]));return(FALSE)
  }
  tmp1 <- getBM(attributes=c(from_type,to_type),values=use_genes,mart=mart,filters=from_type)
  w1 <- apply(tmp1,1,function(x)length(which(is.na(x)==TRUE | x=='')))
  transfer_tab <- tmp1[which(w1==0),]
  return(transfer_tab)
}
## get Id transfer table, for bubble plot !!!
get_IDtransfer2symbol2type <- function(from_type=NULL,use_genes=NULL,dataset=db_info[1],use_level='transcript'){
  message(sprintf('Your setting is at %s level',use_level))
  mart <- useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
  attributes <- listAttributes(mart)
  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,db_info[1]));return(FALSE)
  }
  if(use_level=='gene') tmp1 <- getBM(attributes=c(from_type,c('external_gene_name','gene_biotype')),values=use_genes,mart=mart,filters=from_type)
  if(use_level=='transcript') tmp1 <- getBM(attributes=c(from_type,c('external_gene_name','transcript_biotype')),values=use_genes,mart=mart,filters=from_type)
  w1 <- apply(tmp1,1,function(x)length(which(is.na(x)==TRUE | x=='')))
  transfer_tab <- tmp1[which(w1==0),]
  return(transfer_tab)
}
## get name by transfer tab
get_name_transfertab <- function(x,transfer_tab){
  t1 <- unique(transfer_tab[which(transfer_tab[,1] %in% x),])
  rownames(t1) <- t1[,1]
  x1 <- t1[x,2]
  w1 <- which(is.na(x1)==TRUE)
  x1[w1] <- x[w1]
  x1
}
######## generate project space and working directory structure
## create file directory for NetBID working space, the network part
# project_main_dir: main directory for the project
# project_name: name of the project
# create network.par
NetBID.network.dir.create <- function(project_main_dir=NULL,prject_name=NULL){
  if(exists('network.par')==TRUE){message('network.par is occupied in the current session,please manually run: rm(network.par) and re-try, otherwise will not change !');
    return(network.par)}
  if(is.null(project_main_dir)==TRUE){message('project_main_dir required, please input and re-try!');return(FALSE)}
  if(is.null(prject_name)==TRUE){message('prject_name required, please input and re-try!');return(FALSE)}
  network.par <- list()
  network.par$main.dir <- project_main_dir
  network.par$project.name <- prject_name
  network.par$out.dir <- sprintf('%s/%s/',network.par$main.dir,network.par$project.name)
  # create output directory
  if (!dir.exists(project_main_dir)) {
    dir.create(project_main_dir, recursive = TRUE)
  }
  if (!dir.exists(network.par$out.dir)) {
    dir.create(network.par$out.dir, recursive = TRUE)
  }
  network.par$out.dir.QC <- paste0(network.par$out.dir, '/QC/')
  if (!dir.exists(network.par$out.dir.QC)) {
    dir.create(network.par$out.dir.QC, recursive = TRUE) ## directory for QC
  }
  network.par$out.dir.DATA <- paste0(network.par$out.dir, '/DATA/')
  if (!dir.exists(network.par$out.dir.DATA)) {
    dir.create(network.par$out.dir.DATA, recursive = TRUE) ## directory for DATA
  }
  network.par$out.dir.SJAR <- paste0(network.par$out.dir, '/SJAR/')
  if (!dir.exists(network.par$out.dir.SJAR)) {
    dir.create(network.par$out.dir.SJAR, recursive = TRUE) ## directory for SJARAcne
  }
  message(sprintf('Project space created, please check %s',network.par$out.dir))
  return(network.par)
}
## create file directory for NetBID working space, the analysis part
# project_main_dir: main directory for the project
# project_name: name of the project
# create analysis.par
NetBID.analysis.dir.create <- function(project_main_dir=NULL,prject_name=NULL,network_dir=NULL, network_project_name=NULL){
  if(exists('analysis.par')==TRUE){message('analysis.par is occupied in the current session,please manually run: rm(analysis.par) and re-try, otherwise will not change !');
    return(analysis.par)}
  if(is.null(project_main_dir)==TRUE){message('project_main_dir required, please input and re-try!');return(FALSE)}
  if(is.null(prject_name)==TRUE){message('prject_name required, please input and re-try!');return(FALSE)}
  if(is.null(network_dir)==TRUE){message('network_dir required, please input and re-try!');return(FALSE)}
  if(is.null(network_project_name)==TRUE){message('network_project_name required, please input and re-try!');return(FALSE)}
  analysis.par <- list()
  analysis.par$main.dir <- project_main_dir
  analysis.par$project.name <- prject_name
  analysis.par$out.dir <- sprintf('%s/%s/',analysis.par$main.dir,analysis.par$project.name)
  # create output directory
  if (!dir.exists(analysis.par$out.dir)) {
    dir.create(analysis.par$out.dir, recursive = TRUE)
  }
  analysis.par$out.dir.QC <- paste0(analysis.par$out.dir, '/QC/')
  if (!dir.exists(analysis.par$out.dir.QC)) {
    dir.create(analysis.par$out.dir.QC, recursive = TRUE) ## directory for QC
  }
  analysis.par$out.dir.DATA <- paste0(analysis.par$out.dir, '/DATA/')
  if (!dir.exists(analysis.par$out.dir.DATA)) {
    dir.create(analysis.par$out.dir.DATA, recursive = TRUE) ## directory for DATA
  }
  analysis.par$out.dir.PLOT <- paste0(analysis.par$out.dir, '/PLOT/')
  if (!dir.exists(analysis.par$out.dir.PLOT)) {
    dir.create(analysis.par$out.dir.PLOT, recursive = TRUE) ## directory for Result Plots
  }
  analysis.par$tf.network.file  <- sprintf('%s/SJAR/%s/output_tf_sjaracne_%s_out_.final/consensus_network_ncol_.txt',
                                          network_dir,network_project_name,network_project_name)
  analysis.par$sig.network.file <- sprintf('%s/SJAR/%s/output_sig_sjaracne_%s_out_.final/consensus_network_ncol_.txt',
                                          network_dir,network_project_name,network_project_name)
  if(file.exists(analysis.par$tf.network.file)){
    message(sprintf('TF network file found in %s',analysis.par$tf.network.file))
  }else{
    message(sprintf('TF network file not found in %s, please check and re-try !',analysis.par$tf.network.file))
    return(FALSE)
  }
  if(file.exists(analysis.par$sig.network.file)){
    message(sprintf('SIG network file found in %s',analysis.par$sig.network.file))
  }else{
    message(sprintf('SIG network file not found in %s, please check and re-try ',analysis.par$sig.network.file))
    return(FALSE)
  }
  message(sprintf('Analysis space created, please check %s',analysis.par$out.dir))
  return(analysis.par)
}
###### overall functions
## check network.par
check_network.par <- function(network.par,step='pre-load'){
  if(class(network.par)!='list'){message('Invalid network.par !');return(FALSE)}
  if(step=='pre-load'){
    n1 <- names(network.par)
    n2 <- setdiff(c('main.dir','project.name','out.dir','out.dir.QC','out.dir.DATA','out.dir.SJAR'),n1)
    if(length(n2)>0){
      message(sprintf('Miss %s network.par, please check and re-try !',paste(n2,collapse=';')))
      return(FALSE)
    }
  }
  if(step %in% c('exp-load','exp-QC','exp-cluster')){
    n1 <- names(network.par)
    n2 <- setdiff(c('main.dir','project.name','out.dir','out.dir.QC','out.dir.DATA','out.dir.SJAR','net.eset'),n1)
    if(length(n2)>0){
      message(sprintf('Miss %s network.par, please check and re-try !',paste(n2,collapse=';')))
      return(FALSE)
    }
  }
  return(TRUE)
}
## check analysis.par
check_analysis.par <- function(analysis.par,step='pre-load'){
  if(class(analysis.par)!='list'){message('Invalid network.par !');return(FALSE)}
  if(step=='pre-load'){
    n1 <- names(analysis.par)
    n2 <- setdiff(c('main.dir','project.name','out.dir','out.dir.QC','out.dir.DATA','out.dir.PLOT','tf.network.file','sig.network.file'),n1)
    if(length(n2)>0){
      message(sprintf('Miss %s analysis.par, please check and re-try !',paste(n2,collapse=';')))
      return(FALSE)
    }
  }
  if(step %in% c('exp-load','exp-QC','exp-cluster')){
    n1 <- names(analysis.par)
    n2 <- setdiff(c('main.dir','project.name','out.dir','out.dir.QC','out.dir.DATA','out.dir.PLOT','cal.eset','tf.network.file','sig.network.file'),n1)
    if(length(n2)>0){
      message(sprintf('Miss %s analysis.par, please check and re-try !',paste(n2,collapse=';')))
      return(FALSE)
    }
  }
  return(TRUE)
}
## step save and load
NetBID2.saveRData <- function(network.par=NULL,analysis.par=NULL,step='exp-load'){
  if(is.null(network.par)==FALSE & is.null(analysis.par)==FALSE){
    message('Can not save network.par and analysis.par at once, please only use one !');return(FALSE)
  }
  if(is.null(network.par)==FALSE){
    check_network.par(network.par = network.par,step=step)
    save(network.par,file=sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step))
    message(sprintf('Successful save to %s',sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step)))
  }
  if(is.null(analysis.par)==FALSE){
    check_analysis.par(analysis.par = analysis.par,step=step)
    save(analysis.par,file=sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step))
    message(sprintf('Successful save to %s',sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step)))
  }
}
NetBID2.loadRData <- function(network.par=NULL,analysis.par=NULL,step='exp-load'){
  if(is.null(network.par)==FALSE & is.null(analysis.par)==FALSE){
    message('Can not load network.par and analysis.par at once, please only use one !');return(FALSE)
  }
  if(is.null(network.par)==FALSE){
    check_network.par(network.par = network.par,step=step)
    load(file=sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step),.GlobalEnv)
    message(sprintf('Successful load from %s',sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step)))
  }
  if(is.null(analysis.par)==FALSE){
    check_analysis.par(analysis.par = analysis.par,step=step)
    load(file=sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step),.GlobalEnv)
    message(sprintf('Successful load from %s',sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step)))
  }
}
############################# load functions
## load exp data from GEO
# network.par required
# GSE,GPL required
# update=0, if the RData exists, will directly load
load.exp.GEO <- function(out.dir=network.par$out.dir.DATA,GSE = NULL,GPL = NULL,getGPL=TRUE,update = 0){
  if(is.null(out.dir)==TRUE){
    message('out.dir required, please re-try')
    return(FALSE)
  }
  if (is.null(GSE) | is.null(GPL)){
    message('GSE and GPL required, please re-try')
    return(FALSE)
  }
  expRData_dir <- sprintf('%s/%s_%s.RData', out.dir, GSE,GPL)
  if (file.exists(expRData_dir) & update == 0) {
    message(sprintf('RData exist in %s and update==0, will directly load from RData .',expRData_dir))
    load(expRData_dir)
  } else{
    eset <- getGEO(GSE, GSEMatrix = TRUE, getGPL = getGPL)
    if (length(eset) > 1)
      idx <- grep(GPL, attr(eset, "names"))
    else
      idx <- 1
    eset <- eset[[idx]]
    save(eset, file = expRData_dir)
    message(sprintf('RData for the eset is saved in %s .',expRData_dir))
  }
  return(eset)
}

## load exp from RNASeq, demo functions for salmon
load.exp.RNASeq.demoSalmon <- function(salmon_dir = '',tx2gene=NULL,use_phenotype_info = NULL,use_sample_col=NULL,use_design_col=NULL,return_type='eset',merge_level='gene') {
    files <- file.path(salmon_dir, list.files(salmon_dir), "quant.sf")
    sample_name <- gsub('(.*)_salmon', '\\1', list.files(salmon_dir))
    names(files) <- sample_name
    w1 <- length(files)
    message(sprintf('%d salmon_output/quant.sf found !',w1))
    if(is.null(tx2gene)){
      gene_info <- read.delim(file = files[1], stringsAsFactors = FALSE)[, 1]
      gen1 <- sapply(gene_info, function(x)unlist(strsplit(x, '\\|')))
      gen1 <- t(gen1)
      if(merge_level=='gene'){
        tx2gene <- data.frame('transcript' = gene_info,'gene' = gen1[,2],stringsAsFactors = FALSE)
      }else{
        tx2gene <- data.frame('transcript' = gene_info,'gene' = gen1[,1],stringsAsFactors = FALSE)
      }
    }
    eset <- load.exp.RNASeq.demo(files,type='salmon',tx2gene=tx2gene,use_phenotype_info=use_phenotype_info,use_sample_col=use_sample_col,use_design_col=use_design_col,
                                 return_type=return_type,merge_level=merge_level)
    return(eset)
}
## load exp from RNASeq, demo functions
# return_type: dds, eset, both
load.exp.RNASeq.demo <- function(files,type='salmon',tx2gene=NULL,use_phenotype_info = NULL,use_sample_col=NULL,use_design_col=NULL, return_type='eset',merge_level='gene') {
  n1 <- colnames(use_phenotype_info)
  if(!use_sample_col %in% n1){
    message(sprintf('%s not in the colnames of use_phenotype_info, please check and re-try !',use_sample_col));return(FALSE)
  }
  if(!use_design_col %in% n1){
    message(sprintf('%s not in the colnames of use_phenotype_info, please check and re-try !',use_design_col));return(FALSE)
  }
  # get intersected samples
  rownames(use_phenotype_info) <- use_phenotype_info[,use_sample_col]
  w1 <- intersect(names(files),rownames(use_phenotype_info))
  files <- files[w1]; use_phenotype_info <- use_phenotype_info[w1,]
  if(length(w1)==0){
    message(sprintf('No sample could match the %s in the use_phenotype_info, please check and re-try !',use_sample_col))
    return(FALSE)
  }
  message(sprintf('%d samples could further processed !',length(w1)))
  # import into txi
  txi <- tximport(files, type = type, tx2gene = tx2gene) ## key step one, tximport
  use_phenotype_info <- use_phenotype_info[colnames(txi$abundance), ]
  tmp_phe <- cbind(group=use_phenotype_info[,use_design_col],use_phenotype_info,stringsAsFactors=FALSE)
  # import into deseq2
  dds <- DESeqDataSetFromTximport(txi, colData = tmp_phe, design =  ~ group) ## key step two, DESeqDataSetFromTximport
  dds <- DESeq(dds)
  if(return_type=='dds'){
    return(dds)
  }else{
    vsd <- vst(dds)
    mat <- assay(vsd)
    eset <- generate.eset(exp_mat=mat, phenotype_info = use_phenotype_info, feature_info = NULL, annotation_info='Salmon')
    if(return_type=='eset') return(eset)
    if(return_type=='both') return(list(eset=eset,dds=dds))
  }
}
####################### eset related functions
## generate eset
# phe, gen can be null
generate.eset <- function(exp_mat=NULL, phenotype_info=NULL, feature_info=NULL, annotation_info='') {
  if(is.null(exp_mat)){
    message('exp_mat required, please re-try !');
    return(FALSE)
  }
  if (is.null(phenotype_info)) {
    phenotype_info <- data.frame(group = colnames(exp_mat), stringsAsFactors = FALSE)
    rownames(phenotype_info) <- colnames(exp_mat)
  }
  if (is.null(feature_info)) {
    feature_info <- data.frame(gene = rownames(exp_mat), stringsAsFactors = FALSE)
    rownames(feature_info) <- rownames(exp_mat)
  }
  if(class(phenotype_info)=='character'){
    phenotype_info <- data.frame(group = phenotype_info, stringsAsFactors = FALSE)
    rownames(phenotype_info) <- colnames(exp_mat)
  }
  if(class(feature_info)=='character'){
    feature_info <- data.frame(gene = feature_info, stringsAsFactors = FALSE)
    rownames(feature_info) <- rownames(exp_mat)
  }
  eset <-
    new(
      "ExpressionSet",
      phenoData = new("AnnotatedDataFrame", phenotype_info),
      featureData = new("AnnotatedDataFrame", feature_info),
      annotation = annotation_info,
      exprs = as.matrix(exp_mat)
    )
  return(eset)
}
###
## merge two eset, remove batch effects
merge.eset <- function(eset1,eset2,group1,group2,phe1,phe2,
                       use_col = NULL,
                       remove_batch = TRUE) {
  mat1 <- exprs(eset1)
  mat2 <- exprs(eset2)
  w1 <- intersect(rownames(mat1), rownames(mat2))
  rmat <- cbind(mat1[w1, ], mat2[w1, ])
  choose1 <- apply(rmat <= quantile(rmat, probs = 0.05), 1, sum) <= ncol(rmat) * 0.90 ## low expressed genes
  rmat <- rmat[choose1, ]
  rmat <- normalizeQuantiles(rmat) ## limma quantile normalization
  #rmat <- apply(rmat, 2, function(x) {
  #    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  #})
  if(length(use_col)>1)
    rphe <- rbind(phe1[colnames(mat1), use_col], phe2[colnames(mat2), use_col])
  if(length(use_col)==1){
    rphe <- c(phe1[colnames(mat1), use_col], phe2[colnames(mat2), use_col])
    rphe <- data.frame(rphe,stringsAsFactors=FALSE); colnames(rphe) <- use_col;
    rownames(rphe) <- colnames(rmat)
  }
  rphe$original_group <- c(rep(group1, ncol(mat1)), rep(group2, ncol(mat2)))
  if (remove_batch == TRUE) {
    rmat <- removeBatchEffect(rmat,batch=rphe$original_group)
  }
  reset <- generate.eset(rmat,phenotype_info = rphe, annotation_info = 'combine')
  return(reset)
}

## update eset feature, mainly for gene ID conversion
# ID conversion, or merge transcript level to expression level, use_feature_info can be other dataframe info
# merge_method: mean median max min
# if multi-from --> one to: merge
# if one-from --> multi to: equal distribute ?
update.eset.feature <- function(use_eset=NULL,use_feature_info=NULL,from_feature=NULL,to_feature=NULL,merge_method='median'){
  if(is.null(use_eset)){
    message('use_eset required, please re-try !');
    return(use_eset)
  }
  if(is.null(from_feature) | is.null(to_feature)){
    message('from_feature, to_feature required, please re-try !')
    return(use_eset)
  }
  if(is.null(use_feature_info)) use_feature_info <- fData(use_eset)
  n1 <- colnames(use_feature_info)
  if(!from_feature %in% n1){
    message(sprintf('%s not in in the colnames of use_feature_info, please re-try!',from_feature));return(use_eset)
  }
  if(!to_feature %in% n1){
    message(sprintf('%s not in in the colnames of use_feature_info, please re-try!',to_feature));return(use_eset)
  }
  mat <- exprs(use_eset)
  use_feature_info <- unique(use_feature_info);
  w1 <- which(use_feature_info[,1]!='' & use_feature_info[,2]!='' & is.na(use_feature_info[,1])==FALSE & is.na(use_feature_info[,2])==FALSE)
  use_feature_info <- use_feature_info[w1,]
  g1 <- rownames(mat) ## rownames for the expmat
  f1 <- as.character(use_feature_info[,from_feature]) ## from feature info
  t1 <- as.character(use_feature_info[,to_feature]) ## to feature info
  w1 <- which(f1 %in% g1); f1 <- f1[w1]; t1 <- t1[w1]; ## only consider features in the rownames of expmat
  if(length(w1)==0){
    message(sprintf('Rownames of the expression matrix was not included in the %s column, please check and re-try !',from_feature))
    return(use_eset)
  }
  message(sprintf('%d out of %d rows from original expression matrix will be keeped !',length(w1),length(g1)))
  fc1 <- table(f1); tc1 <- table(t1); fw1 <- which(fc1>1); tw1 <- which(tc1>1); ## check duplicate records
  if(length(fw1)>0){
    message(sprintf('Original feature %s has %d items with duplicate records, will distribute the original values equal to all related items !
            if do not want this, please check and retry !',from_feature,length(fw1)))
    #return(use_eset)
    w2 <- which(f1 %in% names(fw1)) ## need to distribute
    w0 <- setdiff(1:length(f1),w2) ## do not need to distribute
    v1 <- mat[f1[w2],]; rownames(v1) <- paste0(f1[w2],'-',t1[w2]);
    f1[w2] <- paste0(f1[w2],'-',t1[w2]); # update transfer table
    mat <- rbind(v1,mat[f1[w0],]) # update mat table
    fc1 <- table(f1); tc1 <- table(t1); fw1 <- which(fc1>1); tw1 <- which(tc1>1); ## update f1, t1 and related values
  }
  if(length(tw1)>0){
    w2 <- which(t1 %in% names(tw1)) ## need to merge
    w0 <- setdiff(1:length(t1),w2) ## do not need to merge
    mat_new_0 <- mat[f1[w0],]; rownames(mat_new_0) <- t1[w0]  ## mat do not need to merge
    if(merge_method=='mean') tmp1 <- aggregate(mat[f1[w2],],list(t1[w2]),function(x){mean(x,na.rm=TRUE)})
    if(merge_method=='median') tmp1 <- aggregate(mat[f1[w2],],list(t1[w2]),function(x){median(x,na.rm=TRUE)})
    if(merge_method=='max') tmp1 <- aggregate(mat[f1[w2],],list(t1[w2]),function(x){max(x,na.rm=TRUE)})
    if(merge_method=='min') tmp1 <- aggregate(mat[f1[w2],],list(t1[w2]),function(x){min(x,na.rm=TRUE)})
    mat_new_1 <- tmp1[,-1]; rownames(mat_new_1) <- tmp1[,1] ## mat merged
    mat_new <- rbind(mat_new_0,mat_new_1)
  }else{
    mat_new <- mat[f1,]
    rownames(mat_new) <- t1
  }
  new_eset <- generate.eset(exp_mat=mat_new, phenotype_info=pData(use_eset), feature_info=NULL, annotation_info=annotation(use_eset))
  return(new_eset)
}
## update eset phenotype, mainly for input more phenotype info
# use_col: GEO-auto, use :ch1
update.eset.phenotype <- function(use_eset=NULL,use_phenotype_info=NULL,use_sample_col=NULL,use_col='GEO-auto'){
  if(is.null(use_eset)){
    message('use_eset required, please re-try !');
    return(use_eset)
  }
  if(is.null(use_sample_col)){
    message('use_sample_col required, please re-try !')
    return(use_eset)
  }
  if(is.null(use_phenotype_info)) use_phenotype_info <- pData(use_eset)
  if(!use_sample_col %in% colnames(use_phenotype_info)){
    message(sprintf('%s not in the colnames of use_phenotype_info, please re-try!',use_sample_col));return(use_eset)
  }
  if(is.null(use_col)) use_col <- colnames(use_phenotype_info)
  mat <- exprs(use_eset)
  s1 <- colnames(mat) ## all samples
  p1 <- use_phenotype_info[,use_sample_col] ## sample in the phenotype info
  w1 <- which(p1 %in% s1); p1 <- p1[w1]; ## only consider samples in the colnames of expmat
  if(length(w1)==0){
    message(sprintf('Colnames of the expression matrix was not included in the %s column, please check and re-try !',use_sample_col))
    return(use_eset)
  }
  message(sprintf('%d out of %d samples from the expression matrix will be keeped !',length(w1),length(s1)))
  mat_new <- mat[,p1]
  use_phenotype_info <- use_phenotype_info[w1,]
  if(use_col[1] == 'GEO-auto'){
    n1 <- colnames(use_phenotype_info)
    w1 <- c('geo_accession','title','source_name_ch1',n1[grep(':ch1',n1)])
    p1 <- use_phenotype_info[,w1]
    colnames(p1)[4:ncol(p1)] <- gsub('(.*):ch1','\\1',colnames(p1)[4:ncol(p1)])
    colnames(p1)[3] <- gsub('(.*)_ch1','\\1',colnames(p1)[3])
    rownames(p1) <- use_phenotype_info[,use_sample_col]
    p1 <- as.data.frame(apply(p1,2,function(x){if(class(x)=='factor'){as.character(x)}else{x}}),stringsAsFactors=FALSE)
    new_phenotype_info <- p1
  }else{
    if(length(setdiff(use_col,n1))>0){
      message(sprintf('%s not in use_phenotype_info, please re-try!',paste(setdiff(use_col,n1),collapse=';')));return(FALSE)
    }
    p1 <- use_phenotype_info[,use_col]
    p1 <- as.data.frame(apply(p1,2,function(x){if(class(x)=='factor'){as.character(x)}else{x}}),stringsAsFactors=FALSE)
    new_phenotype_info <- p1
  }
  new_eset <- generate.eset(exp_mat=mat_new, phenotype_info=new_phenotype_info, feature_info=fData(use_eset), annotation_info=annotation(use_eset))
  return(new_eset)
}

## IQR fileter;
# loose_gene: genes will use loose threshold
IQR.filter <- function(exp_mat,use_genes=rownames(exp_mat),thre = 0.5,loose_gene=NULL,loose_thre=0.1) {
    use_genes <- intersect(use_genes,rownames(exp_mat))
    use_genes <- setdiff(use_genes,'')
    exp_mat <- exp_mat[use_genes,]
    iqr <- apply(exp_mat, 1, IQR) ## calculate IQR for each gene
    choose0 <- use_genes[iqr > quantile(iqr, loose_thre)] ## for loose_gene
    choose1 <- use_genes[iqr > quantile(iqr, thre)] ## for all genes
    choose2 <- unique(c(intersect(loose_gene, choose0), choose1)) ## union set
    use_vec <- rep(FALSE,length.out=length(use_genes));names(use_vec) <- use_genes
    use_vec[choose2] <- TRUE
    print(table(use_vec))
    return(use_vec)
}
### simple functions
## normalization for RNASeq Count
RNASeqCount.normalize.scale <- function(d,
                                        total = NULL,
                                        pseudoCount = 1) {
  if (!is.data.frame(d))
    d <- data.frame(d)
  if (!all(d > 0))
    d <- d + pseudoCount
  s <- apply(d, 2, sum)
  m <-
    ifelse(is.null(total), as.integer(mean(s)), as.integer(total)) ## total or mean sum
  options(digits = 2 + nchar(m))
  fac <- m / s
  for (i in 1:length(s)) {
    d[, i] <- round(d[, i] * fac[i], 0)
  }
  if (!all(d > 0))
    d <- d + pseudoCount
  d
}
dist2 <- function (x, fun = function(a, b) mean(abs(a - b), na.rm = TRUE),
          diagonal = 0)
{
  if (!(is.numeric(diagonal) && (length(diagonal) == 1)))
    stop("'diagonal' must be a numeric scalar.")
  if (missing(fun)) {
    res = apply(x, 2, function(w) colMeans(abs(x - w), na.rm = TRUE))
  }
  else {
    res = matrix(diagonal, ncol = ncol(x), nrow = ncol(x))
    if (ncol(x) >= 2) {
      for (j in 2:ncol(x)) for (i in 1:(j - 1)) res[i,
                                                    j] = res[j, i] = fun(x[, i], x[, j])
    }
  }
  colnames(res) = rownames(res) = colnames(x)
  return(res)
}
########################### activity-related functions
## functions for activity score calculation, mean, absmean, maxmean, weighted mean ?
es <- function(z, es.method = "mean") {
  if (es.method == "maxmean") {
    n <- length(z)
    m1 <- ifelse(sum(z > 0) > 0, sum(z[z > 0]) / n, 0)
    m2 <- ifelse(sum(z < 0) > 0, sum(z[z < 0]) / n, 0)
    if (m1 > -m2)
      es <- m1
    else
      es <- m2
  }
  else if (es.method == 'absmean') {
    es <- mean(abs(z),na.rm=TRUE)
  }
  else if (es.method == 'mean') {
    es <- mean(z,na.rm=TRUE)
  }
  else if (es.method == 'median') {
    es <- median(z,na.rm=TRUE)
  }
  else if (es.method == 'max') {
    es <- max(z,na.rm=TRUE)
  }
  else if (es.method == 'min') {
    es <- min(z,na.rm=TRUE)
  }
  return(es)
}
std <- function(x) {
  x <- x[!is.na(x)]
  (x - mean(x,na.rm=TRUE)) / sd(x,na.rm=TRUE)
}
## calculate activity for network targets, standard input !!!
# es.method: mean, weighted mean, maxmean, absmean
cal.Activity <- function(all_target=NULL, cal_mat=NULL, es.method = 'mean') {
    ## mean, absmean, maxmean, weightedmean
    use_genes <- row.names(cal_mat)
    #all_target <- all_target[intersect(use_genes, names(all_target))]
    ac.mat <-
      matrix(NA, ncol = ncol(cal_mat), nrow = length(all_target)) ## generate activity matrix, each col for sample, each row for source target
    #z-normalize each sample
    cal_mat <- apply(cal_mat, 2, std)
    for (i in 1:length(all_target)) {
      x <- names(all_target)[i]
      x1 <- all_target[[x]]
      x2 <- unique(intersect(rownames(x1), use_genes)) ## filter target by cal genes
      x1 <- x1[x2, ]
      target_num <- length(x2)
      if (target_num == 0)
        next
      if (target_num == 1){
        ac.mat[i, ] <- cal_mat[x2,]
        next
      }
      if (es.method != 'weightedmean')
        ac.mat[i, ] <- apply(cal_mat[x2,], 2, es, es.method)
      if (es.method == 'weightedmean') {
        weight <- x1$MI * sign(x1$spearman)
        if (length(x2) == 1)
          ac.mat[i, ] <-
            sapply(cal_mat * weight, es, 'mean')[x2]
        else
          ac.mat[i, ] <- apply(cal_mat[x2,] * weight, 2, es, 'mean')
      }
    }
    rownames(ac.mat) <- names(all_target)
    colnames(ac.mat) <- colnames(cal_mat)
    return(ac.mat)
  }
## calculate activity for genesets
# mean, absmean, maxmean
cal.Activity.GS <- function(use_gs2gene=all_gs2gene[c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG')], cal_mat=NULL, es.method = 'mean') {
  while(class(use_gs2gene[[1]])=='list'){
    nn <- unlist(lapply(use_gs2gene,names))
    use_gs2gene <- unlist(use_gs2gene,recursive = FALSE)
    names(use_gs2gene)<-nn
  }
  use_genes <- row.names(cal_mat)
  ac.mat <-
    matrix(NA, ncol = ncol(cal_mat), nrow = length(use_gs2gene)) ## generate activity matrix, each col for sample, each row for source target
  #z-normalize each sample
  cal_mat <- apply(cal_mat, 2, std)
  for (i in 1:length(use_gs2gene)) {
    x <- names(use_gs2gene)[i]
    x1 <- use_gs2gene[[x]]
    x2 <- unique(intersect(x1, use_genes)) ## filter target by cal genes
    target_num <- length(x2)
    if (target_num == 0)
      next
    if (target_num == 1){
      ac.mat[i, ] <- cal_mat[x2,]
      next
    }
    ac.mat[i, ] <- apply(cal_mat[x2,], 2, es, es.method)
  }
  rownames(ac.mat) <- names(use_gs2gene)
  colnames(ac.mat) <- colnames(cal_mat)
  w1 <- apply(ac.mat,1,function(x)length(which(is.na(x)==TRUE)))
  ac.mat <- ac.mat[which(w1==0),]
  return(ac.mat)
}
########################### DE-related functions
#############################################################
#' Network-based Bayesian Inference of Drivers
#'
#' Driver inference by intergrating multiple NetBID results
#'
#' @param netbid.list list of netbid outputs, e.g. list(mRNA=dmrna,wProtein=dwprotein,pProtein=dpprotein)
#' @param signed If TURE, consider directions when combining NetBID results, If FALSE, ignore the directions
#' @return data.frame of case vs. ctrl inference with statistics including cobmbined and individual netSize (network size), z, pval, FDR.BH (FDR by BH method), log2FC (if signed)
#' @examples
#'
#' ###transform network into list of probes
#' dci.kinase.probes<-network.tbl2list(dci.kinase.tbl[,c('source.symbol','target')])
#' dmrna<-netbid(dci.kinase.probes,mrna,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')
#' ##whole proteomcis
#' dci.kinase.wprotein<-transformNetworkListIds(dci.kinase.symbols,ref=fData(wprotein)[,c('geneSymbol','proteinId')])
#' dwprotein<-netbid(dci.kinase.wprotein,wprotein,group.case='CD8pos.DC',group.ctrl='CD8neg.DC')

#' dcomb<-netbidi(list(mRNA=dmrna,wProtein=dwprotein))
#'
#'
#'
#' @export
netbidi<-function(netbid.list,signed=FALSE){
  n<-length(netbid.list)
  if(n<2) stop('At least two netbid outputs are required for netbidi analysis!\n')

  rs<-NULL
  for(i in 1:n){

    if(names(netbid.list)[i]==''|is.null(names(netbid.list)[i]))
      names(netbid.list)[i]<-i

    names(netbid.list[[i]])[-1]<-paste(names(netbid.list[[i]])[-1],names(netbid.list)[i],sep='_')

    netbid.list[[i]]$id<-as.character(netbid.list[[i]]$id)
    if(i==1){
      rs<-netbid.list[[i]]
    }else{
      rs<-dplyr::full_join(rs,netbid.list[[i]],by='id')
    }
  }
  conds<-names(netbid.list)
  if(sum(grepl('log2FC',names(rs)))>0){
    FC.cols<-paste('log2FC',conds,sep='_')
  }else{
    FC.cols<-NULL
  }
  comb<-combinePvalues(rs,pvalue.cols=paste('pval',conds,sep='_'),sign.cols=paste('z',conds,sep='_'),FC.cols =FC.cols,logTransformed=TRUE,log.base=2,method=c('Stouffer'),twosided=TRUE,signed=signed,byRow=TRUE)
  comb<-comb[,-1]
  names(comb)<-paste(names(comb),'comb',sep='_')

  rs<-cbind(comb,rs)

  col.sel<-c('id',
             grep('^netSize_',names(rs),value=TRUE)
             ,
             grep('^z_',names(rs),value=TRUE),
             grep('^log2FC_',names(rs),value=TRUE),
             grep('^pval_',names(rs),value=TRUE),
             grep('^FDR.BH_',names(rs),value=TRUE)
  )

  rs<-rs[,c(col.sel,setdiff(names(rs),col.sel))]
  rs<-dplyr::arrange(rs,pval_comb)
  rs
}
############bid: Bayesian Inference of Drivers by differential analysis of case vs. ctrl
#' @export
#source('pipeline_functions_bid.R')
getDE.BID.2G <-function(eset,G1=NULL, G0=NULL,G1_name=NULL,G0_name=NULL,verbose=TRUE){
  if(verbose==TRUE){
    print(sprintf('G1:%s', paste(G1, collapse = ';')))
    print(sprintf('G0:%s', paste(G0, collapse = ';')))
  }
  exp_mat <- exprs(eset)
  G1 <- intersect(G1,colnames(exp_mat))
  G0 <- intersect(G0,colnames(exp_mat))
  exp_mat <- exp_mat[,c(G1,G0)]
  d<-data.frame(id=rownames(exp_mat),exp_mat,stringsAsFactors=FALSE)
  comp <- factor(c(rep(1,length.out=length(G1)),rep(0,length.out=length(G0))))
  de<-plyr::ddply(d,'id','combRowEvid.2grps',comp=comp,family=gaussian,method='Bayesian',n.iter=5000,nitt=25000,burnin=5000,thin=1,pooling=c('full'),logTransformed=TRUE,restand=FALSE,average.method=c('geometric'))
  names(de)<-gsub('.full','',names(de))
  de$P.Value<-de$pval
  de$adj.P.Val<-p.adjust(de$P.Value,'fdr')
  de$logFC<-sign(de$FC)*log2(abs(de$FC))
  de$AveExpr <- de$AveExp
  de$`Z-statistics`<- de$z
  de$ID <- de$id
  de<-de[,c('ID','logFC','AveExpr','t','P.Value','adj.P.Val','Z-statistics')]
  rownames(de) <- de$ID
  de <- de[order(de$P.Value),]
  tT <- de
  new_mat <- exp_mat
  exp_G1 <- rowMeans(new_mat[rownames(tT),G1])
  exp_G0 <- rowMeans(new_mat[rownames(tT),G0])
  tT <- cbind(tT,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  return(tT)
}
#############################################################
#
# class_label can be obtained by get_class
get_class2design <- function(class_label){
  design <- model.matrix(~0+class_label);colnames(design) <- unique(class_label);
  rownames(design) <- names(class_label)
  return(design)
  #design.mat <-as.data.frame(matrix(0, nrow = length(class_label), ncol = length(unique(class_label))))
  #rownames(design.mat) <- names(class_label) ## sample
  #colnames(design.mat) <- unique(class_label)
  #for(i in 1:length(class_label)){design.mat[names(class_label)[i],class_label[i]]<-1;}
  #return(design.mat)
}
## getDE score by limma
getDE.limma.2G <- function(eset=NULL, G1=NULL, G0=NULL,G1_name=NULL,G0_name=NULL,verbose=TRUE,random_effect=NULL) {
  if(verbose==TRUE){
    print(sprintf('G1:%s', paste(G1, collapse = ';')))
    print(sprintf('G0:%s', paste(G0, collapse = ';')))
  }
  all_samples <- colnames(exprs(eset))
  use_samples <- intersect(c(G0, G1), all_samples)
  phe <- as.data.frame(pData(eset)[use_samples, ]);
  rownames(phe) <- use_samples
  new_eset <- generate.eset(exp_mat=exprs(eset)[, use_samples],phenotype_info=phe, feature_info=fData(eset), 'test')
  new_mat  <- exprs(new_eset)
  ##
  design.mat <-as.data.frame(matrix(NA, nrow = length(use_samples), ncol = 1))
  rownames(design.mat) <-use_samples
  colnames(design.mat) <- 'group'
  design.mat[intersect(G0, use_samples), 'group'] <- 'G0'
  design.mat[intersect(G1, use_samples), 'group'] <- 'G1'
  #  design <- model.matrix( ~ group + 0, design.mat)
  group <- factor(design.mat$group)
  design <- model.matrix(~0+group);
  colnames(design) <- levels(group); rownames(design) <- colnames(new_mat)

  if(is.null(random_effect)==TRUE){
    fit <- lmFit(new_mat,design)
  }else{
    corfit <- duplicateCorrelation(eset,design,block=random_effect)
    fit <- lmFit(new_mat,design,block=random_effect,correlation=corfit$consensus)
  }
  contrasts <- makeContrasts(G1-G0,levels=design)
  fit2 <- contrasts.fit(fit,contrasts=contrasts)
  fit2 <- eBayes(fit2,trend=TRUE)
  #summary(decideTests(fit2, method="global"))
  ##
  tT <- topTable(fit2, adjust = "fdr", number = Inf,coef=1)
  tT <- tT[order(tT$P.Value, decreasing = FALSE), ]
  tT <- cbind(ID=rownames(tT),tT,stringsAsFactors=FALSE)
  exp_G1 <- rowMeans(new_mat[rownames(tT),G1])
  exp_G0 <- rowMeans(new_mat[rownames(tT),G0])
  w1 <- which(tT$P.Value<=0);
  if(length(w1)>0) tT$P.Value[w1] <- .Machine$double.xmin;
  z_val <- sapply(tT$P.Value*sign(tT$logFC),function(x)combinePvalVector(x,twosided = TRUE)[1])
  if(is.null(random_effect)==TRUE){
    tT <- cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  }else{
    tT <- cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1,'Ave.G0_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G0'],
                'Ave.G1_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G1'])
  }
  if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  return(tT)
}
## combine p-values
#combine Pvalues by Fisher or Stouffer's method
#input:
#dat, a dataframe containing pvalue.cols
#sign.cols: column labels indicating sign of stat, if NULL, use absolute values to calculate overall stat for Stouffer's method
#combine one pvalues vector, the signed version
# For Fisher's method: if twosided, pvalues are transformed into single
combinePvalVector <-
  function(pvals,
           method = c('Stouffer', 'Fisher'),
           signed = TRUE,
           twosided = TRUE) {
    #remove NA pvalues
    pvals <- pvals[!is.na(pvals) & !is.null(pvals)]
    pvals[which(pvals==0)] <- .Machine$double.xmin*sign(pvals)
    if (sum(is.na(pvals)) >= 1) {
      stat <- NA
      pval <- NA
    } else{
      if (twosided & (sum(pvals > 1 | pvals < -1) >= 1))
        stop('pvalues must between 0 and 1!\n')
      if (!twosided & (sum(pvals > 0.5 | pvals < -0.5) >= 1))
        stop('One-sided pvalues must between 0 and 0.5!\n')

      if (missing(method))
        method <- 'Stouffer'

      if (!signed) {
        pvals <- abs(pvals)
      }

      signs <- sign(pvals)
      signs[signs == 0] <- 1

      if (grepl('Fisher', method, ignore.case = TRUE)) {
        if (twosided & signed) {
          neg.pvals <- pos.pvals <- abs(pvals) / 2
          pos.pvals[signs < 0] <- 1 - pos.pvals[signs < 0]
          neg.pvals[signs > 0] <- 1 - neg.pvals[signs > 0]
        } else{
          neg.pvals <- pos.pvals <- abs(pvals)
        }
        pvals <-
          c(1, -1) * c(
            pchisq(
              -2 * sum(log(as.numeric(pos.pvals))),
              df = 2 * length(pvals),
              lower.tail = FALSE
            ) / 2,
            pchisq(
              -2 * sum(log(as.numeric(neg.pvals))),
              df = 2 * length(pvals),
              lower.tail = FALSE
            ) / 2
          )
        pval <- min(abs(pvals))[1]
        #if two pvals are equal, pick up the first one
        stat <-
          sign(pvals[abs(pvals) == pval])[1] * qnorm(pval, lower.tail = F)[1]
        pval <- 2 * pval
      }
      else if (grepl('Stou', method, ignore.case = TRUE)) {
        if (twosided) {
          zs <- signs * qnorm(abs(pvals) / 2, lower.tail = FALSE)
          stat <- sum(zs) / sqrt(length(zs))
          pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
        }
        else{
          zs <- signs * qnorm(abs(pvals), lower.tail = FALSE)
          stat <- sum(zs) / sqrt(length(zs))
          pval <- pnorm(abs(stat), lower.tail = FALSE)
        }
      }
      else{
        stop('Only \"Fisher\" or \"Stouffer\" method is supported!!!\n')
      }
    }
    return(c(`Z-statistics` = stat, `P.Value` = pval))
  }
########################### ms table related functions
## merge TF/SIG AC
merge.TF_SIG.AC <- function(TF_AC=NULL,SIG_AC=NULL){
  mat_TF <- exprs(TF_AC)
  mat_SIG <- exprs(SIG_AC)
  funcType <- c(rep('TF',nrow(mat_TF)),rep('SIG',nrow(mat_SIG)))
  rn <- c(rownames(mat_TF),rownames(mat_SIG))
  rn_label <- paste(rn,funcType,sep='_')
  mat_combine <- rbind(mat_TF,mat_SIG[,colnames(mat_TF)])
  rownames(mat_combine) <- rn_label
  eset_combine <- generate.eset(exp_mat=mat_combine,phenotype_info=pData(TF_AC)[colnames(mat_combine),],
                                feature_info=NULL,annotation_info='activity in net-dataset')
  return(eset_combine)
}
## merge TF/SIG network
merge.TF_SIG.network <- function(TF_network=NULL,SIG_network=NULL){
  s_TF  <- names(TF_network$target_list)
  s_SIG <- names(SIG_network$target_list)
  funcType <- c(rep('TF',length(s_TF)),rep('SIG',length(s_SIG)))
  rn <- c(s_TF,s_SIG)
  rn_label <- paste(rn,funcType,sep='_')
  target_list_combine <- c(TF_network$target_list,SIG_network$target_list)
  names(target_list_combine) <- rn_label
  n_TF <- TF_network$network_dat
  n_TF$source <- paste(n_TF$source,'TF',sep='_')
  n_SIG <- SIG_network$network_dat
  n_SIG$source <- paste(n_SIG$source,'SIG',sep='_')
  net_dat <- rbind(n_TF,n_SIG)
  igraph_obj <- graph_from_data_frame(net_dat[,c('source','target','MI')],directed=TRUE)
  return(list(network_dat=net_dat,target_list=target_list_combine,igraph_obj=igraph_obj))
}
## generate ms table
# master table : gene-label, geneSymbol, network id-label,funcType, avgAC,z.DA, avgExp,z.DE, other info
# display_col do not need to include ave related columns
# if not merge before
generate.masterTable.TF_SIG <- function(use_comp=NULL,DE=NULL,DA_TF=NULL,DA_SIG=NULL,TF_network=NULL,SIG_network=NULL,main_id_type=NULL,
                              tf_sigs=tf_sigs,z_col=NULL,display_col=NULL){
  if(is.null(use_comp)){message('No input for use_comp, please check and re-try!');return(FALSE)}
  if(is.null(main_id_type)){message('No input for main_id_type, please check and re-try!');return(FALSE)}
  if(is.null(DE)){message('No input for DE, please check and re-try!');return(FALSE)}
  if(is.null(DA_TF)){message('No input for DA_TF, please check and re-try!');return(FALSE)}
  if(is.null(DA_SIG)){message('No input for DA_SIG, please check and re-try!');return(FALSE)}
  if(is.null(TF_network)){message('No input for TF_network, please check and re-try!');return(FALSE)}
  if(is.null(SIG_network)){message('No input for SIG_network, please check and re-try!');return(FALSE)}
  if(is.null(tf_sigs)){message('tf_sigs not loaded, please check and re-try!');return(FALSE)}
  if(is.null(z_col)){message('No input for z_col, please check and re-try!');return(FALSE)}
  if(length(setdiff(use_comp,names(DE)))>0){message('%s not calculated, please check and re-try!',setdiff(use_comp,names(DE)));return(FALSE)}
  if(length(setdiff(use_comp,names(DA_TF)))>0){message('%s not calculated, please check and re-try!',setdiff(use_comp,names(DA_TF)));return(FALSE)}
  if(length(setdiff(use_comp,names(DA_SIG)))>0){message('%s not calculated, please check and re-try!',setdiff(use_comp,names(DA_SIG)));return(FALSE)}
  funcType <- c(rep('TF',nrow(DA_TF[[1]])),rep('SIG',nrow(DA_SIG[[1]])))
  rn <- c(rownames(DA_TF[[1]]),rownames(DA_SIG[[1]]))
  rn_label <- paste(rn,funcType,sep='_')
  tf_size <- unlist(lapply(TF_network[rownames(DA_TF[[1]])],nrow))
  sig_size <- unlist(lapply(SIG_network[rownames(DA_SIG[[1]])],nrow))
  use_size <- c(tf_size,sig_size)
  # id issue
  current_id <- names(tf_sigs$tf)[-1]
  use_info <- unique(rbind(tf_sigs$tf$info,tf_sigs$sig$info))
  if(main_id_type %in% current_id){
    use_info <- use_info[which(use_info[,main_id_type] %in% rn),]
  }else{
    transfer_tab <- get_IDtransfer(from_type=main_id_type,to_type=current_id[1],use_genes=rn)
    use_info <- merge(use_info,transfer_tab,by.x=current_id[1],by.y=current_id[1])
    use_info <- use_info[which(use_info[,main_id_type] %in% rn),]
  }
  use_info <- unique(use_info)
  # merge info
  tmp1 <- aggregate(use_info,list(use_info[,main_id_type]),function(x){
    x1 <- x[which(x!='')]
    x1 <- x1[which(is.na(x1)==FALSE)]
    paste(sort(unique(x1)),collapse=';')
  })
  tmp1 <- tmp1[,-1]; rownames(tmp1) <- tmp1[,main_id_type]
  geneSymbol <- tmp1[rn,'external_gene_name'] ## this column for function enrichment
  if('external_transcript_name' %in% colnames(tmp1)){ ## this column for display
    gene_label <- paste(tmp1[rn,'external_transcript_name'],funcType,sep = '_')
  }else{
    gene_label <-paste(tmp1[rn,'external_gene_name'],funcType,sep = '_')
  }
  #
  label_info <- data.frame('gene_label'=gene_label,'geneSymbol'=geneSymbol,
                           'originalID'=rn,'originalID_label'=rn_label,'funcType'=funcType,'Size'=use_size,stringsAsFactors=FALSE)
  add_info <- tmp1[rn,]
  #
  combine_info <- lapply(use_comp,function(x){
    rn <- c(rownames(DA_TF[[x]]),rownames(DA_SIG[[x]]))
    rn_label <- paste(rn,funcType,sep='_')
    avg_col <- colnames(DA_TF[[x]])[grep('^Ave',colnames(DA_TF[[x]]))]
    DA_info <- rbind(DA_TF[[x]][,c(z_col,avg_col,setdiff(display_col,c(z_col,avg_col)))],
                     DA_SIG[[x]][,c(z_col,avg_col,setdiff(display_col,c(z_col,avg_col)))])
    avg_col <- colnames(DE[[x]])[grep('^Ave',colnames(DA_TF[[x]]))]
    DE_info <- as.data.frame(DE[[x]])[rn,c(z_col,avg_col,setdiff(display_col,c(z_col,avg_col)))]
    colnames(DA_info) <- paste0(colnames(DA_info),'.',x,'_DA')
    colnames(DE_info) <- paste0(colnames(DE_info),'.',x,'_DE')
    colnames(DA_info)[1] <- paste0('Z.',x,'_DA')
    colnames(DE_info)[1] <- paste0('Z.',x,'_DE')
    out <- cbind(DA_info,DE_info,stringsAsFactors=FALSE)
    rownames(out) <- rn_label
    out
  })
  combine_info_DA <- do.call(cbind,lapply(combine_info,function(x)x[rn_label,grep('_DA$',colnames(x))]))
  combine_info_DE <- do.call(cbind,lapply(combine_info,function(x)x[rn_label,grep('_DE$',colnames(x))]))
  # put them together
  ms_tab <- cbind(label_info,combine_info_DA,combine_info_DE,add_info)
  return(ms_tab)
}
# if the network + DA is merged before
generate.masterTable <- function(use_comp=NULL,DE=NULL,DA=NULL,network=NULL,main_id_type=NULL,
                                        tf_sigs=tf_sigs,z_col=NULL,display_col=NULL){
  if(is.null(use_comp)){message('No input for use_comp, please check and re-try!');return(FALSE)}
  if(is.null(main_id_type)){message('No input for main_id_type, please check and re-try!');return(FALSE)}
  if(is.null(DE)){message('No input for DE, please check and re-try!');return(FALSE)}
  if(is.null(DA)){message('No input for DA, please check and re-try!');return(FALSE)}
  if(is.null(network)){message('No input for network, please check and re-try!');return(FALSE)}
  if(is.null(tf_sigs)){message('tf_sigs not loaded, please check and re-try!');return(FALSE)}
  if(is.null(z_col)){message('No input for z_col, please check and re-try!');return(FALSE)}
  if(length(setdiff(use_comp,names(DE)))>0){message('%s not calculated, please check and re-try!',setdiff(use_comp,names(DE)));return(FALSE)}
  if(length(setdiff(use_comp,names(DA)))>0){message('%s not calculated, please check and re-try!',setdiff(use_comp,names(DA)));return(FALSE)}
  # get original ID
  ori_rn <- rownames(DA[[1]])
  w1 <- grep('(.*)_TF',ori_rn); w2 <- grep('(.*)_SIG',ori_rn)
  funcType <- rep(NA,length.out=length(ori_rn));rn <- funcType
  funcType[w1] <- 'TF'; funcType[w2] <- 'SIG';
  rn[w1] <- gsub('(.*)_TF',"\\1",ori_rn[w1]);rn[w2] <- gsub('(.*)_SIG',"\\1",ori_rn[w2]);
  rn_label <- ori_rn
  use_size <- unlist(lapply(network[rownames(DA[[1]])],nrow))
  # id issue
  current_id <- names(tf_sigs$tf)[-1]
  use_info <- unique(rbind(tf_sigs$tf$info,tf_sigs$sig$info))
  if(main_id_type %in% current_id){
    use_info <- use_info[which(use_info[,main_id_type] %in% rn),]
  }else{
    transfer_tab <- get_IDtransfer(from_type=main_id_type,to_type=current_id[1],use_genes=rn)
    use_info <- merge(use_info,transfer_tab,by.x=current_id[1],by.y=current_id[1])
    use_info <- use_info[which(use_info[,main_id_type] %in% rn),]
  }
  use_info <- unique(use_info)
  # merge info
  tmp1 <- aggregate(use_info,list(use_info[,main_id_type]),function(x){
    x1 <- x[which(x!='')]
    x1 <- x1[which(is.na(x1)==FALSE)]
    paste(sort(unique(x1)),collapse=';')
  })
  tmp1 <- tmp1[,-1]; rownames(tmp1) <- tmp1[,main_id_type]
  geneSymbol <- tmp1[rn,'external_gene_name'] ## this column for function enrichment
  if('external_transcript_name' %in% colnames(tmp1)){ ## this column for display
    gene_label <- paste(tmp1[rn,'external_transcript_name'],funcType,sep = '_')
  }else{
    gene_label <-paste(tmp1[rn,'external_gene_name'],funcType,sep = '_')
  }
  #
  label_info <- data.frame('gene_label'=gene_label,'geneSymbol'=geneSymbol,
                           'originalID'=rn,'originalID_label'=rn_label,'funcType'=funcType,'Size'=use_size,stringsAsFactors=FALSE)
  add_info <- tmp1[rn,]
  #
  combine_info <- lapply(use_comp,function(x){
    DA[[x]] <- DA[[x]][rn_label,]
    DE[[x]] <- as.data.frame(DE[[x]])[rn,]
    avg_col <- colnames(DA[[x]])[grep('^Ave',colnames(DA[[x]]))]
    uc <- c(z_col,avg_col,setdiff(display_col,c(z_col,avg_col))); uc <- intersect(uc,colnames(DA[[x]]))
    DA_info <- DA[[x]][rn_label,uc]
    avg_col <- colnames(DE[[x]])[grep('^Ave',colnames(DE[[x]]))]
    uc <- c(z_col,avg_col,setdiff(display_col,c(z_col,avg_col))); uc <- intersect(uc,colnames(DE[[x]]))
    DE_info <- as.data.frame(DE[[x]])[rn,uc]
    colnames(DA_info) <- paste0(colnames(DA_info),'.',x,'_DA')
    colnames(DE_info) <- paste0(colnames(DE_info),'.',x,'_DE')
    colnames(DA_info)[1] <- paste0('Z.',x,'_DA')
    colnames(DE_info)[1] <- paste0('Z.',x,'_DE')
    out <- cbind(DA_info,DE_info,stringsAsFactors=FALSE)
    rownames(out) <- rn_label
    out
  })
  combine_info_DA <- do.call(cbind,lapply(combine_info,function(x)x[rn_label,grep('_DA$',colnames(x))]))
  combine_info_DE <- do.call(cbind,lapply(combine_info,function(x)x[rn_label,grep('_DE$',colnames(x))]))
  # put them together
  ms_tab <- cbind(label_info,combine_info_DA,combine_info_DE,add_info)
  return(ms_tab)
}
## output ms table to excel files
######## output to excel
out2excel <- function(all_ms_tab,out.xls,mark_gene=NULL,mark_col,workbook_name='ms_tab',z_column=NULL,sig_thre=1.64){
  wb <- createWorkbook(workbook_name)
  if(!'list' %in% class(all_ms_tab)){
    all_ms_tab <- list('Sheet1'=as.data.frame(all_ms_tab))
  }
  z_column_index <- 'defined'
  if(is.null(z_column)==TRUE) z_column_index <- 'auto'
  i <- 0
  headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#4F81BD", border="TopBottom", borderColour = "#4F81BD",wrapText=TRUE) ## style for header line
  for(sheetname in names(all_ms_tab)){ ## list, each item create one sheet
    i <- i +1
    d <- as.data.frame(all_ms_tab[[sheetname]])
    if(z_column_index=='auto') use_z_column <- colnames(d)[grep('Z\\.',colnames(d),ignore.case = TRUE)] else use_z_column <- z_column
    addWorksheet(wb,sheetName=sheetname)
    writeData(wb,sheet = i,d)
    addStyle(wb, sheet = i, headerStyle, rows = 1, cols = 1:ncol(d), gridExpand = TRUE) ## add header style
    all_c <- list()
    for(z_col in use_z_column){ ## find colnames with z. (only applied to pipeline excel)
      j <- which(colnames(d)==z_col)
      z1 <- d[,j]; c1 <- z2col(z1,sig_thre=sig_thre)
      all_c[[as.character(j)]] <- c1
    }
    mat_c <- do.call(cbind,all_c)
    uni_c <- unique(unlist(all_c))
    for(r in uni_c){
      w1 <- which(mat_c==r)
      nn <- nrow(mat_c)
      rr <- w1%%nn+1
      cc <- w1%/%nn+1
      cc[which(rr==1)] <- cc[which(rr==1)]-1
      rr[which(rr==1)] <- nn+1
      addStyle(wb, sheet = i, createStyle(fgFill=r), rows =rr, cols = as.numeric(names(all_c)[cc]))
    }
    if(is.null(mark_gene)==FALSE){
      for(k in names(mark_gene)){
        addStyle(wb, sheet = i, createStyle(fgFill=mark_col[k]), rows =which(toupper(d[,which(colnames(d)=='geneSymbol')]) %in% toupper(mark_gene[[k]]))+1,
                 cols = which(colnames(d)=='geneSymbol')) ## find column with geneSymbol and mark color
      }
    }
    w1 <- which(gsub("\\s","",as.matrix(d))=='TRUE')
    nn <- nrow(d)
    rr <- w1%%nn+1
    cc <- w1%/%nn+1
    cc[which(rr==1)] <- cc[which(rr==1)]-1
    rr[which(rr==1)] <- nn+1
    addStyle(wb, sheet = i, createStyle(fontColour='#FF0000'), rows =rr, cols = as.numeric(cc)) ## find column with TRUE/FALSE and mark with color
  }
  saveWorkbook(wb, out.xls, overwrite = TRUE)
  ##
}
########################### gene set related functions
# geneset info
# ## h:hallmark genesets; c1:positional gene sets; c2: curated gene set;c3: motif; c4: computational;c5:GO;c6:oncogenic;c7:immune
# CGP: chemical and genetic perturbations
# CP: Canonical pathways
# CP:BIOCARTA: BioCarta gene sets
# CP:KEGG: KEGG gene sets
# CP:REACTOME: Reactome gene sets
# MIR: microRNA targets
# TFT: transcription factor targets
# CGN: cancer gene neighborhoods
# CM: cancer modules
gs.preload <- function(use_spe='Homo sapiens',update=0){ ## only support geneSymbol (because pipeline-generated master table will contain geneSymbol column)
  if(is.null(main.dir)){message('main.dir required, please check and re-try !');return(FALSE)}
  db.dir <- sprintf("%s/db/",main.dir)
  all_spe <- msigdbr_show_species()
  if(!use_spe %in% all_spe){
    message(sprintf('%s not in %s, please check and re-try !',use_spe,paste(all_spe,collapse=';')))
    return(FALSE)
  }
  use_spe1 <- gsub(' ','_',use_spe)
  out_file <- sprintf('%s/%s_gs2gene.RData',db.dir,use_spe1)
  if(file.exists(out_file)==FALSE | update==1){
    message('Begin generating all_gs2gene !')
    all_gs_info <-  msigdbr(species = use_spe) ## use msigdbr_show_species() to check possible available species
    # for gs_cat
    all_gs_cat <- unique(all_gs_info$gs_cat)
    all_gs2gene_1 <- lapply(all_gs_cat,function(x){
      x1 <- all_gs_info[which(all_gs_info$gs_cat==x),]
      all_gs <- unique(x1$gs_name)
      x2 <- lapply(all_gs, function(y){
        unique(x1$gene_symbol[which(x1$gs_name==y)])
      })
      names(x2) <- all_gs;x2
    })
    names(all_gs2gene_1) <- all_gs_cat
    all_gs_subcat <- setdiff(unique(all_gs_info$gs_subcat),'')
    all_gs2gene_2 <- lapply(all_gs_subcat,function(x){
      x1 <- all_gs_info[which(all_gs_info$gs_subcat==x),]
      all_gs <- unique(x1$gs_name)
      x2 <- lapply(all_gs, function(y){
        unique(x1$gene_symbol[which(x1$gs_name==y)])
      })
      names(x2) <- all_gs;x2
    })
    names(all_gs2gene_2) <- all_gs_subcat
    all_gs2gene <- c(all_gs2gene_1,all_gs2gene_2)
    #
    all_gs2gene <- all_gs2gene[sort(names(all_gs2gene))]
    gs_size <- unlist(lapply(all_gs2gene,length))
    # info for cat
    info_cat <- c('C1'='positional gene sets','C2'='curated gene set','C3'='motif','C4'='computational','C5'='GO','C6'='oncogenic','C7'='immune','H'='hallmark genesets')
    info_subcat <- c('CGP'='chemical and genetic perturbations','CP'='Canonical pathways','CP:BIOCARTA'='BioCarta gene sets','CP:KEGG'='KEGG gene sets',
                     'CP:REACTOME'='Reactome gene sets','MIR'='microRNA targets','TFT'='transcription factor targets','CGN'='cancer gene neighborhoods','CM'='cancer modules',
                     'BP'='Biological Process','MF'='Molecular Function','CC'='Cellular Component')
    cat_rel <- unique(as.data.frame(all_gs_info[,c('gs_cat','gs_subcat')]))
    all_gs2gene_info <- data.frame(cat_rel[,1],info_cat[cat_rel[,1]],gs_size[cat_rel[,1]],cat_rel[,2],info_subcat[cat_rel[,2]],gs_size[cat_rel[,2]],stringsAsFactors = FALSE)
    colnames(all_gs2gene_info) <- c('Category','Category_Info','Category_Size','Sub-Category','Sub-Category_Info','Sub-Category_Size')
    all_gs2gene_info <- all_gs2gene_info[order(all_gs2gene_info[,1]),]
    all_gs2gene_info[,c(1,2,4,5)] <- as.data.frame(apply(all_gs2gene_info[,c(1,2,4,5)],2,function(x){x[which(is.na(x)==TRUE)] <- '';x}),stringsAsFactors=FALSE)
    save(all_gs2gene,all_gs2gene_info,file=out_file)
  }
  load(out_file,.GlobalEnv)
  message('all_gs2gene loaded, you could see all_gs2gene_info to check the details !')
  return(TRUE)
}
##


########################### visualization functions
## simple functions to get info
get_obs_label <- function(phe_info,i){
  obs_label<-phe_info[,i];names(obs_label) <- rownames(phe_info);obs_label
}
get_int_group <- function(eset){
  phe <- pData(eset)
  feature_len <- apply(phe,2,function(x)length(unique(x)))
  intgroup <- colnames(phe)[which(feature_len>1 & feature_len<nrow(phe))]
  return(intgroup)
}
## get jaccard accuracy
get_jac <- function(pred_label, obs_label) {
  jac1 <- c()
  for (i in unique(pred_label)) {
    jac_index <- c()
    x1 <- names(pred_label)[which(pred_label == i)]
    for (j in unique(obs_label)) {
      x2 <- names(obs_label)[which(obs_label == j)]
      jac_index <-
        c(jac_index, length(intersect(x1, x2)) / length(union(x1, x2)))
    }
    jac1 <- c(jac1, max(jac_index) * length(x1))
  }
  jac1 <- sum(jac1) / length(pred_label)
  return(jac1)
}
## transfer z to color
z2col <- function(x,n_len=60,sig_thre=0.01,col_min_thre=0.01,col_max_thre=3,
                  blue_col=brewer.pal(9,'Set1')[2],red_col=brewer.pal(9,'Set1')[1]){ ## create vector for z-score, can change sig threshold
  x[which(is.na(x)==TRUE)] <- 0
  x[which(x==Inf)]<-  max(x[which(x!=Inf)])+1
  x[which(x==-Inf)]<- min(x[which(x!=-Inf)])-1
  if(col_min_thre<0) col_min_thre<-0.01
  if(col_max_thre<0) col_max_thre<-3
  #c1 <- brewer.pal(9,'Set1')
  c2 <- colorRampPalette(c(blue_col,'white',red_col))(n_len)
  r1 <- 1.05*max(abs(x)) ## -r1~r1
  if(r1 < col_max_thre){
    r1 <- col_max_thre
  }
  if(col_min_thre>r1){
    r2 <- seq(-r1,r1,length.out=n_len+1)
  }else{
    r21 <- seq(-r1,-col_min_thre,length.out=n_len/2)
    r22 <- seq(col_min_thre,r1,length.out=n_len/2)
    r2 <- c(r21,r22)
  }
  x1 <- cut(x,r2)
  names(c2) <- levels(x1)
  x2 <- c2[x1]
  x2[which(abs(x)<sig_thre)] <- 'white'
  x2
}
## get color info for class
# x:vector for class
# use_color : use color sets
get.class.color <- function(x,use_color=NULL,pre_define=NULL,use_inner=TRUE) {
  if(is.null(pre_define)==FALSE & is.null(names(pre_define))==TRUE){
    message('No class name for the color vector, please check and re-try !');return(FALSE);
  }
  if(use_inner==TRUE){ ## use inner pre_defined !!!
    pre_define <- c('blue', 'red', 'yellow', 'green','yellow', 'green') ## pre-defined colors for MB
    names(pre_define) <- c('WNT', 'SHH', 'Group3', 'Group4','GroupC', 'GroupD') ##pre-defined color name for MB
  }
  if(is.null(use_color)==TRUE){
    if (length(intersect(x, names(pre_define))) == 0){
      use_color <- brewer.pal(12, 'Set3')
    }else{
      use_color <- brewer.pal(12, 'Set3')[c(1,3,5,7,9,10,11)]
    }
  }
  if (length(intersect(x, names(pre_define))) == 0) {
    cc2 <- colorRampPalette(use_color)(length(unique(x)))
    names(cc2) <- unique(x)
    cc2 <- cc2[x]
  } else{
    x1 <- unique(x)
    x2 <- setdiff(x1, names(pre_define))
    cc1 <- NULL
    if (length(x2) > 0) {
      cc1 <- colorRampPalette(use_color)(length(x2))
      names(cc1) <- x2
    }
    cc2 <- c(pre_define, cc1)
    cc2 <- cc2[x]
  }
  return(cc2)
}
## get color box text ## refer from web
# https://stackoverflow.com/questions/45366243/text-labels-with-background-colour-in-r
boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA,
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5,
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font')){

  ## The Character expansion factro to be used:
  theCex <- graphics::par('cex')*cex
  ## Is y provided:
  if (missing(y)) y <- x
  ## Recycle coords if necessary:
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }
  }
  ## Width and height of text
  textHeight <- graphics::strheight(labels, cex = theCex, font = font)
  textWidth <- graphics::strwidth(labels, cex = theCex, font = font)
  ## Width of one character:
  charWidth <- graphics::strwidth("e", cex = theCex, font = font)
  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)
    }
  } else {
    adj <- c(0.5, 0.5)
  }

  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }
  } else {
    offsetVec <- c(0, 0)
  }
  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }
  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]

  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth
  graphics::rect(xleft = xMid - rectWidth/2,ybottom = yMid - rectHeight/2,
                 xright = xMid + rectWidth/2,ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg,xpd=TRUE)
  ## Place the text:
  graphics::text(xMid, yMid, labels, col = col.text, cex = theCex, font = font,adj = c(0.5, 0.5),xpd=TRUE)
  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,yMid + rectHeight/2))
  }
}
## draw 2-dimension plot
draw.2D <- function(X,Y,class_label,xlab='PC1',ylab='PC2',label_cex=0.8,main='',point_cex=1){
  if(is.null(class_label)==TRUE){
    message('No class_label, please check and re-try !');return(FALSE);
  }
  if(length(X)!=length(Y)){
    message('Input two dimension vector with different length, please check and re-try !');return(FALSE);
  }
  if(length(X)!=length(class_label)){
    message('Input dimension vector has different length with the class_label, please check and re-try !');return(FALSE);
  }
  par(mar = c(4, 4, 4, 15))
  class_label <- as.character(class_label)
  cls_cc <- get.class.color(class_label) ## get color for each label
  plot(Y ~ X,pch = 16,cex = point_cex,col = cls_cc,main=main,xlab=xlab,ylab=ylab)
  legend(par()$usr[2],par()$usr[4],sort(unique(class_label)),fill = cls_cc[sort(unique(class_label))],
         horiz = FALSE,xpd = TRUE,border = NA,bty = 'n',cex=label_cex)
  return(TRUE)
}
## another way to draw 2-dimension plot
# X,Y
# class_label: class label vector
draw.2D.ellipse <- function(X,Y,class_label,xlab='PC1',ylab='PC2',label_cex=0.8,main='',point_cex=1){
  if(is.null(class_label)==TRUE){
    message('No class_label, please check and re-try !');return(FALSE);
  }
  if(length(X)!=length(Y)){
    message('Input two dimension vector with different length, please check and re-try !');return(FALSE);
  }
  if(length(X)!=length(class_label)){
    message('Input dimension vector has different length with the class_label, please check and re-try !');return(FALSE);
  }
  par(mar=c(5,5,5,5))
  get_transparent <- function(x,alpha=0.1){
    rgb(t(col2rgb(x)/255),alpha=alpha)
  }
  class_label <- as.character(class_label)
  cls_cc <- get.class.color(class_label) ## get color for each color
  par(mar=c(5,8,5,8))
  plot(Y~X,pch = 19,cex = point_cex,xlim=c(min(X)-IQR(X)/3,IQR(X)/3+max(X)),
    col = cls_cc,bty='n',bg='white',xlab=xlab,ylab=ylab,main=main)
  ## add ellipse
  for(i in unique(class_label)){
    w0 <- which(class_label==i)
    w1 <- cbind(X[w0],Y[w0])
    if(is.null(dim(w1))){
      w1 <- t(as.matrix(w1))
    }
    c1 <- get_transparent(cls_cc[which(class_label==i)][1])
    m1 <- colMeans(w1)
    d1 <- unlist(lapply(1:nrow(w1),function(x)sqrt(sum((w1[x,]-m1)^2))))
    if(length(d1)==1){
      draw.ellipse(m1[1],m1[2],a=(par()$usr[2]-par()$usr[1])/30,b=(par()$usr[2]-par()$usr[1])/30,col=c1,border=NA,xpd=TRUE)
      text(m1[1],m1[2],i,xpd=TRUE,adj=0,cex=0.8)
    }
    if(length(d1)==2){
      draw.ellipse(m1[1],m1[2],a=d1[1],b=d1[1],col=c1,border=NA,xpd=TRUE)
      text(m1[1],m1[2],i,xpd=TRUE,adj=0,cex=0.8)
    }
    if(length(d1)>=3){
      w2 <- order(d1,decreasing=TRUE)
      d1 <- d1[w2]
      w1 <- w1[w2,]
      d3 <- do.call(rbind,lapply(1:nrow(w1),function(x)w1[x,]-m1))
      a <- sqrt(d3[1,1]^2+d3[1,2]^2)
      angle <- 360/(2*pi)*atan(d3[1,2]/d3[1,1])
      beta  <- 360/(2*pi)*atan(d3[,2]/d3[,1])
      dx <- d1*cos(2*pi*(beta-angle)/360)
      dy <- d1*sin(2*pi*(beta-angle)/360)
      #
      b <- sqrt((dy^2)/(1-dx^2/a^2))
      b <- max(b[-1])
      draw.ellipse(m1[1],m1[2],a=a*1.05,b=b*1.05,col=c1,border=NA,angle=360/(2*pi)*atan(d3[1,2]/d3[1,1]),xpd=TRUE)
      #
      s1 <- sample(1:nrow(w1),1)
      d1 <- (par()$usr[2]-par()$usr[1])/30
      m2 <- w1[s1,]
      if(m2[1]>m1[1]){
        boxtext(m2[1]+d1,m2[2],labels=i,adj=0,cex=label_cex,col.bg=c1)
        segments(x0=w1[s1,1],y0=w1[s1,2],x1=m2[1]+d1,y1=m2[2],col='dark grey')
      } else{
        boxtext(m2[1]-d1,m2[2],labels=i,adj=1,cex=label_cex,col.bg=c1)
        segments(x0=w1[s1,1],y0=w1[s1,2],x1=m2[1]-d1,y1=m2[2],col='dark grey')
      }
    }
  }
  return(TRUE)
}
## QC plot for eset
# eset, input
# outdir: output directory
# intgroup, interested group
plot.eset.QC <- function(eset,outdir = '.',do.logtransform = FALSE,intgroup=NULL,prefix = 'afterQC_',choose_plot=c('heatmap','pca','density','meansd')) {
    if (!file.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
      message(paste0("The output directory: \"", outdir, "\" is created!"))
    }else
      message(paste0("The output will overwrite the files in directory: \"",outdir,"\""))
    if(is.null(intgroup)){
      intgroup <- get_int_group(eset)
    }
    if(length(intgroup)==0){
      message('No intgroup, please check and re-try!');return(FALSE)
    }
    message('Preparing the data...')
    x  <- prepdata(eset, do.logtransform = do.logtransform, intgroup = intgroup)

    ## pca
    if('pca' %in% choose_plot){
      pca <- prcomp(t(na.omit(x$M)))
      x$key$rect$col <- get.class.color(x$key$text[[1]])
      fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'pca'))
      pdf(fp, width = 14, height = 9)
      for (i in 1:length(intgroup)) {
        class_label <- x$pData[[x$intgroup[i]]]
        class_label[which(is.na(class_label)==TRUE)] <- 'NA'
        draw.2D(as.data.frame(pca$x)$PC1,as.data.frame(pca$x)$PC2,class_label=class_label,xlab='PC1',ylab='PC2',
                label_cex=0.8,main=sprintf('PCA/Kmeans plot for %s',intgroup[i]))
        draw.2D.ellipse(as.data.frame(pca$x)$PC1,as.data.frame(pca$x)$PC2,class_label=class_label,xlab='PC1',ylab='PC2',
                        label_cex=0.8,main=sprintf('PCA/Kmeans plot for %s',intgroup[i]))
      }
      dev.off()
      message('Finish PCA plot !')
    }

    ## heatmap
    if('heatmap' %in% choose_plot){
      fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'heatmap'))
      pdf(fp, width = 12, height = 9)
      par(mar = c(6, 6, 6, 6))
      m <- dist2(x$M)
      dend <- as.dendrogram(hclust(as.dist(m), method = "single"))
      ord <- order.dendrogram(dend)
      m <- m[ord, ord]
      for(i in 1:length(intgroup)){
        class_label <- x$pData[rownames(m),x$intgroup[i]]
        class_label[which(is.na(class_label)==TRUE)] <- 'NA'
        cls_cc <- get.class.color(class_label)
        heatmap(
          m,Colv = NA,Rowv = NA,labRow = NA,labCol = NA,margin = c(10, 10),scale = 'none',
          col = colorRampPalette(c('blue', 'yellow'))(100),
          RowSideColors = cls_cc,
          ColSideColors = cls_cc
        )
        legend(0,0,legend=unique(class_label),
               fill = cls_cc[unique(class_label)],
               xpd = TRUE,border = NA,bty = 'n',horiz = TRUE)
      }
      dev.off()
      message('Finish Heatmap plot !')
    }


    ## meansd
    if('meansd' %in% choose_plot){
      fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'meansd'))
      pdf(fp, width = 12, height = 9)
      meanSdPlot(eset)
      dev.off()
      message('Finish MeanSD plot !')
    }

    ## density
    if('density' %in% choose_plot){
      fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'density'))
      pdf(fp, width = 12, height = 9)
      for(i in 1:length(intgroup)){
        all_dens <- list()
        for (j in 1:ncol(x$M)) {
          all_dens[[j]] <- density(x$M[,j],na.rm=TRUE)
        }
        plot(1,col = 'white',xlim=c(min(unlist(lapply(all_dens,function(x)min(x$x)))),max(unlist(lapply(all_dens,function(x)max(x$x))))),
             type = 'l',xlab = '',ylab='Density',main = sprintf('Density plot for %s',intgroup[i]),
             ylim=c(min(unlist(lapply(all_dens,function(x)min(x$y)))),max(unlist(lapply(all_dens,function(x)max(x$y))))))
        class_label <- x$pData[,x$intgroup[i]]
        cls_cc <- get.class.color(class_label)
        for (j in 1:ncol(x$M)) {
          lines(all_dens[[j]], col = cls_cc[j])
        }
        legend('topright',legend=unique(class_label),
               fill = cls_cc[unique(class_label)],
               xpd = TRUE,border = NA,bty = 'n',horiz = FALSE)
      }
      dev.off()
      message('Finish Density plot !')
    }

    return(TRUE)
}
## pca+kmeans in 2D
# plot_type: 2D or 2D.ellipse
# obs_label, the value should be the class label with names equal to sample names
plot.2D.pca.kmeans <- function(mat=NULL,all_k=NULL,obs_label=NULL,legend_pos = 'topleft',legend_cex = 0.8,plot_type='2D.ellipse',point_cex=1){
  if(is.null(mat)==TRUE){
    message('Please input mat, check and re-try !');return(FALSE)
  }
  if(is.null(obs_label)==TRUE){
    message('Please input obs_label, check and re-try !');return(FALSE)
  }
  if(is.null(all_k)==TRUE){
    all_k <- 2:min(length(obs_label)-1,2*length(unique(obs_label)))
  }
  if(length(setdiff(all_k,2:length(obs_label)))>0){
    message('some value in all_k exceed the maximum sample size, check and re-try !');return(FALSE);
  }
  pc <- prcomp(t(mat))$x
  all_jac <- list()
  all_k_res <- list()
  for(k in all_k){
    tmp_k <- list()
    for(i in 1:10){
      tmp_k[[i]] <- kmeans(pc,centers=as.numeric(k))
    }
    pred_label <- lapply(tmp_k,function(x)x$cluster)
    jac <- unlist(lapply(pred_label,function(x){get_jac(x, obs_label)}))
    top_i <- which.max(jac)
    all_k_res[[k]] <- tmp_k[[top_i]]
  }
  for(k in all_k){
    pred_label <- all_k_res[[k]]$cluster
    jac <- get_jac(pred_label, obs_label)
    all_jac[[as.character(k)]] <- signif(jac,4)
  }
  print(all_jac)
  use_k <- all_k[which.max(all_jac)]
  message(sprintf('Best Jaccard Accuracy occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  pred_label <- all_k_res[[use_k]]$cluster
  d1 <- data.frame(id=colnames(mat),X=pc[,1],Y=pc[,2],label=pred_label,stringsAsFactors=FALSE)
  layout(t(matrix(1:2)))
  if(plot_type=='2D.ellipse'){
    draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='PC1',ylab='PC2',label_cex=0.8,point_cex=point_cex)
    draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,xlab='PC1',ylab='PC2',label_cex=0.8,point_cex=point_cex)
  }
  if(plot_type=='2D'){
    draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='PC1',ylab='PC2',label_cex=0.8,point_cex=point_cex)
    draw.2D(d1$X,d1$Y,class_label=d1$label,xlab='PC1',ylab='PC2',label_cex=0.8,point_cex=point_cex)
  }
  layout(1);
  return(pred_label)
}
## pca+kmeans in 3D
plot.3D.pca.kmeans <- function(mat=NULL,all_k=NULL,obs_label=NULL,legend_pos = 'topleft',legend_ncol = 1,legend_cex = 0.8){
  if(is.null(mat)==TRUE){
    message('Please input mat, check and re-try !');return(FALSE)
  }
  if(is.null(obs_label)==TRUE){
    message('Please input obs_label, check and re-try !');return(FALSE)
  }
  if(is.null(all_k)==TRUE){
    all_k <- 2:min(length(obs_label)-1,2*length(unique(obs_label)))
  }
  if(length(setdiff(all_k,2:length(obs_label)))>0){
    message('some value in all_k exceed the maximum sample size, check and re-try !');return(FALSE);
  }
  pc <- prcomp(t(mat))$x
  all_jac <- list()
  all_k_res <- list()
  for(k in all_k){
    tmp_k <- list()
    for(i in 1:10){
      tmp_k[[i]] <- kmeans(pc,centers=as.numeric(k))
    }
    pred_label <- lapply(tmp_k,function(x)x$cluster)
    jac <- unlist(lapply(pred_label,function(x){get_jac(x, obs_label)}))
    top_i <- which.max(jac)
    all_k_res[[k]] <- tmp_k[[top_i]]
    #all_k_res[[k]] <- kmeans(pc,centers=as.numeric(k))
  }
  for(k in all_k){
    pred_label <- all_k_res[[k]]$cluster
    jac <- get_jac(pred_label, obs_label)
    all_jac[[as.character(k)]] <- signif(jac,4)
  }
  print(all_jac)
  use_k <- all_k[which.max(all_jac)]
  message(sprintf('Best Jaccard Accuracy occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  pred_label <- all_k_res[[use_k]]$cluster
  d1 <- data.frame(id=colnames(mat),X=pc[,1],Y=pc[,2],Z=pc[,3],obs=obs_label[colnames(mat)],label=pred_label,stringsAsFactors=FALSE)
  par(mar=c(3,3,10,3))
  layout(t(matrix(1:2)))
  cls_cc <- get.class.color(unique(d1$obs))
  scatter3D(
    d1$X,
    d1$Y,
    d1$Z,
    pch = 16,
    xlab = 'PC1',
    ylab = 'PC2',
    zlab='PC3',bty='g',colvar=as.numeric(factor(d1$obs,levels=unique(d1$obs))),col=cls_cc,colkey = FALSE
  )
  legend(
    legend_pos,
    legend=unique(obs_label[d1$id]),
    fill = cls_cc[unique(obs_label[d1$id])],
    border = NA,
    bty = 'n',
    ncol = legend_ncol,
    cex = legend_cex
  )
  ##
  cls_cc <- get.class.color(as.character(unique(d1$label)))
  scatter3D(
    d1$X,
    d1$Y,
    d1$Z,
    pch = 16,
    xlab = 'PC1',
    ylab = 'PC2',zlab='PC3',bty='g',colvar=as.numeric(factor(d1$label,levels=unique(d1$label))),col=cls_cc,colkey = FALSE
  )
  legend(
    legend_pos,
    legend=as.character(sort(unique(d1$label))),
    fill = cls_cc[as.character(sort(unique(
      d1$label
    )))],
    border = NA,
    bty = 'n',
    ncol = legend_ncol,
    cex = legend_cex
  )
  return(pred_label)
}
## for umap visualization, use kmeans as cluster strategy
plot.2D.umap.kmeans <- function(mat=NULL,all_k=NULL,obs_label=NULL,legend_pos = 'topleft',legend_cex = 0.8,plot_type='2D.ellipse',point_cex=1){
  if(is.null(mat)==TRUE){
    message('Please input mat, check and re-try !');return(FALSE)
  }
  if(is.null(obs_label)==TRUE){
    message('Please input obs_label, check and re-try !');return(FALSE)
  }
  if(is.null(all_k)==TRUE){
    all_k <- 2:min(length(obs_label)-1,2*length(unique(obs_label)))
  }
  if(length(setdiff(all_k,2:length(obs_label)))>0){
    message('some value in all_k exceed the maximum sample size, check and re-try !');return(FALSE);
  }
  #pc <- prcomp(t(mat))$x
  use_mat_umap <- umap(t(mat))
  #
  all_jac <- list()
  all_k_res <- list()
  for(k in all_k){
    tmp_k <- list()
    for(i in 1:10){
      tmp_k[[i]] <- kmeans(use_mat_umap$layout,centers=as.numeric(k))
    }
    pred_label <- lapply(tmp_k,function(x)x$cluster)
    jac <- unlist(lapply(pred_label,function(x){get_jac(x, obs_label)}))
    top_i <- which.max(jac)
    all_k_res[[k]] <- tmp_k[[top_i]]
    #all_k_res[[k]] <- kmeans(use_mat_umap$layout,centers=as.numeric(k))
  }
  for(k in all_k){
    pred_label <- all_k_res[[k]]$cluster
    jac <- get_jac(pred_label, obs_label)
    all_jac[[as.character(k)]] <- signif(jac,4)
  }
  print(all_jac)
  use_k <- all_k[which.max(all_jac)]
  message(sprintf('Best Jaccard Accuracy occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  pred_label <- all_k_res[[use_k]]$cluster

  d1 <- data.frame(id=colnames(mat),X=use_mat_umap$layout[,1],Y=use_mat_umap$layout[,2],label=pred_label,stringsAsFactors=FALSE)
  layout(t(matrix(1:2)))
  if(plot_type=='2D.ellipse'){
    draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='',ylab='',label_cex=0.8,point_cex=point_cex)
    draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,xlab='',ylab='',label_cex=0.8,point_cex=point_cex)
  }
  if(plot_type=='2D'){
    draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='',ylab='',label_cex=0.8,point_cex=point_cex)
    draw.2D(d1$X,d1$Y,class_label=d1$label,xlab='',ylab='',label_cex=0.8,point_cex=point_cex)
  }
  layout(1);
  return(pred_label)
}
## functions for MICA result visualization
# visualization_type: tsne, umap,mds
plot.2D.MICA <- function(outdir=NULL,prjname=NULL,all_k=NULL,obs_label=NULL,legend_pos = 'topleft',legend_cex = 0.8,point_cex=1,plot_type='2D.ellipse',
                         visualization_type='tsne') {
  # choose best k
    all_jac <- get_jac_MICA(outdir=outdir, all_k=all_k, obs_label=obs_label, prjname = prjname)
    use_k <- all_k[which.max(all_jac)]
    message(sprintf('Best Jaccard Accuracy occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  #
    use_file <- sprintf('%s/scMINER_%s/scMINER_%s_MDS_%s/scMINER_MICA_out/%s.ggplot.txt',
                        outdir,prjname,prjname,use_k,prjname)
    d1 <- read.delim(use_file, stringsAsFactors = FALSE) ## get cluster results
    if(visualization_type=='umap' | visualization_type=='mds'){
      use_file <- sprintf('%s/scMINER_%s/scMINER_%s_MDS_%s/scMINER_MICA_out/%s_clust.h5',
                          outdir,prjname,prjname,use_k,prjname)
      fid <- H5Fopen(use_file)
      dist_mat <- fid$`mds`$block0_values
      if(visualization_type=='mds'){
        X <- fid$mds$block0_values[1,];Y <- fid$mds$block0_values[2,]
        d1$X <- X; d1$Y <- Y;
      }else{
        use_mat_umap <- umap(t(dist_mat))
        X <- use_mat_umap$layout[,1];Y <- use_mat_umap$layout[,2]
        d1$X <- X; d1$Y <- Y;
      }
      H5Fclose(fid)
    }
    layout(t(matrix(1:2)))
    if(plot_type=='2D.ellipse'){
      draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',label_cex=legend_cex,point_cex=point_cex)
      draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',label_cex=legend_cex,point_cex=point_cex)
    }
    if(plot_type=='2D'){
      draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',label_cex=legend_cex,point_cex=point_cex)
      draw.2D(d1$X,d1$Y,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',label_cex=legend_cex,point_cex=point_cex)
    }
    ## jaccard accuracy
    rownames(d1) <- d1$id
    pred_label <- d1[names(obs_label), ]$label
    names(pred_label) <- names(obs_label)
    jac <- get_jac(pred_label, obs_label)
    print(sprintf('Jaccard Accuracy:%s', jac))
    return(pred_label)
}
# get all jaccard accuracy for MICA
get_jac_MICA <- function(outdir, all_k, obs_label, prjname = NULL) {
    all_jac <- list()
    for (k in all_k) {
      use_file <-
        sprintf(
          '%s/scMINER_%s/scMINER_%s_MDS_%s/scMINER_MICA_out/%s.ggplot.txt',
          outdir,prjname,prjname,k,prjname
        )
      d1 <- read.delim(use_file, stringsAsFactors = FALSE)
      ## jaccard accuracy
      rownames(d1) <- d1$id
      pred_label <-
        d1[names(obs_label), ]$label
      names(pred_label) <- names(obs_label)
      jac <- get_jac(pred_label, obs_label)
      all_jac[[as.character(k)]] <- jac
      print(sprintf('Jaccard Accuracy for %d:%s', k, jac))
    }
    return(all_jac)
}
## draw volcano plot
# label_type: origin / distribute
draw.volcanoPlot <- function(dat=NULL,label_col=NULL,logFC_col=NULL,Pv_col=NULL,logFC_thre=1.5, Pv_thre=0.01,
                             xlab='log2 Fold Change',ylab='P-value',show_label=FALSE,label_cex=0.5,legend_cex=0.7,
                             label_type='distribute',main='',pdf_file=NULL){
  dat <- unique(dat[,c(label_col,logFC_col,Pv_col)])
  dat <- dat[order(dat[,3],decreasing=TRUE),]
  dat <- dat[which(is.na(dat[,2])==FALSE),]
  x <- as.numeric(dat[,logFC_col])
  y <- as.numeric(dat[,Pv_col]);
  y <- -log10(y)
  s1 <- which(abs(x)>=logFC_thre & y>= -log10(Pv_thre))
  if(is.null(pdf_file)==FALSE){
    geneWidth <- 0
    geneHeight <- 0
    if(length(s1)>0){
      s11 <- s1[which(x[s1]>=0)]
      s12 <- s1[which(x[s1]<0)]
      geneWidth  <- max(strwidth(dat[s1,label_col],'inches',cex=label_cex))
      geneHeight <- max(strwidth(toupper(letters),'inches',cex=label_cex))*max(length(s11),length(s12))*1.2
      if(show_label==TRUE & label_type=='distribute'){
        pdf(pdf_file,width=10+geneWidth*2,height=max(10,geneHeight))
      }else{
        pdf(pdf_file,width=10,height=10)
      }
    }else{
      pdf(pdf_file,width=10,height=10)
    }
  }
  par(mai=c(1.5,2,1.5,1))
  mm <- max(abs(x))
  if(show_label==TRUE & label_type=='distribute'){
    plot(y~x,pch=16,col=get_transparent('grey',0.7),xlab=xlab,ylab='',
         xlim=c(-3*mm/7*geneWidth-1.5*mm,1.5*mm+3*mm/7*geneWidth),ylim=c(0,max(y)*1.5),yaxt='n',main=main,cex.lab=1.2,cex.main=1.6)
  }else{
    plot(y~x,pch=16,col=get_transparent('grey',0.7),xlab=xlab,ylab='',
         xlim=c(-mm*1.5,mm*1.5),ylim=c(0,max(y)*1.5),yaxt='n',main=main,cex.lab=1.2,cex.main=1.6)
  }
  axis(side=2,at=seq(0,round(max(y)*1.5)),labels=c(1,format(10^-seq(1,round(max(y)*1.5)),scientific = TRUE)),las=2)
  mtext(side=2,line = 4,ylab,cex=1.2)
  z_val <- sapply(dat[,Pv_col]*sign(x),combinePvalVector)[1,]
  if(logFC_thre>0){abline(v=logFC_thre,lty=2,lwd=0.5);abline(v=-logFC_thre,lty=2,lwd=0.5)}
  if(Pv_thre<1) abline(h=-log10(Pv_thre),lty=2,lwd=0.5);
  points(y~x,pch=16,col=get_transparent('grey',0.7))
  #s1 <- which(abs(x)>=logFC_thre & y>= -log10(Pv_thre))
  s_col <- z2col(c(10,-10),sig_thre=0); names(s_col) <- as.character(c(1,-1))
  points(y[s1]~x[s1],pch=16,col=s_col[as.character(sign(x[s1]))])
  legend(0,par()$usr[4],c('Down-regulated','Not-Significant','Up-regulated'),
         fill=c(s_col[2],get_transparent('grey',0.7),s_col[1]),border=NA,bty='o',
         bg='white',box.col='black',horiz=TRUE,xjust=0.5,cex=legend_cex)
  if(show_label==TRUE){
    s11 <- s1[which(x[s1]>=0)]
    s12 <- s1[which(x[s1]<0)]
    dd <- (par()$usr[2]-par()$usr[1])/100
    if(label_type == 'origin'){
      if(length(s11)>0) text(x[s11]+dd,y[s11],dat[s11,label_col],cex=label_cex,adj=0)
      if(length(s12)>0) text(x[s12]-dd,y[s12],dat[s12,label_col],cex=label_cex,adj=1)
    }else{
      if(length(s11)>0){
        dd <- (par()$usr[4]-par()$usr[3])/(length(s11)+1)
        rect(xright=par()$usr[2],ybottom=par()$usr[3],xleft=max(abs(x))*1.1,ytop=par()$usr[4],col='white')
        ypos <- seq(from=par()$usr[3],by=dd,length.out=length(s11))+dd
        text(max(abs(x))*1.125,ypos,dat[s11,label_col],cex=label_cex,adj=0)
        segments(x0=max(abs(x))*1.1,x1=x[s11],y0=ypos,y1=y[s11],lwd=0.4,col='orange')
      }
      if(length(s12)>0){
        dd <- (par()$usr[4]-par()$usr[3])/(length(s12)+1)
        rect(xleft=par()$usr[1],ybottom=par()$usr[3],xright=-max(abs(x))*1.1,ytop=par()$usr[4],col='white')
        ypos <- seq(from=par()$usr[3],by=dd,length.out=length(s12))+dd
        text(-max(abs(x))*1.125,ypos,dat[s12,label_col],cex=label_cex,adj=1)
        segments(x0= -max(abs(x))*1.1,x1=x[s12],y0=ypos,y1=y[s12],lwd=0.4,col='green')
      }
    }
  }
  sig_info <- dat[s1,]
  if(is.null(pdf_file)==FALSE) dev.off()
  return(sig_info)
  #return(TRUE)
}
## draw heatmap, use ComplexHeatmap
# check Heatmap()
draw.heatmap <- function(mat=NULL,use_genes=rownames(mat),use_gene_label=rownames(mat),use_samples=colnames(mat),use_sample_label=colnames(mat),
                         phenotype_info=NULL,use_phe=NULL,main='',scale='none',pdf_file=NULL,
                         cluster_rows=TRUE,cluster_columns=TRUE,
                         clustering_distance_rows='pearson',clustering_distance_columns='pearson',
                         row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12),
                         show_row_names=TRUE,show_column_names=TRUE,...){
  names(use_gene_label) <- use_genes
  names(use_sample_label) <- use_samples
  use_genes <- intersect(use_genes,rownames(mat))
  use_samples <- intersect(use_samples,colnames(mat))
  use_mat <- mat[use_genes,use_samples]
  rownames(use_mat) <- use_gene_label[rownames(use_mat)]
  colnames(use_mat) <- use_sample_label[colnames(use_mat)]
  row_names_max_width <- max(strwidth(rownames(use_mat),'inches',cex=row_names_gp[[1]]/7))
  row_names_max_width <- unit(row_names_max_width,'inches')
  column_names_max_height <- max(strwidth(colnames(use_mat),'inches',cex=column_names_gp[[1]]/7))
  column_names_max_height <- unit(column_names_max_height,'inches')
  if(scale=='row'){use_mat <- t(apply(use_mat,1,std))}
  if(scale=='column'){use_mat <- apply(use_mat,2,std)}
  if(length(use_phe)==0){
    if(scale=='none'){
      ht1 <- Heatmap(use_mat, column_title = main,name='Raw value',
                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                     row_names_gp=row_names_gp,column_names_gp=column_names_gp,
                     show_row_names=show_row_names,show_column_names=show_column_names,
                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
    if(scale!='none'){
      ht1 <- Heatmap(use_mat, column_title = main, name='Z value',
                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                     row_names_gp=row_names_gp,column_names_gp=column_names_gp,
                     show_row_names=show_row_names,show_column_names=show_column_names,
                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
  }else{
    if(length(use_phe)==1){
      use_phe_info <- as.data.frame(phenotype_info[,use_phe],stringsAsFactors=FALSE)
      rownames(use_phe_info) <- rownames(phenotype_info)
      colnames(use_phe_info) <- gsub(' ','.',use_phe)
    }else{
      use_phe_info <- phenotype_info[,use_phe]
      colnames(use_phe_info) <- gsub(' ','.',use_phe)
    }
    use_phe <- colnames(use_phe_info)
    l2c <- get.class.color(unique(as.character(as.matrix(use_phe_info))),use_color=brewer.pal(12, 'Paired'))
    use_col <- lapply(use_phe,function(x)l2c[unique(use_phe_info[,x])])
    names(use_col) <- use_phe
    ha_column <- HeatmapAnnotation(df = data.frame(use_phe_info),col = use_col)
    if(scale=='none'){
      ht1 <- Heatmap(use_mat, column_title = main, top_annotation = ha_column,name='Raw value',
                   cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                   show_row_names=show_row_names,show_column_names=show_column_names,
                   clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                   row_names_gp=row_names_gp,column_names_gp=column_names_gp,
                   row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
    if(scale!='none'){
      ht1 <- Heatmap(use_mat, column_title = main, top_annotation = ha_column,name='Z value',
                   cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                   clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                   row_names_gp=row_names_gp,column_names_gp=column_names_gp,
                   show_row_names=show_row_names,show_column_names=show_column_names,
                   row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
  }
  ht_list <- ht1
  if(is.null(pdf_file)==FALSE){
    ww <- 1.5*column_names_gp[[1]]/72*ncol(use_mat)+max(strwidth(rownames(use_mat),'inches'))+5
    hh <- 1.5*row_names_gp[[1]]/72*nrow(use_mat)+max(strwidth(colnames(use_mat),'inches'))+3
    pdf(pdf_file,width=ww,height=hh)
  }
  draw(ht_list,heatmap_legend_side='left',annotation_legend_side='right')
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}
################################ Function enrichment related functions
find.gsByGene <- function(gene=NULL,use_gs=NULL){
  if(is.null(use_gs)==TRUE){
    all <- unlist(all_gs2gene,recursive = FALSE)
  } else {
    if('all' %in% use_gs){
      all <- unlist(all_gs2gene,recursive = FALSE)
    }else{
      all <-  unlist(all_gs2gene[use_gs],recursive = FALSE)
    }
  }
  x1 <- unlist(lapply(all,function(x){
    length(intersect(x,gene))
  }))
  x2 <- names(x1[which(x1==length(gene))])
  return(x2)
}
##
# functions to merge gs
merge_gs <- function(all_gs2gene=NULL,use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5')){
  nn <- unlist(lapply(all_gs2gene[use_gs],names))
  use_gs2gene <- unlist(all_gs2gene[use_gs],recursive = FALSE)
  names(use_gs2gene)<-nn
  use_gs2gene
}
# simple functions
list2mat <- function(input_list){
  all_x <- unique(unlist(input_list))
  all_y <- unique(names(input_list))
  mat1 <- matrix(0,nrow=length(all_x),ncol=length(all_y))
  rownames(mat1) <- all_x; colnames(mat1) <- all_y;
  for(i in names(input_list)){
    mat1[input_list[[i]],i] <- 1
  }
  return(mat1)
}
vec2list <- function(input_v,sep=NULL){
  if(is.null(sep)==TRUE){
    tmp2 <- list()
    input_vn <- names(input_v)
    input_v <- as.character(input_v); names(input_v) <- input_vn
    for(i in 1:length(input_v)){
      if(input_v[i] %in% names(tmp2)){
        tmp2[[input_v[i]]] <- c(tmp2[[input_v[i]]],names(input_v)[i])
      }else{
        tmp2[[input_v[i]]] <- names(input_v)[i]
      }
    }
  }else{
    tmp1 <- aggregate(names(input_v),list(input_v),function(x)paste(x,collapse=sep))
    tmp2 <- tmp1$x; names(tmp2) <- tmp1$Group.1
  }
  tmp2
}
## Function enrichment by fisher
funcEnrich.Fisher <- function(input_list=NULL,bg_list=NULL,gs2gene=NULL,use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG'),
                              min_gs_size=5,max_gs_size=500,Pv_adj='fdr',Pv_thre=0.1){
  if(is.null(gs2gene)==TRUE){
    if(length(setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)))>0){
      message(sprintf('Input %s not in all_gs2gene, please check all_gs2gene_info (items in Category or Sub-Category) and re-try!',
                      paste(setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)),collapse=';')));
      return(FALSE)
    }
    if(length(use_gs)>1){
      gs2gene <- unlist(all_gs2gene[use_gs],recursive = FALSE)
      names(gs2gene)<-unlist(lapply(all_gs2gene[use_gs],names))
    }else{
      gs2gene <- all_gs2gene[[use_gs]]
    }
  }
  all_gs <- names(gs2gene)
  input_list <- unique(input_list)
  bg_list <- unique(bg_list)
  if(!is.null(bg_list)){
    use_gs2gene <- lapply(gs2gene,function(x){intersect(x,bg_list)})
    names(use_gs2gene) <- names(gs2gene)
  }else{
    use_gs2gene <- gs2gene
  }
  bg_list <- unique(unlist(use_gs2gene))
  ## size selection
  s1 <- unlist(lapply(use_gs2gene,length))
  w1 <- which(s1>=min_gs_size & s1<=max_gs_size)
  use_gs2gene <- use_gs2gene[w1]
  all_gs <- names(use_gs2gene) ## all tested gene set number
  ## input filter
  input_list <- intersect(input_list,bg_list)
  bg_or <- length(input_list)/length(bg_list)
  s1 <- unlist(lapply(use_gs2gene,function(x){
    length(intersect(input_list,x))/length(x)
  }))
  w1 <- which(s1>bg_or)
  use_gs2gene <- use_gs2gene[w1]
  empty_vec <- as.data.frame(matrix(NA,ncol=9));colnames(empty_vec) <- c('#Name','Total_item','Num_item','Num_list','Num_list_item','Ori_P','Adj_P','Odds_Ratio','Intersected_items')
  if(length(w1)==0) return(empty_vec)
  ## fisher~
  pv <- lapply(use_gs2gene,function(x){
    n11 <- length(intersect(input_list,x))
    n12 <- length(intersect(input_list,setdiff(bg_list,x)))
    n21 <- length(setdiff(x,input_list))
    n22 <- length(setdiff(bg_list,unique(c(input_list,x))))
    ft <- fisher.test(cbind(c(n11,n12),c(n21,n22)))$p.value
    or <- n11/n12/(n21/n22)
    c(length(bg_list),length(x),length(input_list),n11,ft,or,paste(intersect(input_list,x),collapse=';'))
  })
  pv <- do.call(rbind,pv)
  pv <- as.data.frame(pv,stringsAsFactors=FALSE)
  colnames(pv) <- c('Total_item','Num_item','Num_list','Num_list_item','Ori_P','Odds_Ratio','Intersected_items')
  pv[1:6] <- lapply(pv[1:6],as.numeric)
  pv$Adj_p <- p.adjust(pv$Ori_P,method=Pv_adj,n=length(all_gs))
  pv$`#Name` <- rownames(pv)
  pv <- pv[,c(9,1:5,8,6:7)]
  pv <- pv[order(pv$Ori_P),]
  use_pv <- pv[which(pv$Adj_p<=Pv_thre),]
  return(use_pv)
}
## Barplot for enrichment
draw.funcEnrich.bar <- function(funcEnrich_res=NULL,top_number=30,
                                Pv_col='Ori_P',item_col='Intersected_items',
                                Pv_thre=0.1,display_genes=FALSE,name_col='#Name',
                                gs_cex=0.5,gene_cex=0.5,main='',bar_col=brewer.pal(8,'RdBu')[7],eg_num=5,
                                pdf_file=NULL){
  if(is.null(top_number)==TRUE) top_number <- nrow(funcEnrich_res)
  funcEnrich_res <- funcEnrich_res[which(funcEnrich_res[,Pv_col]<=Pv_thre),]
  if(nrow(funcEnrich_res)>top_number) funcEnrich_res <- funcEnrich_res[1:top_number,]
  pv_val <- funcEnrich_res[,Pv_col]
  s1 <- funcEnrich_res[[item_col]]
  s1 <- sapply(s1,function(x){
    x1 <- unlist(strsplit(x,';'))
    if(length(x1)>eg_num){x2 <- sprintf('Total %d items, e.g: %s',length(x1),paste(x1[1:5],collapse=';'));x<-x2;}
    return(x)
  })
  textWidth  <- max(strwidth(rownames(funcEnrich_res),units='inch',cex=gs_cex))
  textWidth1  <- max(strwidth(s1,units='inch',cex=gene_cex))
  if(is.null(pdf_file)==FALSE){
    if(display_genes==TRUE)
      pdf(pdf_file,width=textWidth+1+textWidth1+4,height=3+nrow(funcEnrich_res)*1.5*max(strheight(funcEnrich_res[,name_col],units='inch',cex=gs_cex)))
    else
      pdf(pdf_file,width=textWidth+1+4,height=3+nrow(funcEnrich_res)*1.5*max(strheight(funcEnrich_res[,name_col],units='inch',cex=gs_cex)))
  }
  if(display_genes==TRUE) par(mai=c(1.5,textWidth+1,1,textWidth1)) else par(mai=c(1.5,textWidth+1,1,1))
  print(par()$mar)
  a<-barplot(rev(-log10(pv_val)),horiz=TRUE,border = NA,col=bar_col,main=main,xaxt='n',
             xlim=c(0,round(max(-log10(pv_val)))));
  mtext(side=1,'P-value',line=3)
  axis(side=1,at=0:round(par()$usr[2]),labels=10^-(0:round(par()$usr[2])),xpd=TRUE)
  text(-0.01,a,rev(funcEnrich_res[,name_col]),adj=1,xpd=TRUE,cex=gs_cex);
  if(display_genes==TRUE) text(rev(-log10(pv_val))+0.1,a,rev(s1),adj=0,xpd=TRUE,cex=gene_cex)
  if(is.null(pdf_file)==FALSE) dev.off()
  return(funcEnrich_res)
}
## cluster plot for enrichment results
draw.funcEnrich.cluster <- function(funcEnrich_res=NULL,top_number=30,Pv_col='Ori_P',name_col='#Name',item_col='Intersected_items',Pv_thre=0.1,
                                    gs_cex=0.7,gene_cex=0.8,pv_cex=0.7,main='',h=0.95,cluster_gs=TRUE,cluster_gene=TRUE,
                                    pdf_file=NULL,use_genes=NULL,return_mat=FALSE){
  if(is.null(top_number)==TRUE) top_number <- nrow(funcEnrich_res)
  funcEnrich_res <- funcEnrich_res[which(funcEnrich_res[,Pv_col]<=Pv_thre),]
  if(nrow(funcEnrich_res)>top_number) funcEnrich_res <- funcEnrich_res[1:top_number,]
  pv_val <- funcEnrich_res[,Pv_col]; names(pv_val) <- rownames(funcEnrich_res)
  all_g2s <- sapply(funcEnrich_res[,item_col],function(x1)unlist(strsplit(x1,';')))
  names(all_g2s) <- funcEnrich_res[,name_col]
  mat1 <- t(list2mat(all_g2s))
  mat1 <- mat1[rev(funcEnrich_res[,name_col]),]
  if(is.null(use_genes)==FALSE) mat1 <- mat1[,intersect(colnames(mat1),use_genes)]
  if(ncol(mat1)==0){message('No genes left, please check and re-try!');return(FALSE)}
  #mat1 <- mat1[,order(colnames(mat1))]
  h_gs <- hclust(dist(mat1,method='binary'))
  h_gene <- hclust(dist(t(mat1),method='binary'))
  gs_cluster <- cutree(h_gs,h=h)
  gene_cluster <- cutree(h_gene,h=h)
  if(cluster_gs==FALSE){gs_cluster <- rep(1,length.out=nrow(mat1));names(gs_cluster)<-rownames(mat1);}
  if(cluster_gene==FALSE){gene_cluster <- rep(1,length.out=ncol(mat1));names(gene_cluster)<-colnames(mat1);}
  cc1 <- colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(gs_cluster)))
  cc2 <- colorRampPalette(brewer.pal(9,'Pastel1'))(length(unique(gene_cluster)))
  cc3 <- colorRampPalette(brewer.pal(9,'Reds')[3:9])(100)
  # get gs order
  if(cluster_gs==TRUE) gs_cluster <- gs_cluster[h_gs$order]
  tmp2 <- vec2list(gs_cluster,sep=NULL)
  tmp2 <- tmp2[rev(order(unlist(lapply(tmp2,function(x)min(funcEnrich_res[x,Pv_col])))))]
  mat1 <- mat1[unlist(tmp2),]
  if(cluster_gene==TRUE) mat1 <- mat1[,h_gene$order]
  gs_cluster <- gs_cluster[rownames(mat1)]
  gene_cluster <- gene_cluster[colnames(mat1)]
  #a <- heatmap(mat1,scale='none',col=c('white','red'),ColSideColors = cc2[gene_cluster[colnames(mat1)]],
  # RowSideColors = cc1[gs_cluster[rownames(mat1)]],distfun = function(x){dist(x,method='binary')},margins=c(5,20),Rowv=NA,Colv=NA)
  pv <- pv_val[rownames(mat1)]
  #####
  gsWidth <- max(strwidth(rownames(mat1),units='inch',cex=gs_cex))
  gsWidth2 <- max(strheight(rownames(mat1),units='inch',cex=gs_cex))*nrow(mat1)*1.5
  geneWidth <- max(strheight(colnames(mat1),units='inch',cex=gene_cex))*ncol(mat1)*1.25
  geneWidth2 <- max(strwidth(colnames(mat1),units='inch',cex=gene_cex))
  pv1 <- format(pv,scientific = TRUE,digits = 3)
  pvWidth   <- max(strwidth(pv1,units='inch',cex=pv_cex))
  ww <- gsWidth+pvWidth+geneWidth+0.5
  hh <- geneWidth2+gsWidth2+0.5
  if(is.null(pdf_file)==FALSE) pdf(file=pdf_file,width=ww,height=hh)
  #####
  par(mai=c(0.5,0.5,geneWidth2,0));
  mr <- 1/pvWidth
  geneWidth <- round(geneWidth*mr)
  pvWidth <- round(pvWidth*mr)
  gsWidth <- round(gsWidth*mr)
  layout(t(matrix(c(rep(1,geneWidth),rep(2,pvWidth),rep(3,gsWidth)),byrow=TRUE)))
  #print(t(matrix(c(rep(1,geneWidth),rep(2,pvWidth),rep(3,gsWidth)))))
  image(t(mat1),col=c('white',cc3[1]),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr;
  gs_cs <- cumsum(table(gs_cluster)[unique(gs_cluster)])
  gene_cs <- cumsum(table(gene_cluster)[unique(gene_cluster)])
  xx <- (pp[2]-pp[1])/length(gene_cluster); yy <- (pp[4]-pp[3])/length(gs_cluster)
  abline(h=gs_cs*yy+pp[3],col='black',lwd=0.25)
  abline(v=gene_cs*xx+pp[1],col='black',lwd=0.25)
  abline(v=pp[1:2],col='black',lwd=0.5);
  abline(h=pp[3:4],col='black',lwd=0.5)
  ## draw gene name
  if(cluster_gene==TRUE){
    text(c(1:length(gene_cluster))*xx+pp[1]-xx/2,pp[4]+0.3*yy,colnames(mat1),xpd=TRUE,adj=0,cex=gene_cex,srt=90)
    rect(xleft=c(1:length(gene_cluster))*xx+pp[1]-xx,xright=c(1:length(gene_cluster))*xx+pp[1],ybottom=pp[4],ytop=pp[4]+yy*0.25,
         col=cc2[gene_cluster[colnames(mat1)]],xpd=TRUE,border=NA)
  }else{
    text(c(1:length(gene_cluster))*xx+pp[1]-xx/2,pp[4]+yy/10,colnames(mat1),xpd=TRUE,adj=0,cex=gene_cex,srt=90)
  }
  # draw p-value
  #pp <- par()$usr;
  par(mai=c(0.5,0,geneWidth2,0));
  plot(1,xaxt='n',yaxt='n',bty='n',xlim=c(pp[1],pp[2]),ylim=c(pp[3],pp[4]),col='white',xlab='',ylab='')
  pp <- par()$usr;
  yy <- (pp[4]-pp[3])/length(gs_cluster)
  pv_c <- z2col(qnorm(1-pv))
  rect(xleft=pp[1],xright=pp[2],ybottom=c(1:length(gs_cluster))*yy+pp[3]-yy,
       ytop=c(1:length(gs_cluster))*yy+pp[3],
       col=pv_c,border = NA)
  text(0.5,c(1:length(gs_cluster))*yy+pp[3]-yy/2,pv1,xpd=TRUE,adj=0.5,cex=pv_cex)
  abline(v=pp[1],col='black',lwd=2);
  # draw gs name
  par(mai=c(0.5,0,geneWidth2,0.5));
  zz <- min(c(xx,yy))
  plot(1,xaxt='n',yaxt='n',bty='n',xlim=c(pp[1],pp[2]),ylim=c(pp[3],pp[4]),col='white',xlab='',ylab='')
  pp <- par()$usr;
  yy <- (pp[4]-pp[3])/length(gs_cluster)
  text(pp[1]+zz*0.2,c(1:length(gs_cluster))*yy+pp[3]-yy/2,rownames(mat1),xpd=TRUE,adj=0,cex=gs_cex)
  # get region for p-value
  abline(v=pp[1:2],col='black',lwd=0.5);
  abline(h=pp[3:4],col='black',lwd=0.5);
  abline(v=pp[1],col='black',lwd=0.5);
  abline(h=gs_cs*yy+pp[3],col='black',lwd=0.25)
  ##
  if(is.null(pdf_file)==FALSE) dev.off()
  layout(1);
  if(return_mat==TRUE){
    return(mat1)
  }else{
    return(TRUE)
  }
}
#### bubble plot
## Pv_thre, threshold for fisher funcEnrich
## driver_type, can add additional driver type info
draw.bubblePlot <- function(driver_list=NULL,show_label=NULL,Z_val=NULL,driver_type=NULL,
                            target_list=NULL,transfer2symbol2type=NULL,
                            bg_list=NULL,min_gs_size=5,max_gs_size=500,gs2gene=NULL,use_gs=NULL,Pv_adj='none',Pv_thre=0.1,
                            top_geneset_number=30,top_driver_number=30,
                            pdf_file=NULL,main='',mark_gene=NULL,driver_cex=1,gs_cex=1){
  ## check NULL

  ## get top driver list
  names(show_label) <- driver_list
  names(Z_val) <- driver_list
  names(driver_type) <- driver_list
  driver_list <- driver_list[order(abs(Z_val),decreasing = TRUE)]
  if(length(driver_list)>top_driver_number){
    driver_list <- driver_list[1:top_driver_number]
  }
  driver_list <- driver_list[order(Z_val[driver_list],decreasing=TRUE)]
  ## get target gene for driver_list
  transfer_tab <- transfer2symbol2type
  rownames(transfer_tab) <- transfer_tab[,1]
  target_gene <- lapply(driver_list,function(x){
    x1 <- target_list[[x]]$target
    x1 <- x1[which(x1 %in% transfer_tab[,1])]
    x2 <- transfer_tab[x1,]
    x3 <- x2[which(x2[,3]=='protein_coding'),]
    target <- unique(x2[,2])
    target_pc <- unique(x3[,2])
    return(list(target,target_pc))
  })
  names(target_gene) <- driver_list
  ##
  f_res <- lapply(target_gene,function(x){
    funcEnrich.Fisher(input_list=x[[1]],bg_list=bg_list,gs2gene=gs2gene,use_gs=use_gs,min_gs_size=min_gs_size,max_gs_size=max_gs_size,
                      Pv_adj='none',Pv_thre=Pv_thre)
  })
  names(f_res) <- names(target_gene)
  ## get display matrix
  all_path <- unique(unlist(lapply(f_res,function(x){x[[1]]})))
  all_path <- all_path[which(is.na(all_path)==FALSE)] ## get all sig path
  f_mat <- lapply(f_res,function(x){
      as.data.frame(x)[all_path,5:6]
  })
  f_mat2 <- do.call(rbind,lapply(f_mat,function(x)qnorm(1-x[[2]])))
  f_mat2[which(is.na(f_mat2)==TRUE | f_mat2==-Inf)] <- 0
  f_mat1 <- do.call(rbind,lapply(f_mat,function(x)x[,1]))
  colnames(f_mat1) <- all_path
  colnames(f_mat2) <- all_path
  f_mat3 <- t(apply(f_mat2,1,function(x){
    o1 <- order(abs(x),decreasing = TRUE)[1:3];
    x1<-x;x1[setdiff(1:length(x1),o1)] <- 0;x1
  }))
  ## use top
  min_path_num <- 5
  max_path_num <- top_geneset_number
  all_path_order <- apply(f_mat3,2,max)
  all_path_order <- sort(all_path_order,decreasing = TRUE)
  if(length(all_path_order) > max_path_num){
    w1 <- 1:length(all_path_order)
    if(length(w1)>=min_path_num & length(w1)<=max_path_num){
      all_path <- names(sort(all_path_order[w1],decreasing=TRUE))
    }else{
      if(length(w1)<min_path_num){
        all_path <- names(sort(all_path_order,decreasing=TRUE)[1:min_path_num])
      }else{
        all_path <- names(sort(all_path_order,decreasing=TRUE)[1:max_path_num])
      }
    }
    f_mat1 <- f_mat1[,all_path]
    f_mat2 <- f_mat2[,all_path]
  }else{
    f_mat1 <- f_mat1[,names(all_path_order)]
    f_mat2 <- f_mat2[,names(all_path_order)]
  }
  nr <- ncol(f_mat1)
  nc <- nrow(f_mat1)
  gsWidth  <- max(strwidth(colnames(f_mat1),'inches',cex=gs_cex))
  gsHeight <- max(strheight(colnames(f_mat1),'inches',cex=gs_cex)*nrow(f_mat1))
  driverWidth  <- max(strwidth(show_label[colnames(f_mat1)],'inches',cex=gs_cex))
  driverHeight <- max(strheight(show_label[colnames(f_mat1)],'inches',cex=gs_cex)*ncol(f_mat1))
  ww <- (1+nc)*0.5+ gsWidth + 1.5 + 2
  hh <- (4+nr)*0.5+ driverWidth + 2+1
  ## output to pdf
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=ww,height=hh)
  layout(1);par(mai=c(driverWidth+2,gsWidth+1.5,1,2))
  plot(1,bty='n',col='white',xlim=c(0,nc+1),ylim=c(-2,nr+1),xaxt='n',yaxt='n',xlab='',ylab='',main=main)
  segments(x0=0,x1=nc,y0=0:nr,y1=0:nr,col='dark grey',xpd=TRUE)
  segments(x0=0:nc,x1=0:nc,y0=0,y1=nr,col='dark grey',xpd=TRUE)
  segments(x0=0:nc,x1=0:nc,y0=0,y1=-3,col='grey',xpd=TRUE)
  #
  text(-0.5,1:nr-0.5,colnames(f_mat1),xpd=TRUE,srt=0,adj=1,cex=gs_cex) ## sig pathways
  if(is.null(mark_gene)==TRUE){
    text(1:nc-0.5,-2.5,show_label[rownames(f_mat1)],xpd=TRUE,srt=90,adj=1,cex=driver_cex) ## sig regulators
  }else{
    bc <- rep('black',length.out=length(rownames(f_mat1)))
    bc[which(rownames(f_mat1) %in% mark_gene)] <- 'red'
    print(table(bc))
    text(1:nc-0.5,-2.5,show_label[rownames(f_mat1)],xpd=TRUE,srt=90,adj=1,col=bc,cex=driver_cex) ## sig regulators
  }
  ## draw circle
  max_size <- max(f_mat1,na.rm=TRUE)
  f_mat1 <- f_mat1/max_size
  cc_r <- matrix(z2col(f_mat2,n_len=30,sig_thre=qnorm(1-0.1)),ncol=ncol(f_mat2),byrow = FALSE)
  for(i in 1:nrow(f_mat1)){
    for(j in 1:ncol(f_mat1)){
      draw.circle(i-0.5,j-0.5,radius=f_mat1[i,j]/2,col=cc_r[i,j])
    }
  }
  ## draw circle legend
  legend_size <- unique(round(seq(1,max_size,length.out=5)))
  for(i in 1:length(legend_size)){
    draw.circle(length(legend_size)-i+1.5,nr+0.5,radius=0.5*legend_size[i]/max_size)
    text(length(legend_size)-i+1.5,nr+1,legend_size[i])
  }
  text(0.5,nr+0.5,'Size ')
  ## draw p-value legend
  p_label <- c(1,0.1,0.05,0.01,0.001,1e-4,1e-10)
  p_thre <- qnorm(1-c(1,0.1,0.05,0.01,0.001,0.0001,1e-10))
  p_col  <- z2col(p_thre,n_len=30,sig_thre=qnorm(1-0.1))
  p_col_m  <- z2col(-p_thre,n_len=30,sig_thre=qnorm(1-0.1))
  ybottom <- seq(nr-6,nr-1,length.out=1+length(p_thre))[1:length(p_thre)]
  ytop <- seq(nr-6,nr-1,length.out=1+length(p_thre))[2:(length(p_thre)+1)]
  for(i in 1:length(p_thre)){
    rect(nc+0.5,ybottom[i],nc+1.5,ytop[i],col=p_col[i],xpd=TRUE)
    text(nc+1.7,(ybottom[i]+ytop[i])/2,adj=0,p_label[i],xpd=TRUE)
  }
  ybottom <- seq(nr-11,nr-6,length.out=1+length(p_thre))[1:length(p_thre)]
  ytop <- seq(nr-11,nr-6,length.out=1+length(p_thre))[2:(length(p_thre)+1)]
  for(i in 2:length(p_thre)){
    rect(nc+0.5,ybottom[i],nc+1.5,ytop[i],col=rev(p_col_m)[i-1],xpd=TRUE)
    text(nc+1.7,(ybottom[i]+ytop[i])/2,adj=0,rev(p_label)[i-1],xpd=TRUE)
  }
  text(nc+1.5,nr-0.5,'P-Value',xpd=TRUE)
  # target size
  max_target_size <- max(unlist(lapply(target_gene,function(x)max(unlist(lapply(x,length))))))
  ori_size <- unlist(lapply(target_gene,function(x)length(x[[1]])))
  pro_size <- unlist(lapply(target_gene,function(x)length(x[[2]])))
  rect(0:(nc-1)+0.15,-2,0:(nc-1)+0.45,1.5*ori_size/max_target_size-2,col='blue')
  rect(0:(nc-1)+0.55,-2,0:(nc-1)+0.85,1.5*pro_size/max_target_size-2,col='green')
  segments(y0=-2,y1=-0.5,x0=0,x1=0)
  segments(y0=seq(-2,-0.5,length.out=3),y1=seq(-2,-0.5,length.out=3),x0=-0.25,x1=0)
  text(-0.3,seq(-2,-0.5,length.out=3),round(seq(0,max_target_size,length.out=3)),adj=1,xpd=TRUE)
  legend(-5,-0.5,c('target_size','target_size\n(protein_coding)'),fill=c('blue','green'),border=NA,bty='n',cex=0.8,xpd=TRUE)
  # add sig color
  sig_col <- z2col(Z_val[driver_list],n_len=30)
  rect(0:(nc-1)+0.35,-0.4,0:(nc-1)+0.65,-0.1,col=sig_col)
  # add driver type !!!
  if(is.null(driver_type)==FALSE){
    cc_tmp <- get.class.color(unique(driver_type[driver_list]))
    points(x=0:(nc-1)+0.5,y=rep(-2.2,length.out=nc),col=cc_tmp[driver_type[driver_list]],pch=16,cex=2.5)
    legend(nc+0.5,0.5,names(cc_tmp),fill=cc_tmp,border=NA,bty='n',cex=0.8,xpd=TRUE)
  }
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}
###
draw.GSEA <- function(rank_profile=NULL,use_genes=NULL,use_direction=NULL,main='',pdf_file=NULL,
                      annotation=NULL,annotation_cex=1.2,left_annotation=NULL,right_annotation=NULL){
  #### start plot
  if(is.null(pdf_file)==FALSE){
    pdf(pdf_file,width=10,height=10)
  }
  if(is.null(use_direction)==FALSE){
    new_rank_profile <- c(rank_profile,-rank_profile)
    names(new_rank_profile) <- c(paste0('POS_',names(rank_profile)),paste0('NEG_',names(rank_profile)))
    rank_profile <- new_rank_profile
    use_genes[which(use_target_direction==1)] <- paste0('POS_',use_genes[which(use_target_direction==1)])
    use_genes[which(use_target_direction==-1)] <- paste0('NEG_',use_genes[which(use_target_direction==-1)])
  }
  rank_profile <- sort(rank_profile,decreasing = TRUE)
  r_len <- length(rank_profile)
  use_pos <- which(names(rank_profile) %in% use_genes)
  layout(matrix(c(rep(4,10),3,2,rep(1,10)),ncol=1))
  # plot
  ## rank for all
  par(mar=c(10,6,0,2))
  mm <- max(abs(rank_profile))
  y1 <- seq(-mm,mm,length.out=7); y1 <- round(y1,1)
  unit <- r_len/10; unit <- round(unit/100)*100
  x1 <- seq(0,r_len,by=unit);x1 <- unique(x1); x1 <- c(x1,max(x1)+unit)
  par(usr=c(0,max(x1),-mm,mm))
  plot(rank_profile,col='grey',pch=16,xaxt='n',yaxt='n',xlab='',ylab='',bty='n',type='n',ylim=c(-mm,mm))
  polygon(x=c(0,1:r_len,r_len),y=c(0,rank_profile,0),col='grey',border=NA)
  if(is.null(left_annotation)==FALSE) text(0+r_len/100,mm,adj=0,left_annotation,col='red',xpd=TRUE)
  if(is.null(right_annotation)==FALSE) text(r_len-r_len/100,-mm,adj=1,right_annotation,col='blue',xpd=TRUE)
  #axis(side=2,at=y1,labels=y1,las=2);
  pp <- par()$usr
  segments(0,0,r_len,0,lwd=0.5)
  segments(0,min(y1),0,max(y1),lwd=1.5)
  text(-(pp[2]-pp[1])/50,y1,y1,adj=1,xpd=TRUE)
  segments(-(pp[2]-pp[1])/100,y1,0,y1,lwd=1.5)
  mtext(side=2,line = 2,'Ranked list metric (PreRanked)',cex=1.2)
  mtext(side=1,line = 3,'Rank in Ordered Dataset',cex=1.2)
  axis(side=1,at=x1,labels=get_label_manual(x1))
  # get zero cross
  w1 <- which.min(abs(rank_profile))
  abline(v=w1,lty=2,col='grey')
  text(w1,-mm/4,sprintf('Zero cross at %d',w1),adj=0.5)
  legend(w1,pp[3]-(pp[4]-pp[3])/4,lty=1,lwd=2,c('Enrichment profile','Hits','Ranking metric scores'),
         col=c('green','black','grey'),xpd=TRUE,horiz=TRUE,xjust=0.5,cex=1.2)
  pm <- par()$usr
  ## get image bar
  par(mar=c(0,6,0,2))
  use_col <- z2col(rank_profile,sig_thre = 0,n_len = 30,blue_col='blue',red_col='red')
  image(x=as.matrix(1:r_len),col=use_col,bty='n',xaxt='n',yaxt='n',xlim=c(pm[1],pm[2])/r_len)
  abline(v=use_pos/r_len,col='grey')
  ## mark gene position
  par(mar=c(0,6,0,2))
  plot(1,col='white',xlab='',ylab='',bty='n',xlim=c(1,r_len),xaxt='n',yaxt='n')
  abline(v=use_pos)
  ## GSEA ES
  par(mar=c(0,6,5,2))
  # get ES score
  es_res <- get_ES(rank_profile,use_genes)
  y2 <- seq(min(es_res$RES),max(es_res$RES),length.out=7); y2 <- round(y2,1)
  plot(es_res$RES,col='green',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',
       xlim=c(1,r_len),type='l',lwd=3,ylim=c(min(es_res$RES),max(y2)),main=main,xpd=TRUE)
  pp <- par()$usr
  #abline(h=0)
  #axis(side=2,at=y2,label=y2,las=2)
  segments(0,0,r_len,0,lwd=1.5)
  segments(0,min(y2),0,max(y2),lwd=1.5)
  text(-(pp[2]-pp[1])/50,y2,y2,adj=1,xpd=TRUE)
  segments(-(pp[2]-pp[1])/100,y2,0,y2,lwd=1.5)
  mtext(side=2,line = 2,'Enrichment score (ES)',cex=1.2)
  # add annotation
  if(is.null(annotation)==TRUE){
    annotation <- sprintf("KS test p-value:%s",format(ks.test(rank_profile,rank_profile[use_genes])$p.value,digits = 3,scientific = TRUE))
  }
  if(es_res$RES[which.max(abs(es_res$RES))]>0)
    text(r_len-r_len/50,max(y2),annotation,adj=1,cex=annotation_cex,xpd=TRUE)
  else
    text(0+r_len/50,min(y2)+(max(y2)-min(y2))/10,annotation,adj=0,cex=annotation_cex,xpd=TRUE)
  if(is.null(pdf_file)==FALSE){
    dev.off()
  }
  layout(1);
  return(TRUE)
}
## get enrichment score
get_ES <- function(rank_profile=NULL,use_genes=NULL,weighted.score.type=1){
  gene.list <- names(rank_profile)
  correl.vector <- rank_profile
  tag.indicator <- sign(match(gene.list, use_genes, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(use_genes)
  Nm <-  N - Nh
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}
##
get_z2p <- function(x){
  if(is.na(x)==TRUE) return('NA')
  x <- abs(x)
  if(x<5) return(format(1-pnorm(x),digits=2,scientific = TRUE))
  low_p <- .Machine$double.xmin
  low_z <- sapply(10^(-(1:(1+-log10(low_p)))),combinePvalVector)
  use_p <- low_z[2,which(low_z[1,]>=x)[1]]
  use_p <- format(use_p, digits=3,scientific = TRUE)
  if(use_p=='NA') use_p <- '<1e-308'
  use_p <- as.character(use_p)
  return(use_p)
}
# profile_trend: pos2neg, neg2pos
# target_nrow:1,2
# target_col:black,or RdBu
# target_col_type: DE,PN
draw.GSEA.NetBID <- function(DE=NULL,profile_col=NULL,profile_trend='pos2neg',
                             driver_list=NULL,show_label=NULL,driver_DA_Z=NULL,driver_DE_Z=NULL,target_list=NULL,
                             top_driver_number=30,target_nrow=2,target_col='RdBu',target_col_type='PN',
                             left_annotation='',right_annotation='',main='',
                             profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  if(!profile_col %in% colnames(DE)){
    message(sprintf('%s not in colnames of DE, please check and re-try!',profile_col))
    return(FALSE)
  }
  if(is.null(names(driver_DA_Z))) names(driver_DA_Z) <- driver_list
  if(is.null(names(driver_DE_Z))) names(driver_DE_Z) <- driver_list
  if(is.null(names(show_label))) names(show_label) <- driver_list
  if(length(driver_list)>top_driver_number){
    driver_DA_Z <- driver_DA_Z[driver_list]
    driver_DE_Z <- driver_DE_Z[driver_list]
    driver_list <- driver_list[order(abs(driver_DA_Z[driver_list]),decreasing = TRUE)][1:top_driver_number]
  }
  if(profile_trend=='pos2neg')
    driver_list <- driver_list[order(driver_DA_Z[driver_list],decreasing = FALSE)]
  else
    driver_list <- driver_list[order(driver_DA_Z[driver_list],decreasing = TRUE)]
  show_label <- show_label[driver_list]
  driver_DA_Z <- driver_DA_Z[driver_list]
  driver_DE_Z <- driver_DE_Z[driver_list]
  ###################
  ## calculate layout

  DE_profile <- DE[,profile_col]
  names(DE_profile) <- rownames(DE)
  if(profile_trend=='pos2neg') DE_profile <- sort(DE_profile,decreasing = TRUE) else DE_profile <- sort(DE_profile)
  n_gene <- length(DE_profile)
  if(target_nrow==2){
    n_driver <- length(driver_list)*2
    ratio1 <- ceiling(n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  } else {
    n_driver <- length(driver_list)
    ratio1 <- ceiling(1.2*n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  }
  #
  if(is.null(pdf_file)==FALSE){
    pdf(pdf_file,width=(rr*2+ratio2)*1.5,height=(ratio1+rr)*1.5)
  }
  # get layout
  layout(matrix(c(rep(0,length.out=rr),rep(1,length.out=ratio2),rep(0,length.out=rr*1),
                  rep(c(rep(4,length.out=rr),rep(2,length.out=ratio2),rep(3,length.out=rr*1)),
                      length.out=ratio1*(ratio2+rr*2))),
                ncol=c(ratio2+rr*2),byrow=TRUE))
  ## plot 1
  par(mar=c(1.5,1.5,4,0))
  mm <- quantile(DE_profile,probs=c(0.0001,0.9999));
  mm <- max(abs(mm)); mm <- c(-mm,mm)
  y1 <- seq(mm[1],mm[2],length.out=5); y1 <- round(y1,1)
  unit <- n_gene/10; unit <- round(unit/100)*100
  x1 <- seq(0,n_gene,by=unit);x1 <- unique(x1); x1 <- c(x1,max(x1)+unit)
  par(usr=c(0,length(DE_profile),mm[1],mm[2]))
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab='',ylab='',bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
  pp <- par()$usr; rr <- (pp[2]-pp[1])/n_gene
  polygon(x=c(pp[1],c(1:n_gene)*rr+pp[1],pp[2]),y=c(0,DE_profile,0),col='grey',border='grey',xpd=TRUE,lwd=0.3)
  if(profile_trend=='pos2neg'){
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[2]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[1]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }else{
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[1]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[2]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }
  #segments(pp[1],pp[1],pp[2],pp[1],lwd=0.5)
  #segments(pp[1],min(y1),pp[1],max(y1),lwd=1.5)
  #text(-(pp[2]-pp[1])/50,y1,y1,adj=1,xpd=TRUE)
  #segments(-(pp[2]-pp[1])/100,y1,0,y1,lwd=1.5)
  axis(side=2,at=y1,labels=y1)
  mtext(side=2,line = 2.5,profile_col,cex=1)
  #axis(side=1,at=x1,labels=get_label_manual(x1))
  segments(pp[1],mm[1],pp[2],mm[1],xpd=TRUE)
  segments(x1*rr,mm[1]-(mm[2]-mm[1])/30,x1*rr,mm[1],xpd=TRUE)
  text(x1*rr,mm[1]-(mm[2]-mm[1])/10,get_label_manual(x1),adj=0.5,xpd=TRUE)
  ## plot2
  par(mar=c(2,1.5,2,0))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,n_gene),xaxt='n',yaxt='n')
  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1)
  segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
  yy2 <- seq(from=pp[3],to=pp[4],length.out=length(driver_list)+1)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey')
  # add columns
  use_target_list <- target_list[driver_list]
  if(target_col_type=='DE'){
    cc <- z2col(DE_profile,sig_thre=profile_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[5:9],blue_col=brewer.pal(9,'Blues')[5:9],
                col_max_thre=max(abs(DE_profile)))
    names(cc) <- names(DE_profile)
    cc[which(cc=='white')] <- 'light grey'
  }
  if(target_nrow==1){
    for(i in 1:length(driver_list)){
      t1 <- use_target_list[[driver_list[[i]]]]
      w0 <- which(names(DE_profile) %in% t1$target)
      w1 <- w0*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],lwd=1.5,
                   col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],lwd=1.5,
                   col=z2col(t1$spearman,sig_thre=0,col_max_thre=1,col_min_thre=0.01,
                             red_col = brewer.pal(9,'Reds')[7:9],blue_col=brewer.pal(9,'Blues')[7:9]))
        }
      }
    }
  }
  if(target_nrow==2){
    for(i in 1:length(driver_list)){
      t1 <- use_target_list[[driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(names(DE_profile) %in% t11)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col=brewer.pal(9,'Reds')[6],lwd=1.5)
          }
        }
      }
      w0 <- which(names(DE_profile) %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col=brewer.pal(9,'Blues')[6],lwd=1.5)
          }
        }
      }
    }
  }
  ## plot 3
  par(mar=c(2,0.5,2,2))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
  yy2 <- seq(from=pp[3],to=pp[4],length.out=length(driver_list)+1)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white',xpd=TRUE)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey',xpd=TRUE)
  abline(v=c(pp[1],(pp[1]+pp[2])/2,pp[2]))
  ## add text
  mm_min <- min(min(abs(driver_DA_Z[driver_list]),na.rm=TRUE)*0.9,min(abs(driver_DE_Z[driver_list]),na.rm=TRUE)*0.9)
  mm_min <- max(mm_min,Z_sig_thre)
  mm_max <- max(max(abs(driver_DA_Z[driver_list]),na.rm=TRUE)*1.1,max(abs(driver_DE_Z[driver_list]),na.rm=TRUE)*1.1)
  c1 <- z2col(driver_DA_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c2 <- z2col(driver_DE_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  for(i in 1:length(driver_list)){
    z1 <- driver_DA_Z[driver_list[i]]
    z2 <- driver_DE_Z[driver_list[i]]
    p1 <- get_z2p(z1)
    p2 <- get_z2p(z2)
    rect(xleft=pp[1],xright=(pp[1]+pp[2])/2,ybottom=yy2[i],ytop=yy2[i+1],col=c1[i],border='dark grey',xpd=TRUE)
    rect(xright=pp[2],xleft=(pp[1]+pp[2])/2,ybottom=yy2[i],ytop=yy2[i+1],col=c2[i],border='dark grey',xpd=TRUE)
    text(x=(pp[1]+(pp[1]+pp[2])/2)/2,y=(yy2[i]+yy2[i+1])/2,p1,adj=0.5)
    text(x=(pp[2]+(pp[1]+pp[2])/2)/2,y=(yy2[i]+yy2[i+1])/2,p2,adj=0.5)
  }
  textheight <- strheight('DA',units='user',cex=1.5)
  text((pp[1]+(pp[1]+pp[2])/2)/2,pp[4]+textheight,'DA',xpd=TRUE,cex=1.5)
  textheight <- strheight('DE',units='user',cex=1.5)
  text((pp[2]+(pp[1]+pp[2])/2)/2,pp[4]+textheight,'DE',xpd=TRUE,cex=1.5)
  ## plot 4
  par(mar=c(2,6,2,0.2))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  #yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1)
  #yy11 <- (yy1[1:(length(yy1)-1)]+yy1[2:length(yy1)])/2
  yy2 <- seq(from=pp[3],to=pp[4],length.out=length(driver_list)+1)
  yy22 <- (yy2[1:(length(yy2)-1)]+yy2[2:length(yy2)])/2
  dyy <- yy22[2]-yy22[1]
  text(show_label,x=(pp[1]+pp[2])/2,y=yy22,xpd=TRUE,adj=1)
  # add target size
  target_size <- do.call(rbind,lapply(use_target_list,function(x){
    x1 <- length(which(x$spearman>=0))
    x2 <- length(which(x$spearman<0))
    c(x1,x2)
  }))
  if(target_nrow==2){
    mm <- max(target_size)
    tt <- pp[2]-(pp[1]+pp[2])*0.55
    for(i in 1:length(driver_list)){
      rect(xleft=(pp[1]+pp[2])*0.55,xright=(pp[1]+pp[2])*0.55+target_size[i,1]/mm*tt,
           ybottom=yy22[i],ytop=yy22[i]+dyy*0.35,col=brewer.pal(9,'Reds')[5],border=NA)
      rect(xleft=(pp[1]+pp[2])*0.55,xright=(pp[1]+pp[2])*0.55+target_size[i,2]/mm*tt,
           ytop=yy22[i],ybottom=yy22[i]-dyy*0.35,col=brewer.pal(9,'Blues')[5],border=NA)
    }
    segments(x0=(pp[1]+pp[2])*0.55,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+(pp[1]+pp[2])*0.55
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  #
  if(target_nrow==1){
    target_size <- rowSums(target_size)
    mm <- max(target_size)
    tt <- pp[2]-(pp[1]+pp[2])*0.55
    for(i in 1:length(driver_list)){
      rect(xleft=(pp[1]+pp[2])*0.55,xright=(pp[1]+pp[2])*0.55+target_size[i]/mm*tt,
           ybottom=yy22[i]-dyy*0.2,ytop=yy22[i]+dyy*0.2,col='dark grey',border=NA)
    }
    segments(x0=(pp[1]+pp[2])*0.55,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+(pp[1]+pp[2])*0.55
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  ##
  if(is.null(pdf_file)==FALSE) dev.off()
  layout(1);
  return(TRUE)
}
###
# profile_trend: pos2neg, neg2pos
# target_col:black,or RdBu
draw.GSEA.NetBID.GS <- function(DE=NULL,name_col=NULL,profile_col=NULL,profile_trend='pos2neg',
                              sig_gs_list=NULL,gs_DA_Z=NULL,use_gs2gene=NULL,
                              top_gs_number=30,target_col='RdBu',
                              left_annotation='',right_annotation='',main='',
                              profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  if(!profile_col %in% colnames(DE)){
    message(sprintf('%s not in colnames of DE, please check and re-try!',profile_col))
    return(FALSE)
  }
  while(class(use_gs2gene[[1]])=='list'){
    nn <- unlist(lapply(use_gs2gene,names))
    use_gs2gene <- unlist(use_gs2gene,recursive = FALSE)
    names(use_gs2gene)<-nn
  }
  use_gs2gene <- use_gs2gene[sig_gs_list]
  if(is.null(name_col)==TRUE){
    DE <- cbind(DE,name=rownames(DE),stringsAsFactors=FALSE)
    name_col <- 'name'
  }
  if(is.null(names(gs_DA_Z))) names(gs_DA_Z) <- sig_gs_list
  if(length(sig_gs_list)>top_gs_number){
    gs_DA_Z <- gs_DA_Z[sig_gs_list]
    sig_gs_list <- sig_gs_list[order(abs(gs_DA_Z[sig_gs_list]),decreasing = TRUE)][1:top_gs_number]
  }
  if(profile_trend=='pos2neg')
    sig_gs_list <- sig_gs_list[order(gs_DA_Z[sig_gs_list],decreasing = FALSE)]
  else
    sig_gs_list <- sig_gs_list[order(gs_DA_Z[sig_gs_list],decreasing = TRUE)]
  show_label <- sig_gs_list
  gs_DA_Z <- gs_DA_Z[sig_gs_list]
  ###################
  ## calculate layout
  if(profile_trend=='pos2neg') DE <- DE[order(DE[,profile_col],decreasing = TRUE),] else DE <- DE[order(DE[,profile_col],decreasing = FALSE),]
  DE_profile <- DE[,profile_col]
  DE_profile_name <- DE[,name_col]
  use_gs2gene <- lapply(use_gs2gene,function(x)intersect(x,DE_profile_name))
  use_gs2gene <- use_gs2gene[sig_gs_list]
  #names(DE_profile) <- rownames(DE)
  n_gene <- length(DE_profile)
  n_driver <- length(sig_gs_list)
  gswidth <- max(strwidth(sig_gs_list,units='inches',cex=1))
  ratio1 <- ceiling(1.2*n_driver/15) ## profile height to rows
  ratio2 <- 4 ## width of profile to DA/DE
  rr1 <- ceiling(gswidth)/2
  rr2 <- 1
  #
  if(is.null(pdf_file)==FALSE){
    pdf(pdf_file,width=(rr1+rr2+ratio2)*1.5,height=(ratio1+rr2)*1.5)
  }
  # get layout
  layout(matrix(c(rep(0,length.out=rr1),rep(1,length.out=ratio2),rep(0,length.out=rr2*1),
                  rep(c(rep(4,length.out=rr1),rep(2,length.out=ratio2),rep(3,length.out=rr2*1)),
                      length.out=ratio1*(ratio2+rr1+rr2))),
                ncol=c(ratio2+rr1+rr2),byrow=TRUE))
  ## plot 1
  par(mar=c(1.5,1.5,4,0))
  mm <- quantile(DE_profile,probs=c(0.0001,0.9999));
  mm <- max(abs(mm)); mm <- c(-mm,mm)
  y1 <- seq(mm[1],mm[2],length.out=5); y1 <- round(y1,1)
  unit <- n_gene/10; unit <- round(unit/100)*100
  x1 <- seq(0,n_gene,by=unit);x1 <- unique(x1); x1 <- c(x1,max(x1)+unit)
  par(usr=c(0,length(DE_profile),mm[1],mm[2]))
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab='',ylab='',bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
  pp <- par()$usr; rr <- (pp[2]-pp[1])/n_gene
  polygon(x=c(pp[1],c(1:n_gene)*rr+pp[1],pp[2]),y=c(0,DE_profile,0),col='grey',border='grey',xpd=TRUE,lwd=0.3)
  if(profile_trend=='pos2neg'){
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[2]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[1]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }else{
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[1]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[2]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }
  axis(side=2,at=y1,labels=y1)
  mtext(side=2,line = 2.5,profile_col,cex=1)
  segments(pp[1],mm[1],pp[2],mm[1],xpd=TRUE)
  segments(x1*rr,mm[1]-(mm[2]-mm[1])/30,x1*rr,mm[1],xpd=TRUE)
  text(x1*rr,mm[1]-(mm[2]-mm[1])/10,get_label_manual(x1),adj=0.5,xpd=TRUE)
  ## plot2
  par(mar=c(2,1.5,2,0))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,n_gene),xaxt='n',yaxt='n')
  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1)
  segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
  yy2 <- seq(from=pp[3],to=pp[4],length.out=length(sig_gs_list)+1)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey')
  # add columns
  use_target_list <- use_gs2gene[sig_gs_list]
  cc <- z2col(DE_profile,sig_thre=profile_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[5:9],blue_col=brewer.pal(9,'Blues')[5:9],
              col_max_thre=max(abs(DE_profile)))
  cc[which(cc=='white')] <- 'light grey'
  for(i in 1:length(sig_gs_list)){
    t1 <- use_target_list[[sig_gs_list[i]]]
    w0 <- which(DE_profile_name %in% t1)
    w1 <- w0*rr+pp[1]
    if(target_col=='black'){
      segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],col='black',lwd=1)
    }else{
      segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],lwd=1.5,col=cc[w0])
    }
  }
  ## plot 3
  par(mar=c(2,0.5,2,2))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,1),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
  yy2 <- seq(from=pp[3],to=pp[4],length.out=length(sig_gs_list)+1)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white',xpd=TRUE)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey',xpd=TRUE)
  ## add text
  mm_min <- min(abs(gs_DA_Z[sig_gs_list]),na.rm=TRUE)*0.9
  mm_min <- max(mm_min,Z_sig_thre)
  mm_max <- max(abs(gs_DA_Z[sig_gs_list]),na.rm=TRUE)*1.1
  c1 <- z2col(gs_DA_Z[sig_gs_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  for(i in 1:length(sig_gs_list)){
    z1 <- gs_DA_Z[sig_gs_list[i]]
    p1 <- get_z2p(z1)
    rect(xleft=pp[1],xright=pp[2],ybottom=yy2[i],ytop=yy2[i+1],col=c1[i],border='dark grey',xpd=TRUE)
    text(x=(pp[1]+pp[2])/2,y=(yy2[i]+yy2[i+1])/2,p1,adj=0.5)
  }
  textheight <- strheight('DA',units='user',cex=1.5)
  text((pp[1]+pp[2])/2,pp[4]+textheight,'DA',xpd=TRUE,cex=1.5)
  ## plot 4
  par(mar=c(2,6,2,0.2))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  yy2 <- seq(from=pp[3],to=pp[4],length.out=length(sig_gs_list)+1)
  yy22 <- (yy2[1:(length(yy2)-1)]+yy2[2:length(yy2)])/2
  dyy <- yy22[2]-yy22[1]
  # add target size
  target_size <- unlist(lapply(use_gs2gene,length))
  mm <- max(target_size)
  rr <- ceiling(gswidth)*1.5
  tt_left <- pp[2]-(pp[2]-pp[1])/(1+rr)
  tt <- (pp[2]-pp[1])/(1+rr)
  text(show_label,x=tt_left-tt/25,y=yy22,xpd=TRUE,adj=1)
  for(i in 1:length(sig_gs_list)){
    rect(xleft=tt_left,xright=tt_left+target_size[i]/mm*tt,
         ybottom=yy22[i]-dyy*0.2,ytop=yy22[i]+dyy*0.2,col='dark grey',border=NA)
  }
  segments(x0=tt_left,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
  sst <- round(seq(0,mm,length.out=3))
  ss <- sst*tt/mm+tt_left
  segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
  text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
  text('Size',x=tt_left-tt/10,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  ##
  if(is.null(pdf_file)==FALSE) dev.off()
  layout(1);
  return(TRUE)
}
###
####### NetBID SINBA plot
# profile_trend: pos2neg, neg2pos
# target_nrow:1,2
# target_col:black,or RdBu
# target_col_type: DE,PN
# top_order: merge, diff
draw.GSEA.NetBID.SINBA <- function(DE=NULL,profile_col=NULL,profile_trend='pos2neg',
                                   seed_driver=NULL,partner_driver_list=NULL,
                                   seed_driver_label=seed_driver,partner_driver_label=partner_driver_list,
                                   driver_DA_Z=NULL,driver_DE_Z=NULL,target_list=NULL,
                                   DA_Z_merge=NULL,target_list_merge=NULL,
                                   top_driver_number=10,top_order='merge',target_nrow=2,target_col='RdBu',target_col_type='PN',
                                   left_annotation='',right_annotation='',main='',
                                   profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  if(!profile_col %in% colnames(DE)){
    message(sprintf('%s not in colnames of DE, please check and re-try!',profile_col))
    return(FALSE)
  }
  driver_list <- c(seed_driver,partner_driver_list)
  driver_list_gene <- gsub('(.*)_.*','\\1',driver_list)
  show_label <- c(seed_driver_label,partner_driver_label)
  # get names
  if(is.null(names(driver_DA_Z))) names(driver_DA_Z) <- driver_list
  if(is.null(names(driver_DE_Z))) names(driver_DE_Z) <- driver_list
  if(is.null(names(show_label))) names(show_label) <- driver_list
  if(is.null(names(DA_Z_merge))) names(DA_Z_merge) <- driver_list
  #
  ori_part_Z <- driver_DA_Z[partner_driver_list]
  ori_seed_Z <- driver_DA_Z[seed_driver]
  diff_Z <- 2*DA_Z_merge[partner_driver_list]-(ori_part_Z+ori_seed_Z)
  names(diff_Z) <- partner_driver_list
  driver_DA_Z <- driver_DA_Z[driver_list]
  driver_DE_Z <- driver_DE_Z[driver_list_gene]; names(driver_DE_Z) <- driver_list
  DA_Z_merge  <- DA_Z_merge[partner_driver_list]
  #
  if(top_order=='merge'){
    if(length(partner_driver_list)>top_driver_number){
      partner_driver_list <- partner_driver_list[order(abs(DA_Z_merge[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
    }
    if(profile_trend=='pos2neg')
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = FALSE)]
    else
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)]
  }else{
    if(length(partner_driver_list)>top_driver_number){
      partner_driver_list <- partner_driver_list[order(abs(diff_Z[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
    }
    if(profile_trend=='pos2neg')
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = FALSE)]
    else
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = TRUE)]
  }
  #
  driver_list <- c(seed_driver,partner_driver_list)
  show_label <- show_label[driver_list]
  driver_DA_Z <- driver_DA_Z[driver_list]
  driver_DE_Z <- driver_DE_Z[driver_list]
  diff_Z <- diff_Z[partner_driver_list]
  DA_Z_merge <- DA_Z_merge[partner_driver_list]
  ##############################################
  ## calculate layout
  DE_profile <- DE[,profile_col]
  names(DE_profile) <- rownames(DE)
  ##
  target_list <- lapply(target_list,function(x)x[which(x$target %in% names(DE_profile)),])
  target_list_merge <- lapply(target_list_merge,function(x)x[which(x$target %in% names(DE_profile)),])

  if(profile_trend=='pos2neg') DE_profile <- sort(DE_profile,decreasing = TRUE) else DE_profile <- sort(DE_profile)
  n_gene <- length(DE_profile)
  if(target_nrow==2){
    n_driver <- length(partner_driver_list)*4+2
    ratio1 <- ceiling(n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  } else {
    n_driver <- length(partner_driver_list)*2+1
    ratio1 <- ceiling(1.5*n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  }
  #
  if(is.null(pdf_file)==FALSE){
    pdf(pdf_file,width=(rr*2+ratio2)*2,height=(ratio1+rr)*2)
  }
  # get layout
  layout(matrix(c(rep(0,length.out=rr),rep(1,length.out=ratio2),rep(0,length.out=rr*1),
                  rep(c(rep(4,length.out=rr),rep(2,length.out=ratio2),rep(3,length.out=rr*1)),
                      length.out=ratio1*(ratio2+rr*2))),
                ncol=c(ratio2+rr*2),byrow=TRUE))
  ## plot 1
  par(mar=c(1.5,1.5,1.5,0))
  mm <- quantile(DE_profile,probs=c(0.0001,0.9999));
  mm <- max(abs(mm)); mm <- c(-mm,mm)
  y1 <- seq(mm[1],mm[2],length.out=5); y1 <- round(y1,1)
  unit <- n_gene/10; unit <- round(unit/100)*100
  x1 <- seq(0,n_gene,by=unit);x1 <- unique(x1); x1 <- c(x1,max(x1)+unit)
  par(usr=c(0,length(DE_profile),mm[1],mm[2]))
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab='',ylab='',bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
  pp <- par()$usr; rr <- (pp[2]-pp[1])/n_gene
  polygon(x=c(pp[1],c(1:n_gene)*rr+pp[1],pp[2]),y=c(0,DE_profile,0),col='grey',border='grey',xpd=TRUE,lwd=0.3)
  if(profile_trend=='pos2neg'){
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[2]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[1]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }else{
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[1]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[2]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }
  axis(side=2,at=y1,labels=y1)
  mtext(side=2,line = 2.5,profile_col,cex=1)
  segments(pp[1],mm[1],pp[2],mm[1],xpd=TRUE)
  segments(x1*rr,mm[1]-(mm[2]-mm[1])/50,x1*rr,mm[1],xpd=TRUE)
  text(x1*rr,mm[1]-(mm[2]-mm[1])/25,get_label_manual(x1),adj=0.5,xpd=TRUE)
  ## plot2
  par(mar=c(2,1.5,2,0))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,n_gene),xaxt='n',yaxt='n')

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow

  rect(xleft = pp[1],xright=pp[2],ybottom = yy1[length(yy1)-target_nrow],ytop=yy1[length(yy1)],border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.3)) ## for seed rows
  rect(xleft = pp[1],xright=pp[2],ybottom = yy2,ytop=yy4,border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.1)) ## for combine rows

  segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
  segments(x0=pp[1],x1=pp[2],y0=yy3,y1=yy3,lwd=1,col='dark grey')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1,col='black')
  segments(x0=pp[1],x1=pp[2],y0=yy1[length(yy1)-target_nrow],y1=yy1[length(yy1)-target_nrow],lwd=1.5,col='black')

  # shorten yy1
  dyy <- yy1[2]-yy1[1]
  yy11 <- yy1-dyy*0.3
  # add columns
  use_target_list <- target_list[driver_list]
  use_merge_target_list <- target_list_merge[partner_driver_list]

  if(target_col_type=='DE'){
    cc <- z2col(DE_profile,sig_thre=profile_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[5:9],blue_col=brewer.pal(9,'Blues')[5:9],
                col_max_thre=max(abs(DE_profile)))
    names(cc) <- names(DE_profile)
    cc[which(cc=='white')] <- 'light grey'
  }
  if(target_nrow==1){
    # for seed driver
    t1 <- use_target_list[[seed_driver]]
    w0 <- which(names(DE_profile) %in% t1$target)
    w1 <- w0*rr+pp[1]
    if(target_col=='black'){
      segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
    }else{
      if(target_col_type=='DE'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=cc[w0])
      }else{
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=c(brewer.pal(9,'Blues')[5],'white',brewer.pal(9,'Reds')[5])[sign(t1$spearman)+2])
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      w0 <- which(names(DE_profile) %in% t1$target)
      w1 <- w0*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=c(brewer.pal(9,'Blues')[5],'white',brewer.pal(9,'Reds')[5])[sign(t1$spearman)+2])
        }
      }
    }
    # for combine
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      w0 <- which(names(DE_profile) %in% t1$target)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(names(DE_profile) %in% t_over)
      w1_over <- w0_over*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=c(brewer.pal(9,'Blues')[5],'white',brewer.pal(9,'Reds')[5])[sign(t1$spearman)+2])
        }
      }
      points(w1_over,rep((yy11[2*i]+yy1[2*i])/2,length.out=length(w1_over)),pch='*',col='black')
    }
  }
  ###################
  if(target_nrow==2){
    # for seed driver
    t1 <- use_target_list[[seed_driver]]
    t11 <- t1[which(t1$spearman>=0),]$target
    t12 <- t1[which(t1$spearman<0),]$target
    w0 <- which(names(DE_profile) %in% t11)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=brewer.pal(9,'Reds')[9],lwd=1.5)
        }
      }
    }
    w0 <- which(names(DE_profile) %in% t12)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=brewer.pal(9,'Blues')[9],lwd=1.5)
        }
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(names(DE_profile) %in% t11)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=brewer.pal(9,'Reds')[5],lwd=1.5)
          }
        }
      }
      w0 <- which(names(DE_profile) %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=brewer.pal(9,'Blues')[5],lwd=1.5)
          }
        }
      }
    }
    # for each partner driver + seed combination
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(names(DE_profile) %in% t11)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(names(DE_profile) %in% t_over) ## setdiff(t_over,names(DE_profile)) !!!
      w1_over <- w0_over*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=brewer.pal(9,'Reds')[5],lwd=1.5)
          }
        }
      }
      points(intersect(w1_over,w1),rep((yy11[4*i-1]+yy1[4*i-1])/2,length.out=length(intersect(w1_over,w1))),pch='*',col='black')
      w0 <- which(names(DE_profile) %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=brewer.pal(9,'Blues')[5],lwd=1.5)
          }
        }
      }
      points(intersect(w1_over,w1),rep((yy11[4*i-2]+yy1[4*i-2])/2,length.out=length(intersect(w1_over,w1))),pch='*',col='black')
    }
    ####
  }
  ###################
  ## plot 3
  par(mar=c(2,0.5,2,2))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,3),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow
  xx1 <- seq(pp[1],pp[2],length.out=4)

  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white',xpd=TRUE)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey',xpd=TRUE)
  abline(v=xx1)

  ## add text
  mm_min <- min(min(abs(driver_DA_Z[partner_driver_list]),na.rm=TRUE)*0.9,min(abs(driver_DE_Z[partner_driver_list]),na.rm=TRUE)*0.9,
                min(abs(diff_Z[partner_driver_list]),na.rm=TRUE)*0.9,min(abs(DA_Z_merge[partner_driver_list]),na.rm=TRUE)*0.9)
  mm_min <- max(mm_min,Z_sig_thre)
  mm_max <- max(max(abs(driver_DA_Z[partner_driver_list]),na.rm=TRUE)*1.1,max(abs(driver_DE_Z[partner_driver_list]),na.rm=TRUE)*1.1,
                max(abs(diff_Z[partner_driver_list]),na.rm=TRUE)*1.1,max(abs(DA_Z_merge[partner_driver_list]),na.rm=TRUE)*1.1)
  c1 <- z2col(driver_DA_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c2 <- z2col(driver_DE_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c3 <- z2col(diff_Z[partner_driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c4 <- z2col(DA_Z_merge[partner_driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)

  # for seed driver
  yy1<-yy3
  z1 <- driver_DA_Z[seed_driver]
  z2 <- driver_DE_Z[seed_driver]
  p1 <- get_z2p(z1)
  p2 <- get_z2p(z2)
  rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[length(yy1)-1],ytop=yy1[length(yy1)],col=c1[1],border='dark grey',xpd=TRUE)
  rect(xright=xx1[3],xleft=xx1[2],ybottom=yy1[length(yy1)-1],ytop=yy1[length(yy1)],col=c2[1],border='dark grey',xpd=TRUE)
  text(x=sum(xx1[1:2])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,p1,adj=0.5)
  text(x=sum(xx1[3:4])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,p2,adj=0.5)
  text(x=sum(xx1[2:3])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,'-',adj=0.5)

  # for partner driver
  for(i in 1:length(partner_driver_list)){
    z1 <- driver_DA_Z[partner_driver_list[i]]
    z2 <- driver_DE_Z[partner_driver_list[i]]
    z3 <- diff_Z[partner_driver_list[i]]
    z4 <- DA_Z_merge[partner_driver_list[i]]
    p1 <- get_z2p(z1)
    p2 <- get_z2p(z2)
    p3 <- format(z3,digits=3)
    p4 <- get_z2p(z4)
    rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[2*i],ytop=yy1[2*i+1],col=c1[i+1],border='dark grey',xpd=TRUE) ## DA_Z
    rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[2*i-1],ytop=yy1[2*i],col=c4[i],border='dark grey',xpd=TRUE) ## merge_Z
    rect(xleft=xx1[3],xright=xx1[4],ybottom=yy2[i],ytop=yy2[i+1],col=c2[i+1],border='dark grey',xpd=TRUE) ## DE
    rect(xleft=xx1[2],xright=xx1[3],ybottom=yy2[i],ytop=yy2[i+1],col=c3[i],border='dark grey',xpd=TRUE) ## delta Z
    text(x=sum(xx1[1:2])/2,y=(yy1[2*i]+yy1[2*i+1])/2,p1,adj=0.5) ## DA_Z
    text(x=sum(xx1[1:2])/2,y=(yy1[2*i-1]+yy1[2*i])/2,p4,adj=0.5) ## merge_Z
    text(x=sum(xx1[3:4])/2,y=(yy2[i]+yy2[i+1])/2,p2,adj=0.5) ## DE
    text(x=sum(xx1[2:3])/2,y=(yy2[i]+yy2[i+1])/2,p3,adj=0.5) ## delta z
  }

  textheight <- strheight('DA',units='user',cex=1.5)
  text(sum(xx1[1:2])/2,pp[4]+textheight,'DA',xpd=TRUE,cex=1.5)
  textheight <- strheight('DE',units='user',cex=1.5)
  text(sum(xx1[3:4])/2,pp[4]+textheight,'DE',xpd=TRUE,cex=1.5)
  textheight <- strheight('deltaZ',units='user',cex=1.5)
  text(sum(xx1[2:3])/2,pp[4]+textheight,'deltaZ',xpd=TRUE,cex=1.5)

  ## plot 4
  par(mar=c(2,6,2,0.2))
  plot(1,col='white',xlab='',ylab='',xlim=c(0,2),xaxt='n',yaxt='n',bty='n')

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow
  xx1 <- seq(pp[1],pp[2],length.out=4)

  yy22 <- (yy2[1:(length(yy2)-1)]+yy2[2:length(yy2)])/2
  dyy22 <- yy22[2]-yy22[1]

  yy33 <- (yy3[1:(length(yy3)-1)]+yy3[2:length(yy3)])/2
  dyy33 <- yy33[2]-yy33[1]

  xleft <- pp[1]+(pp[2]-pp[1])*0.55
  tt <- pp[2]-xleft

  text(show_label[1],x=pp[1]+(pp[1]+pp[2])*0.53,y=(yy3[length(yy3)-1]+yy3[length(yy3)])/2,xpd=TRUE,adj=1,cex=1.2) ## label for seed
  text(show_label[2:length(show_label)],x=pp[1]+(pp[1]+pp[2])*0.52,y=yy22,xpd=TRUE,adj=1) ## label for partner

  # add target size
  use_target_list <- target_list[driver_list]
  use_merge_target_list <- target_list_merge[partner_driver_list]
  target_size <- do.call(rbind,lapply(use_target_list,function(x){
    x1 <- length(which(x$spearman>=0))
    x2 <- length(which(x$spearman<0))
    c(x1,x2)
  }))
  merge_target_size <- do.call(rbind,lapply(use_merge_target_list,function(x){
    x1 <- length(which(x$spearman>=0))
    x2 <- length(which(x$spearman<0))
    c(x1,x2)
  }))
  # for seed driver
  if(target_nrow==2){
    mm <- max(merge_target_size)
    i <- length(yy1)-1
    rect(xleft=xleft,xright=xleft+target_size[1,1]/mm*tt,
         ybottom=yy1[i],ytop=yy1[i]+dyy22/2*0.35,col=brewer.pal(9,'Reds')[5],border=NA)
    rect(xleft=xleft,xright=xleft+target_size[1,2]/mm*tt,
         ytop=yy1[i],ybottom=yy1[i]-dyy22/2*0.35,col=brewer.pal(9,'Blues')[5],border=NA)
  }else{
    target_size <- rowSums(target_size)
    merge_target_size <- rowSums(merge_target_size)
    mm <- max(merge_target_size)
    i <- length(yy1)
    rect(xleft=xleft,xright=xleft+target_size[1]/mm*tt,
         ybottom=(yy1[i]+yy1[i-1])/2-dyy22*0.2/2,ytop=(yy1[i]+yy1[i-1])/2+dyy22*0.2/2,col='dark grey',border=NA)
  }
  # for partner
  if(target_nrow==2){
    mm <- max(merge_target_size)
    for(i in 1:length(partner_driver_list)){
      rect(xleft=xleft,xright=xleft+target_size[i+1,1]/mm*tt,
           ybottom=yy33[2*i],ytop=yy33[2*i]+dyy33*0.35,col=brewer.pal(9,'Reds')[5],border=NA)
      rect(xleft=xleft,xright=xleft+target_size[i+1,2]/mm*tt,
           ytop=yy33[2*i],ybottom=yy33[2*i]-dyy33*0.35,col=brewer.pal(9,'Blues')[5],border=NA)
      # merge
      rect(xleft=xleft,xright=xleft+merge_target_size[i,1]/mm*tt,
           ybottom=yy33[2*i-1],ytop=yy33[2*i-1]+dyy33*0.35,col=brewer.pal(9,'Reds')[5],border=NA)
      rect(xleft=xleft,xright=xleft+merge_target_size[i,2]/mm*tt,
           ytop=yy33[2*i-1],ybottom=yy33[2*i-1]-dyy33*0.35,col=brewer.pal(9,'Blues')[5],border=NA)
    }
    segments(x0=xleft,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+xleft
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  #
  if(target_nrow==1){
    #target_size <- rowSums(target_size)
    #merge_target_size <- rowSums(merge_target_size)
    mm <- max(merge_target_size)
    for(i in 1:length(partner_driver_list)){
      rect(xleft=xleft,xright=xleft+target_size[i+1]/mm*tt,ybottom=yy33[i*2]-dyy33*0.2,ytop=yy33[i*2]+dyy33*0.2,col='dark grey',border=NA) #each
      rect(xleft=xleft,xright=xleft+merge_target_size[i]/mm*tt,ybottom=yy33[i*2-1]-dyy33*0.2,ytop=yy33[i*2-1]+dyy33*0.2,col='dark grey',border=NA) #merge
    }
    segments(x0=xleft,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+xleft
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  ## add lines
  segments(x0=xleft,x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='grey',xpd=TRUE)
  abline(v=xleft,col='grey')
  #abline(v=pp[2],col='grey')
  ## test for significant overlap
  total_possible_target <- unique(unlist(lapply(target_list,function(x)x$target)))
  for(i in 1:length(partner_driver_list)){
    res1 <- test.targetNet.overlap(seed_driver,partner_driver_list[i],
                                   target1=intersect(use_target_list[[seed_driver]]$target,names(DE_profile)),
                                   target2=intersect(use_target_list[[partner_driver_list[i]]]$target,names(DE_profile)),
                                   total_possible_target = total_possible_target)
    pv <- format(res1[1],digits=2,scientific=TRUE)
    ov <- round(res1[3])
    if(res1[1]<0.05){
      text(sprintf('Overlap:%d, P value:%s',ov,pv),x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.8,col='dark red')
    }else{
      if(ov==0){
        text('No overlap',x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.7)
      }else{
        text(sprintf('Overlap:%d, P value:%s',ov,pv),x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.7)
      }
    }
  }
  ##
  if(is.null(pdf_file)==FALSE) dev.off()
  layout(1);
  return(TRUE)
}
############## SINBA related functions
## merge target list
merge_target_list <- function(driver1=NULL,driver2=NULL,target_list=NULL){
  t1 <- target_list[driver1][[1]]
  t2 <- target_list[driver2][[1]]
  ov <- intersect(t1$target,t2$target)
  rownames(t1) <- t1$target
  rownames(t2) <- t2$target
  if(length(ov)==0){
    t_out <- rbind(t1,t2)
  }else{
    t_out1 <- rbind(t1[setdiff(rownames(t1),ov),],t2[setdiff(rownames(t2),ov),])
    t_out2 <- list()
    for(each_ov in ov){
      x1 <- t1[each_ov,2:3]
      x2 <- t2[each_ov,2:3]
      if(x1[1]>x2[1]) t_out2[[each_ov]] <- t1[each_ov,] else t_out2[[each_ov]] <- t2[each_ov,]
    }
    t_out2 <- do.call(rbind,t_out2)
    t_out <- rbind(t_out1,t_out2)
  }
  return(t_out)
}
#
### draw boxplot + stripchart for multiple groups
draw.categoryValue <- function(ac_val=NULL,exp_val=NULL,use_obs_class=NULL,category_color=NULL,
                               stripchart_color='black',strip_cex=1,class_order=NULL,class_srt=90,class_cex=1,pdf_file=NULL,
                               main_ac='',main_exp=''){
  if(is.null(class_order)){
    class_order <- sort(unique(use_obs_class))
  }
  if(is.null(category_color)==TRUE){
    class_col <- get.class.color(class_order)
    class_col1 <- get_transparent(class_col,0.5)
  }else{
    class_col1 <- category_color
  }
  stripchart_color <- get_transparent(stripchart_color,0.7)
  c1 <- 0
  if(is.null(ac_val)==FALSE){c1 <- c1+1}
  if(is.null(exp_val)==FALSE){c1 <- c1+1}
  if(is.null(pdf_file)==FALSE){
    labelWidth <- max(strwidth(class_order,'inches',cex=class_cex)*sin(class_srt*pi/180))
    hh <- 5+1.5+labelWidth
    if(c1==1) pdf(pdf_file,width=1.5+3,height=hh)
    if(c1==2) pdf(pdf_file,width=1.5+3*2,height=hh)
  }
  if(c1>1) layout(t(matrix(1:c1)))
  par(mai=c(labelWidth+0.5,1,1,0.5))
  if(is.null(ac_val)==FALSE){
    ddf <- data.frame(data=ac_val,class=factor(use_obs_class,levels=class_order))
    a <- boxplot(data~class,data=ddf,ylab='Activity Value',col=class_col1,outline=FALSE,border='dark grey',cex.lab=1.2,names=NA,bty='n',
                 ylim=c(min(ddf$data),max(ddf$data)),main=main_ac)
    text(1:length(class_order),par()$usr[3]-(par()$usr[4]-par()$usr[3])/20,adj=0.5+class_srt/180,class_order,srt=class_srt,xpd=TRUE,cex=class_cex)
    stripchart(data~class,data=ddf,add=TRUE,pch=16,method='jitter',vertical=TRUE,col=stripchart_color,cex=strip_cex)
  }
  if(is.null(exp_val)==FALSE){
    ddf <- data.frame(data=exp_val,class=factor(use_obs_class,levels=class_order))
    a <- boxplot(data~class,data=ddf,col=class_col1,ylab='Expression Value',outline=FALSE,border='dark grey',cex.lab=1.2,names=NA,bty='n',
                 ylim=c(min(ddf$data),max(ddf$data)),main=main_exp)
    text(1:length(class_order),par()$usr[3]-(par()$usr[4]-par()$usr[3])/20,adj=0.5+class_srt/180,class_order,srt=class_srt,xpd=TRUE,cex=class_cex)
    stripchart(data~class,data=ddf,add=TRUE,pch=16,method='jitter',vertical=TRUE,col=stripchart_color,cex=strip_cex)
  }
  layout(1)
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}
### simple functions
get_transparent <- function(x,alpha=0.1){
  rgb(t(col2rgb(x)/255),alpha=alpha)
}
get_label_manual <- function(x){
  x1 <- sapply(x,function(x2){
    x3 <- unlist(strsplit(as.character(x2),''))
    x4 <- length(x3)%/%3 ## add number
    if(x4>0){
      pp <- length(x3)-seq(1,x4)*3; x3[pp] <- paste0(x3[pp],','); paste(x3,collapse='')
    }else{
      x2
    }
  })
  unlist(x1)
}
##### network visualization
## draw network structure
draw.targetNet <- function(source_label='',source_z=NULL,edge_score=NULL,label_cex=0.7,init.angle=0,pdf_file=NULL){
  edge_score<- sort(edge_score)
  tmp1 <- sapply(unique(names(edge_score)),function(x){
    x1 <- edge_score[which(names(edge_score)==x)]
    x1[which.max(abs(x1))]
  })
  names(tmp1) <- unique(names(edge_score))
  edge_score <- tmp1
  edge_score<- sort(edge_score)
  g1 <- names(edge_score)
  ec <- z2col(edge_score*100,sig_thre=0,n_len=length(edge_score));names(ec) <- names(edge_score)
  ec <- get_transparent(ec,alpha=0.8)
  ew <- 2*(abs(edge_score)-min(abs(edge_score)))/(max(abs(edge_score))-min(abs(edge_score)))+0.5; names(ew) <- names(edge_score)
  t2xy <- function(tt,radius=1) {
    t2p <- pi*2 * tt + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  geneWidth <- max(strwidth(g1,'inches',cex=label_cex))
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=6+2*geneWidth,height=6+2*geneWidth)
  par(mai=c(1,1,1,1))
  plot(1,xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white',xlab='',ylab='',xaxt='n',yaxt='n')
  pp <- par()$usr
  tt <- seq(-0.5,0.5,length.out=length(g1)+1)[-1];init.angle <- init.angle;p1<-t2xy(tt,radius=0.8);
  for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g1[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
  p1<-t2xy(tt,radius=0.78);#segments(x0=0,y0=0,x1=p1$x,y1=p1$y,col=ec[g1],lwd=1);
  p2<-t2xy(tt,radius=0.77);arrows(x0=0,y0=0,x1=p2$x,y1=p2$y,col=ec,lwd=ew,angle=10,length=0.1);
  points(p1$x,p1$y,pch=16,col='dark grey')
  geneWidth <- strwidth(source_label,'inches',cex=0.8)
  if(is.null(source_z)==TRUE)
    points(0,0,col='light grey',cex=geneWidth*18,pch=16)
  else
    points(0,0,col=z2col(source_z),cex=geneWidth*18,pch=16)
  #points(0,0,col='light grey',cex=14,pch=16)
  text(0,0,source_label,adj=0.5,xpd=TRUE)
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}
##
draw.targetNet.TWO <- function(source1_label='',source2_label='',
                               source1_z=NULL,source2_z=NULL,
                               edge_score1=NULL,edge_score2=NULL,
                               label_cex=0.7,pdf_file=NULL,
                               total_possible_target=NULL,show_test=FALSE){
  tmp1 <- sapply(unique(names(edge_score1)),function(x){
    x1 <- edge_score1[which(names(edge_score1)==x)];x1[which.max(abs(x1))]
  })
  names(tmp1) <- unique(names(edge_score1));edge_score1 <- tmp1
  tmp1 <- sapply(unique(names(edge_score2)),function(x){
    x1 <- edge_score2[which(names(edge_score2)==x)];x1[which.max(abs(x1))]
  })
  names(tmp1) <- unique(names(edge_score2));edge_score2 <- tmp1
  edge_score1<- sort(edge_score1,decreasing = FALSE)
  edge_score2<- sort(edge_score2,decreasing = TRUE)
  g12 <- intersect(names(edge_score1),names(edge_score2))
  g1  <- setdiff(names(edge_score1),names(edge_score2))
  g2  <- setdiff(names(edge_score2),names(edge_score1))
  ec1 <- z2col(edge_score1*100,sig_thre=0,n_len=length(edge_score1));names(ec1) <- names(edge_score1)
  ec2 <- z2col(edge_score2*100,sig_thre=0,n_len=length(edge_score2));names(ec2) <- names(edge_score2)
  ew1 <- 2*(abs(edge_score1)-min(abs(edge_score1)))/(max(abs(edge_score1))-min(abs(edge_score1)))+0.5; names(ew1) <- names(edge_score1)
  ew2 <- 2*(abs(edge_score2)-min(abs(edge_score2)))/(max(abs(edge_score2))-min(abs(edge_score2)))+0.5; names(ew2) <- names(edge_score2)
  t2xy <- function(tt,radius=1) {
    t2p <- pi*2 * tt + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  geneWidth <- max(strwidth(g1,'inches',cex=label_cex))
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=10+4*geneWidth,height=8+2*geneWidth)
  plot(1,xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white',xlab='',ylab='',xaxt='n',yaxt='n')
  par(mai=c(1,1,1,1))
  if(length(g1)>0){
    tt <- seq(-0.225,0.225,length.out=length(g1));init.angle <- -180;p1<-t2xy(tt,radius=0.8);
    for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g1[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    p1<-t2xy(tt,radius=0.78);
    #segments(x0=-0.2,y0=0,x1=p1$x,y1=p1$y,col=ec1[g1],lwd=1)
    p2<-t2xy(tt,radius=0.76);arrows(x0=-0.2,y0=0,x1=p1$x,y1=p1$y,col=ec1[g1],lwd=ew1[g1],angle=10,length=0.1);
    points(p1$x,p1$y,pch=16,col='dark grey')
  }
  if(length(g2)>0){
    tt <- seq(-0.225,0.225,length.out=length(g2));init.angle <- 0;p1<-t2xy(tt,radius=0.8);
    for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g2[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    p1<-t2xy(tt,radius=0.78);
    #segments(x0=0.2,y0=0,x1=p1$x,y1=p1$y,col=ec2[g2],lwd=1)
    p2<-t2xy(tt,radius=0.76);arrows(x0=0.2,y0=0,x1=p1$x,y1=p1$y,col=ec2[g2],lwd=ew2[g2],angle=10,length=0.1)
    points(p1$x,p1$y,pch=16,col='dark grey')
  }
  if(length(g12)>0){
    tt <- seq(min(0.1*length(g12),0.7),-min(0.1*length(g12),0.7),length.out=length(g12));
    #segments(x0=-0.2,y0=0,x1=0,y1=tt,col=ec1[g12],lwd=ew1[g12])
    #segments(x0=0.2,y0=0,x1=0,y1=tt,col=ec2[g12],lwd=ew2[g12])
    arrows(x0=-0.2,y0=0,x1=0,y1=tt,col=ec1[g12],lwd=ew1[g12],angle=10,length=0.1)
    arrows(x0=0.2,y0=0,x1=0,y1=tt,col=ec2[g12],lwd=ew2[g12],angle=10,length=0.1)
    boxtext(0,tt,labels=g12,col.bg=get_transparent('light grey',0.3),cex=label_cex)
    #text(0,tt,g12,adj=0.5,cex=label_cex)
  }
  geneWidth1 <- strwidth(source1_label,'inches',cex=0.8)
  geneWidth2 <- strwidth(source2_label,'inches',cex=0.8)

  if(is.null(source2_z)==TRUE)
    points(0.2,0,col='light grey',cex=geneWidth1*18,pch=16)
  else
    points(0.2,0,col=z2col(source2_z),cex=geneWidth1*18,pch=16)

  text(0.2,0,source2_label,adj=0.5,cex=0.8)
  if(is.null(source1_z)==TRUE)
    points(-0.2,0,col='light grey',cex=geneWidth1*18,pch=16)
  else
    points(-0.2,0,col=z2col(source1_z),cex=geneWidth1*18,pch=16)

  text(-0.2,0,source1_label,adj=0.5,cex=0.8)
  # fisher test for target
  if(is.null(total_possible_target)==FALSE & show_test==TRUE){
    res <- test.targetNet.overlap(source1_label,source2_label,names(edge_score1),names(edge_score2),total_possible_target)
    if(is.null(pdf_file)==FALSE) dev.off()
    return(res)
  }
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}
test.targetNet.overlap <- function(source1_label,source2_label,target1,target2,total_possible_target){
  t1  <- unique(target1)
  t2  <- unique(target2)
  print(sprintf('%s has %d unique targets !',source1_label,length(t1)))
  print(sprintf('%s has %d unique targets !',source2_label,length(t2)))
  n11 <- length(intersect(t1,t2))
  n12 <- length(setdiff(t1,t2))
  n21 <- length(setdiff(t2,t1))
  if(class(total_possible_target) %in% c('integer','numeric')){
    n22 <- total_possible_target-length(union(t1,t2))
  }else{
    n22 <- length(setdiff(total_possible_target,c(t1,t2)))
  }
  mm  <- cbind(c(n11,n21),c(n12,n22))
  ft  <- fisher.test(mm)$p.value
  or  <- n11/n12/(n21/n22)
  rownames(mm) <- c(sprintf('In %s target',source1_label),sprintf('Not in %s target',source1_label))
  colnames(mm) <- c(sprintf('In %s target',source2_label),sprintf('Not in %s target',source2_label))
  print(mm)
  res <- c('P.Value'=ft,'Odds_Ratio'=or,'Intersected_Number'=mm[1,1])
  return(res)
}
##
####
# plot_pattern: no, random_one, all
get.spByGene <- function(igraph_obj=NULL,driver_list=NULL,target_list=NULL,mode='all',transfer_tab=NULL,plot_pattern='random_one',...){
  ori_driver_list <- driver_list
  ori_target_list <- target_list
  if(is.null(driver_list) & length(target_list)==1){
    driver_list <- unlist(lapply(neighborhood(igraph_obj,nodes=target_list,mode='in'),names))
  }
  if(is.null(target_list) & length(driver_list)==1){
    target_list <- unlist(lapply(neighborhood(igraph_obj,nodes=driver_list,mode='out'),names))
  }
  driver_list <- unique(intersect(driver_list,names(V(igraph_obj))))
  target_list <- unique(intersect(target_list,names(V(igraph_obj))))
  if(length(driver_list)==0 & length(target_list)==0){
    message('No gene included, please check and re-try!');return(FALSE)
  }
  #r1 <- all_shortest_paths(igraph_obj,from=use_list,to=use_list,mode=mode)$res
  #r2 <- unique(unlist(lapply(r1,names)))
 print(driver_list)
 print(target_list)

 if(length(driver_list)==0) driver_list <- target_list
 if(length(target_list)==0) target_list <- driver_list

  all_r1 <- list()
  for(d in driver_list){
    for(t in target_list){
      if(d!=t) all_r1[[sprintf('%s %s',d,t)]] <- all_shortest_paths(igraph_obj,from=d,to=t,mode=mode)$res
    }
  }
  #print(all_r1)

  if(plot_pattern!='all') r2 <- unique(unlist(lapply(all_r1,function(x)names(x[[1]]))))
  if(plot_pattern=='all') r2 <- unique(unlist(lapply(all_r1,function(x)names(unlist(x)))))

  all_sp <- induced_subgraph(igraph_obj,r2)

  ## get layout matrix,x,y
  net <- all_sp
  if(plot_pattern!='no'){
    plot.spByGene(net,driver_list=ori_driver_list,target_list=ori_target_list,transfer_tab=transfer_tab,...)
  }
  return(list(all_sp=all_r1,use_sp=net))
}
### plot igraph
plot.spByGene <- function(net,driver_list=NULL,target_list=NULL,transfer_tab=NULL,vertex.size=25,vertex.label.cex=1,
                          vertex.frame.color="white",vertex.label.color="black",edge.color='black',
                          layout='auto',...){
  n1 <- names(V(net))
  cc <- brewer.pal(11,'Set3')
  c1 <- ifelse(n1 %in% driver_list,cc[3],ifelse(n1 %in% target_list,cc[1],cc[9]))
  names(c1) <- n1
  d1 <- degree(net,mode='out')
  n2 <- n1[order(d1,decreasing = TRUE)]
  y <- ifelse(d1[n2]>0,2,1)
  x <- c(seq(0,1,length.out=length(y[which(y==2)])),seq(0,1,length.out=length(y[which(y==1)])))
  l <- cbind(x,y)
  use_label <- names(V(net))
  if(is.null(transfer_tab)==FALSE){
    transfer_tab <- as.matrix(transfer_tab)
    tt1 <- rbind(cbind(paste0(transfer_tab[,1],'_TF'),paste0(transfer_tab[,2],'_TF')),
                 cbind(paste0(transfer_tab[,1],'_SIG'),paste0(transfer_tab[,2],'_SIG')),
                 transfer_tab[,1:2])
    tt1 <- unique(tt1)
    rownames(tt1) <- tt1[,1]
    w1 <- which(use_label %in% rownames(tt1))
    use_label[w1] <- tt1[use_label[w1],2]
  }
  if(class(layout) != 'function'){
    plot.igraph(net,vertex.label=use_label,vertex.color=c1,layout=l,vertex.size=vertex.size,vertex.frame.color=vertex.frame.color,
                vertex.label.color=vertex.label.color,vertex.label.cex=vertex.label.cex,
                edge.color=edge.color,...)
  }else{
    plot.igraph(net,vertex.label=use_label,vertex.color=c1,layout=layout,vertex.size=vertex.size,vertex.frame.color=vertex.frame.color,
                vertex.label.color=vertex.label.color,vertex.label.cex=vertex.label.cex,
                edge.color=edge.color,...)
  }
  return(TRUE)
}
##
########################## analysis part functions
## get_target, save network data into list
## get target genes for each source
get_net2target_list <- function(net_dat=NULL) {
  all_source <- unique(net_dat$source)
  all_target <- lapply(all_source, function(x) {
    n1 <- net_dat[which(net_dat$source == x), c('target', 'MI', 'spearman')]
    rownames(n1) <- n1$target
    return(n1)
  })
  names(all_target) <- all_source
  return(all_target)
}
get.SJAracne.network <- function(network_file=NULL){
  if(is.null(network_file)){
    message('No input network file, please check and re-try !');return(FALSE)
  }
  net_dat      <- read.delim(file=network_file,stringsAsFactors = FALSE)
  target_list  <- get_net2target_list(net_dat)
  igraph_obj <- graph_from_data_frame(net_dat[,c('source','target','MI')],directed=TRUE)
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}
## functions to plot QC for network
plot.network.QC <- function(igraph_obj,outdir=NULL,prefix=''){
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message(paste0("The output directory: \"", outdir, "\" is created!"))
  }else
    message(paste0("The output will overwrite the files in directory: \"",outdir,"\""))
  if(class(igraph_obj)!='igraph'){
    message('Should input igraph object ! ');return(FALSE);
  }
  res_file <- sprintf('%s/%snetwork_info.pdf',outdir,prefix)
  pdf(res_file);
  plot(density(degree(igraph_obj)),xlab='Degree',main='Density plot for degree');
  hist(degree(igraph_obj),xlab='Degree',main='Histogram of degree')
  check_scalefree(igraph_obj);
  dev.off()
}
# functions to check the network
check_scalefree <- function(gr1) {
  fp1 <- degree_distribution(gr1)
  dd <- as.data.frame(cbind(k = 1:max(degree(gr1)), pk = fp1[-1]))
  r2 <-
    lm(log(dd$pk + 1 / length(V(gr1))) ~ log(dd$k))
  r3 <- summary(r2)$adj.r.squared
  plot(pk ~ k,data = dd,log = 'xy',main = sprintf('R2:%f', r3))
  return(r3)
}

########################## server-run functions (need local run version)
get_server_path <- function(x){
  gsub(RP.main_dir, RP.bash.main_dir, x)
}
## run MICA clustering, row(sample), col(gene)
# run MICA with mat variable (first step output into file)
SJ.MICA.prepare <- function(mat,outdir = '.',prjname = NULL,all_k=NULL,retransformation="False",perplexity=30) {
  if(is.null(prjname)){
    message('prjname should be input, please check and re-try !');return(FALSE)
  }
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  if(is.null(all_k)){
    all_k <- 2:12
  }
  input_exp <- file.path(outdir, 'input.exp')
  mat1 <- t(mat)
  mat1 <- cbind(rownames(mat1), mat1)
  colnames(mat1)[1] <- 'id'
  write.table(mat1,file = input_exp,row.names = FALSE,col.names = TRUE,sep = '\t',quote = FALSE)
  message(sprintf('Output exp file into %s',input_exp))
  run_file <- file.path(outdir, 'run_MICA_S1.sh')
  cmd <- sprintf('%s\nsh %s/run_MIE.sh %s %s/ %s',RP.load,RP.MICA.main,prjname,outdir,input_exp)
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  run_file <- file.path(outdir, 'run_MICA_S2.sh')
  cmd <-
    sprintf('%s\nsh %s/run_MICA.sh %s %s/ %s %s %s \"%s\"',RP.load,RP.MICA.main,prjname,
      outdir,input_exp,retransformation,perplexity,paste(all_k,collapse = ' ')
    )
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  return(TRUE)
}
# run MICA with existing file
SJ.MICA.prepare.withfile <- function(input_exp,outdir = '.',prjname = NULL,all_k=NULL,retransformation="False",perplexity=30) {
  if(is.null(prjname)){
    message('prjname should be input, please check and re-try !');return(FALSE)
  }
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  run_file <- file.path(outdir, 'run_MICA_S1.sh')
  cmd <- sprintf('%s\nsh %s/run_MIE.sh %s %s/ %s',RP.load,RP.MICA.main,prjname,outdir,input_exp)
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  run_file1 <- run_file
  run_file <- file.path(outdir, 'run_MICA_S2.sh')
  cmd <- sprintf('%s\nsh %s/run_MICA.sh %s %s/ %s %s %s \"%s\"',RP.load,RP.MICA.main,prjname,outdir,
      input_exp,retransformation,perplexity,paste(all_k,collapse = ' '))
  cmd <- get_server_path(cmd)
  cat(cmd, file = run_file, sep = '\n')
  run_file <- get_server_path(run_file)
  print(sprintf('Check %s', run_file))
  return(paste0('sh ',c(run_file1,run_file)))
}
## prepare for SJaracne run on SJ server
## prepare SJAracne dataset for net-dataset
# project_name expression_matrix hub_genes outdir
SJ.SJAracne.prepare <-
  function(eset,use.samples = rownames(pData(eset)),TF_list=NULL,SIG_list=NULL,
           SJAR.project_name = '',
           SJAR.main_dir = NULL,
           SJAR.bash.main_dir = NULL,
           mem = 40960,IQR.thre=0.5,IQR.loose_thre=0.1) {
    if(is.null(TF_list)==TRUE){
      message('Empty TF_list, please check and re-try !');return(FALSE)
    }
    if(is.null(SIG_list)==TRUE){
      message('Empty SIG_list, please check and re-try !');return(FALSE)
    }
    if(is.null(SJAR.project_name)==TRUE){
      message('Empty SJAR.project_name, please check and re-try !');return(FALSE)
    }
    SJAR.outdir <- file.path(SJAR.main_dir, SJAR.project_name)
    if (!file.exists(SJAR.outdir)) {
      dir.create(SJAR.outdir, recursive = TRUE)
    }
    SJAR.outdir.tf  <-
      file.path(SJAR.main_dir, SJAR.project_name, 'output_tf_')
    SJAR.outdir.sig <-
      file.path(SJAR.main_dir, SJAR.project_name, 'output_sig_')
    SJAR.expression_matrix <-
      file.path(SJAR.main_dir, SJAR.project_name, 'input.exp')
    SJAR.hub_genes.tf <-
      file.path(SJAR.main_dir, SJAR.project_name, 'tf.txt')
    SJAR.hub_genes.sig <-
      file.path(SJAR.main_dir, SJAR.project_name, 'sig.txt')
    SJAR.bash_file.tf <-
      file.path(SJAR.main_dir, SJAR.project_name, 'run_tf.sh')
    SJAR.bash_file.sig <-
      file.path(SJAR.main_dir, SJAR.project_name, 'run_sig.sh')
    ## process exp matrix
    d <- exprs(eset)[, use.samples]
    # filter genes with count=0
    d <- d[!apply(d == 0, 1, all), ]
    # filter genes with IQR
    choose1 <- IQR.filter(d, rownames(d),thre = 0.5,loose_thre=0.1,loose_gene=unique(c(TF_list,SIG_list)))
    d <- d[choose1, ]
    use.genes <- rownames(d)
    # write exp data to exp format
    expdata <- data.frame(cbind(isoformId = fData(eset)[use.genes, ], geneSymbol = fData(eset)[use.genes, ], d))
    #
    write.table(
      expdata,
      file = SJAR.expression_matrix,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    ##
    cat(intersect(use.genes, TF_list),file = SJAR.hub_genes.tf,sep = '\n')
    cat(intersect(use.genes, SIG_list),file = SJAR.hub_genes.sig,sep = '\n')
    # write scripts to bash file for tf
    cmd_tf1 <-
      sprintf(
        '%s %s/generate_pipeline.py %s %s %s %s --run False --host CLUSTER',
        RP.python3.path,
        RP.SJAR.main,
        SJAR.project_name,
        SJAR.expression_matrix,
        SJAR.hub_genes.tf,
        SJAR.outdir.tf
      )
    SJAR.outdir.tf.script <-
      sprintf('%ssjaracne_%s_scripts_/',
              SJAR.outdir.tf,
              SJAR.project_name)
    cmd_tf2 <-
      sprintf(
        'cat %s/02_bootstrap_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.tf.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name
      )
    cmd_tf3 <-
      sprintf(
        'cat %s/03_getconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.tf.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name
      )
    cmd_tf4 <-
      sprintf(
        'cat %s/04_getenhancedconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.tf.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name,
        SJAR.outdir.tf,
        SJAR.project_name
      )
    all_cmd <-
      paste(RP.SJAR.pre_cmd,
            cmd_tf1,
            cmd_tf2,
            cmd_tf3,
            cmd_tf4,
            sep = '\n')
    all_cmd <- get_server_path(all_cmd)
    cat(all_cmd, file = SJAR.bash_file.tf, sep = '\n')
    # write scripts to bash file for sig
    cmd_sig1 <-
      sprintf(
        '%s %s/generate_pipeline.py %s %s %s %s --run False --host CLUSTER',
        RP.python3.path,
        RP.SJAR.main,
        SJAR.project_name,
        SJAR.expression_matrix,
        SJAR.hub_genes.sig,
        SJAR.outdir.sig
      )
    SJAR.outdir.sig.script <-
      sprintf('%ssjaracne_%s_scripts_/',
              SJAR.outdir.sig,
              SJAR.project_name)
    cmd_sig2 <-
      sprintf(
        'cat %s/02_bootstrap_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.sig.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name
      )
    cmd_sig3 <-
      sprintf(
        'cat %s/03_getconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.sig.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name
      )
    cmd_sig4 <-
      sprintf(
        'cat %s/04_getenhancedconsensusnetwork_%s.sh | while read line; do bsub -R "rusage[mem=%d]" -P %s -oo %s%s.log -eo %s%s.err $line;done',
        SJAR.outdir.sig.script,
        SJAR.project_name,
        mem,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name,
        SJAR.outdir.sig,
        SJAR.project_name
      )
    all_cmd <-
      paste(RP.SJAR.pre_cmd,
            cmd_sig1,
            cmd_sig2,
            cmd_sig3,
            cmd_sig4,
            sep = '\n')
    all_cmd <- get_server_path(all_cmd)
    cat(all_cmd, file = SJAR.bash_file.sig, sep = '\n')
    ##
    message(sprintf('Check %s',SJAR.bash_file.tf))
    message(sprintf('Check %s',SJAR.bash_file.sig))
    return(
      list(
        outdir.tf = SJAR.outdir.tf,
        outdir.sig = SJAR.outdir.sig,
        bash.tf = get_server_path(SJAR.bash_file.tf),
        bash.sig = get_server_path(SJAR.bash_file.sig)
      )
    )
}
## auto generate bash files for Step1, Step2, Step3, Step4 for all bash files under one directory
SJ.SJAracne.step <- function(out.dir.SJAR) {
  tf_bash <- list.files(out.dir.SJAR, pattern = 'tf.sh', recursive = TRUE)
  sig_bash <-
    list.files(out.dir.SJAR, pattern = 'sig.sh', recursive = TRUE)
  s1 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 4 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  s2 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 5 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  s3 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 6 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  s4 <-
    lapply(c(tf_bash, sig_bash), function(x) {
      system(paste0('head -n 7 ', out.dir.SJAR, x, '| tail -n 1'), intern = TRUE)
    })
  cat(
    c(RP.SJAR.pre_cmd, unlist(s1)),
    file = sprintf('%s/step1.bash', out.dir.SJAR),
    sep = '\n'
  )
  cat(
    c(RP.SJAR.pre_cmd, unlist(s2)),
    file = sprintf('%s/step2.bash', out.dir.SJAR),
    sep = '\n'
  )
  cat(
    c(RP.SJAR.pre_cmd, unlist(s3)),
    file = sprintf('%s/step3.bash', out.dir.SJAR),
    sep = '\n'
  )
  cat(
    c(RP.SJAR.pre_cmd, unlist(s4)),
    file = sprintf('%s/step4.bash', out.dir.SJAR),
    sep = '\n'
  )
  message(get_server_path(sprintf('Check: %s', sprintf('%s/step1-4.bash', out.dir.SJAR))))
  return(TRUE)
}
####
