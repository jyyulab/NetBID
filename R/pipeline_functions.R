#' @import Biobase limma tximport igraph biomaRt openxlsx msigdbr ConsensusClusterPlus kableExtra
#' @importFrom GEOquery getGEO
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plot3D scatter3D
#' @importFrom plotrix draw.ellipse draw.circle
#' @importFrom impute impute.knn
#' @importFrom umap umap umap.defaults
#' @importFrom rhdf5 H5Fopen H5Fclose
#' @importFrom DESeq2 DESeqDataSetFromTximport DESeq
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom graphics plot
#' @importFrom aricode clustComp
#' @importFrom GSVA gsva
#' @importFrom MCMCglmm MCMCglmm
#' @importFrom arm bayesglm
#' @importFrom reshape melt
#' @importFrom ordinal clm clmm
#' @importFrom rmarkdown render pandoc_available html_document
#' @importFrom Matrix rowSums
#' @importFrom SummarizedExperiment assay
#' @importFrom lme4 lmer
#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf rgb
#' @importFrom graphics abline arrows axis barplot boxplot hist image layout legend lines mtext par points polygon rect segments strheight stripchart strwidth text
#' @importFrom stats IQR aggregate as.dendrogram as.dist cutree density dist fisher.test gaussian glm hclust kmeans ks.test lm median model.matrix na.omit order.dendrogram p.adjust pchisq pnorm prcomp pt qnorm quantile sd splinefun
#' @importFrom utils read.delim write.table

################################
library(Biobase) ## basic functions for bioconductor
library(GEOquery) ## for samples from GEO
library(limma) ## for data normalization for micro-array
library(DESeq2) ## for data normalization for RNASeq
library(tximport) ## for data import from Salmon/sailfish/kallisto/rsem/stringtie output

library(RColorBrewer) ## for color scale
library(colorspace) ## for color scale
library(plot3D) ## for 3D plot
library(igraph) ## for network related functions
library(plotrix) ## for draw.ellipse

library(biomaRt) ## for gene id conversion
library(openxlsx) ## for output into excel
library(impute) ## for impute

library(msigdbr) ## for msigDB gene sets
library(ComplexHeatmap) ## for complex heatmap
library(umap) ## for umap visualization
library(rhdf5) ## for read in MICA results
library(GSVA)

library(MCMCglmm)
library(arm)
library(reshape)
library(ordinal)
library(rmarkdown)

library(aricode)
##
check_para <- function(para_name,envir){
  if(base::exists(para_name,envir=envir)==FALSE){message(sprintf('%s missing !',para_name));return(0)}
  if(is.null(base::get(para_name,envir=envir))==TRUE){message(sprintf('%s is NULL !',para_name));return(0)}
  return(1)
}
check_option <- function(para_name,option_list,envir){
  if(!base::get(para_name,envir=envir) %in% option_list){
    message(sprintf('Only accept %s set at: %s !',para_name,base::paste(option_list,collapse=';')));return(0)
  }
  return(1)
}
clean_charVector <- function(x){
  x1 <- names(x)
  x <- as.character(x);
  x[which(x=='')] <- 'NULL';
  x[which(is.null(x)==TRUE)] <- 'NULL'
  x[which(is.na(x)==TRUE)] <- 'NA'
  names(x) <- x1
  x
}
##
#
#' Preload database files into R workspace for NetBID2
#'
#' \code{db.preload} is a pre-processing function for NetBID2. It preloads needed data into R workspace,
#' and saves it locally under db/ directory with specified species name and analysis level.
#'
#' Users need to set the species name (e.g. human, mouse) and
#' analysis level (transcript or gene level). TF list and SIG list are optional, if not specified, list from package data will be used as default.
#'
#' @param use_level character, users can choose "transcript" or "gene". Default is "gene".
#' @param use_spe character, the name of an interested species (e.g. "human", "mouse", "rat"). Default is "human".
#' @param update logical, if TRUE, previous loaded RData will be updated. Default is FALSE.
#' @param TF_list a character vector, the list of TF (Transcription Factor) names. If NULL, the pre-defined list in the package will be used.
#' Default is NULL.
#' @param SIG_list a character vector, the list of SIG (Signaling Factor) names. If NULL, the pre-defined list in the package will be use.
#' Default is NULL.
#' @param input_attr_type character, the type of the TF_list and SIG_list.
#' Details please check biomaRt, \url{https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html}.
#' If TF_list and SIG_list are not specified, the list in the NetBID2 package will be used.
#' This only support "external_gene_name" and "ensembl_gene_id".
#' Default is "external_gene_name".
#' @param main.dir character, the main directory for NetBID2.
#' If NULL, will be \code{system.file(package = "NetBID2")}. Default is NULL.
#' @param db.dir character, a path for saving the RData.
#' Default is \code{db} directory under the \code{main.dir}, if \code{main.dir} is provided.
#' @param useCache Boolean, parameter pass to \code{getBM()} indicating whether the results cache should be used. Setting to FALSE will disable reading and writing of the cache.
#'
#' @return Return TRUE if loading is successful, otherwise return FALSE. Two variables will be loaded into R workspace, \code{tf_sigs} and \code{db_info}.
#' @examples
#' db.preload(use_level='gene',use_spe='human',update=FALSE)
#'
#' \dontrun{
#' db.preload(use_level='transcript',use_spe='human',update=FALSE)
#' db.preload(use_level='gene',use_spe='mouse',update=FALSE)
#' }
#' @export
db.preload <- function(use_level='transcript',use_spe='human',update = FALSE,
                       TF_list=NULL,SIG_list=NULL,input_attr_type='external_gene_name',
                       main.dir=NULL,
                       db.dir=sprintf("%s/db/",main.dir),useCache = TRUE){
  all_input_para <- c('use_level','use_spe','update')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('use_level',c('transcript','gene'),envir=environment()),
                 check_option('update',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  ## load annotation info, including: TF/Sig list, gene info
  if(is.null(main.dir)==TRUE){
    main.dir <- system.file(package = "NetBID2")
    message(sprintf('main.dir not set, will use package directory: %s',main.dir))
  }
  if(is.null(db.dir)==TRUE){
    db.dir <- sprintf("%s/db/",main.dir)
  }
  message(sprintf('Will use directory %s as the db.dir',db.dir))
  message(sprintf('Your setting for species is %s, with level at %s',use_spe,use_level))
  use_spe <- toupper(use_spe)
  output.db.dir <- sprintf('%s/%s',db.dir,use_spe)
  RData.file <- sprintf('%s/%s_%s.RData', output.db.dir,use_spe,use_level)
  if (!file.exists(RData.file) | update==TRUE) { ## not exist or need to update
    ## get info from use_spe
    ensembl <- biomaRt::useMart("ensembl")
    all_ds  <- biomaRt::listDatasets(ensembl)
    w1 <- grep(sprintf("^%s GENES",use_spe),toupper(all_ds$description))
    if(base::length(w1)==0){
      tmp_use_spe <- unlist(strsplit(use_spe,' ')); tmp_use_spe <- tmp_use_spe[base::length(tmp_use_spe)]
      w1 <- grep(sprintf(".*%s_GENE_ENSEMBL",toupper(tmp_use_spe)),toupper(all_ds$dataset))
      if(base::length(w1)==1){use_spe <- toupper(strsplit(all_ds[w1,2],' ')[[1]][1]); output.db.dir <- sprintf('%s/%s',db.dir,use_spe)}
    }
    if(base::length(w1)==1){
      w2 <- all_ds[w1,1]
      mart <- biomaRt::useMart(biomart="ensembl", dataset=w2) ## get id for input spe
      message(sprintf('Read in ensembl annotation file for %s and output all db files in %s/%s !',use_spe,db.dir,use_spe))
    }
    if(base::length(w1)==0){
      message(sprintf('Check input use_spe parameter: %s, not included in the ensembl database',use_spe))
      return(FALSE)
    }
    if(base::length(w1)>1){
      w2 <- base::paste(all_ds[w1,2],collapse=';')
      message(sprintf('Check input use_spe parameter: %s, more than one species match in ensembl database : %s,
                      please check and re-try',use_spe,w2))
      return(FALSE)
    }
  }
  RData.file <- sprintf('%s/%s_%s.RData', output.db.dir,use_spe,use_level)
  if (update == TRUE | !file.exists(RData.file)) {
    if(!file.exists(output.db.dir)){
      dir.create(output.db.dir)
    }
    ## get attributes for mart
    filters <- biomaRt::listFilters(mart)
    attributes <- biomaRt::listAttributes(mart)
    ensembl.attr.transcript <- c('ensembl_transcript_id','ensembl_gene_id',
                                 'external_transcript_name','external_gene_name',
                                 'gene_biotype','gene_biotype',
                                 'chromosome_name','strand','start_position','end_position','band','transcript_start','transcript_end',
                                 'description','phenotype_description','refseq_mrna')
    ensembl.attr.gene <- c('ensembl_gene_id','external_gene_name',
                           'gene_biotype',
                           'chromosome_name','strand','start_position','end_position','band',
                           'description','phenotype_description','refseq_mrna')
    if(use_spe=='HUMAN'){
      ensembl.attr.transcript <- c(ensembl.attr.transcript,'hgnc_symbol','entrezgene_id')
      ensembl.attr.gene <- c(ensembl.attr.gene,'hgnc_symbol','entrezgene_id')
    }
    ## do not output hgnc in non-human species
    if(use_spe != 'HUMAN')
      ensembl.attr.transcript <- base::setdiff(ensembl.attr.transcript,'hgnc_symbol')
    if(use_spe != 'HUMAN')
      ensembl.attr.gene <- base::setdiff(ensembl.attr.gene,'hgnc_symbol')
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
      if(use_spe != 'HUMAN' & use_spe != 'MOUSE'){ ## if spe not human/mouse
        TF_f <- sprintf('%s/%s_TF_%s.txt',db.dir,'HUMAN',filter_attr)
        message(sprintf('Will use %s file as the input TF_list!',TF_f))
        TF_list <- read.delim(TF_f,stringsAsFactors=FALSE,header=F)$V1
        filter_attr <- 'external_gene_name'
        tmp1 <- biomaRt::getBM(attributes=c('hsapiens_homolog_associated_gene_name','external_gene_name'),values=TRUE,mart=mart,filters='with_hsapiens_homolog',useCache = useCache)
        TF_list <- base::unique(tmp1[which(tmp1[,1] %in% TF_list),2])
      }else{
        TF_f <- sprintf('%s/%s_TF_%s.txt',db.dir,use_spe,filter_attr)
        message(sprintf('Will use %s file as the input TF_list!',TF_f))
        TF_list <- read.delim(TF_f,stringsAsFactors=FALSE,header=F)$V1
      }
    }
    if(use_level=='transcript'){
      message(sprintf('Begin read TF list information from ensembl for %s !',use_spe))
      TF_info  <- biomaRt::getBM(attributes = ensembl.attr.transcript,values=TF_list, mart=mart, filters=filter_attr,useCache = useCache)
    }
    if(use_level=='gene'){
      message(sprintf('Begin read TF list information from ensembl for %s !',use_spe))
      TF_info  <- biomaRt::getBM(attributes = ensembl.attr.gene,values=TF_list, mart=mart, filters=filter_attr,useCache = useCache)
    }
    # for SIG list
    if(is.null(SIG_list)){
      if(use_spe != 'HUMAN' & use_spe != 'MOUSE'){ ## if spe not human/mouse
        SIG_f <- sprintf('%s/%s_SIG_%s.txt',db.dir,'HUMAN',filter_attr)
        message(sprintf('Will use %s file as the input SIG_list!',SIG_f))
        SIG_list <- read.delim(SIG_f,stringsAsFactors=FALSE,header=F)$V1
        filter_attr <- 'external_gene_name'
        tmp1 <- biomaRt::getBM(attributes=c('hsapiens_homolog_associated_gene_name','external_gene_name'),values=TRUE,mart=mart,filters='with_hsapiens_homolog',useCache = useCache)
        SIG_list <- base::unique(tmp1[which(tmp1[,1] %in% SIG_list),2])
      }else{
        SIG_f <- sprintf('%s/%s_SIG_%s.txt',db.dir,use_spe,filter_attr)
        message(sprintf('Will use %s file as the input SIG_list!',SIG_f))
        SIG_list <- read.delim(SIG_f,stringsAsFactors=FALSE,header=F)$V1
      }
    }
    if(use_level=='transcript'){
      message(sprintf('Begin read SIG list information from ensembl for %s !',use_spe))
      SIG_info  <- biomaRt::getBM(attributes = ensembl.attr.transcript,values=SIG_list, mart=mart, filters=filter_attr,useCache = useCache)
    }
    if(use_level=='gene'){
      message(sprintf('Begin read SIG list information from ensembl for %s !',use_spe))
      SIG_info  <- biomaRt::getBM(attributes = ensembl.attr.gene,values=SIG_list, mart=mart, filters=filter_attr,useCache = useCache)
    }
    # check input not in the list
    miss_TF  <- base::unique(base::setdiff(TF_list,TF_info[[filter_attr]]))
    miss_SIG <- base::unique(base::setdiff(SIG_list,SIG_info[[filter_attr]]))
    if(base::length(miss_TF)>0){message(sprintf("%d TFs could not match,please check and choose to re-try : %s",
                                          base::length(miss_TF),base::paste(sort(miss_TF),collapse=';')))}
    if(base::length(miss_SIG)>0){message(sprintf("%d SIGs could not match,please check and choose to re-try : %s",
                                           base::length(miss_SIG),base::paste(sort(miss_SIG),collapse=';')))}
    ####### output full info
    tf_sigs <- list();tf_sigs$tf <- list();tf_sigs$sig <- list();
    tf_sigs$tf$info  <- TF_info; tf_sigs$sig$info <- SIG_info;
    for(each_id_type in base::intersect(c('ensembl_transcript_id','ensembl_gene_id',
                                    'external_transcript_name','external_gene_name','hgnc_symbol',
                                    'entrezgene_id','refseq_mrna'),colnames(TF_info))){
      tf_sigs$tf[[each_id_type]] <- base::setdiff(base::unique(TF_info[[each_id_type]]),"")
      tf_sigs$sig[[each_id_type]] <- base::setdiff(base::unique(SIG_info[[each_id_type]]),"")
    }
    db_info <- all_ds[w1,]
    save(tf_sigs,db_info=db_info,file = RData.file)
  }
  load(RData.file,.GlobalEnv)
  return(TRUE)
  }

#' Get Transcription Factor (TF) and Signaling Factor (SIG) List
#'
#' \code{get.TF_SIG.list} is a function converts gene ID into the corresponding TF/SIG list,
#' with selected gene/transcript type.
#'
#' @param use_genes a vector of characters, genes will be used in the network construction.
#' If NULL, no filter will be performed to the TF/SIG list. Default is NULL.
#' @param use_gene_type character, the attribute name inherited from the biomaRt package.
#' Some options are, "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" and "refseq_mrna".
#' All options can be accessed by calling \code{biomaRt::useMart} (e.g. mart <- biomaRt::useMart('ensembl',db_info[1]); biomaRt::listAttributes(mart)$name).
#'
#' The type must match the gene type from the input \code{use_genes}. Default is "external_gene_name".
#' @param ignore_version logical, if TRUE, the version "ensembl_gene_id_version" or "ensembl_transcript_id_version" will be ignored.
#' Default is FALSE.
#' @param dataset character, the dataset used for ID conversion (e.g. "hsapiens_gene_ensembl").
#' If NULL, use \code{db_info[1]} from \code{db.preload}. Default is NULL.
#'
#'
#' @return Return a list containing two elements. \code{tf} is the TF list, \code{sig} is the SIG list.
#'
#' @examples
#' db.preload(use_level='transcript',use_spe='human',update=FALSE)
#' use_genes <- c("ENST00000210187","ENST00000216083","ENST00000216127",
#'              "ENST00000216416","ENST00000217233","ENST00000221418",
#'              "ENST00000504956","ENST00000507468")
#' res_list  <- get.TF_SIG.list(use_gene_type = 'ensembl_transcript_id',
#'                                use_genes=use_genes,
#'                                dataset='hsapiens_gene_ensembl')
#' print(res_list)
#'
#' \dontrun{
#' }
#'
#' @export
get.TF_SIG.list <- function(use_genes=NULL,
                            use_gene_type='external_gene_name',ignore_version=FALSE,
                            dataset=NULL){
  #
  all_input_para <- c('use_genes')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('ignore_version',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(dataset)==TRUE){
    check_res <- check_para('db_info',envir=environment())
    if(base::min(check_res)==0){message('Please use db.preload() to get db_info or set dataset, check and re-try!');return(FALSE)}
    dataset <- db_info[1]
  }
  if(is.null(tf_sigs)==TRUE){
    message('tf_sigs not loaded yet, please run db.preload() before processing !');return(FALSE);
  }
  n1 <- names(tf_sigs$tf)[-1]
  if(use_gene_type %in% n1){
    if(is.null(use_genes)==TRUE){
      TF_list <- base::unique(tf_sigs$tf[[use_gene_type]])
      SIG_list <- base::unique(tf_sigs$sig[[use_gene_type]])
    }else{
      TF_list <- base::unique(base::intersect(use_genes,tf_sigs$tf[[use_gene_type]]))
      SIG_list <- base::unique(base::intersect(use_genes,tf_sigs$sig[[use_gene_type]]))
    }
  }else{
    if(grepl('version$',use_gene_type)==TRUE & ignore_version==TRUE){
      use_genes_no_v <- gsub('(.*)\\..*','\\1',use_genes)
      transfer_tab <- data.frame(to_type=use_genes,from_type=use_genes_no_v,stringsAsFactors = FALSE)
      print(str(transfer_tab))
      TF_list <- transfer_tab[which(transfer_tab$from_type %in% tf_sigs$tf[[gsub('(.*)_version','\\1',use_gene_type)]]),'to_type']
      SIG_list <- transfer_tab[which(transfer_tab$from_type %in% tf_sigs$sig[[gsub('(.*)_version','\\1',use_gene_type)]]),'to_type']
    }else{
      mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
      #filters <- biomaRt::listFilters(mart)
      attributes <- biomaRt::listAttributes(mart)
      if(!use_gene_type %in% attributes$name){
        message(sprintf('%s not in the attributes for %s, please check and re-try !',use_gene_type,dataset));return(FALSE)
      }
      transfer_tab <- get_IDtransfer(from_type=use_gene_type,to_type=n1[1],ignore_version = ignore_version)
      print(str(transfer_tab))
      TF_list <- get_name_transfertab(use_genes=tf_sigs$tf[n1[1]][[1]],transfer_tab = transfer_tab,ignore_version = ignore_version,ignore_order=TRUE,from_type=n1[1],to_type=use_gene_type)
      SIG_list <- get_name_transfertab(use_genes=tf_sigs$sig[n1[1]][[1]],transfer_tab = transfer_tab,ignore_version = ignore_version,ignore_order=TRUE,from_type=n1[1],to_type=use_gene_type)
    }
    TF_list <- base::unique(base::intersect(use_genes,TF_list))
    SIG_list <- base::unique(base::intersect(use_genes,SIG_list))
  }
  message(sprintf('%d TFs and %s SIGs are included in the expression matrix !',base::length(TF_list),base::length(SIG_list)))
  return(list(tf=TF_list,sig=SIG_list))
}

#' Creates Data Frame for ID Conversion
#'
#' \code{get_IDtransfer} creates a data frame for ID conversion using biomaRt. For example, to convert Ensembl ID into gene symbol.
#'
#' @param from_type character, the attribute name match the current ID type (the type of \code{use_genes}).
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' The "attribute" is inherited from the biomaRt package. For details, user can call \code{biomaRt::listAttributes()} to display all available attributes in the selected dataset.
#' @param to_type character, the attribute name to convert into.
#' @param add_type character, the additional attribute name to add into the conversion data frame.
#' @param use_genes a vector of characters, the genes for ID conversion.
#' If NULL, all genes will be selected.
#' @param dataset character, name of the dataset used for ID conversion. For example, "hsapiens_gene_ensembl".
#' If NULL, \code{db_info[1]} will be used. \code{db_info} requires the calling of \code{db.preload} in the previous steps.
#' Default is NULL.
#' @param ignore_version logical, if it is set to TRUE and \code{from_type} is "ensembl_gene_id_version" or "ensembl_transcript_id_version",
#' the version of the original ID will be ignored in ID mapping.
#' Default is FALSE.
#' @param useCache Boolean, parameter pass to \code{getBM()} indicating whether the results cache should be used. Setting to FALSE will disable reading and writing of the cache.
#'
#' @return
#' Return a data frame for ID conversion.
#'
#' @examples
#' use_genes <- c("ENST00000210187","ENST00000216083","ENST00000216127",
#'              "ENST00000216416","ENST00000217233","ENST00000221418")
#' transfer_tab <- get_IDtransfer(from_type = 'ensembl_transcript_id',
#'                                to_type='external_gene_name',
#'                                use_genes=use_genes,
#'                                dataset='hsapiens_gene_ensembl')
#'                                ## get transfer table !!!
#' res1 <- get_name_transfertab(use_genes,transfer_tab=transfer_tab)
#' transfer_tab_withtype <- get_IDtransfer2symbol2type(from_type = 'ensembl_transcript_id',
#'                                                    use_genes=use_genes,
#'                                                    dataset='hsapiens_gene_ensembl')
#'                                                    ## get transfer table !!!
#' \dontrun{
#' }
#' @export
get_IDtransfer <- function(from_type=NULL,to_type=NULL,add_type=NULL,use_genes=NULL,dataset=NULL,ignore_version=FALSE,useCache = TRUE){
  #
  all_input_para <- c('from_type','to_type')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('ignore_version',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(dataset)==TRUE){
    check_res <- check_para('db_info',envir=environment())
    if(base::min(check_res)==0){message('Please use db.preload() to get db_info or set dataset, check and re-try!');return(FALSE)}
    dataset <- db_info[1]
  }
  mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
  attributes <- biomaRt::listAttributes(mart)

  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,dataset));return(FALSE)
  }
  if(!to_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',to_type,dataset));return(FALSE)
  }
  ori_from_type <- from_type
  ori_use_genes <- use_genes
  if(from_type %in% c('ensembl_gene_id_version','ensembl_transcript_id_version')){
    if(ignore_version==FALSE){
      message(sprintf('Attention: %s in %s will be updated with new version number, please check the output.
                      If lots of missing, try to set ignore_version=TRUE and try again !',from_type,dataset));
    }else{
      from_type <- gsub('(.*)_version','\\1',from_type)
      if(is.null(use_genes)==FALSE) use_genes <- gsub('(.*)\\..*','\\1',use_genes)
    }
  }
  if(is.null(use_genes)==TRUE | base::length(use_genes)>100){
    tmp1 <- biomaRt::getBM(attributes=c(from_type,to_type,add_type),values=1,mart=mart,filters='strand',useCache = useCache)
    tmp2 <- biomaRt::getBM(attributes=c(from_type,to_type,add_type),values=-1,mart=mart,filters='strand',useCache = useCache)
    tmp1 <- base::rbind(tmp1,tmp2)
    if(is.null(use_genes)==FALSE){
      tmp1 <- tmp1[which(tmp1[,1] %in% use_genes),]
    }
  }else{
    tmp1 <- biomaRt::getBM(attributes=c(from_type,to_type,add_type),values=use_genes,mart=mart,filters=from_type,useCache = useCache)
  }
  if(ori_from_type %in% c('ensembl_gene_id_version','ensembl_transcript_id_version') & is.null(use_genes)==FALSE & ignore_version==TRUE){
    tmp2 <- data.frame(ori_from_type=ori_use_genes,from_type=use_genes,stringsAsFactors=FALSE)
    names(tmp2) <- c(ori_from_type,from_type)
    tmp1 <- base::merge(tmp2,tmp1,by.y=from_type,by.x=from_type)[c(from_type,to_type,add_type,ori_from_type)]
  }
  w1 <- apply(tmp1,1,function(x)base::length(which(is.na(x)==TRUE | x=="")))
  transfer_tab <- tmp1[which(w1==0),]
  for(i in 1:ncol(transfer_tab)){
    transfer_tab[,i] <- as.character(transfer_tab[,i])
  }
  return(transfer_tab)
}

#' Create Data Frame for ID Conversion Between Species
#'
#' \code{get_IDtransfer_betweenSpecies} creates a data frame to convert ID between species.
#'
#' @param from_spe character, name of the original species (e.g. "human", "mouse", "rat") that \code{use_genes} belongs to. Default is "human".
#' @param to_spe character, name of the target species (e.g. "human", "mouse", "rat"). Default is "mouse".
#' @param from_type character, the attribute name match the current ID type (the type of \code{use_genes}).
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' The "attribute" is inherited from the \code{biomaRt} package. For details, user can call \code{biomaRt::listAttributes()} function to display all available attributes in the selected dataset.
#' @param to_type character, the attribute name match the target ID type.
#' @param use_genes a vector of characters, the genes for ID conversion. Must be the genes with ID type of \code{from_type}, and from species \code{from_spe}.
#' If NULL, all the possible genes will be shown in the conversion table. Default is NULL.
#' @param useCache Boolean, parameter pass to \code{getBM()} indicating whether the results cache should be used. Setting to FALSE will disable reading and writing of the cache.
#'
#' @return Return a data frame for ID conversion, from one species to another.
#'
#' @examples
#' use_genes <- c("ENST00000210187","ENST00000216083","ENST00000216127",
#'              "ENST00000216416","ENST00000217233","ENST00000221418")
#' transfer_tab <- get_IDtransfer_betweenSpecies(from_spe='human',
#'                                to_spe='mouse',
#'                                from_type = 'ensembl_transcript_id',
#'                                to_type='external_gene_name',
#'                                use_genes=use_genes)
#'                                ## get transfer table !!!
#' transfer_tab <- get_IDtransfer_betweenSpecies(from_spe='human',
#'                                to_spe='mouse',
#'                                from_type = 'ensembl_transcript_id',
#'                                to_type='ensembl_transcript_id_version',
#'                                use_genes=use_genes)
#'                                ## get transfer table !!!
#' \dontrun{
#' transfer_tab <- get_IDtransfer_betweenSpecies(from_spe='human',
#'                                to_spe='mouse',
#'                                from_type='refseq_mrna',
#'                                to_type='refseq_mrna')
#' }
#' @export
get_IDtransfer_betweenSpecies <- function(from_spe='human',to_spe='mouse',
                                          from_type=NULL,to_type=NULL,
                                          use_genes=NULL,useCache = TRUE){
  #
  all_input_para <- c('from_spe','to_spe','from_type','to_type')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  from_spe <- toupper(from_spe)
  to_spe <- toupper(to_spe)
  ensembl <- biomaRt::useMart("ensembl")
  all_ds  <- biomaRt::listDatasets(ensembl)
  w1 <- grep(sprintf("^%s GENES",from_spe),toupper(all_ds$description))
  if(base::length(w1)==1){
    from_spe_ds <- all_ds[w1,1]
    mart1 <- biomaRt::useMart(biomart="ensembl", dataset=from_spe_ds) ## get id for input spe
  }
  if(base::length(w1)==0){
    message(sprintf('Check input from_spe parameter: %s, not included in the ensembl database',from_spe))
    return(FALSE)
  }
  if(base::length(w1)>1){
    w2 <- base::paste(all_ds[w1,2],collapse=';')
    message(sprintf('Check input from_spe parameter: %s, more than one species match in ensembl database : %s,
                    please check and re-try',from_spe,w2))
    return(FALSE)
  }
  w1 <- grep(sprintf("^%s GENES",to_spe),toupper(all_ds$description))
  if(base::length(w1)==1){
    to_spe_ds <- all_ds[w1,1]
    mart2 <- biomaRt::useMart(biomart="ensembl", dataset=to_spe_ds) ## get id for input spe
  }
  if(base::length(w1)==0){
    message(sprintf('Check input to_spe parameter: %s, not included in the ensembl database',to_spe))
    return(FALSE)
  }
  if(base::length(w1)>1){
    w2 <- base::paste(all_ds[w1,2],collapse=';')
    message(sprintf('Check input to_spe parameter: %s, more than one species match in ensembl database : %s,
                    please check and re-try',to_spe,w2))
    return(FALSE)
  }
  #### mart1 mart2
  attributes <- biomaRt::listAttributes(mart1)
  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,from_spe));return(FALSE)
  }
  attributes <- biomaRt::listAttributes(mart2)
  if(!to_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',to_type,to_spe));return(FALSE)
  }
  ## get homolog between from_spe to to_spe
  cn1 <- gsub('(.*)_gene_ensembl','\\1',from_spe_ds)
  cn2 <- attributes$name ## attribute names in mart2
  cn3 <- cn2[grep(sprintf('%s_homolog_associated_gene_name',cn1),cn2)]
  if(base::length(cn3)!=1){
    message('No homolog info found in Biomart, sorry !');return(FALSE)
  }
  tmp1 <- get_IDtransfer(from_type=from_type,to_type='external_gene_name',use_genes=use_genes,dataset=from_spe_ds)
  tmp2 <- biomaRt::getBM(attributes=c(cn3,'external_gene_name'),values=TRUE,
                mart=mart2,filters=sprintf('with_%s_homolog',cn1),useCache = useCache)
  colnames(tmp1) <- sprintf('%s_%s',colnames(tmp1),from_spe)
  colnames(tmp2) <- sprintf('%s_%s',colnames(tmp2),to_spe)
  tmp3 <- base::merge(tmp1,tmp2,by.x=sprintf('external_gene_name_%s',from_spe),by.y=sprintf('%s_%s',cn3,to_spe))
  transfer_tab <- tmp3[,c(2,3,1)]
  if(to_type != 'external_gene_name'){
    tmp4 <- get_IDtransfer(from_type='external_gene_name',to_type=to_type,use_genes=tmp3[,3],dataset=to_spe_ds)
    colnames(tmp4) <- sprintf('%s_%s',colnames(tmp4),to_spe)
    tmp5 <- base::merge(tmp3,tmp4)
    transfer_tab <- tmp5[,c(3,4,2,1)]
  }
  return(transfer_tab)
}


#' Create Data Frame for ID Conversion With Biotype Information
#'
#' \code{get_IDtransfer2symbol2type} creates a data frame to convert original ID into gene symbol and gene biotype (gene level),
#' or into transcript symbol and transcript biotype (transcript level).
#'
#' @param from_type character, the attribute name matches the current ID type (the type of use_genes).
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' The "attribute" is inherited from the biomaRt package.
#' For details, user can call \code{biomaRt::listAttributes()} function to display all available attributes in the selected dataset.
#' @param use_genes a vector of characters, the genes for ID conversion.
#' If NULL, all genes will be selected.
#' @param dataset character, name of the dataset used for ID conversion.
#' For example, "hsapiens_gene_ensembl".
#' If NULL, \code{db_info[1]} will be used. \code{db_info} requires the calling of \code{db.preload} in the previous steps.
#' Default is NULL.
#' @param use_level character, users can chose between "transcript" and "gene". Default is "gene".
#' @param ignore_version logical, if it is set to TRUE and \code{from_type} is "ensembl_gene_id_version" or "ensembl_transcript_id_version",
#' the version of the original ID will be ignored in ID mapping.
#'
#' @return
#' Return a data frame for ID conversion, from ID to gene symbol and gene biotype.
#'
#' @examples
#' use_genes <- c("ENST00000210187","ENST00000216083","ENST00000216127",
#'              "ENST00000216416","ENST00000217233","ENST00000221418")
#' transfer_tab <- get_IDtransfer(from_type = 'ensembl_transcript_id',
#'                                to_type='external_gene_name',use_genes=use_genes,
#'                                dataset='hsapiens_gene_ensembl')
#'                                ## get transfer table !!!
#' res1 <- get_name_transfertab(use_genes,transfer_tab=transfer_tab)
#' transfer_tab_withtype <- get_IDtransfer2symbol2type(from_type = 'ensembl_transcript_id',
#'                                                    use_genes=use_genes,
#'                                                    dataset='hsapiens_gene_ensembl',
#'                                                    use_level='transcript')
#'                                                    ## get transfer table !!!
#' \dontrun{
#' }
#' @export
get_IDtransfer2symbol2type <- function(from_type=NULL,use_genes=NULL,dataset=NULL,use_level='gene',ignore_version=FALSE){
  #
  all_input_para <- c('from_type')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('ignore_version',c(TRUE,FALSE),envir=environment()),
                 check_option('use_level',c('gene','transcript'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(dataset)==TRUE){
    check_res <- check_para('db_info',envir=environment())
    if(base::min(check_res)==0){message('Please use db.preload() to get db_info or set dataset, check and re-try!');return(FALSE)}
    dataset <- db_info[1]
  }
  message(sprintf('Your setting is at %s level',use_level))
  mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
  attributes <- biomaRt::listAttributes(mart)
  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,dataset));return(FALSE)
  }
  if(use_level=='gene') tmp1 <- get_IDtransfer(from_type=from_type,to_type='external_gene_name',add_type='gene_biotype',use_genes=use_genes,dataset=dataset,ignore_version=ignore_version)
  if(use_level=='transcript')   tmp1 <- get_IDtransfer(from_type=from_type,to_type='external_transcript_name',add_type='transcript_biotype',use_genes=use_genes,dataset=dataset,ignore_version=ignore_version)
  transfer_tab <- tmp1
  return(transfer_tab)
}

#' Convert Original Gene ID into Target Gene ID
#'
#' \code{get_name_transfertab} converts the original gene IDs into target gene IDs, with conversion table provided.
#'
#' @param use_genes a vector of characters, the genes for ID conversion.
#' @param transfer_tab data.frame, the conversion table. Users can create it by calling \code{get_IDtransfer}.
#' @param from_type character, the attribute name match the current ID type (the type of use_genes).
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' The "attribute" is inherited from the biomaRt package. For details, user can call \code{biomaRt::listAttributes()} to see all available attributes in the selected dataset.
#' If NULL, will use the first column of \code{transfer_tab}.
#' @param to_type character, the attribute name to convert into. If NULL, will use the second column of \code{transfer_tab}.
#' @param ignore_version logical, if TRUE and \code{from_type} is "ensembl_gene_id_version" or "ensembl_transcript_id_version", the version will be ignored. Default is FALSE.
#' @param ignore_order logical, whether need to ignore the output order to match the input list of \code{use_genes}. Default is FALSE.
#' @return Return a vector of converted gene IDs.
#'
#' @examples
#' use_genes <- c("ENST00000210187","ENST00000216083","ENST00000216127",
#'              "ENST00000216416","ENST00000217233","ENST00000221418")
#' transfer_tab <- get_IDtransfer(from_type = 'ensembl_transcript_id',
#'                                to_type='external_gene_name',use_genes=use_genes,
#'                                dataset='hsapiens_gene_ensembl')
#'                                ## get transfer table !!!
#' res1 <- get_name_transfertab(use_genes=use_genes,transfer_tab=transfer_tab)
#' transfer_tab_withtype <- get_IDtransfer2symbol2type(from_type = 'ensembl_transcript_id',
#'                                                    use_genes=use_genes,
#'                                                    dataset='hsapiens_gene_ensembl')
#'                                                    ## get transfer table !!!
#' \dontrun{
#' }
#' @export
get_name_transfertab <- function(use_genes=NULL,transfer_tab=NULL,from_type=NULL,to_type=NULL,ignore_version=FALSE,ignore_order=FALSE){
  #
  all_input_para <- c('use_genes','transfer_tab')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('ignore_version',c(TRUE,FALSE),envir=environment()),
                 check_option('ignore_order',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(from_type)==TRUE){from_type=colnames(transfer_tab)[1];}
  if(is.null(to_type)==TRUE){to_type=colnames(transfer_tab)[2];}
  if(ignore_version==TRUE){
    w1 <- which(colnames(transfer_tab)==from_type)
    transfer_tab[,w1] <- gsub('(.*)\\..*','\\1',transfer_tab[,w1])
    from_type <- gsub('(.*)_version','\\1',from_type)
    colnames(transfer_tab)[w1] <- from_type
    use_genes <- gsub('(.*)\\..*','\\1',use_genes)
  }
  transfer_tab <- base::unique(transfer_tab[,c(from_type,to_type)])
  x <- use_genes;
  t1 <- base::unique(transfer_tab[which(transfer_tab[,from_type] %in% x),])
  c1 <- base::unique(t1[,from_type])
  if(base::length(c1)<nrow(t1) & ignore_order==FALSE){
    message('Gene ID in from type contain multiple items!');return(FALSE)
  }
  if(ignore_order==TRUE){
    x1 <- t1[,to_type]
  }else{
    rownames(t1) <- t1[,from_type]
    x1 <- t1[x,to_type]
    w1 <- which(is.na(x1)==TRUE)
    x1[w1] <- x[w1]
  }
  return(x1)
}

#' Manipulation of Working Directories for NetBID2 Network Construction Step
#'
#' \code{NetBID.network.dir.create} is used to help users create an organized working directory for the network construction step in NetBID2 analysis.
#' However, it is not essential for the analysis.
#' It creates a hierarchcial working directory and returns a list contains this directory information.
#'
#' This function needs users to define the main working directory and the project's name.
#' It creates a main working directory with a subdirectory of the project.
#' It also automatically creates three subfolders (QC, DATA and SJAR) within the project folder. QC/,
#' storing Quality Control related plots; DATA/, saving data in RData format;
#' SJAR/, storing files needed for running SJAracne command.
#' This function also returns a list object (example, \code{network.par} in the demo) with directory information wrapped inside.
#' This list is an essential for
#' network construction step, all the important intermediate data generated later will be wrapped inside.

#' @param project_main_dir character, name or absolute path of the main working directory.
#' @param project_name character, name of the project folder.
#'
#' @return \code{NetBID.network.dir.create} returns a list object, containing main.dir (path of the main working directory),
#' project.name (project name), out.dir (path of the project folder, which contains three subfolders), out.dir.QC,
#' out.dir.DATA and out.dir.SJAR.

#' @examples
#'
#' \dontrun{
#' # Creating a main working directory under the current working directory by folder name
#' network.par <- NetBID.network.dir.create("MyMainDir","MyProject")
#' # Or creating a main working directory under the current working directory by relative path
#' network.par <- NetBID.network.dir.create("./MyMainDir","MyProject")
#' # Or creating a main working directory to a specific path by absolute path
#' network.par <- NetBID.network.dir.create("~/Desktop/MyMainDir","MyProject")
#' }
#' @export
NetBID.network.dir.create <- function(project_main_dir=NULL,project_name=NULL){
  #
  if(base::exists('network.par')==TRUE){
    stop('network.par is occupied in the current session,please manually run: rm(network.par) and re-try, otherwise will not change !');
  }
  #
  all_input_para <- c('project_main_dir','project_name')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  network.par <- list()
  network.par$main.dir <- project_main_dir
  network.par$project.name <- project_name
  network.par$out.dir <- sprintf('%s/%s',network.par$main.dir,network.par$project.name)
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

#' Manipulation of Working Directories for NetBID2 Driver Estimation Step
#'
#' \code{NetBID.analysis.dir.create} is used to help users create an organized working directory
#' for the driver estimation step in NetBID2 analysis.
#' However, it is not essential for the analysis.
#' It creates a hierarchcial working directory and returns a list contains this directory information.
#'
#' This function requires user to define the main working directory and the project’s name.
#' It creates a main working directory with a subdirectory of the project.
#' It also automatically creates three subfolders (QC, DATA and PLOT) within the project folder.
#' QC/, storing Quality Control related plots; DATA/, saving data in RData format; PLOT/, storing output plots.
#' This function also returns a list object (e.g. \code{analysis.par} in the demo) with directory information wrapped inside.
#' This list is an essential for driver construction step, all the important intermediate data generated later will be wrapped inside.
#'
#' @param project_main_dir character, name or absolute path of the main working directory for driver analysis.
#' @param project_name character, name of the project folder.
#' @param network_dir character, name or absolute path of the main working directory for network construction.
#' @param network_project_name character, the project name of network construction. Or use the project name of SJARACNe.
#' This parameter is optional. If one didn't run NetBID2 network construction part in the pipeline, he could set it to NULL.
#' If one like to follow the NetBID2 pipeline, he should set it to the path of the TF network file and the SIG network file.
#' @param tf.network.file character, the path of the TF network file (e.g. "XXX/consensus_network_ncol_.txt").
#' Default is the path of network_project_name.
#' @param sig.network.file character, the path of the SIG network file (e.g. "XXX/consensus_network_ncol_.txt").
#' Default is the path of network_project_name.
#'
#' @return Returns a list object, containing main.dir (path of the main working directory), project.name (project name),
#' out.dir (path of the project folder, which contains three subfolders), out.dir.QC, out.dir.DATA and out.dir.PLOT.
#'
#' @examples
#'
#' \dontrun{
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' #
#' project_main_dir <- 'demo1/'
#' project_name <- 'driver_test'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             project_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' }
#' @export
NetBID.analysis.dir.create <- function(project_main_dir=NULL,project_name=NULL,
                                       network_dir=NULL,
                                       network_project_name=NULL,
                                       tf.network.file=NULL,
                                       sig.network.file=NULL){
  #
  if(base::exists('analysis.par')==TRUE){
    stop('analysis.par is occupied in the current session,please manually run: rm(analysis.par) and re-try, otherwise will not change !');
  }
  #
  all_input_para <- c('project_main_dir','project_name')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  if((is.null(network_dir)==TRUE | is.null(network_project_name)==TRUE) & (is.null(tf.network.file)==TRUE | is.null(sig.network.file)==TRUE)){
    message('Either network_dir,network_project_name or tf.network.file,sig.network.file is required, please check and re-try !')
    return(FALSE);
  }
  #
  analysis.par <- list()
  analysis.par$main.dir <- project_main_dir
  analysis.par$project.name <- project_name
  analysis.par$out.dir <- sprintf('%s/%s/',analysis.par$main.dir,analysis.par$project.name)
  analysis.par$tf.network.file <- ''
  analysis.par$sig.network.file <- ''
  if(is.null(tf.network.file)==FALSE){
    analysis.par$tf.network.file  <- tf.network.file
  }else{
    tf_net1 <- sprintf('%s/SJAR/%s/output_tf_sjaracne_%s_out_.final/consensus_network_ncol_.txt',
                       network_dir,network_project_name,network_project_name) ## old version of sjaracne
    tf_net2 <- sprintf('%s/SJAR/SJARACNE_%s_TF/consensus_network_ncol_.txt',
                       network_dir,network_project_name) ## new version of sjaracne
    if(file.exists(tf_net2)) analysis.par$tf.network.file <- tf_net2 else analysis.par$tf.network.file <- tf_net1
  }
  if(is.null(tf.network.file)==FALSE){
    analysis.par$sig.network.file <- sig.network.file
  }else{
    sig_net1 <- sprintf('%s/SJAR/%s/output_sig_sjaracne_%s_out_.final/consensus_network_ncol_.txt',
                        network_dir,network_project_name,network_project_name) ## old version of sjaracne
    sig_net2 <- sprintf('%s/SJAR/SJARACNE_%s_SIG/consensus_network_ncol_.txt',
                        network_dir,network_project_name) ## new version of sjaracne
    if(file.exists(sig_net2)) analysis.par$sig.network.file <- sig_net2  else analysis.par$sig.network.file <- sig_net1
  }
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
  #
  message(sprintf('Analysis space created, please check %s',analysis.par$out.dir))
  return(analysis.par)
}

#' Save Data Produced by Corresponding NetBID2 Pipeline Step.
#'
#' \code{NetBID.saveRData} is a function to save complicated list object generated by certain steps of NetBID2's pipeline
#' (e.g. load gene expression file from GEO, 'exp-load').
#' This function is not essential, but it is highly suggested for easier pipeline step checkout and reference.
#'
#' There are two important steps in the NetBID2 pipeline, network construction and driver analysis.
#' User can save these two complicated list objects, network.par and analysis.par.
#' Assigning the \code{step} name to save the RData for easier reference.
#' Calling \code{NetBID.loadRData} to load the corresponding step RData, users can avoid repeating the former steps.
#'
#' @param network.par list, stores all related datasets from network construction pipeline step.
#' @param analysis.par list, stores all related datasets from driver analysis pipeline step.
#' @param step character, name of the pipeline step decided by user for easier reference.
#'
#' @examples
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' NetBID.saveRData(analysis.par=analysis.par,step='ms-tab_test')
#' }
#' @export
NetBID.saveRData <- function(network.par=NULL,analysis.par=NULL,step='exp-load'){
  #
  all_input_para <- c('step')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(network.par)==FALSE & is.null(analysis.par)==FALSE){
    message('Can not save network.par and analysis.par at once, please only use one !');return(FALSE)
  }
  if(is.null(network.par)==FALSE){
    save(network.par,file=sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step))
    message(sprintf('Successful save to %s',sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step)))
  }
  if(is.null(analysis.par)==FALSE){
    save(analysis.par,file=sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step))
    message(sprintf('Successful save to %s',sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step)))
  }
}

#' Reload Saved RData Created by \code{NetBID.saveRData}.
#'
#' \code{NetBID.loadRData} is a function loads RData saved by \code{NetBID.saveRData} function.
#' It prevents user from repeating former pipeline steps.
#'
#' @param network.par list, stores all related datasets from network construction step.
#' @param analysis.par list, stores all related datasets from driver analysis step.
#' @param step character, name of the pipeline step. It should be previously assigned by user when calling \code{NetBID.saveRData} function.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#'
#' @export
NetBID.loadRData <- function(network.par=NULL,analysis.par=NULL,step='exp-load'){
  #
  all_input_para <- c('step')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(network.par)==FALSE & is.null(analysis.par)==FALSE){
    message('Can not load network.par and analysis.par at once, please only use one !');return(FALSE)
  }
  if(is.null(network.par)==FALSE){
    load(file=sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step),.GlobalEnv)
    message(sprintf('Successful load from %s',sprintf('%s/network.par.Step.%s.RData',network.par$out.dir.DATA,step)))
  }
  if(is.null(analysis.par)==FALSE){
    load(file=sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step),.GlobalEnv)
    message(sprintf('Successful load from %s',sprintf('%s/analysis.par.Step.%s.RData',analysis.par$out.dir.DATA,step)))
  }
}

#' Download Gene Expression Series From GEO Database with Platform Specified
#'
#' \code{load.exp.GEO} downloads user assigned Gene Expression Series (GSE file) along with its Platform from GEO dataset.
#' It returns an ExpressionSet class object and saves it as RData. If the GSE RData already exists, it will be loaded directly.
#' It also allows users to update the Gene Expression Series RData saved before.
#'
#' @param out.dir character, the file path used to save the GSE RData. If the data already exsits, it will be loaded from this path.
#' @param GSE character, the GEO Series Accession ID.
#' @param GPL character, the GEO Platform Accession ID.
#' @param getGPL logical, if TRUE, the corresponding GPL file will be downloaded. Default is TRUE.
#' @param update logical, if TRUE, the previous stored Gene ExpressionSet RData will be updated. Default is FALSE
#'
#' @return Return an ExpressionSet class object.
#' @examples
#'
#' \dontrun{
#' # Download the GSE116028 which performed on GPL6480 platform
#' # from GEO and save it to the current directory
#' # Assign this ExpressionSet object to net_eset
#' net_eset <- load.exp.GEO(out.dir='./',
#'                          GSE='GSE116028',
#'                          GPL='GPL6480',
#'                          getGPL=TRUE,
#'                          update=FALSE)
#' }
#' @export
load.exp.GEO <- function(out.dir = NULL,GSE = NULL,GPL = NULL,getGPL=TRUE,update = FALSE){
  #
  all_input_para <- c('out.dir','GSE','GPL')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('getGPL',c(TRUE,FALSE),envir=environment()),
                 check_option('update',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(!grepl('^GSE',GSE)){
    message('Only support GSE ID')
    return(FALSE)
  }
  expRData_dir <- sprintf('%s/%s_%s.RData', out.dir, GSE,GPL)
  if (file.exists(expRData_dir) & update == FALSE) {
    message(sprintf('RData exist in %s and update==TRUE, will directly load from RData .',expRData_dir))
    load(expRData_dir)
  } else{
    eset <- GEOquery::getGEO(GSE, GSEMatrix = TRUE, getGPL = getGPL)
    if (base::length(eset) > 1)
      idx <- grep(GPL, attr(eset, "names"))
    else
      idx <- 1
    eset <- eset[[idx]]
    if(GPL!=annotation(eset)) {GPL <- annotation(eset); message(sprintf('Real GPL:%s',GPL))}
    expRData_dir <- sprintf('%s/%s_%s.RData', out.dir, GSE,GPL)
    save(eset, file = expRData_dir)
    message(sprintf('RData for the eset is saved in %s .',expRData_dir))
  }
  return(eset)
}

#' Load Gene Expression Set from Salmon Output (demo version)
#'
#' \code{load.exp.RNASeq.demoSalmon} is a function to read in Salmon results and convert it to eSet/DESeqDataSet class object.
#'
#' This function helps users to read in results created by Salmon.
#' Due to the complicated manipulations (e.g. reference sequence) in processing Salmon, this demo function may not be suitable for all scenarios.
#'
#' @param salmon_dir character, the directory to save the results created by Salmon.
#' @param tx2gene data.frame or NULL, this parameter will be passed to \code{tximport}. For details, please check \code{tximport}.
#' If NULL, will read in one of the transcript names from Salmon's results. Note, it works when using e.g. "gencode.v29.transcripts.fa" from GENCODE as reference.
#' @param use_phenotype_info data.frame, the data frame contains phenotype information. It must have the columns \code{use_sample_col} and \code{use_design_col}.
#' @param use_sample_col character, the column name, indicating which column in \code{use_phenotype_info} should be used as the sample name.
#' @param use_design_col character, the column name, indicating which column in \code{use_phenotype_info} should be used as the design feature for samples.
#' @param return_type character, the class of the return object.
#' "txi" is the output of tximport. It is a list containing three matrices, abundance, counts and length.
#' "counts" is the matrix of raw count.
#' "tpm" is the raw tpm.
#' "fpm", "cpm" is the fragments/counts per million mapped fragments (fpm/cpm).
#' "raw-dds" is the DESeqDataSet class object, which is the original one without processing.
#' "dds" is the DESeqDataSet class object, which is processed by DESeq.
#' "eset" is the ExpressionSet class object, which is processed by DESeq and vst.
#' Default is "tpm".
#' @param merge_level character, users can choose between "gene" and "transcript".
#' "gene", the original salmon results will be mapped to the transcriptome and the expression matrix will be merged to the gene level.
#' This only works when using e.g. "gencode.vXX.transcripts.fa" from GENCODE as the reference.
#' @export
load.exp.RNASeq.demoSalmon <- function(salmon_dir = NULL,tx2gene=NULL,
                                       use_phenotype_info = NULL,
                                       use_sample_col=NULL,
                                       use_design_col=NULL,
                                       return_type='tpm',
                                       merge_level='gene') {
  #
  all_input_para <- c('salmon_dir','use_phenotype_info','use_sample_col','use_design_col','return_type','merge_level')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('return_type',c('txi','counts','tpm','fpm','cpm','raw-dds','dds','eset'),envir=environment()),
                 check_option('merge_level',c('gene','transcript'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  files <- file.path(salmon_dir, list.files(salmon_dir), "quant.sf")
  if(!grepl('_salmon',use_phenotype_info[1,use_sample_col])) sample_name <- gsub('(.*)_salmon', '\\1', list.files(salmon_dir)) ## if no _salmon is also ok
  names(files) <- sample_name
  w1 <- base::length(files)
  message(sprintf('%d %s/quant.sf found !',w1,salmon_dir))
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
  eset <- load.exp.RNASeq.demo(files,type='salmon',
                               tx2gene=tx2gene,
                               use_phenotype_info=use_phenotype_info,
                               use_sample_col=use_sample_col,
                               use_design_col=use_design_col,
                               return_type=return_type,
                               merge_level=merge_level)
  return(eset)
}

#' Load Gene Expression Set from RNA-Seq Results (demo version)
#'
#' \code{load.exp.RNASeq.demo} is a function to read in RNA-Seq results and convert it to \code{eSet/DESeqDataSet} class object.
#'
#' This function helps users to read in RNA-Seq results from various sources.
#' Due to the complicated manipulations (e.g. reference sequence) in processing RNA-Seq, this demo function may not be suitable for all scenarios.
#'
#' @param files a vector of characters, the filenames for the transcript-level abundances. It will be passed to \code{tximport}.
#' For details, please check \code{tximport}.
#' @param type character, the type of software used to generate the abundances. It will be passed to \code{tximport}.
#' For details, please check \code{tximport}.
#' @param tx2gene data.frame or NULL, this parameter will be passed to \code{tximport}. For details, please check \code{tximport}.
#' @param use_phenotype_info data.frame, the data frame contains phenotype information. It must have the columns \code{use_sample_col} and \code{use_design_col}.
#' @param use_sample_col character, the column name, indicating which column in \code{use_phenotype_info} should be used as the sample name.
#' @param use_design_col character, the column name, indicating which column in \code{use_phenotype_info} should be used as the design feature for samples.
#' @param return_type character, the class of the return object.
#' "txi" is the output of \code{tximport}. It is a list containing three matrices, abundance, counts and length.
#' "counts" is the matrix of raw count.
#' "tpm" is the raw tpm.
#' "fpm", "cpm" is the fragments/counts per million mapped fragments.
#' "raw-dds" is the DESeqDataSet class object, which is the original one without processing.
#' "dds" is the DESeqDataSet class object, which is processed by \code{DESeq}.
#' "eset" is the ExpressionSet class object, which is processed by \code{DESeq} and \code{vst}.
#' Default is "tpm".
#' @param merge_level character, users can choose between "gene" and "transcript".
#' "gene", the original salmon results will be mapped to the transcriptome and the expression matrix will be merged to the gene level.
#' This only works when using e.g. "gencode.vXX.transcripts.fa" from GENCODE as the reference.
#' @export
load.exp.RNASeq.demo <- function(files,type='salmon',
                                 tx2gene=NULL,
                                 use_phenotype_info = NULL,
                                 use_sample_col=NULL,
                                 use_design_col=NULL,
                                 return_type='tpm',
                                 merge_level='gene') {
  #
  all_input_para <- c('files','type','tx2gene','use_phenotype_info','use_sample_col','use_design_col','return_type','merge_level')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('return_type',c('txi','counts','tpm','fpm','cpm','raw-dds','dds','eset'),envir=environment()),
                 check_option('merge_level',c('gene','transcript'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  n1 <- colnames(use_phenotype_info)
  if(!use_sample_col %in% n1){
    message(sprintf('%s not in the colnames of use_phenotype_info,
                    please check and re-try !',use_sample_col));return(FALSE)
  }
  if(!use_design_col %in% n1){
    message(sprintf('%s not in the colnames of use_phenotype_info,
                    please check and re-try !',use_design_col));return(FALSE)
  }
  # get intersected samples
  rownames(use_phenotype_info) <- use_phenotype_info[,use_sample_col]
  w1 <- base::intersect(names(files),rownames(use_phenotype_info))
  files <- files[w1]; use_phenotype_info <- use_phenotype_info[w1,]
  if(base::length(w1)==0){
    message(sprintf('No sample could match the %s in the use_phenotype_info, please check and re-try !',use_sample_col))
    return(FALSE)
  }
  message(sprintf('%d samples could further processed !',base::length(w1)))
  # import into txi
  txi <- tximport::tximport(files, type = type, tx2gene = tx2gene) ## key step one, tximport
  if(return_type=='counts'){
    return(txi$counts)
  }
  if(return_type=='tpm'){
    return(txi$abundance)
  }
  if(return_type=='txi'){
    return(txi)
  }
  use_phenotype_info <- use_phenotype_info[colnames(txi$abundance), ]
  tmp_phe <- base::cbind(group=use_phenotype_info[,use_design_col],use_phenotype_info,stringsAsFactors=FALSE)
  # import into deseq2
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = tmp_phe, design =  ~ group) ## key step two, DESeqDataSetFromTximport
  if(return_type=='raw-dds'){
    return(dds)
  }
  if(return_type=='fpm' | return_type=='cpm'){
    return(DESeq2::fpm(dds))
  }
  dds <- DESeq2::DESeq(dds)
  if(return_type=='dds'){
    return(dds)
  }else{
    vsd <- DESeq2::vst(dds)
    mat <- SummarizedExperiment::assay(vsd)
    eset <- generate.eset(exp_mat=mat, phenotype_info = use_phenotype_info, feature_info = NULL, annotation_info='Salmon')
    if(return_type=='eset') return(eset)
    if(return_type=='both') return(list(eset=eset,dds=dds))
  }
  }

#' Generate ExpressionSet Object
#'
#' \code{generate.eset} generates ExpressionSet class object to contain and describe the high-throughput assays.
#' Users need to define its slots, which are expression matrix (required),
#' phenotype information and feature information (optional).
#' It is very useful when only expression matrix is available.
#'
#' @param exp_mat matrix, the expression data matrix. Each row represents a gene/transcript/probe, each column represents a sample.
#' @param phenotype_info data.frame, the phenotype information for all the samples in \code{exp_mat}.
#' In the phenotype data frame, each row represents a sample, each column represents a phenotype feature.
#' The row names must match the column names of \code{exp_mat}. If NULL, it will generate a single-column data frame.
#' Default is NULL.
#' @param feature_info data.frame, the feature information for all the genes/transcripts/probes in \code{exp_mat}.
#' In the feature data frame, each row represents a gene/transcript/probe and each column represents an annotation of the feature.
#' The row names must match the row names of \code{exp_mat}. If NULL, it will generate a single-column data frame.
#' Default is NULL.
#' @param annotation_info character, the annotation set by users for easier reference. Default is "".
#'
#' @return Return an ExressionSet object.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000),nrow=1000,ncol=10)
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' eset <- generate.eset(exp_mat=mat1)
#' @export
generate.eset <- function(exp_mat=NULL, phenotype_info=NULL, feature_info=NULL, annotation_info="") {
  #
  all_input_para <- c('exp_mat')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(dim(exp_mat))==TRUE){
    exp_mat <- t(as.matrix(exp_mat));rownames(exp_mat) <- 'g1'
  }
  if (is.null(phenotype_info)) {
    phenotype_info <- data.frame(group = colnames(exp_mat), stringsAsFactors = FALSE)
    rownames(phenotype_info) <- colnames(exp_mat)
  }
  if (is.null(feature_info)) {
    feature_info <- data.frame(gene = rownames(exp_mat), stringsAsFactors = FALSE)
    rownames(feature_info) <- rownames(exp_mat)
  }
  if((class(phenotype_info)=='character' | is.null(dim(phenotype_info))==TRUE) & is.null(names(phenotype_info))==TRUE){
    phenotype_info <- data.frame(group = phenotype_info, stringsAsFactors = FALSE)
    rownames(phenotype_info) <- colnames(exp_mat)
  }
  if((class(feature_info)=='character' | is.null(dim(feature_info))==TRUE) & is.null(names(feature_info))==TRUE){
    feature_info <- data.frame(gene = feature_info, stringsAsFactors = FALSE)
    rownames(feature_info) <- rownames(exp_mat)
  }
  #
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

#' Merge Two ExpressionSet Class Objects into One
#'
#' \code{merge_eset} merges two ExpressionSet class objects and returns one ExpresssionSet object.
#' If genes in the two ExpressionSet objects are identical, the expression matrix will be combined directly.
#' Otherwise, Z-transformation is strongly suggested to be performed before combination (set std=TRUE).
#'
#' @param eset1 ExpressionSet class, the first ExpressionSet.
#' @param eset2 ExpressionSet class, the second ExpressionSet.
#' @param group1 character, name of the first ExpressionSet.
#' @param group2 character, name of the second ExpressionSet.
#' @param use_col a vector of characters, the column names in the phenotype information to be kept.
#' If NULL, shared column names of \code{eset1} and \code{eset2} will be used. Default is NULL.
#' @param group_col_name character, name of the column which contains the names defined in \code{group1} and \code{group2}.
#' This column is designed to show which original ExpressionSet each sample comes from before combination.
#' Default name of this column is "original_group".
#' @param remove_batch logical, if TRUE, remove the batch effects from these two expression datasets. Default is FALSE.
#' @param std logical, whether to perform std to the original expression matrix. Default is FALSE.
#'
#' @return Return an ExressionSet class object.
#' @examples
#' mat1 <- matrix(rnorm(10000),nrow=1000,ncol=10)
#' colnames(mat1) <- paste0('Sample1_',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' eset1 <- generate.eset(exp_mat=mat1)
#' mat2 <- matrix(rnorm(10000),nrow=1000,ncol=10)
#' colnames(mat2) <- paste0('Sample2_',1:ncol(mat1))
#' rownames(mat2) <- paste0('Gene',1:nrow(mat1))
#' eset2 <- generate.eset(exp_mat=mat2)
#' new_eset <- merge_eset(eset1,eset2)
#' @export
merge_eset <- function(eset1,eset2,
                       group1=NULL,group2=NULL,
                       group_col_name='original_group',
                       use_col = NULL,
                       remove_batch = FALSE,std=FALSE) {
  #
  all_input_para <- c('eset1','eset2','group_col_name')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('remove_batch',c(TRUE,FALSE),envir=environment()),
                 check_option('std',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  mat1 <- Biobase::exprs(eset1)
  mat2 <- Biobase::exprs(eset2)
  w1 <- base::intersect(rownames(mat1), rownames(mat2))
  if(base::length(w1)==0){
    message('No overlap genes between two eSet, please check and re-try!');return(FALSE);
  }
  if((base::length(w1)<nrow(mat1) | base::length(w1)<nrow(mat2)) & std==FALSE){
    message('Original two esets contain different gene list, strongly suggest to do z transformation (set std=TRUE) across all samples before merge!');
  }
  if(std==TRUE){
    ## z-transformation
    mat1 <- apply(mat1,2,do.std) # std to samples
    mat2 <- apply(mat2,2,do.std)
  }
  rmat <- base::cbind(as.data.frame(mat1)[w1, ], as.data.frame(mat2)[w1,])
  rmat <- as.matrix(rmat)
  #choose1 <- apply(rmat <= quantile(rmat, probs = 0.05), 1, sum) <= ncol(rmat) * 0.90 ## low expressed genes
  #rmat <- rmat[choose1, ]
  phe1 <- Biobase::pData(eset1)
  phe2 <- Biobase::pData(eset2)
  phe1 <- as.data.frame(apply(phe1,2,clean_charVector),stringsAsFactors=F)
  phe2 <- as.data.frame(apply(phe2,2,clean_charVector),stringsAsFactors=F)
  if(base::length(use_col)==0){
    use_col <- base::intersect(colnames(phe1),colnames(phe2))
  }
  rphe <- list();
  if(base::length(use_col)>1)
    rphe <- base::rbind(phe1[colnames(mat1), use_col], phe2[colnames(mat2), use_col])
  if(base::length(use_col)==1){
    rphe <- c(phe1[colnames(mat1), use_col], phe2[colnames(mat2), use_col])
    rphe <- data.frame(rphe,stringsAsFactors=FALSE); colnames(rphe) <- use_col;
    rownames(rphe) <- colnames(rmat)
  }
  if(base::length(use_col)==0){message('Warning: no intersected phenotype column!');}
  if(is.null(group1)==TRUE) group1 <- 'group1'
  if(is.null(group2)==TRUE) group2 <- 'group2'
  rphe[[group_col_name]]<- c(rep(group1, ncol(mat1)), rep(group2, ncol(mat2)))
  if (remove_batch == TRUE) {
    rmat <- limma::removeBatchEffect(rmat,batch=rphe[[group_col_name]])
  }
  if(class(rphe)=='list'){rphe <- as.data.frame(rphe,stringsAsFactors=FALSE); rownames(rphe) <- colnames(rmat)}
  reset <- generate.eset(rmat,phenotype_info = rphe, annotation_info = 'combine')
  return(reset)
}

#' Reassign featureData slot of ExpressionSet and Update feature information
#'
#' \code{update_eset.feature} reassigns the featureData slot of ExpressionSet object based on user's demand. It is mainly used for gene ID conversion.
#'
#' User can pass a conversion table to \code{use_eset} for the ID conversion. A conversion table can be obtained from the original featureDta slot
#' (if one called the \code{load.exp.GEO} function and set getGPL==TRUE) or by running the \code{get_IDtransfer} function.
#' The mapping between original ID and target ID can be summerised into 4 categories.
#' 1) One-to-one, simply replaces the original ID with target ID;
#' 2) Many-to-one, the expression value for the target ID will be merged from its original ID;
#' 3) One-to-many, the expression value for the original ID will be distributed to the matched target IDs;
#' 4) Many-to-many, apply part 3) first, then part 2).
#'
#' @param use_eset ExpressionSet class object.
#' @param use_feature_info data.frame, a data frame contains feature information, it can be obtained by calling \code{fData} function.
#' @param from_feature character, orginal ID. Must be one of the column names in \code{use_feature_info} and correctly characterize the \code{use_eset}'s row names.
#' @param to_feature character, target ID. Must be one of the column names in \code{use_feature_info}.
#' @param merge_method character, the agglomeration method to be used for merging gene expression value.
#' This should be one of, "median", "mean", "max" or "min". Default is "median".
#' @param distribute_method character, the agglomeration method to be used for distributing the gene expression value.
#' This should be one of, "mean" or "equal". Default is "equal".
#'
#' @return Return an ExressionSet object with updated feature information.
#' @examples
#' mat1 <- matrix(rnorm(10000),nrow=1000,ncol=10)
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' eset <- generate.eset(exp_mat=mat1)
#' test_transfer_table <- data.frame(
#'                        'Gene'=c('Gene1','Gene1','Gene2','Gene3','Gene4'),
#'                        'Transcript'=c('T11','T12','T2','T3','T3'))
#' new_eset <- update_eset.feature(use_eset=eset,
#'                                 use_feature_info=test_transfer_table,
#'                                 from_feature='Gene',
#'                                 to_feature='Transcript',
#'                                 merge_method='median',
#'                                 distribute_method='equal'
#'                                 )
#' print(Biobase::exprs(eset)[test_transfer_table$Gene,])
#' print(Biobase::exprs(new_eset))
#'
#' @export update_eset.feature
update_eset.feature <- function(use_eset=NULL,use_feature_info=NULL,from_feature=NULL,to_feature=NULL,
                                merge_method='median',distribute_method='equal'){
  #
  all_input_para <- c('use_eset','from_feature','to_feature')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('merge_method',c("median","mean","max","min"),envir=environment()),
                 check_option('distribute_method',c('mean','equal'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(use_feature_info)) use_feature_info <- Biobase::fData(use_eset)
  n1 <- colnames(use_feature_info)
  if(!from_feature %in% n1){
    message(sprintf('%s not in in the colnames of use_feature_info, please re-try!',from_feature));return(use_eset)
  }
  if(!to_feature %in% n1){
    message(sprintf('%s not in in the colnames of use_feature_info, please re-try!',to_feature));return(use_eset)
  }
  mat <- Biobase::exprs(use_eset)
  use_feature_info <- base::unique(use_feature_info);
  w1 <- which(use_feature_info[,1]!="" & use_feature_info[,2]!="" & is.na(use_feature_info[,1])==FALSE & is.na(use_feature_info[,2])==FALSE)
  use_feature_info <- use_feature_info[w1,]
  g1 <- rownames(mat) ## rownames for the expmat
  f1 <- clean_charVector(use_feature_info[,from_feature]) ## from feature info
  t1 <- clean_charVector(use_feature_info[,to_feature]) ## to feature info
  w1 <- which(f1 %in% g1); f1 <- f1[w1]; t1 <- t1[w1]; ## only consider features in the rownames of expmat
  if(base::length(w1)==0){
    message(sprintf('Rownames of the expression matrix was not included in the %s column, please check and re-try !',from_feature))
    return(use_eset)
  }
  message(sprintf('%d transfer pairs related with %d rows from original expression matrix will be keeped !',base::length(w1),base::length(g1)))
  fc1 <- base::table(f1); tc1 <- base::table(t1); fw1 <- which(fc1>1); tw1 <- which(tc1>1); ## check duplicate records
  if(base::length(fw1)>0){
    message(sprintf('Original feature %s has %d items with duplicate records, will distribute the original values equal to all related items !
                    if do not want this, please check and retry !',from_feature,base::length(fw1)))
    #return(use_eset)
    w2 <- which(f1 %in% names(fw1)) ## need to distribute
    w0 <- base::setdiff(1:base::length(f1),w2) ## do not need to distribute
    if(distribute_method=='equal'){
      v1 <- mat[f1[w2],]; ## distribute equal
    }
    if(distribute_method=='mean'){
      v1 <- mat[f1[w2],]; ## distribute mean
      tt <- as.numeric(base::table(f1[w2])[f1[w2]])
      v1 <- v1/tt;
    }
    rownames(v1) <- paste0(f1[w2],'-',t1[w2]);
    f1[w2] <- paste0(f1[w2],'-',t1[w2]); # update transfer table
    mat <- base::rbind(v1,mat[f1[w0],]) # update mat table
    fc1 <- base::table(f1); tc1 <- base::table(t1); fw1 <- which(fc1>1); tw1 <- which(tc1>1); ## update f1, t1 and related values
  }
  if(base::length(tw1)>0){
    w2 <- which(t1 %in% names(tw1)) ## need to merge
    w0 <- base::setdiff(1:base::length(t1),w2) ## do not need to merge
    mat_new_0 <- mat[f1[w0],]; rownames(mat_new_0) <- t1[w0]  ## mat do not need to merge
    if(merge_method=='mean') tmp1 <- stats::aggregate(mat[f1[w2],,drop=FALSE],list(t1[w2]),function(x){base::mean(x,na.rm=TRUE)})
    if(merge_method=='median') tmp1 <- stats::aggregate(mat[f1[w2],,drop=FALSE],list(t1[w2]),function(x){stats::median(x,na.rm=TRUE)})
    if(merge_method=='max') tmp1 <- stats::aggregate(mat[f1[w2],,drop=FALSE],list(t1[w2]),function(x){base::max(x,na.rm=TRUE)})
    if(merge_method=='min') tmp1 <- stats::aggregate(mat[f1[w2],,drop=FALSE],list(t1[w2]),function(x){base::min(x,na.rm=TRUE)})
    mat_new_1 <- tmp1[,-1]; rownames(mat_new_1) <- tmp1[,1] ## mat merged
    mat_new <- base::rbind(mat_new_0,mat_new_1)
  }else{
    mat_new <- mat[f1,]
    rownames(mat_new) <- t1
  }
  new_eset <- generate.eset(exp_mat=mat_new, phenotype_info=Biobase::pData(use_eset), feature_info=NULL, annotation_info=annotation(use_eset))
  return(new_eset)
}

#' Reassign the phenoData slot of ExpressionSet and Update phenotype information
#'
#' \code{update_eset.phenotype} reassigns the phenoData slot of ExpressionSet based on user's demand.
#' It is mainly used to modify sample names and extract interested phenotype information for further sample clustering.
#'
#' @param use_eset ExpressionSet class object.
#' @param use_phenotype_info data.frame, a dataframe contains phenotype information, can be obtained by calling \code{pData} function.
#' @param use_sample_col character, must be one of the column names in \code{use_phenotype_info}.
#' @param use_col character, the columns will be kept from \code{use_phenotype_info}.
#' 'auto', only extracting 'cluster-meaningful' sample features (e.g. it is meaningless to use 'gender' as clustering feature, if all samples are female).
#' 'GEO-auto' means it will extract the following selected columns,
#' "geo_accession", "title", "source_name_ch1", and columns ended with ":ch1". Default is "auto".
#' @return Return an ExressionSet object with updated phenotype information.
#' @examples
#' \dontrun{
#' net_eset <- load.exp.GEO(out.dir='./test',
#'                          GSE='GSE116028',
#'                          GPL='GPL6480',
#'                          getGPL=TRUE,
#'                          update=FALSE)
#' net_eset <- update_eset.phenotype(use_eset=net_eset,
#'                                   use_phenotype_info=Biobase::pData(net_eset),
#'                                   use_sample_col='geo_accession',
#'                                   use_col='GEO-auto')
#' }
#' @export update_eset.phenotype
update_eset.phenotype <- function(use_eset=NULL,use_phenotype_info=NULL,use_sample_col=NULL,use_col='auto'){
  #
  all_input_para <- c('use_eset')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(use_eset)){
    message('use_eset required, please re-try !');
    return(use_eset)
  }
  if(is.null(use_phenotype_info)) use_phenotype_info <- Biobase::pData(use_eset)
  if(is.null(use_sample_col)==FALSE){
    if(!use_sample_col %in% colnames(use_phenotype_info)){
      stop(sprintf('%s not in the colnames of use_phenotype_info, please re-try!',use_sample_col));#return(use_eset)
    }
  }
  if(is.null(use_col)) use_col <- colnames(use_phenotype_info)
  mat <- Biobase::exprs(use_eset)
  s1 <- colnames(mat) ## all samples
  if(is.null(use_sample_col)==TRUE){
    p1 <- rownames(use_phenotype_info)
  }else{
    p1 <- use_phenotype_info[,use_sample_col] ## sample in the phenotype info
  }
  w1 <- which(p1 %in% s1); p1 <- p1[w1]; ## only consider samples in the colnames of expmat
  if(base::length(w1)==0){
    if(is.null(use_sample_col)==TRUE){
      message('Colnames of the expression matrix was not included in the rownames of use_phenotype_info, please check and re-try !')
    }else{
      message(sprintf('Colnames of the expression matrix was not included in the %s column, please check and re-try !',use_sample_col))
    }
    return(use_eset)
  }
  message(sprintf('%d out of %d samples from the expression matrix will be keeped !',base::length(w1),base::length(s1)))
  mat_new <- mat[,p1]
  use_phenotype_info <- use_phenotype_info[w1,]
  n1 <- colnames(use_phenotype_info)
  if(use_col[1] == 'GEO-auto'){
    w1 <- c('geo_accession','title','source_name_ch1',n1[grep(':ch1',n1)])
    p1 <- use_phenotype_info[,w1]
    colnames(p1)[4:ncol(p1)] <- gsub('(.*):ch1','\\1',colnames(p1)[4:ncol(p1)])
    colnames(p1)[3] <- gsub('(.*)_ch1','\\1',colnames(p1)[3])
    if(is.null(use_sample_col)==FALSE) rownames(p1) <- use_phenotype_info[,use_sample_col]
    if(base::length(w1)>1) p1 <- as.data.frame(apply(p1,2,clean_charVector),stringsAsFactors=FALSE)
    if(base::length(w1)==1) p1 <- as.data.frame(clean_charVector(p1),stringsAsFactors=FALSE)
    new_phenotype_info <- p1;
  }else{
    if(use_col[1] == 'auto'){
      u1 <- apply(use_phenotype_info,2,function(x)base::length(base::unique(x)))
      w1 <- which(u1>=2 & u1<=nrow(use_phenotype_info)-1)
      if(base::length(w1)==0){
        message('No column could match the auto criteria, please check and re-try!');return(FALSE)
      }
      p1 <- use_phenotype_info[,w1]
      if(base::length(w1)>1) p1 <- as.data.frame(apply(p1,2,clean_charVector),stringsAsFactors=FALSE)
      if(base::length(w1)==1) p1 <- as.data.frame(clean_charVector(p1),stringsAsFactors=FALSE)
      new_phenotype_info <- use_phenotype_info[,w1];names(new_phenotype_info) <- names(use_phenotype_info)[w1];
    }else{
      if(base::length(base::setdiff(use_col,n1))>0){
        message(sprintf('%s not in use_phenotype_info, please re-try!',base::paste(base::setdiff(use_col,n1),collapse=';')));return(FALSE)
      }
      p1 <- use_phenotype_info[,use_col]
      if(base::length(use_col)>1) p1 <- as.data.frame(apply(p1,2,clean_charVector),stringsAsFactors=FALSE)
      if(base::length(use_col)==1) p1 <- as.data.frame(clean_charVector(p1),stringsAsFactors=FALSE)
      new_phenotype_info <- p1; names(new_phenotype_info) <- use_col;
    }
  }
  rownames(new_phenotype_info) <- rownames(use_phenotype_info)
  #print(new_phenotype_info)
  message(sprintf('%d out of %d sample features will be keeped !',ncol(new_phenotype_info),ncol(use_phenotype_info)))
  new_eset <- generate.eset(exp_mat=mat_new, phenotype_info=new_phenotype_info, feature_info=Biobase::fData(use_eset), annotation_info=annotation(use_eset))
  return(new_eset)
}

#' IQR (interquartile range) filter to extract genes from expression matrix
#'
#' \code{IQR.filter} is a function to extract genes from the expression matrix by setting threshold to their IQR value.
#' IQR (interquartile range) is a measure of statistical dispersion. It is calculated for each gene across all the samples.
#' By setting threshold value, genes with certain statistical dispersion across samples will be filtered out.
#' This step is mainly used to perform sample cluster and to prepare the input for SJAracne.
#'
#' @param exp_mat matrix, the gene expression matrix. Each row represents a gene/transcript/probe, each column represents a sample.
#' @param use_genes a vector of characters, the gene list needed to be filtered. Default is the row names of \code{exp_mat}.
#' @param thre numeric, the threshold for IQR of the genes in \code{use_genes}. Default is 0.5.
#' @param loose_gene a vector of characters, the gene list that only need to pass the \code{loose_thre}.
#' This parameter is designed for the input of possible drivers used in SJAracne. Default is NULL.
#' @param loose_thre numeric, the threshold for IQR of the genes in \code{loose_gene}. Default is 0.1.
#' @return Return a vector with logical values indicate which genes should be kept.
#' @examples
#' mat1 <- matrix(rnorm(15000),nrow=1500,ncol=10)
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' choose1 <- IQR.filter(mat1,thre=0.5,
#'                      loose_gene=paste0('Gene',1:100))
#' @export
IQR.filter <- function(exp_mat,use_genes=rownames(exp_mat),thre = 0.5,loose_gene=NULL,loose_thre=0.1) {
  #
  all_input_para <- c('exp_mat','use_genes','thre')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  use_genes <- base::intersect(use_genes,rownames(exp_mat))
  use_genes <- base::setdiff(use_genes,"")
  exp_mat <- exp_mat[use_genes,]
  iqr <- apply(exp_mat, 1, stats::IQR) ## calculate IQR for each gene
  choose0 <- use_genes[iqr > quantile(iqr, loose_thre)] ## for loose_gene
  choose1 <- use_genes[iqr > quantile(iqr, thre)] ## for all genes
  choose2 <- base::unique(c(base::intersect(loose_gene, choose0), choose1)) ## union set
  use_vec <- rep(FALSE,length.out=base::length(use_genes));names(use_vec) <- use_genes
  use_vec[choose2] <- TRUE
  print(base::table(use_vec))
  return(use_vec)
}

#' Normalization of RNA-Seq Reads Count
#'
#' \code{RNASeqCount.normalize.scale} is a simple version to normalize the RNASeq reads count data.
#'
#' Users can also load \code{load.exp.RNASeq.demo}, and follow the \code{DESeq2} pipeline for RNASeq data processing.
#' Warning, \code{load.exp.RNASeq.demo} and \code{load.exp.RNASeq.demoSalmon} in NetBID2 may not cover all the possible scenarios.
#'
#' @param mat matrix, matrix of RNA-Seq reads data. Each row is a gene/transcript, each column is a sample.
#' @param total integer, total RNA-Seq reads count. If NULL, will use the mean of each column's summation. Default is NULL.
#' @param pseudoCount integer, the integer added to avoid "-Inf" showing up during log transformation. Default is 1.
#'
#' @return Return a numeric matrix, containing the normalized RNA-Seq reads count.
#'
#' @examples
#' mat1 <- matrix(rnbinom(10000, mu = 10, size = 1),nrow=1000,ncol=10)
#' colnames(mat1) <- paste0('Sample1',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' norm_mat1 <- RNASeqCount.normalize.scale(mat1)
#' @export
RNASeqCount.normalize.scale <- function(mat,
                                        total = NULL,
                                        pseudoCount = 1) {
  #
  all_input_para <- c('mat','pseudoCount')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  d <- mat
  if (!is.data.frame(d))
    d <- data.frame(d)
  if (!all(d > 0))
    d <- d + pseudoCount
  s <- apply(d, 2, sum)
  m <-
    ifelse(is.null(total), as.integer(base::mean(s)), as.integer(total)) ## total or mean sum
  options(digits = 2 + nchar(m))
  fac <- m / s
  for (i in 1:base::length(s)) {
    d[, i] <- d[, i] * fac[i]
    #d[, i] <- round(d[, i] * fac[i], 0)
  }
  if (!all(d > 0))
    d <- d + pseudoCount
  d
}

## inner function for dist2
dist2.mod <- function (x, fun = function(a, b) base::mean(abs(a - b), na.rm = TRUE),
                       diagonal = 0)
{
  if (!(is.numeric(diagonal) && (base::length(diagonal) == 1)))
    stop("'diagonal' must be a numeric scalar.")
  if (missing(fun)) {
    res = apply(x, 2, function(w) base::colMeans(abs(x - w), na.rm = TRUE))
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
    n <- base::length(z)
    m1 <- ifelse(sum(z > 0) > 0, sum(z[z > 0]) / n, 0)
    m2 <- ifelse(sum(z < 0) > 0, sum(z[z < 0]) / n, 0)
    if (m1 > -m2)
      es <- m1
    else
      es <- m2
  }
  else if (es.method == 'absmean') {
    es <- base::mean(abs(z),na.rm=TRUE)
  }
  else if (es.method == 'mean') {
    es <- base::mean(z,na.rm=TRUE)
  }
  else if (es.method == 'median') {
    es <- stats::median(z,na.rm=TRUE)
  }
  else if (es.method == 'max') {
    es <- base::max(z,na.rm=TRUE)
  }
  else if (es.method == 'min') {
    es <- base::min(z,na.rm=TRUE)
  }
  return(es)
}
do.std <- function(x) {
  x <- x[!is.na(x)]
  (x - base::mean(x,na.rm=TRUE)) / sd(x,na.rm=TRUE)
}

#' Calculate Activity Value for Each Driver
#'
#' \code{cal.Activity} calculates the activity value for each driver.
#' This function requires two inputs, the driver-to-target list object \code{target_list} and the expression matrix.
#'
#' @param target_list list, the driver-to-target list object. Either igraph_obj or target_list is necessary for this function.
#' The names of the list elements are drivers.
#' Each element is a data frame, usually contains at least three columns.
#' "target", target gene names;
#' "MI", mutual information;
#' "spearman", spearman correlation coefficient.
#' "MI" and "spearman" is necessary if es.method="weightedmean".
#' Users can call \code{get.SJAracne.network} to get this list from the network file generated by SJAracne (the second element of the return list)
#' or prepare the list object by hand but should match the data format described above.
#' @param igraph_obj igraph object, optional. Either igraph_obj or target_list is necessary for this function.
#' Users can call \code{get.SJAracne.network} to get this object from the network file generated by SJAracne (the third element of the return list),
#' or prepare the igraph network object by hand (Directed network and the edge attributes should include "weight" and "sign" if es.method="weightedmean").
#' @param cal_mat numeric matrix, the expression matrix of genes/transcripts.
#' @param es.method character, method applied to calculate the activity value. User can choose from "mean", "weightedmean", "maxmean" and "absmean".
#' Default is "weightedmean".
#' @param std logical, if TRUE, the expression matrix will be normalized by column. Default is TRUE.
#' @param memory_constrain logical, if TRUE, the calculation strategy will not use Matrix Cross Products, which is memory consuming.
#' Default is FALSE.
#' @return Return a matrix of activity values. Rows are drivers, columns are samples.
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,
#'                        cal_mat=Biobase::exprs(analysis.par$cal.eset),
#'                        es.method='weightedmean')
#' ac_mat <- cal.Activity(igraph_obj=analysis.par$merge.network$igraph_obj,
#'                        cal_mat=Biobase::exprs(analysis.par$cal.eset),
#'                        es.method='maxmean')
#' @export
cal.Activity <- function(target_list=NULL, igraph_obj = NULL, cal_mat=NULL, es.method = 'weightedmean',std=TRUE,memory_constrain=FALSE) {
  #
  all_input_para <- c('cal_mat','es.method','std','memory_constrain')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('memory_constrain',c(TRUE,FALSE),envir=environment()),
                 check_option('std',c(TRUE,FALSE),envir=environment()),
                 check_option('es.method',c('mean','weightedmean','maxmean','absmean'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(target_list)==TRUE & is.null(igraph_obj)==TRUE){
    message('Either target_list or igraph_obj is required, please check and re-try!');return(FALSE);
  }
  if(is.null(target_list)==FALSE & memory_constrain==TRUE){
    ac.mat <- cal.Activity.old(target_list=target_list, cal_mat=cal_mat, es.method = es.method,std=std)
    return(ac.mat)
  }
  if(is.null(target_list)==TRUE & memory_constrain==TRUE){
    message('Only accepts target_list input when memory_constrain=TRUE, please check and re-try!');return(FALSE);
  }
  if(nrow(cal_mat)==0){
    message('No genes in the cal_mat, please check and re-try!');return(FALSE);
  }
  if(std==TRUE) cal_mat <- apply(cal_mat, 2, do.std)
  if(is.null(igraph_obj)==FALSE){
    gr <- igraph_obj
    mat1 <- get_igraph2matrix(gr,es.method=es.method)
    mat2 <- get_igraph2matrix(gr,es.method='mean')
    all_source <- get_gr2driver(gr)
  }else{
    if(is.null(target_list)==FALSE){
      mat1 <- get_target_list2matrix(target_list,es.method=es.method)
      mat2 <- get_target_list2matrix(target_list,es.method='mean')
      all_source <- names(target_list)
    }
  }
  ##
  mat1_source <- mat1[all_source,,drop=FALSE]
  w1 <- base::intersect(rownames(cal_mat),colnames(mat1_source))
  if(base::length(w1)==0){
    message('No intersected genes found for the cal_mat and target in the network, please check and re-try!');
    return(FALSE)
  }
  use_mat1_source <- mat1_source[,w1,drop=FALSE] ## network info
  use_mat2_source <- mat2[all_source,w1,drop=FALSE] ## network binary info

  ## weighted mean + mean
  if(es.method %in% c('weightedmean','mean')){
    use_cal_mat <- cal_mat[w1,,drop=FALSE] ## expression info
    out_mat <- use_mat1_source %*% use_cal_mat
    out_mat <- out_mat/Matrix::rowSums(use_mat2_source) ## get mean
  }
  ## absmean
  if(es.method == 'absmean'){
    use_cal_mat <- cal_mat[w1,,drop=FALSE] ## expression info
    out_mat <- use_mat1_source %*% abs(use_cal_mat)
    out_mat <- out_mat/Matrix::rowSums(use_mat2_source) ## get mean
  }
  ## maxmean
  if(es.method == 'maxmean'){
    use_cal_mat <- cal_mat[w1,,drop=FALSE] ## expression info
    use_cal_mat_pos <- use_cal_mat;use_cal_mat_pos[which(use_cal_mat_pos<0)] <- 0;
    use_cal_mat_neg <- use_cal_mat;use_cal_mat_neg[which(use_cal_mat_neg>0)] <- 0;
    out_mat_pos  <- use_mat1_source %*% use_cal_mat_pos
    out_mat_pos  <- out_mat_pos/Matrix::rowSums(use_mat2_source) ## get mean
    out_mat_neg  <- use_mat1_source %*% use_cal_mat_neg
    out_mat_neg  <- out_mat_neg/Matrix::rowSums(use_mat2_source) ## get mean
    out_mat_sign <- sign(abs(out_mat_pos)-abs(out_mat_neg))
    out_mat_sign_pos <- out_mat_sign; out_mat_sign_pos[out_mat_sign_pos!=1] <-0;
    out_mat_sign_neg <- out_mat_sign; out_mat_sign_neg[out_mat_sign_neg!= -1] <-0;
    out_mat <- out_mat_pos*out_mat_sign_pos-out_mat_neg*out_mat_sign_neg
  }
  ## median, min, max , not supported
  # output
  ac.mat <- as.matrix(out_mat)
  w1 <- which(is.na(ac.mat[,1])==FALSE)
  if(base::length(w1)==0){
    message('Fail in calculating activity, please check the ID type in cal_mat and target_list and try again !')
  }
  ac.mat <- ac.mat[w1,,drop=FALSE]
  return(ac.mat)
}

## inner functions
cal.Activity.old <- function(target_list=NULL, cal_mat=NULL, es.method = 'weightedmean',std=TRUE) {
  ## mean, absmean, maxmean, weightedmean
  use_genes <- row.names(cal_mat)
  if(base::length(use_genes)==0){
    message('No genes in the cal_mat, please check and re-try!');return(FALSE);
  }
  all_target <- target_list
  #all_target <- all_target[base::intersect(use_genes, names(all_target))] ## if the driver is not included in cal_mat but its target genes are included, will also calculate activity
  ac.mat <-
    matrix(NA, ncol = ncol(cal_mat), nrow = base::length(all_target)) ## generate activity matrix, each col for sample, each row for source target
  #z-normalize each sample
  if(std==TRUE) cal_mat <- apply(cal_mat, 2, do.std)
  for (i in 1:base::length(all_target)) {
    x <- names(all_target)[i]
    x1 <- all_target[[x]]
    x2 <- base::unique(base::intersect(rownames(x1), use_genes)) ## filter target by cal genes
    x1 <- x1[x2, ] ## target info
    target_num <- base::length(x2)
    if (target_num == 0)
      next
    if (target_num == 1){
      if (es.method != 'weightedmean') ac.mat[i, ] <- cal_mat[x2,] # 20230228
      if (es.method == 'weightedmean') ac.mat[i, ] <- cal_mat[x2,]*x1$MI * sign(x1$spearman) # 20230228
      next
    }
    if (es.method != 'weightedmean')
      ac.mat[i, ] <- apply(cal_mat[x2,,drop=FALSE], 2, es, es.method)
    if (es.method == 'weightedmean') {
      weight <- x1$MI * sign(x1$spearman)
      ac.mat[i, ] <- apply(cal_mat[x2,,drop=FALSE] * weight, 2, es, 'mean')
    }
  }
  rownames(ac.mat) <- names(all_target)
  colnames(ac.mat) <- colnames(cal_mat)
  w1 <- which(is.na(ac.mat[,1])==FALSE)
  if(base::length(w1)==0){
    message('Fail in calculating activity, please check the ID type in cal_mat and target_list and try again !')
  }
  ac.mat <- ac.mat[w1,]
  return(ac.mat)
}
get_target_list2matrix <- function(target_list=NULL,es.method = 'weightedmean') {
  all_source <- names(target_list)
  all_target <- base::unique(unlist(lapply(target_list,function(x)rownames(x))))
  mat1 <- matrix(0,nrow=base::length(all_source),ncol=base::length(all_target))
  rownames(mat1) <- all_source;
  colnames(mat1) <- all_target;
  for(i in all_source){
    if(es.method!='weightedmean') mat1[i,rownames(target_list[[i]])] <- 1;
    if(es.method=='weightedmean') mat1[i,rownames(target_list[[i]])] <- target_list[[i]]$MI*sign(target_list[[i]]$spearman);
  }
  return(mat1)
}
get_igraph2matrix <- function(gr=NULL,es.method = 'weightedmean'){
  if(es.method=='weightedmean'){
    if('weight' %in% igraph::edge_attr_names(gr) & 'sign' %in% igraph::edge_attr_names(gr)){
      mat1 <- igraph::as_adjacency_matrix(gr,type='both',attr = 'weight')
      mat2 <- igraph::as_adjacency_matrix(gr,type='both',attr = 'sign')
      mat1 <- mat1*mat2
    }else{
      message('weight, sign attributes are not included in the input igraph object, weightedmean can not be used !')
      return(FALSE)
    }
  }else{
    mat1 <- igraph::as_adjacency_matrix(gr,type='both')
  }
  return(mat1)
}
get_gr2driver <- function(gr,mode='out'){
  d1 <- igraph::degree(gr,mode=mode)
  names(d1[which(d1>0)])
}

#' Clean Activity-based profile
#'
#' \code{processDriverProfile} is a helper function to pre-process Activity-based profile.
#'
#' @param Driver_profile a numeric vector, contain statistics for drivers
#' (e.g driver's target size, driver's Z-statistics)
#' @param Driver_name a character vector, contain name for drivers.
#' The length of `Driver_profile` and `Driver_name` must be equal and the order of item must match.
#' @param choose_strategy character, strategy of selection if duplicate driver name (e.g TP53_TF, TP53_SIG).
#' Choose from "min","max","absmin","absmax". Default is "min".
#' @param return_type character, strategy of return type.
#' Choose from "driver_name","gene_name", "driver_statistics", "gene_statistics".
#' If choose "*_name", only the name vector is returned.
#' If choose "*_statistics", the statistics vector is returned with character name.
#' Default is "driver_name".
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' Driver_profile <- ms_tab$P.Value.G4.Vs.WNT_DA
#' Driver_name <- ms_tab$gene_label
#' res1 <- processDriverProfile(Driver_profile=Driver_profile,
#'                               Driver_name=Driver_name,
#'                               choose_strategy='min')
#' res2 <- processDriverProfile(Driver_profile=Driver_profile,
#'                               Driver_name=Driver_name,
#'                               return_type = 'gene_name',
#'                               choose_strategy='min')
#' Driver_profile <- ms_tab$Z.G4.Vs.WNT_DA
#' res3 <- processDriverProfile(Driver_profile=Driver_profile,
#'                               Driver_name=Driver_name,
#'                               choose_strategy='absmax',
#'                               return_type = 'driver_statistics')
#' res4 <- processDriverProfile(Driver_profile=Driver_profile,
#'                               Driver_name=Driver_name,
#'                               choose_strategy='absmax',
#'                               return_type = 'gene_statistics')
#' driver_size <- ms_tab$Size
#' res5 <- processDriverProfile(Driver_profile=driver_size,
#'                               Driver_name=Driver_name,
#'                               choose_strategy='max')
#' @export
processDriverProfile <- function(Driver_profile,Driver_name,
                                  choose_strategy='min',
                                  return_type='driver_name'){
  all_input_para <- c('Driver_profile','Driver_name','choose_strategy')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('choose_strategy',c('min','max','absmin','absmax'),envir=environment()),
                 check_option('return_type',c('driver_name','gene_name','driver_statistics','gene_statistics'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}

  ##
  l1 <- length(Driver_profile)
  l2 <- length(Driver_name)
  if(l1!=l2 | l1==0 | l2==0){
    message('Driver_profile and Driver_name need same length.
            Please check Driver_profile, Driver_name, and re-try!');return(FALSE)
  }
  gene_name <- gsub('(.*)_(TF|SIG)',"\\1",Driver_name)
  uni_gene_name <- unique(gene_name)
  tmp1 <- lapply(uni_gene_name,function(x){
    w1 <- which(gene_name == x)
    x1 <- Driver_profile[w1]
    x2 <- rank(x1)
    x3 <- rank(abs(x1))
    if(choose_strategy=='min'){w2 <- w1[which.min(x2)]}
    if(choose_strategy=='max'){w2 <- w1[which.max(x2)]}
    if(choose_strategy=='absmin'){w2 <- w1[which.min(x3)]}
    if(choose_strategy=='absmax'){w2 <- w1[which.max(x3)]}
    w2
  })
  remain_item <- unlist(tmp1)
  if(return_type=='driver_name') new_vec <- Driver_name[remain_item]
  if(return_type=='gene_name') new_vec <- gene_name[remain_item]
  if(return_type=='driver_statistics'){new_vec <- Driver_profile[remain_item]; names(new_vec) <- Driver_name[remain_item]}
  if(return_type=='gene_statistics'){new_vec <-  Driver_profile[remain_item]; names(new_vec) <- gene_name[remain_item]}
  return(new_vec)
}

#' Calculate Activity Value for Gene Sets
#'
#' \code{cal.Activity.GS} calculates activity value for each gene set, and return a numeric matrix with rows of gene sets and columns of samples.
#'
#' @param use_gs2gene list, contains elements of gene sets. Element name is gene set name, each element contains a vector of genes belong to that gene set.
#' Default is using \code{all_gs2gene[c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG')]}, which is loaded from \code{gs.preload}.
#' @param cal_mat numeric matrix, gene/transcript expression matrix.
#' If want to input activity matrix, need to use `processDriverProfile()` to pre-process the dataset.
#' Detailed could see demo.
#' @param es.method character, method to calculate the activity value. Users can choose from "mean", "absmean", "maxmean", "gsva", "ssgsea", "zscore" and "plage".
#' The details for using the last four options, users can check \code{gsva}. Default is "mean".
#' @param std logical, if TRUE, the expression matrix will be normalized by column. Default is TRUE.
#' @return Return an activity matrix with rows of gene sets and columns of samples.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,
#'                         use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
#' exp_mat_gene <- Biobase::exprs(analysis.par$cal.eset)
#' ## each row is a gene symbol, if not, must convert ID first
#' ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,
#'                         cal_mat = exp_mat_gene)
#' ## if want to input activity-matrix
#' ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,
#'           cal_mat=Biobase::exprs(analysis.par$cal.eset),
#'           es.method='weightedmean')
#' # pre-process the activity matrix by selecting the one
#' # with larger target size for duplicate drivers
#' Driver_name <- rownames(ac_mat)
#' ms_tab <- analysis.par$final_ms_tab
#' driver_size <- ms_tab[Driver_name,]$Size
#' use_driver <- processDriverProfile(Driver_profile=driver_size,
#'                               Driver_name=Driver_name,
#'                               choose_strategy='max',
#'                               return_type='driver_name')
#' use_driver_gene_name <-
#'               processDriverProfile(Driver_profile=driver_size,
#'                               Driver_name=Driver_name,
#'                               choose_strategy='max',
#'                               return_type='gene_name')
#' ac_mat_gene <- ac_mat[use_driver,]
#' rownames(ac_mat_gene) <- use_driver_gene_name
#' driver_ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,
#'                         cal_mat = ac_mat_gene)
#' @export
cal.Activity.GS <- function(use_gs2gene=all_gs2gene[c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG')], cal_mat=NULL, es.method = 'mean',std=TRUE) {
  #
  all_input_para <- c('use_gs2gene','cal_mat','es.method','std')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('std',c(TRUE,FALSE),envir=environment()),
                 check_option('es.method',c("mean","absmean","maxmean","gsva","ssgsea","zscore","plage"),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  while(class(use_gs2gene[[1]])=='list'){
    nn <- unlist(lapply(use_gs2gene,names))
    use_gs2gene <- unlist(use_gs2gene,recursive = FALSE)
    names(use_gs2gene)<-nn
  }
  use_genes <- row.names(cal_mat)
  if(es.method %in% c('gsva','ssgsea','zscore','plage')){
    ac.mat <- GSVA::gsva(generate.eset(cal_mat),use_gs2gene,method=es.method)
    ac.mat <- Biobase::exprs(ac.mat)
    return(ac.mat)
  }
  ac.mat <-
    matrix(NA, ncol = ncol(cal_mat), nrow = base::length(use_gs2gene)) ## generate activity matrix, each col for sample, each row for source target
  #z-normalize each sample
  if(std==TRUE) cal_mat <- apply(cal_mat, 2, do.std)
  for (i in 1:base::length(use_gs2gene)) {
    x <- names(use_gs2gene)[i]
    x1 <- use_gs2gene[[x]]
    x2 <- base::unique(base::intersect(x1, use_genes)) ## filter target by cal genes
    target_num <- base::length(x2)
    if (target_num == 0)
      next
    if (target_num == 1){
      ac.mat[i, ] <- cal_mat[x2,]
      next
    }
    ac.mat[i, ] <- apply(cal_mat[x2,,drop=FALSE], 2, es, es.method)
  }
  rownames(ac.mat) <- names(use_gs2gene)
  colnames(ac.mat) <- colnames(cal_mat)
  w1 <- apply(ac.mat,1,function(x)base::length(which(is.na(x)==TRUE)))
  ac.mat <- ac.mat[which(w1==0),,drop=FALSE]
  return(ac.mat)
}

#' Differential Expression Analysis and Differential Activity Analysis Between 2 Sample Groups Using Bayesian Inference
#'
#' \code{getDE.BID.2G} is a function performs differential gene expression analysis and differential driver activity analysis between
#' control group (parameter G0) and experimental group (parameter G1).
#'
#' @param eset ExpressionSet class object, contains gene expression data or driver activity data.
#' @param output_id_column character, the column names of Biobase::fData(eset).
#' This option is useful when the \code{eset} expression matrix is at transcript-level, and user is expecting to see the gene-level comparison statistics.
#' If NULL, rownames of the Biobase::fData(eset) will be used.
#' Default is NULL.
#' @param G1 a vector of characters, the sample names of experimental group.
#' @param G0 a vecotr of characters, the sample names of control group.
#' @param G1_name character, the name of experimental group (e.g. "Male"). Default is "G1".
#' @param G0_name character, the name of control group (e.g. "Female"). Default is "G0".
#' @param logTransformed logical, if TRUE, log tranformation of the expression value will be performed.
#' @param method character, users can choose between "MLE" and "Bayesian".
#' "MLE", the maximum likelihood estimation, will call generalized linear model(glm/glmer) to perform data regression.
#' "Bayesian", will call Bayesian generalized linear model (bayesglm) or multivariate generalized linear mixed model (MCMCglmm) to perform data regression.
#' Default is "Bayesian".
#' @param family character or family function or the result of a call to a family function.
#' This parameter is used to define the model's error distribution. See \code{?family} for details.
#' Currently, options are gaussian, poisson, binomial(for two-group sample classes)/category(for multi-group sample classes)/ordinal(for multi-group sample classes with class_ordered=TRUE).
#' If set with gaussian or poission, the response variable in the regression model will be the expression level, and the independent variable will be the sample's phenotype.
#' If set with binomial, the response variable in the regression model will be the sample phenotype, and the independent variable will be the expression level.
#' For binomial, category and ordinal input, the family will be automatically reset, based on the sample's class level and the setting of \code{class_ordered}.
#' Default is gaussian.
#' @param pooling character, users can choose from "full","no" and "partial".
#' "full", use probes as independent observations.
#' "no", use probes as independent variables in the regression model.
#' "partial", use probes as random effect in the regression model.
#' Default is "full".
#' @param verbose logical, if TRUE, sample names of both groups will be printed. Default is TRUE.
#'
#' @return
#' Return a data frame. Rows are genes/drivers, columns are "ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "Z-statistics", "Ave.G1" and "Ave.G0".
#' Names of the columns may vary from different group names. Sorted by P.Value.
#'
#' @examples
#' mat <- matrix(c(0.50099,-1.2108,-1.0524,
#'                 0.34881,-0.13441,0.87112,
#'                 1.84579,-2.0356,-2.6025,
#'                 1.62954,1.88281,1.29604),nrow=2,byrow=TRUE)
#' rownames(mat) <- c('A1','A2')
#' colnames(mat) <-  c('Case-rep1','Case-rep2','Case-rep3',
#'                 'Control-rep1','Control-rep2','Control-rep3')
#' tmp_eset <- generate.eset(mat,feature_info=data.frame(row.names=rownames(mat),
#'             probe=rownames(mat),gene=rep('GeneX',2),
#'             stringsAsFactors = FALSE))
#' res1 <- getDE.BID.2G(tmp_eset,output_id_column='probe',
#'         G1=c('Case-rep1','Case-rep2','Case-rep3'),
#'         G0=c('Control-rep1','Control-rep2','Control-rep3'))
#' res2 <- getDE.BID.2G(tmp_eset,output_id_column='gene',
#'         G1=c('Case-rep1','Case-rep2','Case-rep3'),
#'         G0=c('Control-rep1','Control-rep2','Control-rep3'))
#' res3 <- getDE.BID.2G(tmp_eset,output_id_column='gene',
#'         G1=c('Case-rep1','Case-rep2','Case-rep3'),
#'         G0=c('Control-rep1','Control-rep2','Control-rep3'),
#'         pooling='partial')
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' each_subtype <- 'G4'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_BID <- getDE.BID.2G(eset=analysis.par$cal.eset,
#'                                 G1=G1,G0=G0,
#'                                 G1_name=each_subtype,
#'                                 G0_name='other')
#' DA_driver_BID <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,
#'                                 G1=G1,G0=G0,
#'                                 G1_name=each_subtype,
#'                                 G0_name='other')
#' }
#' @export
getDE.BID.2G <-function(eset,output_id_column=NULL,G1=NULL, G0=NULL,G1_name=NULL,G0_name=NULL,method='Bayesian',family=gaussian,pooling='full',logTransformed=TRUE,verbose=TRUE){
  #
  all_input_para <- c('eset','G1','G0','method','family','pooling','logTransformed','verbose')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('logTransformed',c(TRUE,FALSE),envir=environment()),
                 check_option('verbose',c(TRUE,FALSE),envir=environment()),
                 check_option('method',c('Bayesian','MLE'),envir=environment()),
                 check_option('pooling',c('full','no','partial'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  exp_mat <- Biobase::exprs(eset)
  G1 <- base::intersect(G1,colnames(exp_mat))
  G0 <- base::intersect(G0,colnames(exp_mat))
  if(verbose==TRUE){
    print(sprintf('G1:%s', base::paste(G1, collapse = ';')))
    print(sprintf('G0:%s', base::paste(G0, collapse = ';')))
  }
  if(base::length(G1)==0 | base::length(G0)==0){
    message('Too few samples, please check the sample name of G1, G0 and samples in eset !');return(FALSE);
  }
  #
  rn <- rownames(exp_mat)
  exp_mat <- exp_mat[,c(G1,G0)]
  if(is.null(dim(exp_mat))==TRUE){exp_mat <- t(as.matrix(exp_mat));rownames(exp_mat)<-rn}
  if(is.null(output_id_column)==TRUE) use_id <- rownames(Biobase::fData(eset)) else use_id <- Biobase::fData(eset)[,output_id_column]
  comp <- c(rep(1,length.out=base::length(G1)),rep(0,length.out=base::length(G0)))
  all_id <- base::unique(use_id)
  de <- lapply(all_id,function(x){
    w1 <- which(use_id==x)
    x1 <- exp_mat[w1,,drop=FALSE]
    bid(mat=x1,use_obs_class=comp,class_order=c(0,1),family=family,method=method,
        nitt=13000,burnin=3000,thin=1,pooling=pooling,class_ordered=FALSE,
        logTransformed=logTransformed,std=FALSE,average.method=c('geometric'),verbose=FALSE)
  })
  de <- as.data.frame(do.call(base::rbind,de))
  de$adj.P.Val<-p.adjust(de$P.Value,'fdr')
  de$logFC<-sign(de$FC)*log2(abs(de$FC))
  de$ID <- all_id
  de<-de[,c('ID','logFC','AveExpr','t','P.Value','adj.P.Val','Z-statistics')]
  rownames(de) <- de$ID
  tT <- de;
  tmp1 <- stats::aggregate(exp_mat,list(use_id),mean)
  new_mat <- tmp1[,-1];rownames(new_mat) <- tmp1[,1]
  tT <- tT[rownames(new_mat),,drop=FALSE]
  exp_G1 <- base::rowMeans(new_mat[,G1,drop=FALSE]);
  exp_G0 <- base::rowMeans(new_mat[,G0,drop=FALSE]);
  tT <- base::cbind(tT,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  tT <- tT[order(tT$P.Value),]
  return(tT)
}

#' Combine Multiple Comparison Results from Differential Expression (DE) or Differential Activity (DA) Analysis
#'
#' \code{combineDE} combines multiple comparisons of DE or DA analysis.
#' Can combine DE with DE, DA with DA and also DE with DA if proper transfer table prepared.
#'
#' For example, there are 4 subgroups in the phenotype, G1, G2, G3 and G4. One DE analysis was performed on G1 vs. G2, and another DE was performed on G1 vs. G3.
#' If user is interested in the DE analysis between G1 vs. (G2 and G3), he can call this function to combine the two comparison results above toghether.
#' The combined P values will be taken care by \code{combinePvalVector}.
#'
#' @param DE_list list, each element in the list is one DE/DA comparison need to be combined.
#' @param DE_name a vector of characters, the DE/DA comparison names.
#' If not NULL, it must match the names of DE_list in correct order.
#' If NULL, names of the DE_list will be used.
#' Default is NULL.
#' @param transfer_tab data.frame, the ID conversion table. Users can call \code{get_IDtransfer} to get this table.
#' The purpose is to correctly mapping ID for \code{DE_list}. The column names must match \code{DE_name}.
#' If NULL, ID column of each DE comparison will be considered as the same type.
#' Default is NULL.
#' @param main_id character, a name of the element in \code{DE_list}. The ID column of that comparison will be used as the ID of the final combination.
#' If NULL, the first element name from \code{DE_list} will be used. Default is NULL.
#' @param method character, users can choose between "Stouffer" and "Fisher". Default is "Stouffer".
#' @param twosided logical, if TRUE, a two-tailed test will be performed.
#' If FALSE, a one-tailed test will be performed, and P value falls within the range of 0 to 0.5. Default is TRUE.
#' @param signed logical, if TRUE, give a sign to the P value, which indicating the direction of testing.
#' Default is TRUE.
#' @return Return a list contains the combined DE/DA analysis. Each single comparison result before combination is wrapped inside
#' (may have with some IDs filtered out, due to the combination). A data frame named "combine" inside the list is the combined analysis.
#' Rows are genes/drivers, columns are combined statistics (e.g. "logFC", "AveExpr", "t", "P.Value" etc.).
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' each_subtype <- 'G4'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,
#'                                G1=G1,G0=G0,
#'                                G1_name=each_subtype,
#'                                G0_name='other')
#' DA_driver_limma <- getDE.limma.2G(eset=analysis.par$merge.ac.eset,
#'                                  G1=G1,G0=G0,
#'                                  G1_name=each_subtype,
#'                                  G0_name='other')
#' DE_list <- list(DE=DE_gene_limma,DA=DA_driver_limma)
#' g1 <- gsub('(.*)_.*','\\1',DE_list$DA$ID)
#' transfer_tab <- data.frame(DE=g1,DA=DE_list$DA$ID,stringsAsFactors = FALSE)
#' res1 <- combineDE(DE_list,transfer_tab=transfer_tab,main_id='DA')
#'
#' \dontrun{
#' each_subtype <- 'G4'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_limma_G4 <- getDE.limma.2G(eset=analysis.par$cal.eset,
#'                                    G1=G1,G0=G0,
#'                                    G1_name=each_subtype,
#'                                    G0_name='other')
#' each_subtype <- 'SHH'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_limma_SHH <- getDE.limma.2G(eset=analysis.par$cal.eset,
#'                                     G1=G1,G0=G0,
#'                                     G1_name=each_subtype,
#'                                     G0_name='other')
#' DE_list <- list(G4=DE_gene_limma_G4,SHH=DE_gene_limma_SHH)
#' res2 <- combineDE(DE_list,transfer_tab=NULL)
#' }
#' @export
combineDE<-function(DE_list,DE_name=NULL,transfer_tab=NULL,main_id=NULL,method='Stouffer',twosided=TRUE,signed=TRUE){
  #
  all_input_para <- c('DE_list','method','twosided','signed')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('twosided',c(TRUE,FALSE),envir=environment()),
                 check_option('signed',c(TRUE,FALSE),envir=environment()),
                 check_option('method',c('Stouffer','Fisher'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  nDE<-base::length(DE_list)
  if(nDE<2) stop('At least two DE outputs are required for combineDE analysis!\n')
  if(is.null(DE_name)==TRUE){
    DE_name <- names(DE_list)
  }
  if(is.null(transfer_tab)==TRUE){
    transfer_tab <- lapply(DE_list,function(x)x$ID)
    names(transfer_tab) <- DE_name
    w1 <- names(which(base::table(unlist(transfer_tab))==nDE))
    if(base::length(w1)==0){message('No intersected IDs found between DE list, please check and re-try!');return(FALSE);}
    transfer_tab <- do.call(base::cbind,lapply(transfer_tab,function(x)w1))
    transfer_tab <- as.data.frame(transfer_tab,stringsAsFactors=FALSE)
  }
  w1 <- lapply(DE_name,function(x){
    which(transfer_tab[[x]] %in% DE_list[[x]]$ID)
  })
  w11 <- which(base::table(unlist(w1))==nDE)
  if(base::length(w11)==0){
    message('No intersected IDs found between DE list, please check and re-try!');return(FALSE);
  }
  w2 <- as.numeric(names(w11))
  transfer_tab <- transfer_tab[w2,]
  combine_info <- lapply(DE_name,function(x1){
    DE_list[[x1]][transfer_tab[[x1]],]
  })
  names(combine_info) <- DE_name
  dd <- do.call(base::cbind,lapply(combine_info,function(x1){
    x1$P.Value*sign(x1$logFC)
  }))
  res1 <- t(apply(dd,1,function(x){
    combinePvalVector(x,method=method,signed=signed,twosided=twosided)
  }))
  res1 <- as.data.frame(res1)
  if(is.null(main_id)==TRUE) main_id <- DE_name[1]
  res1$adj.P.Val <- p.adjust(res1$P.Value,'fdr')
  res1$logFC <- base::rowMeans(do.call(base::cbind,lapply(combine_info,function(x)x$logFC)))
  res1$AveExpr <- base::rowMeans(do.call(base::cbind,lapply(combine_info,function(x)x$AveExpr)))
  res1 <- base::cbind(ID=transfer_tab[,main_id],transfer_tab,res1,stringsAsFactors=FALSE)
  rownames(res1) <- res1$ID
  combine_info$combine <- res1
  return(combine_info)
}


# inner function: class_label can be obtained by get_class
get_class2design <- function(class_label){
  design <- model.matrix(~0+class_label);colnames(design) <- base::unique(class_label);
  rownames(design) <- names(class_label)
  return(design)
  #design.mat <-as.data.frame(matrix(0, nrow = base::length(class_label), ncol = base::length(base::unique(class_label))))
  #rownames(design.mat) <- names(class_label) ## sample
  #colnames(design.mat) <- base::unique(class_label)
  #for(i in 1:base::length(class_label)){design.mat[names(class_label)[i],class_label[i]]<-1;}
  #return(design.mat)
}

#' Differential Expression Analysis and Differential Activity Analysis Between 2 Sample Groups Using Limma
#'
#' \code{getDE.limma.2G} is a function performs differential gene expression analysis and differential driver activity analysis
#' between control group (parameter G0) and experimental group (parameter G1), using limma related functions.
#'
#' @param eset ExpressionSet class object, contains gene expression data or driver activity data.
#' @param G1 a vector of characters, the sample names of experimental group.
#' @param G0 a vecotr of characters, the sample names of control group.
#' @param G1_name character, the name of experimental group (e.g. "Male"). Default is "G1".
#' @param G0_name character, the name of control group (e.g. "Female"). Default is "G0".
#' @param verbose logical, if TRUE, sample names of both groups will be printed. Default is TRUE.
#' @param random_effect a vector of characters, vector or factor specifying a blocking variable.
#' Default is NULL, no random effect will be considered.
#'
#' @return
#' Return a data frame. Rows are genes/drivers, columns are "ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Z-statistics", "Ave.G1" and "Ave.G0".
#' Names of the columns may vary from different group names. Sorted by P-values.
#'
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' each_subtype <- 'G4'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,
#'                                 G1=G1,G0=G0,
#'                                 G1_name=each_subtype,
#'                                 G0_name='other')
#' DA_driver_limma <- getDE.limma.2G(eset=analysis.par$merge.ac.eset,
#'                                 G1=G1,G0=G0,
#'                                 G1_name=each_subtype,
#'                                 G0_name='other')
#' @export
getDE.limma.2G <- function(eset=NULL, G1=NULL, G0=NULL,G1_name=NULL,G0_name=NULL,verbose=TRUE,random_effect=NULL) {
  #
  all_input_para <- c('eset','G1','G0','verbose')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('verbose',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  exp_mat <- Biobase::exprs(eset)
  G1 <- base::intersect(G1,colnames(exp_mat))
  G0 <- base::intersect(G0,colnames(exp_mat))
  if(verbose==TRUE){
    print(sprintf('G1:%s', base::paste(G1, collapse = ';')))
    print(sprintf('G0:%s', base::paste(G0, collapse = ';')))
  }
  if(base::length(G1)==0 | base::length(G0)==0){
    message('Too few samples, please check the sample name of G1, G0 and samples in eset !');return(FALSE);
  }
  #
  all_samples <- colnames(Biobase::exprs(eset))
  use_samples <- c(G0, G1)
  phe <- as.data.frame(Biobase::pData(eset)[use_samples, ,drop=FALSE]);
  rownames(phe) <- use_samples
  new_eset <- generate.eset(exp_mat=Biobase::exprs(eset)[, use_samples,drop=F],phenotype_info=phe, feature_info=Biobase::fData(eset))
  new_mat  <- Biobase::exprs(new_eset)
  ##
  design.mat <-as.data.frame(matrix(NA, nrow = base::length(use_samples), ncol = 1))
  rownames(design.mat) <-use_samples
  colnames(design.mat) <- 'group'
  design.mat[base::intersect(G0, use_samples), 'group'] <- 'G0'
  design.mat[base::intersect(G1, use_samples), 'group'] <- 'G1'
  #  design <- model.matrix( ~ group + 0, design.mat)
  group <- factor(design.mat$group)
  design <- model.matrix(~0+group);
  colnames(design) <- levels(group); rownames(design) <- colnames(new_mat)

  if(is.null(random_effect)==TRUE){
    fit <- limma::lmFit(new_mat,design)
  }else{
    random_effect <- random_effect[colnames(Biobase::exprs(new_eset))]
    corfit <- limma::duplicateCorrelation(new_eset,design,block=random_effect)
    fit <- limma::lmFit(new_mat,design,block=random_effect,correlation=corfit$consensus)
  }
  contrasts <- limma::makeContrasts(G1-G0,levels=design)
  fit2 <- limma::contrasts.fit(fit,contrasts=contrasts)
  fit2 <- limma::eBayes(fit2,trend=TRUE)
  #summary(decideTests(fit2, method="global"))
  ##
  tT <- limma::topTable(fit2, adjust.method = "fdr", number = Inf,coef=1)
  if(nrow(tT)==1){
    rownames(tT) <- rownames(new_mat)
  }
  tT <- base::cbind(ID=rownames(tT),tT,stringsAsFactors=FALSE)
  tT <- tT[rownames(new_mat),,drop=FALSE]
  exp_G1 <- base::rowMeans(new_mat[,G1,drop=FALSE]);
  exp_G0 <- base::rowMeans(new_mat[,G0,drop=FALSE]);
  w1 <- which(tT$P.Value<=0);
  if(base::length(w1)>0) tT$P.Value[w1] <- .Machine$double.xmin;
  #z_val <- sapply(tT$P.Value*sign(tT$logFC),function(x)combinePvalVector(x,twosided = TRUE)[1])
  z_val <- sapply(tT$P.Value*sign(tT$logFC),function(x)ifelse(x ==0, 0, combinePvalVector(x,twosided = TRUE)[1])) ## remove zero
  if(is.null(random_effect)==TRUE){
    tT <- base::cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  }else{
    tT <- base::cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1,
                'Ave.G0_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G0'],
                'Ave.G1_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G1'])
  }
  if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  tT <- tT[order(tT$P.Value, decreasing = FALSE), ]
  return(tT)
}

#' Combine P Values Using Fisher's Method or Stouffer's Method
#'
#' \code{combinePvalVector} is a function to combine multiple comparison's P values using Fisher's method or Stouffer's method.
#'
#' @param pvals a vector of numerics, the P values from multiple comparison need to be combined.
#' @param method character, users can choose between "Stouffer" and "Fisher". Default is "Stouffer".
#' @param signed logical, if TRUE, will give a sign to the P value to indicate the direction of testing.
#' Default is TRUE.
#' @param twosided logical, if TRUE, P value is calculated in a one-tailed test.
#' If FALSE, P value is calculated in a two-tailed test, and it falls within the range 0 to 0.5.
#' Default is TRUE.
#' @return Return a vector contains the "Z-statistics" and "P.Value".
#' @examples
#' combinePvalVector(c(0.1,1e-3,1e-5))
#' combinePvalVector(c(0.1,1e-3,-1e-5))
#' @export
combinePvalVector <-
  function(pvals,
           method = 'Stouffer',
           signed = TRUE,
           twosided = TRUE) {
    #
    all_input_para <- c('pvals','method','signed','twosided')
    check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
    if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
    check_res <- c(check_option('signed',c(TRUE,FALSE),envir=environment()),
                   check_option('twosided',c(TRUE,FALSE),envir=environment()),
                   check_option('method',c('Stouffer','Fisher'),envir=environment()))
    if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
    #
    #remove NA pvalues
    pvals <- pvals[!is.na(pvals) & !is.null(pvals)]
    pvals[which(abs(pvals)<=0)] <- .Machine$double.xmin
    if (sum(is.na(pvals)) >= 1) {
      stat <- NA
      pval <- NA
    } else{
      if (twosided & (sum(pvals > 1 | pvals < -1) >= 1))
        stop('pvalues must between 0 and 1!\n')
      if (!twosided & (sum(pvals > 0.5 | pvals < -0.5) >= 1))
        stop('One-sided pvalues must between 0 and 0.5!\n')

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
              df = 2 * base::length(pvals),
              lower.tail = FALSE
            ) / 2,
            pchisq(
              -2 * sum(log(as.numeric(neg.pvals))),
              df = 2 * base::length(pvals),
              lower.tail = FALSE
            ) / 2
          )
        pval <- base::min(abs(pvals))[1]
        #if two pvals are equal, pick up the first one
        stat <-
          sign(pvals[abs(pvals) == pval])[1] * qnorm(pval, lower.tail = F)[1]
        pval <- 2 * pval
      }
      else if (grepl('Stou', method, ignore.case = TRUE)) {
        if (twosided) {
          zs <- signs * qnorm(abs(pvals) / 2, lower.tail = FALSE)
          stat <- sum(zs) / sqrt(base::length(zs))
          pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
        }
        else{
          zs <- signs * qnorm(abs(pvals), lower.tail = FALSE)
          stat <- sum(zs) / sqrt(base::length(zs))
          pval <- pnorm(abs(stat), lower.tail = FALSE)
        }
      }
      else{
        stop('Only \"Fisher\" or \"Stouffer\" method is supported!!!\n')
      }
    }
    return(c(`Z-statistics` = stat, `P.Value` = pval))
  }


#' Merge Activity Values from TF (transcription factors) ExpressionSet Object and Sig (signaling factors) ExpressionSet Object
#'
#' \code{merge_TF_SIG.AC} combines the activity value from TF (transcription factors) and Sig (signaling factors) ExpressionSet objects together,
#' and adds "_TF" or "_SIG" suffix to drivers for easier distinction.
#'
#'
#' @param TF_AC ExpressionSet object, containing the activity values for all TFs.
#' @param SIG_AC ExpressionSet object, containing the activity values for all SIGs.
#'
#' @return Return an ExpressionSet object.
#' @examples
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             project_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
#' ## get eset (here for demo, use network.par$net.eset)
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' analysis.par$cal.eset <- network.par$net.eset
#' ac_mat_TF <- cal.Activity(target_list=analysis.par$tf.network$target_list,
#'                        cal_mat=Biobase::exprs(analysis.par$cal.eset),
#'                        es.method='weightedmean')
#' ac_mat_SIG <- cal.Activity(target_list=analysis.par$tf.network$target_list,
#'                        cal_mat=Biobase::exprs(analysis.par$cal.eset),
#'                        es.method='weightedmean')
#' analysis.par$ac.tf.eset  <- generate.eset(exp_mat=ac_mat_TF,
#'                                           phenotype_info=Biobase::pData(analysis.par$cal.eset))
#' analysis.par$ac.sig.eset <- generate.eset(exp_mat=ac_mat_SIG,
#'                                           phenotype_info=Biobase::pData(analysis.par$cal.eset))
#' analysis.par$merge.ac.eset <- merge_TF_SIG.AC(TF_AC=analysis.par$ac.tf.eset,
#'                                           SIG_AC=analysis.par$ac.sig.eset)
#' @export
merge_TF_SIG.AC <- function(TF_AC=NULL,SIG_AC=NULL){
  #
  all_input_para <- c('TF_AC','SIG_AC')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  mat_TF <- Biobase::exprs(TF_AC)
  mat_SIG <- Biobase::exprs(SIG_AC)
  funcType <- c(rep('TF',nrow(mat_TF)),rep('SIG',nrow(mat_SIG)))
  rn <- c(rownames(mat_TF),rownames(mat_SIG))
  rn_label <- base::paste(rn,funcType,sep='_')
  mat_combine <- base::rbind(mat_TF,mat_SIG[,colnames(mat_TF)])
  rownames(mat_combine) <- rn_label
  eset_combine <- generate.eset(exp_mat=mat_combine,phenotype_info=Biobase::pData(TF_AC)[colnames(mat_combine),],
                                feature_info=NULL,annotation_info='activity in dataset')
  return(eset_combine)
}

#' Merge TF (transcription factor) Network and Sig (signaling factor) Network
#'
#' \code{merge_TF_SIG.network} takes TF network and Sig network and combine them together.
#' The merged list object contains three elements, a data.frame contains all the combined network information \code{network_dat},
#' a driver-to-target list object \code{target_list}, and an igraph object of the network \code{igraph_obj}.
#'
#' @param TF_network list, the TF network created by \code{get.SJAracne.network} function.
#' @param SIG_network list, the SIG network created by \code{get.SJAracne.network} function.
#' @return
#' Return the a list containing three elements, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#' @examples
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             project_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
#' analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,
#'                                                    SIG_network=analysis.par$sig.network)
#' @export
merge_TF_SIG.network <- function(TF_network=NULL,SIG_network=NULL){
  #
  all_input_para <- c('TF_network','SIG_network')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  s_TF  <- names(TF_network$target_list)
  s_SIG <- names(SIG_network$target_list)
  funcType <- c(rep('TF',base::length(s_TF)),rep('SIG',base::length(s_SIG)))
  rn <- c(s_TF,s_SIG)
  rn_label <- base::paste(rn,funcType,sep='_')
  target_list_combine <- c(TF_network$target_list,SIG_network$target_list)
  names(target_list_combine) <- rn_label
  n_TF <- TF_network$network_dat
  if(nrow(n_TF)>0) n_TF$source <- base::paste(n_TF$source,'TF',sep='_')
  n_SIG <- SIG_network$network_dat
  if(nrow(n_SIG)>0) n_SIG$source <- base::paste(n_SIG$source,'SIG',sep='_')
  net_dat <- base::rbind(n_TF,n_SIG)
  igraph_obj <- graph_from_data_frame(net_dat[,c('source','target')],directed=TRUE)
  if('MI' %in% colnames(net_dat)) igraph_obj <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  if('spearman' %in% colnames(net_dat)) igraph_obj <- set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
  return(list(network_dat=net_dat,target_list=target_list_combine,igraph_obj=igraph_obj))
}

#' Generate the Master Table for Drivers
#'
#' \code{generate.masterTable} generates a master table to show the mega information of all tested drivers.
#'
#' The master table gathers TF (transcription factor) information, Sig (signaling factor) information, all the DE (differential expression analysis)
#' and DA (differential activity analysis) from multiple comparisons. It also shows each driver's target gene size and other additional information
#' (e.g. gene biotype, chromosome name, position etc.).
#'
#' @param use_comp a vector of characters, the name of multiple comparisons. It will be used to name the columns of master table.
#' @param DE list, a list of DE comparisons, each comparison is a data.frame. The element name in the list must contain the name in \code{use_comp}.
#' @param DA list, a list of DA comparisons, each comparison is a data.frame. The element name in the list must contain the name in \code{use_comp}.
#' @param target_list list, a driver-to-target list. The names of the list elements are drivers. Each element is a data frame, usually contains three columns.
#' "target", target gene names; "MI", mutual information; "spearman", spearman correlation coefficient.
#' It is highly suggested to follow the NetBID2 pipeline, and the \code{TF_network} could be generated by \code{get_net2target_list} and \code{get.SJAracne.network}.
#' @param main_id_type character, the type of driver's ID. It comes from the attribute name in biomaRt package.
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' For details, user can call \code{biomaRt::listAttributes()} to display all available attributes in the selected dataset.
#' @param transfer_tab data.frame, the data frame for ID conversion. This can be obtained by calling \code{get_IDtransfer}.
#' If NULL and \code{main_id_type} is not in the column names of \code{tf_sigs}, it will use the conversion table within the function.
#' Default is NULL.
#' @param tf_sigs list, contains all the detailed information of TF and Sig. Users can call \code{db.preload} for access.
#' @param z_col character, name of the column in \code{DE} and \code{DA} contains the Z statistics. Default is "Z-statistics".
#' @param display_col character, name of the column in \code{DE} and \code{DA} need to be kept in the master table. Default is c("logFC","P.Value").
#' @param column_order_strategy character, users can choose between "type" and "comp". Default is "type".
#' If set as type, the columns will be ordered by column type; If set as comp, the columns will be ordered by comparison.
#' @return Return a data frame contains the mega information of all tested drivers.
#' The column "originalID" and "originalID_label" is the same ID as from the original dataset.
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' #analysis.par$final_ms_tab ## this is master table generated before
#' ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,
#'                        cal_mat=Biobase::exprs(analysis.par$cal.eset),es.method='weightedmean')
#' analysis.par$ac.merge.eset  <- generate.eset(exp_mat=ac_mat,
#'                                              phenotype_info=Biobase::pData(analysis.par$cal.eset))
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' all_subgroup <- base::unique(phe_info$subgroup) ##
#' for(each_subtype in all_subgroup){
#'  comp_name <- sprintf('%s.Vs.others',each_subtype) ## each comparison must give a name !!!
#'  G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#'  G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#'  DE_gene_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,
#'                                  G1_name=each_subtype,G0_name='other')
#'  analysis.par$DE[[comp_name]] <- DE_gene_limma
#'  DA_driver_limma <- getDE.limma.2G(eset=analysis.par$ac.merge.eset,G1=G1,G0=G0,
#'                                    G1_name=each_subtype,G0_name='other')
#'  analysis.par$DA[[comp_name]] <- DA_driver_limma
#' }
#' all_comp <- names(analysis.par$DE) ## get all comparison name for output
#' db.preload(use_level='gene',use_spe='human',update=FALSE);
#' test_ms_tab <- generate.masterTable(use_comp=all_comp,
#'                                            DE=analysis.par$DE,
#'                                            DA=analysis.par$DA,
#'                                            target_list=analysis.par$merge.network$target_list,
#'                                            tf_sigs=tf_sigs,
#'                                            z_col='Z-statistics',
#'                                            display_col=c('logFC','P.Value'),
#'                                            main_id_type='external_gene_name')
#' @export
generate.masterTable <- function(use_comp=NULL,DE=NULL,DA=NULL,
                                 target_list=NULL,main_id_type=NULL,transfer_tab=NULL,
                                 tf_sigs=NULL,
                                 z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                 column_order_strategy='type'){
  #
  all_input_para <- c('use_comp','DE','DA','target_list','main_id_type','tf_sigs','z_col','display_col','column_order_strategy')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=base::environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('column_order_strategy',c('type','comp'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(base::length(base::setdiff(use_comp,names(DE)))>0){message(sprintf('%s not included in DE, please check and re-try!',base::setdiff(use_comp,names(DE))));return(FALSE)}
  if(base::length(base::setdiff(use_comp,names(DA)))>0){message(sprintf('%s not included in DA, please check and re-try!',base::setdiff(use_comp,names(DA))));return(FALSE)}
  # get original ID
  ori_rn <- rownames(DA[[1]]) ## with label
  w1 <- grep('(.*)_TF',ori_rn); w2 <- grep('(.*)_SIG',ori_rn)
  funcType <- rep(NA,length.out=base::length(ori_rn));rn <- funcType
  funcType[w1] <- 'TF'; funcType[w2] <- 'SIG';
  rn[w1] <- gsub('(.*)_TF',"\\1",ori_rn[w1]);rn[w2] <- gsub('(.*)_SIG',"\\1",ori_rn[w2]);
  rn_label <- ori_rn
  use_size <- unlist(lapply(target_list[rownames(DA[[1]])],nrow))
  # id issue
  current_id <- names(tf_sigs$tf)[-1]
  use_info <- base::unique(base::rbind(tf_sigs$tf$info,tf_sigs$sig$info))
  if(main_id_type %in% current_id){
    use_info <- use_info[which(use_info[,main_id_type] %in% rn),]
  }else{
    if(is.null(transfer_tab)==TRUE){
      transfer_tab <- get_IDtransfer(from_type=main_id_type,to_type=current_id[1],use_genes=rn,ignore_version = TRUE)
      use_info <- base::merge(use_info,transfer_tab,by.x=current_id[1],by.y=current_id[1])
    }else{
      transfer_tab <- transfer_tab[which(transfer_tab[,main_id_type] %in% rn),]
      uid <- base::intersect(colnames(transfer_tab),colnames(use_info))
      if(base::length(uid)==0){message('No ID type in the transfer_tab could match ID type in tf_sigs, please check and re-try!');return(FALSE);}
      uid <- uid[1]
      use_info <- base::merge(use_info,transfer_tab,by.x=uid,by.y=uid)
    }
    use_info <- use_info[which(use_info[,main_id_type] %in% rn),]
  }
  use_info <- base::unique(use_info)
  if(nrow(use_info)==0){message('ID issue error, please check main_id_type setting!');return(FALSE);}
  # merge info
  tmp1 <- stats::aggregate(use_info,list(use_info[,main_id_type]),function(x){
    x1 <- x[which(x!="")]
    x1 <- x1[which(is.na(x1)==FALSE)]
    base::paste(sort(base::unique(x1)),collapse=';')
  })
  tmp1 <- tmp1[,-1]; rownames(tmp1) <- tmp1[,main_id_type]
  geneSymbol <- tmp1[rn,'external_gene_name'] ## this column for function enrichment
  if('external_transcript_name' %in% colnames(tmp1)){ ## this column for display
    gene_label <- base::paste(tmp1[rn,'external_transcript_name'],funcType,sep = '_')
  }else{
    gene_label <-base::paste(tmp1[rn,'external_gene_name'],funcType,sep = '_')
  }
  #
  #label_info <- data.frame('gene_label'=gene_label,'geneSymbol'=geneSymbol,
  #                         'originalID'=rn,'originalID_label'=rn_label,'funcType'=funcType,'Size'=use_size,stringsAsFactors=FALSE)
  label_info <- data.frame('originalID_label'=rn_label,'originalID'=rn,'gene_label'=gene_label,'geneSymbol'=geneSymbol,
                           'funcType'=funcType,'Size'=use_size,stringsAsFactors=FALSE)
  w1 <- which(is.na(geneSymbol)==TRUE)
  label_info[w1,'geneSymbol'] <- label_info[w1,'originalID']
  label_info[w1,'gene_label'] <- label_info[w1,'originalID_label']
  add_info <- tmp1[rn,]
  #
  combine_info <- lapply(use_comp,function(x){
    DA[[x]] <- DA[[x]][rn_label,,drop=F]
    DE[[x]] <- as.data.frame(DE[[x]])[rn,]
    avg_col <- colnames(DA[[x]])[grep('^Ave',colnames(DA[[x]]))]
    uc <- c(z_col,avg_col,base::setdiff(display_col,c(z_col,avg_col))); uc <- base::intersect(uc,colnames(DA[[x]]))
    DA_info <- DA[[x]][rn_label,uc,drop=F]
    avg_col <- colnames(DE[[x]])[grep('^Ave',colnames(DE[[x]]))]
    uc <- c(z_col,avg_col,base::setdiff(display_col,c(z_col,avg_col))); uc <- base::intersect(uc,colnames(DE[[x]]))
    DE_info <- as.data.frame(DE[[x]])[rn,uc,drop=F]
    colnames(DA_info) <- paste0(colnames(DA_info),'.',x,'_DA')
    colnames(DE_info) <- paste0(colnames(DE_info),'.',x,'_DE')
    colnames(DA_info)[1] <- paste0('Z.',x,'_DA')
    colnames(DE_info)[1] <- paste0('Z.',x,'_DE')
    out <- base::cbind(DA_info,DE_info,stringsAsFactors=FALSE)
    rownames(out) <- rn_label
    out
  })
  combine_info_DA <- do.call(base::cbind,lapply(combine_info,function(x)x[rn_label,grep('_DA$',colnames(x)),drop=T]))
  combine_info_DE <- do.call(base::cbind,lapply(combine_info,function(x)x[rn_label,grep('_DE$',colnames(x)),drop=T]))
  # re-organize the columns for combine info
  if(column_order_strategy=='type' & length(use_comp)>1){
    col_ord <- c('Z','AveExpr',display_col)
    tmp1 <- lapply(col_ord,function(x){
      x1 <- grep(sprintf('^%s\\.',x),colnames(combine_info_DA))
      if(length(x1)>0) combine_info_DA[,x1] else return(NULL)
    })
    combine_info_DA <- do.call(base::cbind,tmp1)
    tmp1 <- lapply(col_ord,function(x){
      combine_info_DE[,grep(sprintf('^%s\\.',x),colnames(combine_info_DE))]
    })
    combine_info_DE <- do.call(base::cbind,tmp1)
  }
  # put them together
  ms_tab <- base::cbind(label_info,combine_info_DA,combine_info_DE,add_info)
  rownames(ms_tab) <- ms_tab$originalID_label
  return(ms_tab)
}

#' Save the Master Table into Excel File
#'
#' \code{out2excel} is a function can save data frame as Excel File. This is mainly for the output of master table generated by \code{generate.masterTable}.
#'
#' @param all_ms_tab list or data.frame, if data.frame, it is generated by \code{generate.masterTable}.
#' If list, each list element is data.frame/master table.
#' The name of the list element will be the sheet name in the excel file.
#' @param out.xlsx character, path and file name of the output Excel file.
#' @param mark_gene list, list of marker genes. The name of the list element is the marked group name. Each element is a vector of marker genes.
#' This is optional, just to add additional information to the file.
#' @param mark_col character, the color to mark the marker genes. If NULL, will use \code{get.class.color} to get the colors.
#' @param mark_strategy character, users can choose between "color" and "add_column".
#' "Color" means the mark_gene will be marked by filling its background color;
#' "add_column" means the mark_gene will be displayed in a separate column with TRUE/FALSE, indicating whether the gene belongs to a mark group or not.
#' @param workbook_name character, name of the workbook for the output Excel. Default is "ms_tab".
#' @param only_z_sheet logical, if TRUE, will create a separate sheet only contains Z-statistics related columns from DA/DE analysis.
#' Default is FALSE.
#' @param z_column character, name of the columns contain Z-statistics. If NULL, find column names start with "Z.".
#' Default is NULL.
#' @param sig_thre numeric, threshold for the Z-statistics. Z values passed the threshold will be colored. The color scale is defined by \code{z2col}.
#' Default is 1.64.
#' @return Return a logical value. If TRUE, the Excel file has been generated successfully.
#' @examples
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab ## this is master table generated before
#' mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),
#'                  SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
#'                  Group3=c('IMPG2','GABRA5','EGFL11','NRL','MAB21L2','NPR3','MYC'),
#'                  Group4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
#' mark_col <- get.class.color(names(mark_gene),
#'             pre_define=c('WNT'='blue','SHH'='red',
#'             'Group3'='yellow','Group4'='green'))
#' outfile <- 'test_out.xlsx'
#' out2excel(ms_tab,out.xlsx = outfile,mark_gene,mark_col)
#' }
#' @export
out2excel <- function(all_ms_tab,out.xlsx,
                      mark_gene=NULL,
                      mark_col=NULL,
                      mark_strategy='color',
                      workbook_name='ms_tab',
                      only_z_sheet=FALSE,
                      z_column=NULL,sig_thre=1.64){
  #
  all_input_para <- c('all_ms_tab','out.xlsx','mark_strategy','workbook_name','only_z_sheet','sig_thre')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('only_z_sheet',c(TRUE,FALSE),envir=environment()),
                 check_option('mark_strategy',c('color','add_column'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  wb <- openxlsx::createWorkbook(workbook_name)
  if(!'list' %in% class(all_ms_tab)){
    all_ms_tab <- list('Sheet1'=as.data.frame(all_ms_tab))
  }
  if(only_z_sheet==TRUE){
    nn <- names(all_ms_tab)
    all_ms_tab <- lapply(all_ms_tab,function(x){
      w1 <- grep('.*_[DE|DA]',colnames(x))
      w2 <- base::setdiff(colnames(x)[w1],colnames(x)[w1][grep('^Z.*',colnames(x)[w1])])
      w3 <- base::setdiff(colnames(x),w2)
      list(x[,w3],x)
    })
    all_ms_tab <- unlist(all_ms_tab,recursive = FALSE)
    nn1 <- lapply(nn,function(x){
      if(x!='Sheet1'){
        c(sprintf('Only_Z_%s',x),sprintf('Full_info_%s',x))
      }else{
        c('Only_Z','Full_info')
      }
    })
    nn1 <- unlist(nn1)
    names(all_ms_tab) <- nn1
  }
  if(is.null(mark_gene)==FALSE){
    if(!'list' %in% class(mark_gene)){
      mark_gene <- list('mark_gene'=mark_gene)
    }
    if(is.null(names(mark_gene))==TRUE){
      message('Must give name to the mark_gene list');return(FALSE)
    }
    if(mark_strategy=='color' & is.null(mark_col)==TRUE){
      mark_col <- get.class.color(names(mark_gene))
    }
    if(mark_strategy=='add_column'){
      new_col_name <- paste0('is',names(mark_gene))
      all_ms_tab <- lapply(all_ms_tab,function(x){
        g1 <- x$geneSymbol
        r1 <- do.call(base::cbind,lapply(mark_gene,function(x1){
          ifelse(g1 %in% x1,'TRUE','FALSE')
        }))
        new_x <- base::cbind(x,r1)
        colnames(new_x) <- c(colnames(x),new_col_name)
        new_x
      })
    }
  }
  z_column_index <- 'defined'
  if(is.null(z_column)==TRUE) z_column_index <- 'auto'
  i <- 0
  headerStyle <- openxlsx::createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#4F81BD",
                             border="TopBottom", borderColour = "#4F81BD",wrapText=TRUE) ## style for header line
  for(sheetname in names(all_ms_tab)){ ## list, each item create one sheet
    i <- i +1
    d <- as.data.frame(all_ms_tab[[sheetname]])
    if(z_column_index=='auto') use_z_column <- colnames(d)[grep('^Z\\.',colnames(d),ignore.case = TRUE)] else use_z_column <- base::intersect(colnames(d),z_column)
    openxlsx::addWorksheet(wb,sheetName=sheetname)
    openxlsx::writeData(wb,sheet = i,d)
    openxlsx::addStyle(wb, sheet = i, headerStyle, rows = 1, cols = 1:ncol(d), gridExpand = TRUE) ## add header style
    all_c <- list()
    for(z_col in use_z_column){ ## find colnames with z. (only applied to pipeline excel)
      j <- which(colnames(d)==z_col)
      z1 <- d[,j]; c1 <- z2col(z1,sig_thre=sig_thre)
      all_c[[as.character(j)]] <- c1
    }
    mat_c <- do.call(base::cbind,all_c)
    uni_c <- base::unique(unlist(all_c))
    for(r in uni_c){
      w1 <- which(mat_c==r)
      nn <- nrow(mat_c)
      rr <- w1%%nn+1
      cc <- w1%/%nn+1
      cc[which(rr==1)] <- cc[which(rr==1)]-1
      rr[which(rr==1)] <- nn+1
      openxlsx::addStyle(wb, sheet = i, createStyle(fgFill=r), rows =rr, cols = as.numeric(names(all_c)[cc]))
    }
    if(is.null(mark_gene)==FALSE){
      for(k in names(mark_gene)){
        openxlsx::addStyle(wb, sheet = i, openxlsx::createStyle(fgFill=mark_col[k]),
                           rows =which(toupper(d[,which(colnames(d)=='geneSymbol')]) %in% toupper(mark_gene[[k]]))+1,
                 cols = which(colnames(d)=='geneSymbol')) ## find column with geneSymbol and mark color
      }
    }
    w1 <- which(gsub("\\s","",as.matrix(d))=='TRUE')
    nn <- nrow(d)
    rr <- w1%%nn+1
    cc <- w1%/%nn+1
    cc[which(rr==1)] <- cc[which(rr==1)]-1
    rr[which(rr==1)] <- nn+1
    openxlsx::addStyle(wb, sheet = i, openxlsx::createStyle(fontColour='#FF0000'), rows =rr, cols = as.numeric(cc)) ## find column with TRUE/FALSE and mark with color
  }
  openxlsx::saveWorkbook(wb, out.xlsx, overwrite = TRUE)
  return(TRUE)
  ##
}

#' Load MSigDB Database into R Workspace
#'
#' \code{gs.preload} downloads data from MSigDB and stores it into two variables in R workspace, \code{all_gs2gene} and \code{all_gs2gene_info}.
#' \code{all_gs2gene} is a list object with elements of gene sets collections.
#' \code{all_gs2gene_info} is a data.frame contains the description of each gene sets.
#'
#' This is a pre-processing function for NetBID2 advanced analysis. User only need to input the species name (e.g. "Homo sapiens", "Mus musculus").
#' It will call \code{msigdbr} to download data from MSigDB and save it as RData under the \code{db/} directory with species name.
#'
#' @param use_spe character, name of interested species (e.g. "Homo sapiens", "Mus musculus").
#' Users can call \code{msigdbr_species()} to access the full list of available species names.
#' Default is "Homo sapiens".
#' @param update logical, if TRUE, the previous loaded RData will be updated. Default is FALSE.
#' @param main.dir character, the main file path of user's NetBID2 project.
#' If NULL, will be set to \code{system.file(package = "NetBID2")}. Default is NULL.
#' @param db.dir character, the file path to save the RData. Default is \code{db} directory under the \code{main.dir}, if one has a \code{main.dir}.
#'
#' @return Reture a logical value. If TRUE, MsigDB database is loaded successfully, with \code{all_gs2gene} and \code{all_gs2gene_info} created
#' in the workspace.
#'
#' @examples
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' gs.preload(use_spe='Mus musculus',update=FALSE)
#' print(all_gs2gene_info)
#' # contain the information for all gene set category, category info, category size,
#' ## sub category,sub category info, sub category size
#' print(names(all_gs2gene)) # the first level of the list is the category and sub-category IDs
#' print(str(all_gs2gene$`CP:KEGG`))
#'
#' \dontrun{
#' gs.preload(use_spe='Homo sapiens',update=TRUE)
#' }
#' @export
gs.preload <- function(use_spe='Homo sapiens',update=FALSE,
                       main.dir=NULL,
                       db.dir=sprintf("%s/db/",main.dir)){
  #
  all_input_para <- c('use_spe','update')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  all_spe <- msigdbr::msigdbr_species()[["species_name"]]
  check_res <- c(check_option('update',c(TRUE,FALSE),envir=environment()),
                 check_option('use_spe',all_spe,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  ## only support geneSymbol (because pipeline-generated master table will contain geneSymbol column)
  if(is.null(main.dir)==TRUE){
    main.dir <- system.file(package = "NetBID2")
    message(sprintf('main.dir not set, will use package directory: %s',main.dir))
  }
  if(is.null(db.dir)==TRUE){
    db.dir <- sprintf("%s/db/",main.dir)
  }
  message(sprintf('Will use directory %s as the db.dir',db.dir))
  use_spe1 <- gsub(' ','_',use_spe)
  out_file <- sprintf('%s/%s_gs2gene.RData',db.dir,use_spe1)
  if(file.exists(out_file)==FALSE | update==TRUE){
    message('Begin generating all_gs2gene !')
    all_gs_info <-  msigdbr::msigdbr(species = use_spe) ## use msigdbr_species() to check possible available species
    # for gs_collection
    all_gs_cat <- base::unique(all_gs_info$gs_collection)
    all_gs2gene_1 <- lapply(all_gs_cat,function(x){
      x1 <- all_gs_info[which(all_gs_info$gs_collection==x),]
      all_gs <- base::unique(x1$gs_name)
      x2 <- lapply(all_gs, function(y){
        base::unique(x1$gene_symbol[which(x1$gs_name==y)])
      })
      names(x2) <- all_gs;x2
    })
    names(all_gs2gene_1) <- all_gs_cat
    all_gs_subcat <- base::setdiff(base::unique(all_gs_info$gs_subcollection),"")
    all_gs2gene_2 <- lapply(all_gs_subcat,function(x){
      x1 <- all_gs_info[which(all_gs_info$gs_subcollection==x),]
      all_gs <- base::unique(x1$gs_name)
      x2 <- lapply(all_gs, function(y){
        base::unique(x1$gene_symbol[which(x1$gs_name==y)])
      })
      names(x2) <- all_gs;x2
    })
    names(all_gs2gene_2) <- all_gs_subcat
    all_gs2gene <- c(all_gs2gene_1,all_gs2gene_2)
    #
    all_gs2gene <- all_gs2gene[sort(names(all_gs2gene))]
    gs_size <- unlist(lapply(all_gs2gene,length))
    # info for cat
    info_cat <- c('C1'='positional gene sets','C2'='curated gene set','C3'='motif','C4'='computational','C5'='GO','C6'='oncogenic','C7'='immune','C8'='cell type signature','H'='hallmark genesets')
    info_subcat <- c('CGP'='chemical and genetic perturbations','CP'='Canonical pathways','CP:BIOCARTA'='BioCarta gene sets','CP:KEGG'='KEGG gene sets',
                     'CP:REACTOME'='Reactome gene sets','MIR'='microRNA targets','TFT'='transcription factor targets','CGN'='cancer gene neighborhoods','CM'='cancer modules',
                     'BP'='Biological Process','MF'='Molecular Function','CC'='Cellular Component')
    cat_rel <- base::unique(as.data.frame(all_gs_info[,c('gs_collection','gs_subcollection')]))
    all_gs2gene_info <- data.frame(cat_rel[,1],info_cat[cat_rel[,1]],gs_size[cat_rel[,1]],cat_rel[,2],info_subcat[cat_rel[,2]],gs_size[cat_rel[,2]],stringsAsFactors = FALSE)
    colnames(all_gs2gene_info) <- c('Category','Category_Info','Category_Size','Sub-Category','Sub-Category_Info','Sub-Category_Size')
    all_gs2gene_info <- all_gs2gene_info[order(all_gs2gene_info[,1]),]
    all_gs2gene_info[,c(1,2,4,5)] <- as.data.frame(apply(all_gs2gene_info[,c(1,2,4,5)],2,function(x){x[which(is.na(x)==TRUE)] <- "";x}),stringsAsFactors=FALSE)
    save(all_gs2gene,all_gs2gene_info,file=out_file)
  }
  load(out_file,.GlobalEnv)
  message('all_gs2gene loaded, you could see all_gs2gene_info to check the details !')
  return(TRUE)
}

######################################################### visualization functions
## simple functions to get info
#' Create a vector of each sample's selected phenotye descriptive information.
#'
#' \code{get_obs_label} creates a vector of each sample's selected phenotype descriptive information.
#' This is a helper function for data visualization.
#'
#' @param phe_info data.frame, the phenotype data of the samples.
#' It is a data frame that can store any number of descriptive columns (covariates) for each sample row.
#' To get the phenotype data, using the accessor function \code{pData}.
#' @param use_col a vector of numerics or characters.
#' Users can select the interested descriptive column(s) by calling index or name of the column(s).
#' @param collapse character, an optional character string to separate the results when the length
#' of \code{use_col} is more than 1. Not NA_character. Default is "|".
#'
#' @return
#' Return a vector of selected phenotype descriptive information (covariates) for each sample.
#' Vector name is the sample name.

#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
#' print(use_obs_class)
#' \dontrun{
#'}
#' @export
get_obs_label <- function(phe_info,use_col,collapse='|'){
  #
  all_input_para <- c('phe_info','use_col','collapse')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  w1 <- base::setdiff(use_col,colnames(phe_info))
  if(length(w1)>0){
    message(sprintf('%s not in the colnames(phe_info),please check and re-try!',paste(w1,collapse=';')));return(FALSE);
  }
  obs_label<-phe_info[,use_col];
  if(base::length(use_col)>1){
    obs_label<-apply(obs_label,1,function(x)base::paste(x,collapse=collapse))
  }
  names(obs_label) <- rownames(phe_info);
  obs_label
}

#' Get interested phenotype groups from pData slot of the ExpressionSet object.
#'
#' \code{get_int_group} is a function to extract interested phenotype groups from the ExpressionSet object
#' with 'cluster-meaningful' sample features.
#'
#' @param eset an ExpressionSet object.
#' @return Return a vector of phenotype groups which could be used for sample cluster analysis.
#'
#' @examples
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' intgroups <- get_int_group(network.par$net.eset)
#' @export
get_int_group <- function(eset){
  all_input_para <- c('eset')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  phe <- Biobase::pData(eset)
  feature_len <- apply(phe,2,function(x)base::length(base::unique(x)))
  intgroup <- colnames(phe)[which(feature_len>1 & feature_len<nrow(phe))]
  return(intgroup)
}

#' Get Score to Measure Similarity Between Observed Classification and Predicted Classification.
#'
#' \code{get_clustComp} calculates a score to measure the similarity between two classifications.
#'
#' @param pred_label a vector of characters, the predicted classification labels.
#' @param obs_label a vector of characters, the observed classification labels.
#' @param strategy character, the method applied to calculate the score.
#' Users can choose "ARI (adjusted rand index)", "NMI (normalized mutual information)" or "Jaccard".
#' Default is "ARI".
#' @return Return a score for the measurement of similarity.
#' @examples
#' obs_label <- c('A','A','A','B','B','C','D')
#' pred_label  <- c(1,1,1,1,2,2,2)
#' get_clustComp(pred_label,obs_label)
#' @export
get_clustComp <- function(pred_label, obs_label,strategy='ARI') {
  #
  all_input_para <- c('pred_label','obs_label','strategy')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('strategy',c('ARI','NMI','Jaccard'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(names(pred_label))==TRUE & is.null(names(obs_label))==FALSE){
    message('The names of pred_label will use the order of obs_label')
    names(pred_label) <- names(obs_label)
  }
  if(is.null(names(obs_label))==TRUE & is.null(names(pred_label))==FALSE){
    message('The names of obs_label will use the order of pred_label')
    names(obs_label) <- names(pred_label)
  }
  if(is.null(names(obs_label))==TRUE & is.null(names(pred_label))==TRUE){
    message('Assume pred_label and obs_label have the same order!')
    names(obs_label) <- as.character(1:base::length(obs_label))
    names(pred_label) <- names(obs_label)
  }
  if(strategy=='Jaccard') res1 <- get_jac(pred_label, obs_label) else res1 <- clustComp(pred_label, obs_label)[[strategy]]
  return(res1)
}

# get jaccard accuracy
get_jac <- function(pred_label, obs_label) {
  jac1 <- c()
  for (i in base::unique(pred_label)) {
    jac_index <- c()
    x1 <- names(pred_label)[which(pred_label == i)]
    for (j in base::unique(obs_label)) {
      x2 <- names(obs_label)[which(obs_label == j)]
      jac_index <-
        c(jac_index, base::length(base::intersect(x1, x2)) / base::length(union(x1, x2)))
    }
    jac1 <- c(jac1, base::max(jac_index) * base::length(x1))
  }
  jac1 <- sum(jac1) / base::length(pred_label)
  return(jac1)
}

#' Visualize Each Sample's Observed Label vs. Predicted Label in Table
#'
#' \code{draw.clustComp} draws a table to show each sample's observed label vs. its predicted label.
#' Each row represents an observed label (e.g. one subgroup of disease), each column represents the predicted label created by classification algorithm (e.g K-means).
#'
#' The table provides more details about the side-by-side PCA biplot created by \code{draw.emb.kmeans}.
#' The purpose is to find if any abnormal sample (outlier) exists. The darker the table cell is,
#' the more samples are gathered in the corresponding label.
#'
#' @param pred_label a vector of characters, the predicted labels created by classification (e.g K-means).
#' @param obs_label a vector of characters, the observed labels annotated by phenotype data.
#' @param strategy character, method to quantify the similarity between predicted labels vs. observed labels.
#' Users can choose from "ARI (adjusted rand index)", "NMI (normalized mutual information)" and "Jaccard".
#' Default is "ARI".
#' @param use_col logical, If TRUE, the table will be colored. The more sample gathered in one table cell, the darker shade it has.
#' Default is TRUE.
#' @param low_K integer, a threshold of sample number to be shown in a single cell.
#' If too many samples gathered in a single table cell, it will be challenging for eyes.
#' By setting the value of this threshold, if the number of samples gathered in one table cell exceeded the threshold,
#' only the number will be shown. Otherwise, all samples' names will be listed.
#' Default is 5.
#' @param highlight_clust a vector of characters, the predicted label need to be highlighted in the figure.
#' @param main character, an overall title for the plot.
#' @param clust_cex numeric, text size for the predicted label (column names). Default is 1.
#' @param outlier_cex numeric, text size for the observed label (row names). Default is 0.3.
#' @return Return a matrix of integers and a table for visualization. Rows are predicted label, columns are observed label.
#' Integer is the number of samples gathered in the corresponding label.
#' @examples
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' mat <- Biobase::exprs(network.par$net.eset)
#' phe <- Biobase::pData(network.par$net.eset)
#' intgroup <- 'subgroup'
#' pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,
#'                              obs_label=get_obs_label(phe,intgroup),
#'                              kmeans_strategy='consensus')
#' draw.clustComp(pred_label,get_obs_label(phe,intgroup),outlier_cex=1,low_K=2,use_col=TRUE)
#' draw.clustComp(pred_label,get_obs_label(phe,intgroup),outlier_cex=1,low_K=2,use_col=FALSE)
#' @export
draw.clustComp <- function(pred_label, obs_label,strategy='ARI',
                           use_col=TRUE,low_K=5,
                           highlight_clust=NULL,
                           main=NULL,clust_cex=1,outlier_cex=0.3) {
  #
  all_input_para <- c('pred_label','obs_label','strategy','use_col','low_K','clust_cex','outlier_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('use_col',c(TRUE,FALSE),envir=environment()),
                 check_option('strategy',c('ARI','NMI','Jaccard'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(names(pred_label))==TRUE & is.null(names(obs_label))==FALSE){
    message('The names of pred_label will use the order of obs_label')
    names(pred_label) <- names(obs_label)
  }
  if(is.null(names(obs_label))==TRUE & is.null(names(pred_label))==FALSE){
    message('The names of obs_label will use the order of pred_label')
    names(obs_label) <- names(pred_label)
  }
  if(is.null(names(obs_label))==TRUE & is.null(names(pred_label))==TRUE){
    message('Assume pred_label and obs_label have the same order!')
    names(obs_label) <- as.character(1:base::length(obs_label))
    names(pred_label) <- names(obs_label)
  }
  nn <- names(pred_label)
  k1 <- get_clustComp(pred_label,obs_label,strategy=strategy)
  if(is.null(main)==TRUE){
    mm <- sprintf('%s:%s',strategy,format(k1,digits=4),
                  format(k1,digits=4))
  }else{
    mm <- main
  }
  t1 <- base::table(list(pred_label[nn],obs_label[nn]))
  graphics::layout(1)
  textWidth <- base::max(strwidthMod(colnames(t1),units='inch',cex=clust_cex))+par.char2inch()[1]*1.5
  par(mai=c(0.5,textWidth,1,1))
  if(use_col==TRUE){
    graphics::image(t1,col=c('white',grDevices::colorRampPalette(brewer.pal(8,'Reds'))(base::length(base::unique(as.numeric(t1)))-1)),bty='n',xaxt='n',yaxt='n',
          main=mm)
  }else{
    graphics::image(t1,bty='n',xaxt='n',yaxt='n',
          main=mm,col='white')
  }

  pp <- par()$usr
  graphics::rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
  xx <- base::seq(pp[1],pp[2],length.out = nrow(t1)+1)
  yy <- base::seq(pp[3],pp[4],length.out = ncol(t1)+1)
  xxx <- (xx[1:(base::length(xx)-1)]+xx[2:base::length(xx)])/2
  yyy <- (yy[1:(base::length(yy)-1)]+yy[2:base::length(yy)])/2

  graphics::abline(h=yy);graphics::abline(v=xx)
  graphics::text(pp[1]-par.char2pos()[1],yyy,colnames(t1),adj=1,xpd=TRUE,col=ifelse(colnames(t1) %in% highlight_clust,2,1),cex=clust_cex)
  graphics::text(xxx,pp[4],rownames(t1),srt=0,xpd=TRUE,
       col=ifelse(rownames(t1) %in% highlight_clust,2,1),cex=clust_cex,pos=3)

  for(i in 1:nrow(t1)){
    for(j in 1:ncol(t1)){
      v1 <- t1[i,j]
      if(v1==0) next
      if(v1>low_K){
        graphics::text(xxx[i],yyy[j],v1,cex=clust_cex)
      }else{
        v2 <- names(obs_label)[which(pred_label==rownames(t1)[i] & obs_label==colnames(t1)[j])]
        v2 <- base::paste(v2,collapse='\n')
        graphics::text(xxx[i],yyy[j],v2,cex=outlier_cex)
      }
    }
  }
  return(t1)
}


#' Set Color Scale for Z Statistics Value
#'
#' \code{z2col} is a helper function in \code{out2excel}. It defines the color scale of the Z statistics value.
#'
#' @param x a vector of numerics, a vector of Z statistics.
#' @param n_len integer, number of unique colors. Default is 60.
#' @param sig_thre numeric, the threshold for significance (absolute value of Z statistics). Z values failed to pass the threshold will be colored "white".
#' @param col_min_thre numeric, the lower threshold for the color bar value. Default is 0.01.
#' @param col_max_thre numeric, the upper threshold for the color bar value. Default is 3.
#' @param blue_col a vector of characters, the blue colors used to show the negative Z values. Default is brewer.pal(9,'Set1')[2].
#' @param red_col a vector of characters, the red colors used to show positive Z values. Default is brewer.pal(9,'Set1')[1].
#' @return Return a vector of color codes.
#' @examples
#' t1 <- sort(rnorm(mean=0,sd=2,n=100))
#' graphics::image(as.matrix(t1),col=z2col(t1))
#' @export
z2col <- function(x,n_len=60,sig_thre=0.01,col_min_thre=0.01,col_max_thre=3,
                  blue_col=brewer.pal(9,'Set1')[2],
                  red_col=brewer.pal(9,'Set1')[1]){
  #
  tmp_x <- setdiff(x,c(Inf,-Inf))
  if(length(tmp_x)==0) return(ifelse(x>0,'red','blue'))
  all_input_para <- c('x','n_len','sig_thre','col_min_thre','col_max_thre','blue_col','red_col')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  ## create vector for z-score, can change sig threshold
  x[which(is.na(x)==TRUE)] <- 0
  x[which(x==Inf)]<-  base::max(x[which(x!=Inf)])+1
  x[which(x==-Inf)]<- base::min(x[which(x!=-Inf)])-1
  if(col_min_thre<0) col_min_thre<-0.01
  if(col_max_thre<0) col_max_thre<-3
  c2 <- grDevices::colorRampPalette(c(blue_col,'white',red_col))(n_len)
  r1 <- 1.05*base::max(abs(x)) ## -r1~r1
  if(r1 < col_max_thre){
    r1 <- col_max_thre
  }
  if(col_min_thre>r1){
    r2 <- seq(-r1,r1,length.out=n_len+1)
  }else{
    r21 <- seq(-r1,-col_min_thre,length.out=n_len/2)
    r22 <- base::seq(col_min_thre,r1,length.out=n_len/2)
    r2 <- c(r21,r22)
  }
  x1 <- cut(x,r2)
  names(c2) <- levels(x1)
  x2 <- c2[x1]
  x2[which(abs(x)<sig_thre)] <- 'white'
  x2
}

#' Create Color Codes for a Vector of Characters
#'
#' \code{get.class.color} creates a vector of color codes for the input character vector. This is a helper function to assign nice looking colors for better visualization.
#'
#' @param x a vector of characters, names or labels.
#' @param use_color a vector of color codes, colors to be assigned to each member of \code{x}. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return a vector of color codes, with input character vector as names.

#' @examples
#' get.class.color(c('ClassA','ClassB','ClassC','ClassA','ClassC','ClassC'))
#' get.class.color(c('ClassA','ClassB','ClassC','SHH','WNT','Group3','Group4'))
#' get.class.color(c('ClassA','ClassB','ClassC','SHH','WNT','Group3','Group4'),
#'                  use_color=brewer.pal(8, 'Set1'))
#'
#' pre_define <- c('blue', 'red', 'yellow', 'green','yellow', 'green')
#'                 ## pre-defined colors for MB
#' names(pre_define) <- c('WNT', 'SHH', 'Group3', 'Group4','GroupC', 'GroupD')
#'                 ##pre-defined color name for MB
#' get.class.color(c('ClassA','ClassB','ClassC','SHH','WNT','Group3','Group4'),
#'                 pre_define=pre_define)
#'
#' \dontrun{
#'}
#' @export
get.class.color <- function(x,use_color=NULL,pre_define=NULL) {
  #
  all_input_para <- c('x')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  x <- clean_charVector(x);
  #
  if(is.null(pre_define)==FALSE & is.null(names(pre_define))==TRUE){
    message('No class name for the color vector, please check and re-try !');return(FALSE);
  }
  x <- clean_charVector(x);
  if(is.null(use_color)==TRUE){
    use_color <- brewer.pal(9, 'Set1')
  }
  if (base::length(base::intersect(x, names(pre_define))) == 0) {
    w1 <- base::length(base::unique(x))
    if(w1 < length(use_color)){
      cc2 <- use_color[1:w1]
    }else{
      cc2 <- grDevices::colorRampPalette(use_color)(base::length(base::unique(x)))
    }
    names(cc2) <- base::unique(x)
    cc2 <- cc2[x]
  } else{
    x1 <- base::unique(x)
    x2 <- base::setdiff(x1, names(pre_define))
    cc1 <- NULL
    w1 <- base::length(x2)
    if (w1 > 0) {
      if(w1 < length(use_color)){
        cc1 <- use_color[1:w1]
      }else{
        cc1 <- grDevices::colorRampPalette(use_color)(w1)
      }
      names(cc1) <- x2
    }
    cc2 <- c(pre_define, cc1)
    cc2 <- cc2[x]
  }
  return(cc2)
}

## get color box text,inner function ## refer from web
# https://stackoverflow.com/questions/45366243/text-labels-with-background-colour-in-r
boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA,
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5,
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font')){

  ## The Character expansion factro to be used:
  theCex <- graphics::par('cex')*cex
  ## Is y provided:
  if (missing(y)) y <- x
  ## Recycle coords if necessary:
  if (base::length(x) != base::length(y)){
    lx <- base::length(x)
    ly <- base::length(y)
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
    if (base::length(adj == 1)){
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
  if (base::length(padding) == 1){
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
  if (base::length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,yMid + rectHeight/2))
  } else {
    invisible(base::cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,yMid + rectHeight/2))
  }
}

#' Visualize Sample Clustering Result in 2D Plot
#'
#' \code{draw.2D} creats a 2D plot to visualize the sample clustering result.
#'
#' @param X a vector of numerics, the x coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the first component.
#' @param Y a vector of numerics, the y coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the second component.
#' @param class_label a vector of characters, labels or categories of samples. The vector name should be sample names.
#' @param xlab character, the label for x-axis. Default is "PC1".
#' @param ylab character, the label for y-axis. Default is "PC2".
#' @param legend_cex numeric, giving the amount by which the text of legend should be magnified relative to the default. Default is 0.8.
#' @param main character, an overall title for the plot. Default is "".
#' @param point_cex numeric, giving the amount by which the size of the data points should be magnified relative to the default. Default is 1.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- stats::prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.2D(X=pc[,1],Y=pc[,2],class_label=pred_label)
#' @export
draw.2D <- function(X,Y,class_label,xlab='PC1',ylab='PC2',legend_cex=0.8,main="",point_cex=1,use_color=NULL,pre_define=NULL){
  #
  all_input_para <- c('X','Y','class_label','xlab','ylab','legend_cex','main','point_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  class_label <- clean_charVector(class_label)
  #
  if(base::length(X)!=base::length(Y)){
    message('Input two dimension vector with different length, please check and re-try !');return(FALSE);
  }
  if(base::length(X)!=base::length(class_label)){
    message('Input dimension vector has different length with the class_label, please check and re-try !');return(FALSE);
  }
  par(mai = c(1, 1, 1, 0.5+base::max(strwidthMod(class_label,units='inch',cex=legend_cex,ori=FALSE,mod=FALSE))))
  cls_cc <- get.class.color(class_label,use_color=use_color,pre_define=pre_define) ## get color for each label
  graphics::plot(Y ~ X,pch = 16,cex = point_cex,col = cls_cc,main=main,xlab=xlab,ylab=ylab)
  graphics::legend(par()$usr[2],par()$usr[4],sort(base::unique(class_label)),fill = cls_cc[sort(base::unique(class_label))],
         horiz = FALSE,xpd = TRUE,border = NA,bty = 'n',cex=legend_cex)
  return(TRUE)
}

#' Visualize Sample Clustering Result in 2D Plot with interactive mode
#'
#' \code{draw.2D.interactive} creats a 2D plot to visualize the sample clustering result with interactive mode realized by plotly.
#'
#' @param X a vector of numerics, the x coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the first component.
#' @param Y a vector of numerics, the y coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the second component.
#' @param sample_label a vector of characters, name of samples to be displayed on the figure.
#' @param color_label a vector of characters, labels used to define the point color.
#' @param shape_label a vector of characters, labels used to define the point shape.
#' @param xlab character, the label for x-axis. Default is "PC1".
#' @param ylab character, the label for y-axis. Default is "PC2".
#' @param main character, an overall title for the plot. Default is "".
#' @param point_cex numeric, giving the amount by which the size of the data points should be magnified relative to the default. Default is 1.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return the plotly class object for interactive visualization.
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- stats::prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.2D.interactive(X=pc[,1],Y=pc[,2],
#'                     sample_label=rownames(pc),
#'                     color_label=pred_label,
#'                     pre_define = c('1'='blue','2'='red','3'='yellow','4'='green'))
#' draw.2D.interactive(X=pc[,1],Y=pc[,2],
#'                     sample_label=rownames(pc),
#'                     shape_label=pred_label)
#' @export
draw.2D.interactive <- function(X,Y,sample_label=NULL,color_label=NULL,shape_label=NULL,
                                xlab='PC1',ylab='PC2',main="",point_cex=1,
                                use_color=NULL,pre_define=NULL){
  if(!'plotly' %in% rownames(installed.packages())){
    message('plotly not installed!');return(FALSE);
  }
  #
  all_input_para <- c('X','Y','sample_label','xlab','ylab','main','point_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  sample_label <- clean_charVector(sample_label)
  if(is.null(color_label)==FALSE){
    color_label <- clean_charVector(color_label)
  }
  if(is.null(shape_label)==FALSE){
    shape_label <- clean_charVector(shape_label)
  }
  #
  if(base::length(X)!=base::length(Y)){
    message('Input two dimension vector with different length, please check and re-try !');return(FALSE);
  }
  if(base::length(X)!=base::length(sample_label)){
    message('Input dimension vector has different length with the sample_label, please check and re-try !');return(FALSE);
  }
  if(is.null(shape_label)==TRUE & is.null(color_label)==TRUE){
    message('Either color_label or shape_label is required!');return(FALSE);
  }
  if(is.null(shape_label)==FALSE & is.null(color_label)==FALSE){

    if(base::length(X)!=base::length(color_label)){
      message('Input dimension vector has different length with the color_label, please check and re-try !');return(FALSE);
    }
    color_label.factor <- as.factor(color_label)
    cls_cc <- get.class.color(levels(color_label.factor),use_color=use_color,pre_define=pre_define) ## get color for each label
    if(base::length(X)!=base::length(shape_label)){
      message('Input dimension vector has different length with the shape_label, please check and re-try !');return(FALSE);
    }
    shape_label.factor <- as.factor(shape_label)
    data <- data.frame(X=X,Y=Y,color_label=color_label.factor,shape_label=shape_label.factor);
    display_text <- paste0(sample_label,':',color_label,':',shape_label)
    p <- plotly::plot_ly(data = data, x = ~X, y = ~Y,
                         marker = list(size = point_cex*12),type='scatter',color=~color_label,colors=cls_cc,
                         hoverinfo='text',text=display_text,
                         mode='markers',symbol=~shape_label) %>%
      plotly::layout(title = main,
                     xaxis = list(zeroline = FALSE,title=list(text=xlab)),#20240327
                     yaxis = list(zeroline = FALSE,title=list(text=ylab),
                                  showlegend=TRUE)
      )
  }
  if(is.null(shape_label)==TRUE & is.null(color_label)==FALSE){
    if(base::length(X)!=base::length(color_label)){
      message('Input dimension vector has different length with the color_label, please check and re-try !');return(FALSE);
    }
    color_label.factor <- as.factor(color_label)
    cls_cc <- get.class.color(levels(color_label.factor),use_color=use_color,pre_define=pre_define) ## get color for each label
    data <- data.frame(X=X,Y=Y,color_label=color_label.factor);
    display_text <- paste0(sample_label,':',color_label)
    p <- plotly::plot_ly(data = data, x = ~X, y = ~Y,
                         marker = list(size = point_cex*12),type='scatter',color=~color_label,colors=cls_cc,
                         hoverinfo='text',text=display_text,
                         mode='markers') %>%
      plotly::layout(title = main,
                     yaxis = list(zeroline = FALSE,title=list(text=xlab)),
                     xaxis = list(zeroline = FALSE,title=list(text=ylab),showlegend=TRUE)
      )
  }
  if(is.null(shape_label)==FALSE & is.null(color_label)==TRUE){
    if(base::length(X)!=base::length(shape_label)){
      message('Input dimension vector has different length with the shape_label, please check and re-try !');return(FALSE);
    }
    shape_label.factor <- as.factor(shape_label)
    data <- data.frame(X=X,Y=Y,shape_label=shape_label.factor);
    display_text <- paste0(sample_label,':',shape_label)
    p <- plotly::plot_ly(data = data, x = ~X, y = ~Y,
                         marker = list(size = point_cex*12),type='scatter',color = I('black'),
                         hoverinfo='text',text=display_text,
                         mode='markers',symbol=~shape_label) %>%
      plotly::layout(title = main,
                     yaxis = list(zeroline = FALSE,title=list(text=xlab)),
                     xaxis = list(zeroline = FALSE,title=list(text=ylab))
      )
  }
  return(p)
}


#' Visualize Sample Clustering Result in 2D Plot with Sample Names
#'
#' \code{draw.2D.text} creates a 2D plot with sample names labeled, to visualize the sample clustering result.
#'
#' @param X a vector of numerics, the x coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the first component.
#' @param Y a vector of numerics, the y coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the second component.
#' @param class_label a vector of characters, labels or categories of samples. The vector name should be sample names.
#' @param xlab character, the label for x-axis. Default is "PC1".
#' @param ylab character, the label for y-axis. Default is "PC2".
#' @param legend_cex numeric, giving the amount by which the text of legend should be magnified relative to the default. Default is 0.8.
#' @param main character, an overall title for the plot. Default is "".
#' @param point_cex numeric, giving the amount by which the size of the data points should be magnified relative to the default. Default is 1.
#' @param class_text a vector of characters, the user-defined sample names to label each data points in the plot.
#' If NULL, will use the names of \code{class_label}. Default is NULL.
#' @param text_cex numeric, giving the amount by which the text of \code{class_text} should be magnified relative to the default. Default is NULL.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- stats::prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.2D.text(X=pc[,1],Y=pc[,2],class_label=pred_label,
#'              point_cex=5,text_cex=0.5)
#' @export
draw.2D.text <- function(X,Y,class_label,class_text=NULL,xlab='PC1',ylab='PC2',legend_cex=0.8,main="",
                         point_cex=1,text_cex=NULL,use_color=NULL,pre_define=NULL){
  #
  all_input_para <- c('X','Y','class_label','xlab','ylab','legend_cex','main','point_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  class_label <- clean_charVector(class_label)
  #
  if(base::length(X)!=base::length(Y)){
    message('Input two dimension vector with different length, please check and re-try !');return(FALSE);
  }
  if(base::length(X)!=base::length(class_label)){
    message('Input dimension vector has different length with the class_label, please check and re-try !');return(FALSE);
  }
  par(mai = c(1, 1, 1, 0.5+base::max(strwidthMod(class_label,units='inch',cex=legend_cex,ori=FALSE,mod=FALSE))))
  if(is.null(class_text)==TRUE){
    class_text <- names(class_label)
  }
  cc <- 10/base::length(class_label)
  if(cc<0.05) cc<-0.05
  if(cc>1) cc<-1
  if(is.null(text_cex)==FALSE) cc <- text_cex
  cls_cc <- get.class.color(class_label,use_color=use_color,pre_define=pre_define) ## get color for each label
  graphics::plot(Y ~ X,pch = 16,cex = point_cex,col = cls_cc,main=main,xlab=xlab,ylab=ylab)
  graphics::text(x=X,y=Y,labels=class_text,cex=cc,xpd=TRUE,adj=0.5)
  #print(cc);print(str(nn))
  graphics::legend(par()$usr[2],par()$usr[4],sort(base::unique(class_label)),fill = cls_cc[sort(base::unique(class_label))],
         horiz = FALSE,xpd = TRUE,border = NA,bty = 'n',cex=legend_cex)
  return(TRUE)
}

#' Visualize Sample Clustering Result in 3D Plot
#'
#' \code{draw.3D} creates a 3D plot to visualize the sample clustering result.

#' @param X a vector of numerics, the x coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the first component.
#' @param Y a vector of numerics, the y coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the second component.
#' @param Z a vector of numerics, the z coordinates of points in the plot. If user would like to create a PCA biplot, this parameter should be the third component.
#' @param class_label a vector of characters, labels or categories of samples. The vector name should be sample names.
#' @param xlab character, the label for x-axis. Default is "PC1".
#' @param ylab character, the label for y-axis. Default is "PC2".
#' @param zlab character, the label for z-axis. Default is "PC3".
#' @param legend_cex numeric, giving the amount by which the text of legend should be magnified relative to the default. Default is 0.8.
#' @param main character, an overall title for the plot. Default is "".
#' @param point_cex numeric, giving the amount by which the size of the data points should be magnified relative to the default. Default is 1.
#' @param legend_pos character, the position of legend. Default is "topright".
#' @param legend_ncol integer, number of columns of legend. Default is 1.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#' @param ... other paramters used in \code{scatter3D}.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- stats::prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.3D(X=pc[,1],Y=pc[,2],Z=pc[,3],class_label=pred_label)
#' @export
draw.3D <- function(X,Y,Z,class_label,xlab='PC1',ylab='PC2',zlab='PC3',
                    legend_cex=0.8,main="",point_cex=1,legend_pos='topright',legend_ncol=1,use_color=NULL,pre_define=NULL,...){

  #
  all_input_para <- c('X','Y','Z','class_label','xlab','ylab','zlab','legend_cex','main','point_cex','legend_pos','legend_ncol')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  class_label <- clean_charVector(class_label)
  #
  if(base::length(X)!=base::length(Y) | base::length(X)!=base::length(Z) | base::length(Z)!=base::length(Y)){
    message('Input three dimension vectors with different length, please check and re-try !');return(FALSE);
  }
  if(base::length(X)!=base::length(class_label)){
    message('Input dimension vector has different length with the class_label, please check and re-try !');return(FALSE);
  }
  c1 <- factor(class_label,levels=sort(base::unique(class_label)))
  cls_cc <- get.class.color(class_label,use_color=use_color,pre_define=pre_define) ## get color for each label
  #print(str(class_label))
  #print(str(sort(base::unique(class_label))))
  #print(str(cls_cc))
  #print(cls_cc[sort(base::unique(class_label))])
  plot3D::scatter3D(X,Y,Z,pch = 16,xlab = xlab,ylab = ylab,zlab=zlab,bty='g',
            colvar=1:base::length(c1),
            col=cls_cc,colkey = FALSE,cex=point_cex,main=main,...)
  graphics::legend(legend_pos,sort(base::unique(class_label)),fill = cls_cc[sort(base::unique(class_label))],
         border = NA,bty = 'n',ncol = legend_ncol,cex = legend_cex)
  return(TRUE)
}

#' Visualize Sample Clustering Result in 2D Plot with Ellipse
#'
#' \code{draw.2D.ellipse} creates a 2D plot with an ellipse drawn around each cluster to visualize the sample clustering result.
#'
#' @param X a vector of numerics, the x coordinates of points in the plot. If user would like to creat a PCA biplot, this parameter should be the first component.
#' @param Y a vector of numerics, the y coordinates of points in the plot. If user would like to creat a PCA biplot, this parameter should be the second component.
#' @param class_label a vector of characters, labels or categories of samples. The vector name should be sample names.
#' @param xlab character, the label for x-axis. Default is "PC1".
#' @param ylab character, the label for y-axis. Default is "PC2".
#' @param legend_cex numeric, giving the amount by which the text of legend should be magnified relative to the default. Default is 0.8.
#' @param main character, an overall title for the plot. Default is "".
#' @param point_cex numeric, giving the amount by which the size of the data points should be magnified relative to the default. Default is 1.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- stats::prcomp(t(mat1))$x
#' pred_label <- stats::kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.2D.ellipse(X=pc[,1],Y=pc[,2],class_label=pred_label)
#' @export
draw.2D.ellipse <- function(X,Y,class_label,xlab='PC1',ylab='PC2',legend_cex=0.8,main="",point_cex=1,use_color=NULL,pre_define=NULL){
  #
  all_input_para <- c('X','Y','class_label','xlab','ylab','legend_cex','main','point_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  class_label <- clean_charVector(class_label)
  #
  if(base::length(X)!=base::length(Y)){
    message('Input two dimension vector with different length, please check and re-try !');return(FALSE);
  }
  if(base::length(X)!=base::length(class_label)){
    message('Input dimension vector has different length with the class_label, please check and re-try !');return(FALSE);
  }
  par(mar=c(5,5,5,5))
  get_transparent <- function(x,alpha=0.1){
    grDevices::rgb(t(grDevices::col2rgb(x)/255),alpha=alpha)
  }
  if(is.null(names(class_label))==TRUE) names(class_label) <- paste0('S_',1:base::length(class_label))
  if(is.null(names(X))==TRUE) names(X) <- names(class_label)
  if(is.null(names(Y))==TRUE) names(Y) <- names(class_label)
  X <- X[names(class_label)]
  Y <- Y[names(class_label)]
  cls_cc <- get.class.color(class_label,use_color=use_color,pre_define=pre_define) ## get color for each label
  par(mar=c(5,8,5,8))
  graphics::plot(Y~X,pch = 19,cex = point_cex,xlim=c(base::min(X)-IQR(X)/3,IQR(X)/3+base::max(X)),
       col = cls_cc,bty='n',bg='white',xlab=xlab,ylab=ylab,main=main)
  S_dist <- stats::dist(base::cbind(X,Y)); S_dist<-as.matrix(S_dist);#print(str(S_dist))
  ## add ellipse
  for(i in base::unique(class_label)){
    w0 <- which(class_label==i);w00 <- which(class_label!=i);
    w1 <- base::cbind(X[w0],Y[w0])
    if(is.null(dim(w1))){
      w1 <- t(as.matrix(w1))
    }
    c1 <- get_transparent(cls_cc[which(class_label==i)][1])
    m1 <- base::colMeans(w1)
    d1 <- unlist(lapply(1:nrow(w1),function(x)sqrt(sum((w1[x,]-m1)^2))))
    if(base::length(d1)==1){
      plotrix::draw.ellipse(m1[1],m1[2],a=(par()$usr[2]-par()$usr[1])/30,b=(par()$usr[2]-par()$usr[1])/30,col=c1,border=NA,xpd=TRUE)
      graphics::text(m1[1],m1[2],i,xpd=TRUE,adj=0,cex=legend_cex)
    }
    if(base::length(d1)==2){
      plotrix::draw.ellipse(m1[1],m1[2],a=d1[1],b=d1[1],col=c1,border=NA,xpd=TRUE)
      graphics::text(m1[1],m1[2],i,xpd=TRUE,adj=0,cex=legend_cex)
    }
    if(base::length(d1)>=3){
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
      b <- base::max(b[-1])
      plotrix::draw.ellipse(m1[1],m1[2],a=a*1.05,b=b*1.05,col=c1,border=NA,angle=360/(2*pi)*atan(d3[1,2]/d3[1,1]),xpd=TRUE)
      #
      #s1 <- sample(1:nrow(w1),1) ## random select one to mark label
      s1 <- S_dist[names(class_label)[w0],names(class_label)[w00]] ## select one point with largest distance to nodes outside the cluster
      s1 <- names(class_label)[w0][base::which.max(apply(s1,1,min))];
      d1 <- (par()$usr[2]-par()$usr[1])/30
      m2 <- w1[s1,]
      if(m2[1]>m1[1]){
        boxtext(m2[1]+d1,m2[2],labels=i,adj=0,cex=legend_cex,col.bg=c1)
        graphics::segments(x0=w1[s1,1],y0=w1[s1,2],x1=m2[1]+d1,y1=m2[2],col='dark grey')
      } else{
        boxtext(m2[1]-d1,m2[2],labels=i,adj=1,cex=legend_cex,col.bg=c1)
        graphics::segments(x0=w1[s1,1],y0=w1[s1,2],x1=m2[1]-d1,y1=m2[2],col='dark grey')
      }
    }
  }
  return(TRUE)
}

#' QC plots for ExpressionSet class object.
#'
#' \code{draw.eset.QC} is a function to draw a set of plots for quality control analysis.
#' 6 types of plots will be created, including heatmap, dimension reduction plots (pca, mds, umap),
#' boxplot, density, correlation and meansd.
#'
#' @param eset ExpressionSet class, quality control analysis target.
#' @param outdir character, the directory path for saving output files.
#' @param do.logtransform logical, if TRUE, the log transformation will be performed on the gene expression value. Default is FALSE.
#' @param intgroup a vector of characters, the interested phenotype groups from the ExpressionSet.
#' If NULL, it will automatcially extract all possible groups by \code{get_int_group}.
#' Default is NULL.
#' @param prefix character, the prefix for the QC figures' name. Default is "".
#' @param choose_plot a vector of characters,
#' choose one or many from 'heatmap', 'pca','mds','umap','boxplot', 'density','correlation' and 'meansd' plots.
#' Default is 'heatmap', 'pca', 'boxplot','density' and 'correlation'
#' @param generate_html logical, if TRUE, it will generate a html file by R Markdown. Otherwise, it will generate separate PDF files.
#' Default is TRUE.
#' @param correlation_strategy character, the strategy to calculate the sample correlation,
#' choose from 'pearson' and 'spearman'. Default is 'pearson'.
#' @param plot_all_point logical, if TRUE, the scatterplot will plot all points in the correlation,
#' otherwise will plot the main trend for reducing the figure size. Default is FALSE.
#' @param emb_plot_type character, plot type for dimension reduction methods (pca, mds, umap).
#' Users can choose from "2D","2D.interactive", "2D.ellipse", "2D.text" and "3D".
#' Default is "2D.ellipse".
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#' @examples
#' \dontrun{
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' intgroups <- get_int_group(network.par$net.eset)
#' network.par$out.dir.QC <- getwd() ## set the output directory
#' draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=intgroups,
#'              pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
#' draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=intgroups,
#'              pre_define=c('WNT'='blue','SHH'='red','G4'='green'),
#'              emb_plot_type = '2D.interactive',choose_plot = 'pca')
#' }
#' @export
draw.eset.QC <- function(eset,outdir = '.',do.logtransform = FALSE,intgroup=NULL,prefix = '',
                         choose_plot=c('heatmap','pca','boxplot','density','correlation'),
                         generate_html=TRUE,
                         correlation_strategy='pearson',plot_all_point=FALSE,
                         emb_plot_type='2D.ellipse',
                         use_color=NULL,pre_define=NULL) {
  #
  all_input_para <- c('eset','outdir','do.logtransform','prefix','choose_plot',
                      'generate_html','correlation_strategy','plot_all_point')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('do.logtransform',c(TRUE,FALSE),envir=environment()),
                 check_option('plot_all_point',c(TRUE,FALSE),envir=environment()),
                 check_option('correlation_strategy',c('pearson','spearman'),envir=environment()),
                 check_option('emb_plot_type',
                              c('2D','2D.ellipse','3D','2D.text','2D.interactive'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  w1 <- base::setdiff(choose_plot,c('heatmap','pca','mds','umap','boxplot','density','correlation','meansd'))
  if(base::length(w1)>0){
    message(sprintf('Wrong input for choose_plot, %s not included (Only accept "heatmap","pca","mds","umap","boxplot","density","correlation","meansd").
                    Please check and re-try!',
                    w1));return(FALSE)
  }
  #
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message(paste0("The output directory: \"", outdir, "\" is created!"))
  }else
    message(paste0("The output will overwrite the files in directory: \"",outdir,"\""))
  if(is.null(intgroup)){
    intgroup <- get_int_group(eset)
    if(base::length(intgroup)>6){
      intgroup <- intgroup[1:6]
      message('Warning, too many meaningful sample phenotype columns, will use the first 6 columns for visualization !')
    }
  }
  if(base::length(intgroup)==0){
    message('No intgroup, please check and re-try!');return(FALSE)
  }
  message('Preparing the data...')
  #x  <- prepdata(eset, do.logtransform = do.logtransform, intgroup = intgroup)
  use_mat  <- Biobase::exprs(eset)
  if(nrow(use_mat)<=3){
    message('Too small gene number for plot (<=3), please check and re-try!');return(FALSE);
  }
  if(do.logtransform==TRUE){
    if(base::min(as.numeric(use_mat))<=0){
      message('Warning, the original expression matrix has values not larger than 0, the log-transformation may introduce NA, please manually modifiy and re-try !')
    }
    use_mat <- log2(use_mat)
  }
  if(generate_html==TRUE){
    if(rmarkdown::pandoc_available()==FALSE){
      stop('pandoc not available, please set Sys.setenv(RSTUDIO_PANDOC=$pandoc_installed_path), or set generate_html=FALSE')
    }
    output_rmd_file <- sprintf('%s/%sQC.Rmd',outdir,prefix)
    file.copy(from=system.file('Rmd/eset_QC.Rmd',package = "NetBID2"),to=output_rmd_file)
    rmarkdown::render(output_rmd_file, rmarkdown::html_document(toc = TRUE))
    return(TRUE)
  }
  ## embedding plot
  emb_plot <- intersect(choose_plot,c('pca','mds','umap'))
  if(length(emb_plot)>0){
    if(emb_plot_type=='2D.interactive'){
      message('Please set generate_html=TRUE to get interactive plot!')
      emb_plot_type <- '2D'
    }
    for(each_emb_method in emb_plot){
      fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, each_emb_method))
      pdf(fp, width = 12, height = 6)
      pp <- 0.3+0.7-(ncol(use_mat)-10)*(0.7/900)
      if(ncol(use_mat)<=10) pp <- 1
      if(ncol(use_mat)>1000) pp <- 0.3
      for (i in 1:base::length(intgroup)){
        tmp1 <- draw.emb.kmeans(use_mat,embedding_method=each_emb_method,obs_label=get_obs_label(Biobase::pData(eset),intgroup[i]),
                                verbose=FALSE,point_cex=pp,main=intgroup[i],
                                use_color=use_color,pre_define=pre_define,plot_type = emb_plot_type)
      }
      dev.off()
      message(sprintf('Finish %s plot !',toupper(each_emb_method)))
    }
  }
  ## heatmap
  if('heatmap' %in% choose_plot){
    fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'heatmap'))
    pdf(fp, width = 12, height = 9)
    par(mar = c(6, 6, 6, 6))
    m <- dist2.mod(use_mat)
    dend <- as.dendrogram(hclust(as.dist(m), method = "single"))
    ord <- order.dendrogram(dend)
    m <- m[ord, ord]
    draw.heatmap(mat=m,phenotype_info = Biobase::pData(eset),use_phe=intgroup,use_color=use_color,pre_define=pre_define)
    dev.off()
    message('Finish Heatmap plot !')
  }

  ## meansd
  if('meansd' %in% choose_plot){
    fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'meansd'))
    pdf(fp, width = 12, height = 9)
    draw.meanSdPlot(eset)
    dev.off()
    message('Finish MeanSD plot !')
  }

  ## boxplot
  if('boxplot' %in% choose_plot){
    fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'boxplot'))
    pdf(fp, width = 14, height = 4+ncol(exprs(eset))/4)
    par(mar=c(4,10,4,10))
    for(i in 1:base::length(intgroup)){
      class_label <- get_obs_label(Biobase::pData(eset),intgroup[i])
      cls_cc <- get.class.color(class_label,use_color=use_color,pre_define=pre_define) ## get color for each label
      graphics::boxplot(use_mat,col = cls_cc,ylab = "",xlab='Value',main = sprintf('Boxplot for %s',intgroup[i]),
                     ylim=c(base::min(use_mat),
                            base::max(use_mat)),horizontal=TRUE,las=2)
      pp <- par()$usr
      legend(pp[2],pp[4],legend=base::unique(class_label),
             fill = cls_cc[base::unique(class_label)],
             xpd = TRUE,border = NA,bty = 'n',horiz = FALSE)
    }
    dev.off()
    message('Finish Boxplot !')
  }

  ## density
  if('density' %in% choose_plot){
    fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'density'))
    pdf(fp, width = 12, height = 9)
    for(i in 1:base::length(intgroup)){
      all_dens <- list()
      for (j in 1:ncol(use_mat)) {
        all_dens[[j]] <- stats::density(use_mat[,j],na.rm=TRUE)
      }
      graphics::plot(1,col = 'white',xlim=c(base::min(unlist(lapply(all_dens,function(x)base::min(x$x)))),base::max(unlist(lapply(all_dens,function(x)base::max(x$x))))),
           type = 'l',xlab = "",ylab='Density',main = sprintf('Density plot for %s',intgroup[i]),
           ylim=c(base::min(unlist(lapply(all_dens,function(x)base::min(x$y)))),base::max(unlist(lapply(all_dens,function(x)base::max(x$y))))))
      class_label <- get_obs_label(Biobase::pData(eset),intgroup[i])
      cls_cc <- get.class.color(class_label,use_color=use_color,pre_define=pre_define) ## get color for each label
      for (j in 1:ncol(use_mat)) {
        lines(all_dens[[j]], col = cls_cc[j])
      }
      legend('topright',legend=base::unique(class_label),
             fill = cls_cc[base::unique(class_label)],
             xpd = TRUE,border = NA,bty = 'n',horiz = FALSE)
    }
    dev.off()
    message('Finish Density plot !')
  }

  ## correlation
  if('correlation' %in% choose_plot){
    fp <- file.path(outdir, sprintf("%s%s.pdf", prefix, 'correlation'))
    pdf(fp, width = 12, height = 12);
    par(mar=c(3,3,3,3))
    for(i in 1:base::length(intgroup)){
      class_label <- get_obs_label(Biobase::pData(eset),intgroup[i])
      draw.correlation(use_mat,class_label,main=intgroup[i],correlation_strategy=correlation_strategy,plot_all_point=plot_all_point,
                       use_color=use_color,pre_define=pre_define)
    }
    dev.off()
    message('Finish Density plot !')
  }
  return(TRUE)
}
# two inner functions
draw.meanSdPlot <- function(eset){
  exp_mat <- Biobase::exprs(eset)
  mean_g <- apply(exp_mat,1,mean)
  mean_g <- rank(mean_g)
  sd_g <- apply(exp_mat,1,sd)
  n <- 50
  mean_g_c <- cut(mean_g,breaks = n)
  sd_g_c <- cut(sd_g,breaks = n)
  mat <- base::table(base::cbind(list(mean_g_c,sd_g_c)))
  tmp1 <- stats::aggregate(sd_g,list(mean_g_c),mean); rownames(tmp1) <- tmp1$Group.1
  sd_g_v <- tmp1[levels(mean_g_c),'x']
  #
  mean_g_l <- t(sapply(levels(mean_g_c),function(x)as.numeric(strsplit(gsub('\\(|\\]','\\1',x),',')[[1]])))
  sd_g_l <- t(sapply(levels(sd_g_c),function(x)as.numeric(strsplit(gsub('\\(|\\]','\\1',x),',')[[1]])))
  mean_g_l_m <- base::rowMeans(mean_g_l)
  sd_g_l_m <- base::rowMeans(sd_g_l)
  # col=get_transparent(brewer.pal(8,'Set1')[2],alpha=0.1)
  par(mai=c(1,1,1,2))
  dat <- stats::spline(x=mean_g_l_m,y=sd_g_v,n=n*10)
  graphics::plot(y~x,data=dat,type='l',col=brewer.pal(8,'Set1')[1],lwd=2,xlab='rank(mean)',ylab='sd',
       ylim=c(0,base::max(sd_g)),xlim=c(0,base::max(mean_g)),cex.lab=1.2)
  pp <- par()$usr
  ag <- 30
  r <- base::max(mean_g)/(2*nrow(mat))
  rx <- r/cos(ag*pi/180)*par.pos2inch()[1] # x-inch
  ry <- rx
  x1 <- c(rx*cos(ag*pi/180),rx*cos(ag*pi/180),0,-rx*cos(ag*pi/180),-rx*cos(ag*pi/180),0,rx*cos(ag*pi/180))
  y1 <- c(ry*sin(ag*pi/180),-ry*sin(ag*pi/180),-ry,-ry*sin(ag*pi/180),+ry*sin(ag*pi/180),ry,ry*sin(ag*pi/180))
  #print(x1/par.pos2inch()[1]);print(y1/par.pos2inch()[2])
  simplehexbin <- function(x,y,r,col){
    graphics::polygon(x=x+x1/par.pos2inch()[1],y=y+y1/par.pos2inch()[2],col=col,border = 'white',lwd=0.1)
  }
  mm <- base::max(as.numeric(mat))
  cc <- grDevices::colorRampPalette(brewer.pal(9,'Blues')[c(2:9)])(mm)
  cc <- rev(cc)
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      if(mat[i,j]>0){
        if(j%%2==0) simplehexbin(x=mean_g_l_m[i],y=sd_g_l_m[j],col=cc[mat[i,j]])
        if(j%%2==1) simplehexbin(x=mean_g_l_m[i]-r,y=sd_g_l_m[j],col=cc[mat[i,j]])
      }
    }
  }
  lines(y~x,data=dat,col=brewer.pal(8,'Set1')[1],lwd=2)
  ypos <- base::seq(pp[4]-par.char2pos()[2]*2,pp[3]+(pp[4]-pp[3])/2,length.out=base::length(cc)+1)
  graphics::text(x=pp[2]+par.char2pos()[1]*1.15,y=pp[4]-par.char2pos()[2],pos=4,'Count',xpd=TRUE)
  graphics::rect(xleft=pp[2]+par.char2pos()[1],xright=pp[2]+par.char2pos()[1]*1.5,ybottom=ypos[1:(base::length(ypos)-1)],ytop=ypos[2:base::length(ypos)],col=rev(cc),border=NA,xpd=TRUE)
  graphics::text(x=pp[2]+par.char2pos()[1]*1.5,y=quantile(ypos,probs=c(0,0.5,1)),c(1,round(mm/2),mm),xpd=TRUE,pos=4)
  return()
}
draw.correlation <- function(use_mat,class_label,main='',correlation_strategy='pearson',plot_all_point=FALSE,use_color=NULL,pre_define=NULL){
  # plot part
  n_sample <- ncol(use_mat);
  if(n_sample>30){
    message('Too many samples for drawing the correlation plot, will select the first 30 samples !
            If want to specify the sample list, please choose the subset of eset as input !')
    use_mat <- use_mat[,1:30]
    class_label <- class_label[1:30]
    n_sample <- 30
    #return(FALSE)
  }
  all_s <- colnames(use_mat)
  uni_class <- base::unique(class_label)
  class_label <- class_label[order(factor(class_label,levels=uni_class))]
  use_mat <- use_mat[,names(class_label)]
  cor_res <- stats::cor(use_mat,method=correlation_strategy)
  cc1 <- get.class.color(class_label,use_color=use_color,pre_define=pre_define);
  cc1 <- get_transparent(cc1,0.5);
  p1 <- cumsum(base::table(class_label)[uni_class]); p1 <- c(0,p1)
  p2 <- p1[1:(base::length(p1)-1)]/2+p1[2:base::length(p1)]/2
   ####
  graphics::plot(1,xlim=c(-2,n_sample+1),ylim=c(-2,n_sample+1),xaxt='n',yaxt='n',xlab='',ylab='',col='white',bty='n',xaxs='i',yaxs='i',main=main)
  pp <- 15/n_sample
  if(pp>1) pp <- 1
  if(pp<0.2) pp <- 0.2
  graphics::abline(v=1:n_sample,lwd=0.5,col='grey'); graphics::abline(h=1:n_sample,lwd=0.5,col='grey');
  graphics::rect(xleft=0:(n_sample-1),xright=1:n_sample,ybottom=n_sample,ytop=n_sample+1,col=cc1,border=NA,xpd=TRUE)
  graphics::rect(ybottom=0:(n_sample-1),ytop=1:n_sample,xleft=n_sample,xright=n_sample+1,col=cc1,border=NA,xpd=TRUE)
  graphics::abline(v=p1);graphics::abline(h=p1);
  graphics::text(p2,n_sample+1/2,uni_class,adj=0.5,cex=pp,xpd=TRUE);
  graphics::text(n_sample+1/2,p2,uni_class,adj=0.5,srt=270,cex=pp,xpd=TRUE);
  graphics::text(x=1:n_sample-0.5,y=0,all_s,adj=1,xpd=TRUE,srt=90,cex=pp,xpd=TRUE)
  graphics::text(y=1:n_sample-0.5,x=0,all_s,pos=2,xpd=TRUE,srt=0,cex=pp,xpd=TRUE)
  for(i in 1:(n_sample-1)){
    for(j in (i+1):n_sample){
      graphics::text(y=i-0.5,x=j-0.5,format(cor_res[i,j],digits=2),cex=pp,xpd=TRUE)
    }
  }
  n_box <- round(500/n_sample)
  bb <- base::seq(0,1,length.out=n_box)
  bb1 <- cut(bb,breaks=bb)
  for(j in 1:(n_sample-1)){
    for(i in (j+1):n_sample){
      #graphics::text(y=i-0.5,x=j-0.5,format(cor_res[i,j],digits=2),cex=pp,col=2)
       ei <- use_mat[,i]; ej <- use_mat[,j]
       mmin <- base::min(c(ei,ej)); mmax <- base::max(c(ei,ej))
       ei_m <- (ei-mmin)/(mmax-mmin)
       ej_m <- (ej-mmin)/(mmax-mmin)
       if(plot_all_point==TRUE){
         graphics::points(y=ei_m+i-1,x=ej_m+j-1,cex=0.1,pch=16,col=get_transparent('dark grey',0.6)) ## memory consuming !!!
       }else{
         ei_mc <- cut(ei_m,breaks = bb)
         ej_mc <- cut(ej_m,breaks = bb)
         tt <- base::table(list(ei_mc,ej_mc))
         mm_t <- quantile(tt,probs=0.99)
         tt_thre <- mm_t/100
         for(ii in 1:(n_box-1)){
           for(jj in 1:(n_box-1)){
             ap <- (tt[ii,jj]/mm_t)^(1/5); if(ap>0.9) ap <- 0.9
             if(tt[ii,jj]>tt_thre) graphics::points(y=bb[ii]+i-1,x=bb[jj]+j-1,col=get_transparent('black',ap),pch=16,cex=5/n_box)
           }
         }
       }
    }
  }
  for(i in 1:n_sample){
    ei <- use_mat[,i]
    dei <- stats::density(ei)
    dei$x <- (dei$x-base::min(dei$x))/(base::max(dei$x)-base::min(dei$x))
    dei$y <- (dei$y-base::min(dei$y))/(base::max(dei$y)-base::min(dei$y))
    lines(x=dei$x+i-1,y=dei$y+i-1,col='black',lwd=1)
  }
##
}

#' Visualize the K-means Clustering Result by several dimension reduction and embedding methods.
#'
#' \code{draw.emb.kmeans} is a data visualization function to show the K-means clustering result of a data matrix.
#' A PCA/MDS/UMAP biplot is generated to visualize the clustering. Two biplots side-by-side will show the comparison
#' between real observation labels (left) and the K-means predicted labels (right).
#'
#' This function is mainly used to check the sample clustering result, in aim to detect if any abnormal (outlier) sample(s) exsist.
#' The input is a high-throughput expression matrix.
#' Each row is a gene/transcript/probe and each column is a sample.
#' Users need to provide the real observation label for each sample.
#' A K-value yielding the optimal classification result will be used to generate the predicted labels.
#' A comparision score (choose from ARI, NMI, Jaccard) will be calculated and shown in the figure.
#'
#' @param mat a numeric data matrix, the columns (e.g. sample) will be clustered using the feature (e.g. genes) rows.
#' @param embedding_method character, embedding method, choose from pca, mds and umap. Default is pca.
#' @param all_k a vector of integers, a pre-defined K value. K is the number of final clusters.
#' If NULL, the function will try all possible K values. Default is NULL.
#' @param obs_label a vector of characters, a vector describes each sample's selected phenotype information,
#' using sample name as vector name. Can be obtained by calling \code{get_obs_label}.
#' @param legend_pos character, position of the plot legend. Default is 'topleft'.
#' @param legend_cex numeric, text size of the plot legend. Default is 0.8.
#' @param plot_type character, plot type. Users can choose from "2D", "2D.ellipse", "2D.interactive","2D.text" and "3D". Default is "2D.ellipse".
#' @param point_cex numeric, size of the point in the plot. Default is 1.
#' @param kmeans_strategy character, K-means clustering algorithm.
#' Users can choose "basic" or "consensus". "consensus" is performed by \code{ConsensusClusterPlus}.
#' Default is "basic".
#' @param choose_k_strategy character, method to choose the K-value.
#' Users can choose from "ARI (adjusted rand index)", "NMI (normalized mutual information)" and "Jaccard".
#' Default is "ARI".
#' @param return_type character, the type of result returned.
#' Users can choose "optimal" or "all". "all", all the K-values in \code{all_k} will be returned.
#' "optimal", only the K-value yielding the optimal classification result will be returned.
#' Default is "optimal".
#' @param main character, title for the plot.
#' @param verbose logical, if TRUE, print out detailed information during calculation. Default is TRUE.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return a vector of predicted labels, if \code{return_type} is set to "optimal".
#' Or a list of all possible K-values, if \code{return_type} is set to be "all".
#' If plot_type='2D.interactive', will return a plotly class object for interactive display.
#' @examples
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' mat <- Biobase::exprs(network.par$net.eset)
#' phe <- Biobase::pData(network.par$net.eset)
#' intgroup <- get_int_group(network.par$net.eset)
#' for(i in 1:base::length(intgroup)){
#'  print(intgroup[i])
#'  pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]))
#'  print(base::table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
#' }
#' pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,
#'                              obs_label=get_obs_label(phe,'subgroup'),
#'                              kmeans_strategy='consensus')
#' ## interactive display
#' draw.emb.kmeans(mat=mat,all_k = NULL,
#'                obs_label=get_obs_label(phe,'subgroup'),
#'                plot_type='2D.interactive',
#'                pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
#' @export
draw.emb.kmeans <- function(mat=NULL,embedding_method='pca',all_k=NULL,obs_label=NULL,legend_pos = 'topleft',legend_cex = 0.8,
                            plot_type='2D.ellipse',point_cex=1,
                            kmeans_strategy='basic',choose_k_strategy='ARI',
                            return_type='optimal',main='',verbose=TRUE,
                            use_color=NULL,pre_define=NULL){
  #
  all_input_para <- c('mat','embedding_method','obs_label','legend_pos','legend_cex','plot_type','point_cex','kmeans_strategy','choose_k_strategy',
                      'return_type','main','verbose')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('verbose',c(TRUE,FALSE),envir=environment()),
                 check_option('plot_type',c("2D","2D.interactive", "2D.ellipse","2D.text","3D"),envir=environment()),
                 check_option('kmeans_strategy',c('basic','consensus'),envir=environment()),
                 check_option('embedding_method',c('pca','mds','umap'),envir=environment()),
                 check_option('choose_k_strategy',c('ARI','NMI','Jaccard'),envir=environment()),
                 check_option('return_type',c('optimal','all'),envir=environment()),
                 check_option('legend_pos',c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right","center"),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  obs_label <- clean_charVector(obs_label)
  #
  if(is.null(all_k)==TRUE){
    all_k <- 2:base::min(base::length(obs_label)-1,2*base::length(base::unique(obs_label)))
  }
  if(base::length(base::setdiff(all_k,2:base::length(obs_label)))>0){
    message('some value in all_k exceed the maximum sample size, check and re-try !');return(FALSE);
  }
  if(embedding_method=='pca'){
    pca <- stats::prcomp(t(mat))
    cluster_mat <- pca$x
  }
  if(embedding_method=='mds'){
    d <- dist(t(mat))
    if(plot_type=='3D'){
      fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim
      cluster_mat <- fit$points
    }else{
      fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
      cluster_mat <- fit$points
    }
  }
  if(embedding_method=='umap'){
    ori_cc <- umap::umap.defaults;
    ori_cc$n_epochs <- 2000;
    ori_cc$n_neighbors <- base::min(15,round(ncol(mat)/2));
    if(plot_type=='3D') ori_cc$n_components <- 3
    cluster_mat <- umap::umap(as.matrix(t(mat)),config=ori_cc)$layout
  }
  all_jac <- list()
  all_k_res <- list()
  if(kmeans_strategy=='basic'){
    for(k in all_k){
      tmp_k <- list()
      for(i in 1:10){
        tmp_k[[i]] <- stats::kmeans(cluster_mat,centers=as.numeric(k))$cluster
      }
      pred_label <- tmp_k
      jac <- unlist(lapply(pred_label,function(x){get_clustComp(x, obs_label,strategy=choose_k_strategy)}))
      top_i <- base::which.max(jac)
      all_k_res[[as.character(k)]] <- tmp_k[[top_i]]
    }
  }else{
    all_k_res <- get_consensus_cluster(mat=mat,all_k=all_k)
  }
  for(k in all_k){
    pred_label <- all_k_res[[as.character(k)]]
    jac <- get_clustComp(pred_label, obs_label,strategy = choose_k_strategy)
    all_jac[[as.character(k)]] <- signif(jac,4)
  }
  if(verbose==TRUE) message('Optimal k is chosen by Score between predicted and observed label')
  if(verbose==TRUE) print(all_jac)
  use_k <- all_k[base::which.max(all_jac)]
  pred_label <- all_k_res[[as.character(use_k)]]
  if(verbose==TRUE) message(sprintf('Best Score occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  ##
  if(plot_type=='3D') d1 <- data.frame(id=colnames(mat),X=cluster_mat[,1],Y=cluster_mat[,2],Z=cluster_mat[,3],label=pred_label,stringsAsFactors=FALSE)
  if(plot_type!='3D') d1 <- data.frame(id=colnames(mat),X=cluster_mat[,1],Y=cluster_mat[,2],label=pred_label,stringsAsFactors=FALSE)
  xlab <- sprintf('Coordinate 1 (%s)',embedding_method);
  ylab <- sprintf('Coordinate 2 (%s)',embedding_method);
  zlab <- sprintf('Coordinate 3 (%s)',embedding_method);
  if(embedding_method=='pca'){
    xlab <- sprintf('PC1(%s%s variance)',format(summary(pca)$importance[2,'PC1']*100,digits=3),'%')
    ylab <- sprintf('PC2(%s%s variance)',format(summary(pca)$importance[2,'PC2']*100,digits=3),'%')
    if(plot_type=='3D') zlab <- sprintf('PC3(%s%s variance)',format(summary(pca)$importance[2,'PC3']*100,digits=3),'%')
  }
  if(plot_type=='2D.interactive'){
    p <- draw.2D.interactive(d1$X,d1$Y,
                             sample_label=names(obs_label),
                             color_label=obs_label[d1$id],
                             shape_label=d1$label,
                             xlab=xlab,
                             ylab=ylab,
                             point_cex=point_cex,main=main,
                             use_color=use_color,pre_define=pre_define)
    return(p)
  }
  graphics::layout(t(matrix(1:2)))
  if(plot_type=='2D.ellipse'){
    draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],
                    xlab=xlab,
                    ylab=ylab,
                    legend_cex=legend_cex,point_cex=point_cex,main=main,
                    use_color=use_color,pre_define=pre_define)
    draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,
                    xlab=xlab,
                    ylab=ylab,
                    legend_cex=legend_cex,point_cex=point_cex,
                    main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
                    use_color=use_color,pre_define=pre_define)
  }
  if(plot_type=='2D'){
    draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],
            xlab=xlab,
            ylab=ylab,
            legend_cex=legend_cex,point_cex=point_cex,main=main,
            use_color=use_color,pre_define=pre_define)
    draw.2D(d1$X,d1$Y,class_label=d1$label,
            xlab=xlab,
            ylab=ylab,
            legend_cex=legend_cex,point_cex=point_cex,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
            use_color=use_color,pre_define=pre_define)
  }

  if(plot_type=='2D.text'){
    draw.2D.text(d1$X,d1$Y,class_label=obs_label[d1$id],class_text=d1$id,
                 xlab=xlab,
                 ylab=ylab,
                 legend_cex=legend_cex,point_cex=point_cex,main=main,
                 use_color=use_color,pre_define=pre_define)
    draw.2D.text(d1$X,d1$Y,class_label=d1$label,class_text=d1$id,
                 xlab=xlab,
                 ylab=ylab,
                 legend_cex=legend_cex,point_cex=point_cex,
                 main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
                 use_color=use_color,pre_define=pre_define)
  }
  if(plot_type=='3D'){
    draw.3D(d1$X,d1$Y,d1$Z,class_label=obs_label[d1$id],
            xlab=xlab,
            ylab=ylab,
            zlab=zlab,
            legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos,main=main,
            use_color=use_color,pre_define=pre_define)
    draw.3D(d1$X,d1$Y,d1$Z,class_label=d1$label,
            xlab=xlab,
            ylab=ylab,
            zlab=zlab,
            legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
            use_color=use_color,pre_define=pre_define)
  }
  graphics::layout(1);
  if(return_type=='optimal') return(pred_label)
  if(return_type=='all') return(all_k_res)
}

## get consensus Kmeans results, using M3C(toooo slow),ConsensusClusterPlus
# return cluster results, RCSI score, P-value
# plot==TRUE, plot RSCI+BETA_P
# get_consensus_cluster
get_consensus_cluster <-function(mat,all_k=2:12,clusterAlg="km",plot='png',...)
{
  maxK <- base::max(all_k)
  res1 <- ConsensusClusterPlus::ConsensusClusterPlus(mat,maxK=maxK,clusterAlg=clusterAlg,plot=plot,...)
  cls_res <- lapply(all_k,function(x){
    x1 <- res1[[x]]$consensusClass
    x1
  })
  names(cls_res) <- as.character(all_k)
  cluster_res <- cls_res
  return(cluster_res)
}


#' Draw Cluster Plot Using MICA (cluster algorithm)
#'
#' \code{draw.MICA} is a function to visualize the cluster result for the samples using MICA (mutual information based clustering analysis) algorithm.
#' Users need to give the MICA project information (directory and name), and the samples real labels.
#' MICA returns the K-value that yields the best clustering performance. Users can pick one comparison score to show in the plot, "ARI", "NMI" or "Jaccard".
#' MICA is not suggested, when sample size is small.
#' This function may not work well for the updated version of MICA.
#'
#' @param outdir character, the output directory for running MICA.
#' @param prjname charater, the project name for running MICA.
#' @param all_k a vector of integers, the pre-defined K-values.
#' If NULL, will use all possible K. Default is NULL.
#' @param obs_label a vector of characters, the observed sample labels or categories.
#' @param legend_pos character, position of the legend in the plot. Default is "topleft".
#' @param legend_cex numeric, giving the amount by which the text of legend should be magnified relative to the default. Default is 0.8.
#' @param plot_type character, type of the plot. Users can choose from "2D", "2D.ellipse", "2D.text" and "3D". Default is "2D.ellipse".
#' @param point_cex numeric, giving the amount by which the size of the data points should be magnified relative to the default. Default is 1.
#' @param choose_k_strategy character, method to choose the K-value.
#' Users can choose from "ARI (adjusted rand index)", "NMI (normalized mutual information)" and "Jaccard". Default is "ARI".
#' @param visualization_type character, users can choose from "tsne", "umap" and "mds". Default is "tsne".
#' @param return_type character, the type of result returned. Users can choose "optimal" or "all".
#' "all", all the K-values in all_k will be returned.
#' "optimal", only the K-value yielding the optimal classification result will be returned. Default is "optimal".
#' @param main character, title for the plot.
#' @param verbose logical, if TRUE, print out detailed information during calculation. Default is TRUE.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return a vector of the predicted label (if return_type is "optimal") and a list of all possible K- values (if return_type is "all").
#' @export
draw.MICA <- function(outdir=NULL,prjname=NULL,all_k=NULL,obs_label=NULL,
                      legend_pos = 'topleft',legend_cex = 0.8,
                      point_cex=1,plot_type='2D.ellipse',
                      choose_k_strategy='ARI',
                      visualization_type='tsne',return_type='optimal',
                      main='',verbose=TRUE,
                      use_color=NULL,pre_define=NULL) {
  #
  all_input_para <- c('outdir','prjname','obs_label','legend_pos','legend_cex','plot_type','point_cex','choose_k_strategy','visualization_type',
                      'return_type','main','verbose','all_k')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('verbose',c(TRUE,FALSE),envir=environment()),
                 check_option('plot_type',c("2D", "2D.ellipse","2D.text","3D"),envir=environment()),
                 check_option('choose_k_strategy',c('ARI','NMI','Jaccard'),envir=environment()),
                 check_option('visualization_type',c('tsne','umap','mds'),envir=environment()),
                 check_option('return_type',c('optimal','all'),envir=environment()),
                 check_option('legend_pos',c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right","center"),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(plot_type=='3D' & visualization_type=='tsne'){
    message('Current tsne not support for 3D');return(FALSE)
  }
  # choose best k
  res1 <- get_clustComp_MICA(outdir=outdir, all_k=all_k, obs_label=obs_label, prjname = prjname,strategy = choose_k_strategy)
  all_k_res <- res1$all_k_res
  all_jac <- res1$all_jac
  use_k <- all_k[base::which.max(all_jac)]
  message(sprintf('Best Score occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  #
  use_file <- sprintf('%s/%s_k%s_ClusterMem.txt',
                      outdir,prjname,use_k,prjname)
  d1 <- read.delim(use_file, stringsAsFactors = FALSE) ## get cluster results
  if(visualization_type=='umap' | visualization_type=='mds'){
    use_file <- sprintf('%s/%s_reduced.h5',
                        outdir,prjname)
    fid <- rhdf5::H5Fopen(use_file)
    dist_mat <- fid$`mds`$block0_values[1:19,] ###
    if(visualization_type=='mds'){
      X <- fid$mds$block0_values[1,];
      Y <- fid$mds$block0_values[2,];
      Z <- fid$mds$block0_values[3,];
      d1$X <- X; d1$Y <- Y;d1$Z <- Z
    }else{
      ori_cc <- umap.defaults;
      ori_cc$n_epochs <- 2000;
      ori_cc$n_neighbors <- round(ncol(dist_mat)/use_k)
      ori_cc$min_dist <- 0.01 # 0.1
      if(plot_type=='3D'){
        ori_cc$n_components <- 3
        use_mat_umap <- umap::umap(t(dist_mat),config=ori_cc)
        X <- use_mat_umap$layout[,1];Y <- use_mat_umap$layout[,2]; Z <- use_mat_umap$layout[,3];
        d1$X <- X; d1$Y <- Y;d1$Z <- Z
      }else{
        use_mat_umap <- umap::umap(t(dist_mat),config=ori_cc)
        X <- use_mat_umap$layout[,1];Y <- use_mat_umap$layout[,2]
        d1$X <- X; d1$Y <- Y;
      }
    }
    rhdf5::H5Fclose(fid)
  }
  graphics::layout(t(matrix(1:2)))
  if(plot_type=='2D.ellipse'){
    draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
                    use_color=use_color,pre_define=pre_define,main=main)
    draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
                    main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
                    use_color=use_color,pre_define=pre_define)
  }
  if(plot_type=='2D'){
    draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
            use_color=use_color,pre_define=pre_define,main=main)
    draw.2D(d1$X,d1$Y,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
            use_color=use_color,pre_define=pre_define)
  }
  if(plot_type=='2D.text'){
    draw.2D.text(d1$X,d1$Y,class_label=obs_label[d1$id],class_text=d1$id,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
                 use_color=use_color,pre_define=pre_define,main=main)
    draw.2D.text(d1$X,d1$Y,class_label=d1$label,class_text=d1$id,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
                 main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
                 use_color=use_color,pre_define=pre_define)
  }
  if(plot_type=='3D'){
    draw.3D(d1$X,d1$Y,d1$Z,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',zlab='MICA-3',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos,
            use_color=use_color,pre_define=pre_define,main=main)
    draw.3D(d1$X,d1$Y,d1$Z,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',zlab='MICA-3',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)),
            use_color=use_color,pre_define=pre_define)
  }
  ##Score
  rownames(d1) <- d1$id
  pred_label <- d1[names(obs_label), ]$label
  names(pred_label) <- names(obs_label)
  jac <- get_clustComp(pred_label, obs_label,strategy = choose_k_strategy)
  print(sprintf('Best Score:%s', jac))
  if(return_type=='optimal') return(pred_label)
  if(return_type=='all') return(all_k_res)
}


# get allScore for MICA
get_clustComp_MICA <- function(outdir, all_k, obs_label, prjname = NULL,strategy = 'ARI') {
  all_jac <- list()
  all_k_res <- list()
  for (k in all_k) {
    use_file <- sprintf('%s/%s_k%s_ClusterMem.txt',
                        outdir,prjname,k,prjname)
    d1 <- read.delim(use_file, stringsAsFactors = FALSE)
    ##Score
    rownames(d1) <- d1$id
    pred_label <-
      d1[names(obs_label), ]$label
    names(pred_label) <- names(obs_label)
    jac <- get_clustComp(pred_label, obs_label,strategy=strategy)
    all_jac[[as.character(k)]] <- jac
    all_k_res[[as.character(k)]] <- pred_label
    print(sprintf('Best Score for %d:%s', k, jac))
  }
  return(list(all_k_res=all_k_res,all_jac=all_jac))
}

#' Draw Volcano Plot for Top DE (differentiated expressed) Genes or DA (differentiated activity) Drivers
#'
#' \code{draw.volcanoPlot} draws the volcano plot to identify and visualize DE genes or DA drivers with fold change threshold and significant P-value from the input dataset.
#' The function will return a data.frame of these highlighted genes/drivers.
#'
#' Top genes or drivers will be colored (blue for down-regulated and red for up-regulated) and labeled with their names.
#' This function requires the input of master table and two thresholds of logFC and P-value.
#'
#' @param dat data.frame, the master table created by function \code{generate.masterTable}. Or a table with columns of the following parameters.
#' @param label_col character, the name of the column in \code{dat} contains gene/driver names.
#' @param logFC_col character, the name of the column in \code{dat} contains logFC values.
#' @param Pv_col character, the name of the column in \code{dat} contains P-values.
#' @param logFC_thre numeric, the threshold of logFC. Genes or drivers with absolute logFC value higher than the threshold will be kept.
#' Default is 0.1.
#' @param Pv_thre numeric, the threshold of P-values. Genes or drivers with P-values lower than the threshold will be kept.
#' Default is 0.01.
#' @param show_plot logical, if TRUE, the plot will be shown in the plot pane. Default is TRUE.
#' @param xlab character, a title for the X axis.
#' @param ylab character, a title for the Y axis.
#' @param show_label logical, if TRUE, labels of selected genes or drivers will be displayed in the plot. Default is FALSE.
#' @param label_cex numeric, giving the amount by which the text of genes/drivers label should be magnified relative to the default. Default is 0.5.
#' @param legend_cex numeric, giving the amount by which the text of legend should be magnified relative to the default. Default is 0.7.
#' @param label_type character, users can choose between "origin" and "distribute". If "origin", all the labels will be displayed without location modification.
#' If "distribute", location of labels will be rearranged to avoid overlap. Default is "distribute".
#' @param main character, an overall title for the plot.
#' @param pdf_file character, the path to save the plot as PDF file. If NULL, no PDF file will be created. Default is NULL.
#'
#' @return Return a data.frame of selected significant genes or drivers, with columns contain \code{label_col}, \code{logFC_col} and \code{Pv_col}.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#'
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',
#'                                logFC_col='logFC.G4.Vs.others_DE',
#'                                Pv_col='P.Value.G4.Vs.others_DE',
#'                                logFC_thre=2,Pv_thre=1e-3,
#'                                main='Volcano Plot for G4.Vs.others_DE',
#'                                show_label=FALSE,
#'                                pdf_file=sprintf('%s/vocalno_nolabel_DE.pdf',
#'                                analysis.par$out.dir.PLOT))
#'}
#' @export
draw.volcanoPlot <- function(dat=NULL,label_col=NULL,logFC_col=NULL,Pv_col=NULL,logFC_thre=0.1, Pv_thre=0.01,
                             show_plot=TRUE,
                             xlab='log2 Fold Change',ylab='P-value',show_label=FALSE,label_cex=0.5,legend_cex=0.8,
                             label_type='distribute',main="",pdf_file=NULL){
  #
  all_input_para <- c('dat','label_col','logFC_col','Pv_col','logFC_thre','Pv_thre','show_plot','xlab','ylab',
                      'show_label','label_cex','legend_cex','label_type','main')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('show_plot',c(TRUE,FALSE),envir=environment()),
                 check_option('show_label',c(TRUE,FALSE),envir=environment()),
                 check_option('label_type',c("origin", "distribute"),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  dat <- base::unique(dat[,c(label_col,logFC_col,Pv_col)])
  dat[,label_col] <- as.character(dat[,label_col])
  dat <- dat[order(dat[,3],decreasing=FALSE),]
  dat <- dat[which(is.na(dat[,2])==FALSE),]
  x <- as.numeric(dat[,logFC_col])
  y <- as.numeric(dat[,Pv_col]);
  y <- -log10(y)
  s1 <- which(abs(x)>=logFC_thre & y>= -log10(Pv_thre))
  if(show_plot==FALSE){
    sig_info <- dat[s1,,drop=FALSE]
    return(sig_info)
  }
  plot_part <- function(ori=FALSE,before_off=FALSE){
    geneWidth <- 0
    geneHeight <- 0
    if(base::length(s1)>0){
      s11 <- s1[which(x[s1]>=0)]
      s12 <- s1[which(x[s1]<0)]
      geneWidth  <- base::max(strwidthMod(dat[s1,label_col],'inches',cex=label_cex,ori=ori))*1.05
      geneHeight <- base::max(strheightMod(dat[s1,label_col],'inches',cex=label_cex))*base::max(base::length(s11),base::length(s12))*1.25+par.char2inch()[2]
    }

    if(before_off==TRUE) dev.off()
    # width: 2|genewidth|5|genewidth|1
    # height: 1.5|geneHeight|1.5
    if(is.null(pdf_file)==FALSE){
      if(show_label==TRUE & label_type=='distribute' & base::length(s1)>0){
        pdf(pdf_file,width=8+geneWidth*2,height=3+base::max(5,geneHeight))
      }else{
        pdf(pdf_file,width=8,height=8)
      }
    }
    par(mai=c(1.5,2,1.5,1))
    mm <- base::max(abs(x))
    xr <- 1.15*mm/2.5*(2.5+geneWidth) ## not not consider 4% overflow by setting xaxs='i'
    if(show_label==TRUE & label_type=='distribute'){
      graphics::plot(y~x,pch=16,col=get_transparent('grey',0.7),xlab=xlab,ylab="",
           xlim=c(-xr,xr),ylim=c(0,base::max(y)*1.5),yaxt='n',main=main,cex.lab=1.2,cex.main=1.6,xaxs='i')
    }else{
      graphics::plot(y~x,pch=16,col=get_transparent('grey',0.7),xlab=xlab,ylab="",
           xlim=c(-mm*1.5,mm*1.5),ylim=c(0,base::max(y)*1.5),yaxt='n',main=main,cex.lab=1.2,cex.main=1.6,xaxs='i')
    }
    yyy <- c(1,round(base::seq(1,base::max(y)*1.5,length.out=base::min(base::length(y),5)))) ## max:5
    graphics::axis(side=2,at=c(0,yyy),labels=c(1,format(10^-yyy,scientific = TRUE)),las=2)
    graphics::mtext(side=2,line = 4,ylab,cex=1.2)
    w1 <- which(dat[,Pv_col]==0) ## 2022-05-13
    if(base::length(w1)>0){dat[w1,Pv_col] <- .Machine$double.xmin;} ## remove zero for logFC
    #z_val <- sapply(dat[,Pv_col]*sign(x),combinePvalVector)[1,]
    z_val <- sapply(dat[,Pv_col]*sign(x),function(xx)ifelse(xx ==0, 0, combinePvalVector(xx,twosided = TRUE)[1]))
    if(logFC_thre>0){graphics::abline(v=logFC_thre,lty=2,lwd=0.5);graphics::abline(v=-logFC_thre,lty=2,lwd=0.5)}
    if(Pv_thre<1) graphics::abline(h=-log10(Pv_thre),lty=2,lwd=0.5);
    graphics::points(y~x,pch=16,col=get_transparent('grey',0.7))
    #s1 <- which(abs(x)>=logFC_thre & y>= -log10(Pv_thre))
    s_col <- z2col(c(10,-10),sig_thre=0); names(s_col) <- as.character(c(1,-1))
    graphics::points(y[s1]~x[s1],pch=16,col=s_col[as.character(sign(x[s1]))])
    graphics::legend(0,par()$usr[4],c('Down-regulated','Not-Significant','Up-regulated'),
           fill=c(s_col[2],get_transparent('grey',0.7),s_col[1]),border=NA,bty='o',
           bg='white',box.col='white',horiz=TRUE,xjust=0.5,cex=legend_cex)
    if(show_label==TRUE){
      s11 <- s1[which(x[s1]>=0)]
      s12 <- s1[which(x[s1]<0)]
      dd <- (par()$usr[2]-par()$usr[1])/100
      if(label_type == 'origin'){
        if(base::length(s11)>0) graphics::text(x[s11]+dd,y[s11],dat[s11,label_col],cex=label_cex,adj=0)
        if(base::length(s12)>0) graphics::text(x[s12]-dd,y[s12],dat[s12,label_col],cex=label_cex,adj=1)
      }else{
        if(base::length(s11)>0){
          dd <- (par()$usr[4]-par()$usr[3]-par.char2pos()[2])/(base::length(s11)+1)
          graphics::rect(xright=par()$usr[2],ybottom=par()$usr[3],xleft=base::max(abs(x))*1.1,ytop=par()$usr[4],col='white',border='white')
          ypos <- base::seq(from=par()$usr[3],by=dd,length.out=base::length(s11))+dd; ypos <- rev(ypos);
          graphics::text(base::max(abs(x))*1.15,ypos,dat[s11,label_col],cex=label_cex,adj=0)
          graphics::segments(x0=base::max(abs(x))*1.1,x1=x[s11],y0=ypos,y1=y[s11],lwd=0.4,col=get_transparent(s_col[1],alpha=0.5))
          graphics::segments(x0=base::max(abs(x))*1.1,x1=base::max(abs(x))*1.14,y0=ypos,y1=ypos,lwd=0.4,col=get_transparent(s_col[1],alpha=0.5))
        }
        if(base::length(s12)>0){
          dd <- (par()$usr[4]-par()$usr[3]-par.char2pos()[2])/(base::length(s12)+1)
          graphics::rect(xleft=par()$usr[1],ybottom=par()$usr[3],xright=-base::max(abs(x))*1.1,ytop=par()$usr[4],col='white',border='white')
          ypos <- base::seq(from=par()$usr[3],by=dd,length.out=base::length(s12))+dd; ypos <- rev(ypos);
          text(-base::max(abs(x))*1.15,ypos,dat[s12,label_col],cex=label_cex,adj=1)
          graphics::segments(x0= -base::max(abs(x))*1.1,x1=x[s12],y0=ypos,y1=y[s12],lwd=0.4,col=get_transparent(s_col[2],alpha=0.5))
          graphics::segments(x0= -base::max(abs(x))*1.1,x1=-base::max(abs(x))*1.14,y0=ypos,y1=ypos,lwd=0.4,col=get_transparent(s_col[2],alpha=0.5))
        }
      }
    }
    graphics::rect(xleft=par()$usr[1],xright=par()$usr[2],ybottom=par()$usr[3],ytop=par()$usr[4])
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}
  sig_info <- dat[s1,,drop=FALSE]
  return(sig_info)
  #return(TRUE)
}
###

#' Plot for combined DE (differentiated expressed)/DA (differentiated activity) Vs. original DE/DA
#'
#' \code{draw.combineDE} draw the image for the combined DE/DA Vs. original DE/DA.
#'
#' This plot function need to input the output from \code{combineDE}.
#'
#' @param DE_list list, a list of DE/DA results, with one more component named "combine" that include the combined results.
#' Strongly suggest to use the output from \code{combineDE}.
#' @param main_id character, the main id for display in the figure, must be one of the name in DE_list.
#' If NULL, will use the first name. Default is NULL.
#' @param top_number number for the top significant genes/drivers in the combine results to be displayed on the plot.
#' Default is 30.
#' @param display_col character, column names used to display. Default is 'P.Value'.
#' @param z_col character, column names for Z statistics used for background color bar. Default is 'Z-statistics'.
#' @param digit_num integer, number of digits to display on the plot. Default is 2.
#' @param row_cex numeric, \code{cex} for the row labels displayed on the plot. Default is 1
#' @param column_cex numeric, \code{cex} for the col labels displayed on the plot. Default is 1
#' @param text_cex numeric, \code{cex} for the text displayed on the plot. Default is 1
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return logical value indicating whether the plot has been successfully generated
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' each_subtype <- 'G4'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_limma_G4 <- getDE.limma.2G(eset=analysis.par$cal.eset,
#'                                    G1=G1,G0=G0,
#'                                    G1_name=each_subtype,
#'                                    G0_name='other')
#' each_subtype <- 'SHH'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#' DE_gene_limma_SHH <- getDE.limma.2G(eset=analysis.par$cal.eset,
#'                                     G1=G1,G0=G0,
#'                                     G1_name=each_subtype,
#'                                     G0_name='other')
#' DE_list <- list(G4=DE_gene_limma_G4,SHH=DE_gene_limma_SHH)
#' res2 <- combineDE(DE_list,transfer_tab=NULL)
#' draw.combineDE(res2)
#' @export
draw.combineDE <- function(DE_list=NULL,main_id=NULL,top_number=30,
                           display_col='P.Value',z_col='Z-statistics',
                           digit_num=2,row_cex=1,column_cex=1,text_cex=1,pdf_file=NULL){
  #
  all_input_para <- c('DE_list','top_number','display_col','z_col',
                      'digit_num','row_cex','column_cex','text_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  DE_name <- base::setdiff(names(DE_list),'combine')
  if(top_number > nrow(DE_list$combine)) top_number <- nrow(DE_list$combine)
  w1 <- DE_list$combine[order(DE_list$combine$P.Value)[1:top_number],]
  dd_Z <- do.call(base::cbind,lapply(DE_name,function(x){
    x1 <- DE_list[[x]]
    x2 <- x1[w1[,x],]
    x2[,z_col]
  }))
  colnames(dd_Z) <- DE_name
  dd <- do.call(base::cbind,lapply(DE_name,function(x){
    x1 <- DE_list[[x]]
    x2 <- x1[w1[,x],]
    x2[,display_col]
  }))
  colnames(dd) <- DE_name
  if(is.null(main_id)==TRUE) main_id <- DE_name[1]
  dat <- data.frame(main_id=w1[,main_id],dd_Z,combine=w1[,z_col],stringsAsFactors=FALSE)
  mat <- as.matrix(dat[,-1])
  rownames(mat) <- dat$main_id
  dd <- base::cbind(dd,combine=w1[,display_col])
  mat1 <- signif(dd,digits=digit_num)
  plot_part <- function(ori=FALSE,before_off=FALSE){
    geneWidth <- base::max(strwidthMod(rownames(mat),units='inch',cex=row_cex,ori=ori))*1.05+par.char2inch()[1]
    textWidth <- base::max(strwidthMod(as.character(mat1),units='inch',cex=text_cex,ori=ori))*ncol(mat1)*1.5
    geneHeight <- base::max(strheightMod(rownames(mat),units='inch',cex=row_cex,ori=ori))*nrow(mat)*1.75
    textHeight <- base::max(strheightMod(as.character(mat1),units='inch',cex=text_cex,ori=ori))*nrow(mat)*1.75
    geneHeight <- base::max(geneHeight,textHeight)
    if(before_off==TRUE) dev.off()
    # geneWidth|textWidth|par.char2inch()[1]*8
    if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=geneWidth+textWidth+par.char2inch()[1]*8,height=1.5+geneHeight)
    par(mai=c(0.5,geneWidth,1,par.char2inch()[1]*8));graphics::layout(1)
    draw.heatmap.local(mat,inner_line=TRUE,out_line=TRUE,col_srt=0,display_text_mat=mat1,row_cex=row_cex,column_cex=column_cex,text_cex=text_cex,inner_col='white')
    #legend
    pp <- par()$usr
    cc1 <- grDevices::colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(5)
    draw.colorbar(col=rev(cc1),min_val=-base::max(abs(mat)),max_val=base::max(abs(mat)),n=5,xleft=pp[2]+par.char2pos()[1]*2,
                  xright=pp[2]+par.char2pos()[1]*3,ytop=pp[3]+(pp[4]-pp[3])/2,ybottom=pp[3]+(pp[4]-pp[3])/2-par.char2pos()[2]*5*text_cex*0.8,cex=text_cex*0.8)
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}
  return(TRUE)
}

#' Create Plot to show Differential Expression and Differential Activity Analysis of Top Drivers
#'
#' \code{draw.NetBID} creates two side-by-side heatmaps to show the result of NetBID analysis.
#' Both the differential expression analysis (the right heatmap) and differential activity analysis (the left heatmap) of top drivers are shown in the plot.
#'
#' @param DA_list list, contains the differential activity (DA) analysis result.
#' @param DE_list list, contains the differential expression (DE) analysis result.
#' @param main_id character, the top genes/drivers from which DA comparison group. Must be one of the names in \code{DA_list}.
#' If NULL, the first element name of \code{DA_list} will be used. Default is NULL.
#' @param top_number integer, the number of the top significant genes/drivers to be displayed in the plot.
#' Default is 30.
#' @param DA_display_col character, which statistic column from NetBID analysis is to be used as the column of the left heatmap.
#' Default is "P.Value".
#' @param DE_display_col character, which statistic column from NetBID analysis is to be used as the column of the right heatmap.
#' Default is "logFC".
#' @param z_col character, which statistic column from NetBID analysis is to be used as the color scale of the heatmap.
#' Default is "Z-statistics".
#' @param digit_num integer, indicating the number of decimal places (round) or significant digits (signif) to be used.
#' Default is 2.
#' @param row_cex numeric, giving the amount by which the text of row names should be magnified relative to the default.
#' Default is 1.
#' @param column_cex numeric, giving the amount by which the text of column names should be maginified relative to the default.
#' Default is 1.
#' @param text_cex numeric, giving the amount by which the text of in the table cell should be maginified relative to the default.
#' Default is 1.
#' @param col_srt numeric, rotation angle of the column labels at the bottom of heatmap.
#' Default is 60.
#' @param pdf_file character, the path to save the plot as PDF file.
#' If NULL, PDF file will not be generated. Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot is created successfully.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE)
#' @export
draw.NetBID <- function(DA_list=NULL,DE_list=NULL,main_id=NULL,top_number=30,
                        DA_display_col='P.Value',DE_display_col='logFC',z_col='Z-statistics',digit_num=2,
                        row_cex=1,column_cex=1,text_cex=1,col_srt=60,pdf_file=NULL){
  #
  all_input_para <- c('DA_list','DE_list','top_number','DA_display_col','DE_display_col','z_col',
                      'digit_num','row_cex','column_cex','text_cex','col_srt')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(top_number<1){message('top_number must be larger than 1!');return(FALSE)}
  DA_name <- names(DA_list)
  DE_name <- names(DE_list)
  if(is.null(main_id)==TRUE) main_id <- DA_name[1]
  if(!main_id %in% DA_name){message('main id not in DA list!');return(FALSE)}
  if(top_number > nrow(DA_list[[main_id]])) top_number <- nrow(DA_list[[main_id]])
  w1 <- rownames(DA_list[[main_id]][order(DA_list[[main_id]]$P.Value)[1:top_number],]) ## display ID
  w2 <- gsub('_TF','',w1);  w2 <- gsub('_SIG','',w2)
  # get display
  dd_DA <- do.call(base::cbind,lapply(DA_name,function(x){
    x1 <- DA_list[[x]]
    x2 <- x1[w1,]
    x2[,DA_display_col]
  }))
  dd_DE <- do.call(base::cbind,lapply(DE_name,function(x){
    x1 <- DE_list[[x]]
    x2 <- x1[w2,]
    x2[,DE_display_col]
  }))
  colnames(dd_DA) <- DA_name; rownames(dd_DA) <- w1;
  colnames(dd_DE) <- DE_name; rownames(dd_DE) <- w1;
  # get Z
  dd_DA_Z <- do.call(base::cbind,lapply(DA_name,function(x){
    x1 <- DA_list[[x]]
    x2 <- x1[w1,]
    x2[,z_col]
  }))
  dd_DE_Z <- do.call(base::cbind,lapply(DE_name,function(x){
    x1 <- DE_list[[x]]
    x2 <- x1[w2,]
    x2[,z_col]
  }))
  colnames(dd_DA_Z) <- DA_name; rownames(dd_DA_Z) <- w1;
  colnames(dd_DE_Z) <- DE_name; #rownames(dd_DE_Z) <- w1;
  #
  mat1 <- signif(dd_DA,digits=digit_num); mat2 <- signif(dd_DE,digits=digit_num);
  if(base::length(DA_list)==1) {mat1 <- as.matrix(mat1); dd_DA_Z<-as.matrix(dd_DA_Z);}
  if(base::length(DE_list)==1) {mat2 <- as.matrix(mat2); dd_DE_Z<-as.matrix(dd_DE_Z);}
  plot_part <- function(ori=FALSE,before_off=FALSE){
    geneWidth <- base::max(strwidthMod(rownames(dd_DA_Z),units='inch',cex=row_cex,ori=ori))*1.05+par.char2inch()[1]
    textWidth1 <- base::max(strwidthMod(as.character(mat1),units='inch',cex=text_cex,ori=ori))*ncol(mat1)*1.5
    geneHeight1 <- base::max(strheightMod(rownames(dd_DA_Z),units='inch',cex=row_cex,ori=ori))*nrow(mat1)*1.75
    textHeight1 <- base::max(strheightMod(as.character(mat1),units='inch',cex=text_cex,ori=ori))*nrow(mat1)*1.75

    textWidth2 <- base::max(strwidthMod(as.character(mat2),units='inch',cex=text_cex,ori=ori))*ncol(mat2)*1.5
    textHeight2 <- base::max(strheightMod(as.character(mat2),units='inch',cex=text_cex,ori=ori))*nrow(mat2)*1.75

    geneHeight <- base::max(c(geneHeight1,textHeight1,textHeight2))


    if(before_off==TRUE) dev.off()
    # geneWidth|textWidth1|0.3|textWidth2|par.char2inch()[1]*8
    # 1|geneHeight|0.5
    if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=geneWidth+textWidth1+0.3+textWidth2+par.char2inch()[1]*8,height=2+geneHeight)
    graphics::layout(t(as.matrix(c(rep(1,ncol(mat1)),rep(2,ncol(mat2))))))
    par(mai=c(1,geneWidth,1,0))
    mm <- max(abs(c(dd_DE_Z,dd_DA_Z)),na.rm=TRUE)

    draw.heatmap.local(dd_DA_Z,inner_line=TRUE,out_line=TRUE,col_srt=col_srt,display_text_mat=mat1,
                       row_cex=row_cex,column_cex=column_cex,text_cex=text_cex,text_col='down',inner_col='white',bb_max=mm)
    pp <- par()$usr; graphics::rect(xleft=pp[1],xright=pp[2],ybottom=pp[4],ytop=pp[4]+2*(pp[4]-pp[3])/nrow(mat1),border='black',col='light grey',xpd=TRUE);
    DA_width <- base::max(strwidthMod(sprintf('(Differential activity:%s)',DA_display_col),units='inch',cex=1,ori=ori))*1.05
    if(textWidth1>DA_width){
      graphics::text(x=pp[1]/2+pp[2]/2,y=pp[4]+1*(pp[4]-pp[3])/nrow(mat1),labels=sprintf('NetBID\n(Differential activity:%s)',DA_display_col),xpd=TRUE,cex=column_cex)
    }else{
      graphics::text(x=pp[1]/2+pp[2]/2,y=pp[4]+1*(pp[4]-pp[3])/nrow(mat1),labels=sprintf('NetBID\n(DA:%s)',DA_display_col),xpd=TRUE,cex=column_cex)
    }
    DE_width <- base::max(strwidthMod('Differential expression',units='inch',cex=1,ori=ori))*1.05

    par(mai=c(1,0.3,1,par.char2inch()[1]*8))
    draw.heatmap.local(dd_DE_Z,inner_line=TRUE,out_line=TRUE,col_srt=col_srt,display_text_mat=mat2,
                       row_cex=row_cex,column_cex=column_cex,text_cex=text_cex,text_col='down',inner_col='white',bb_max=mm)
    pp <- par()$usr; graphics::rect(xleft=pp[1],xright=pp[2],ybottom=pp[4],ytop=pp[4]+2*(pp[4]-pp[3])/nrow(mat1),border='black',col='light grey',xpd=TRUE);
    if(textWidth2>DA_width){
      graphics::text(x=pp[1]/2+pp[2]/2,y=pp[4]+1*(pp[4]-pp[3])/nrow(mat1),labels=sprintf('Differential expression:\n%s',DE_display_col),xpd=TRUE,cex=column_cex)
    }else{
      graphics::text(x=pp[1]/2+pp[2]/2,y=pp[4]+1*(pp[4]-pp[3])/nrow(mat1),labels=sprintf('DE:\n%s',DE_display_col),xpd=TRUE,cex=column_cex)
    }

    #legend
    cc1 <- grDevices::colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(5)
    draw.colorbar(col=rev(cc1),min_val=-mm,max_val=mm,
                  n=5,xleft=pp[2]+par.char2pos()[1]*2,
                  xright=pp[2]+par.char2pos()[1]*3,ytop=pp[3]+(pp[4]-pp[3])/2,ybottom=pp[3]+(pp[4]-pp[3])/2-par.char2pos()[2]*5*text_cex*0.8,cex=text_cex*0.8)
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}
  return(TRUE)
}

### local functions to draw heatmap
draw.heatmap.local <- function(mat,inner_line=FALSE,out_line=TRUE,inner_col='black',
                               n=20,col_srt=60,display_text_mat=NULL,row_cex=1,column_cex=1,text_cex=1,text_col='up',
                               bb_max=NULL){
  if(ncol(mat)>1) mat1 <- mat[nrow(mat):1,] else mat1 <- as.matrix(mat[nrow(mat):1,])
  if(is.null(display_text_mat)==FALSE){
    if(ncol(display_text_mat)>1) display_text_mat <- display_text_mat[nrow(display_text_mat):1,] else display_text_mat <- as.matrix(display_text_mat[nrow(display_text_mat):1,])
  }
  colnames(mat1) <- colnames(mat)
  cc1 <- grDevices::colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(n*2)
  if(is.null(bb_max)==TRUE) bb_max <- base::max(abs(mat1),na.rm=TRUE)
  bb1 <- base::seq(0,bb_max,length.out=n)
  #mat1[which(is.na(mat1)==TRUE)] <- 0;
  if(out_line==TRUE) graphics::image(t(mat1),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n')
  if(out_line==FALSE) graphics::image(t(mat1),col=cc1,breaks=sort(c(-bb1,0,bb1)),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  xx <- base::seq(pp[1],pp[2],length.out=1+ncol(mat1))
  yy <- base::seq(pp[3],pp[4],length.out=1+nrow(mat1))
  xxx <- xx[1:(base::length(xx)-1)]/2+xx[2:base::length(xx)]/2
  yyy <- yy[1:(base::length(yy)-1)]/2+yy[2:base::length(yy)]/2
  if(inner_line==TRUE){
    graphics::abline(v=xx,col=inner_col)
    graphics::abline(h=yy,col=inner_col)
  }
  graphics::text(pp[1]-par.char2pos()[1],yyy,rownames(mat1),adj=1,xpd=TRUE,cex=row_cex) ## distance for one character
  if(text_col=='up'){
    if(col_srt==0){
      graphics::text(xxx,pp[4],colnames(mat1),pos=3,xpd=TRUE,srt=col_srt,cex=column_cex) ## do not use pos, for srt?
    }else{
      graphics::text(xxx,pp[4]+base::max(strheightMod(colnames(mat1))/2)+base::max(strheightMod(colnames(mat1))/10),colnames(mat1),adj=0.5-col_srt/90,xpd=TRUE,srt=col_srt,cex=column_cex) ## do not use pos, for srt?
    }
  }
  if(text_col=='down'){
    if(col_srt!=0) graphics::text(xxx,pp[3]-0.1*(pp[4]-pp[3])/nrow(mat1),colnames(mat1),adj=1,xpd=TRUE,srt=col_srt,cex=column_cex)
    if(col_srt==0) graphics::text(xxx,pp[3]-0.1*(pp[4]-pp[3])/nrow(mat1),colnames(mat1),adj=0.5,xpd=TRUE,srt=col_srt,cex=column_cex)
  }
  if(is.null(display_text_mat)==FALSE){
    for(i in 1:nrow(display_text_mat)){
      for(j in 1:ncol(display_text_mat)){
        graphics::text(xxx[j],yyy[i],display_text_mat[i,j],cex=text_cex)
      }
    }
  }
  return(TRUE)
}
draw.colorbar <- function(col=NULL,min_val=NULL,max_val=NULL,n=5,digit_num=2,direction='vertical',xleft=0,xright=1,ytop=1,ybottom=0,cex=1){
  val_bar<-base::seq(max_val,min_val,length.out = n)
  if(is.null(col)==TRUE) col <- z2col(val_bar)
  if(direction=='vertical'){
    y_pos <- base::seq(ytop,ybottom,length.out=n+1)
    yy <- y_pos[2:base::length(y_pos)]/2+y_pos[1:(base::length(y_pos)-1)]/2
    graphics::rect(xleft=xleft,xright=xright,ytop=y_pos[2:base::length(y_pos)],ybottom=y_pos[1:(base::length(y_pos)-1)],col=col,border='light grey',xpd=TRUE)
    graphics::text(xright,yy,signif(val_bar,digits = digit_num),xpd=TRUE,cex=cex,pos=4)
    graphics::text(xleft/2+xright/2,ytop,'Z value',cex=cex,xpd=TRUE,pos=3)
  }
}

#' Draw Heatmap Plot to Display the Expression Level or Activity Level of Genes and Drivers
#'
#' \code{draw.heatmap} plots the heatmap to see the expression level or activity level of genes and drivers across selected samples.
#'
#' @param mat numeric matrix, the expression/activity matrix. Rows are genes or drivers, columns are selected samples.
#' @param use_genes a vector of characters, selected genes (e.g. "originID"). Default is row names of \code{mat}.
#' @param use_gene_label a vector of characters, a vector of labels for \code{use_genes} (e.g. "geneSymbol" or "gene_label"). Default is \code{use_genes}.
#' @param use_samples a vector of characters, selected samples. Default is column names of \code{mat}.
#' @param use_sample_label a vector of characters, a vector of labels for \code{use_samples}. Default is \code{use_samples}.
#' @param phenotype_info data.frame, phenotype of samples. Users can call \code{Biobase::pData(eset)} to create.
#' The row names should match the column names in \code{mat}. Default is NULL.
#' @param use_phe a list of characters, selected phenotype of samples. A subset of columns from \code{phenotype_info}.
#' Default is NULL.
#' @param main character, an overall title for the plot. Default is "".
#' @param scale character, users can choose from "row", "column" and "none". Indicating if the values should be
#' centered and scaled in either the row direction or the column direction, or none. Default is "none".
#' @param pdf_file character, the file path to save plot as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#' @param cluster_rows,cluster_columns logical, the same parameters in \code{Heatmap}.
#' Please check \code{?Heatmap} for more details. Default is TRUE.
#' @param clustering_distance_rows,clustering_distance_columns character, the same parameters used in \code{Heatmap}.
#' Please check \code{?Heatmap} for more details. Default is "pearson".
#' @param show_row_names,show_column_names logical, the same parameters used in \code{Heatmap}.
#' Please check \code{?Heatmap} for more details. Default is TRUE.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#' @param ..., for more options, please check \code{?Heatmap} for more details.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1/driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' exp_mat <- Biobase::exprs(analysis.par$cal.eset) ## rownames matches originalID
#' ac_mat <- Biobase::exprs(analysis.par$merge.ac.eset) ## rownames matches originalID_label
#' phe_info <- Biobase::pData(analysis.par$cal.eset) ## phenotype information
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' sig_driver <- sig_driver[1:10,]
#' draw.heatmap(mat=exp_mat,use_genes=ms_tab[rownames(sig_driver),'originalID'],
#'             use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
#'             use_samples=colnames(exp_mat),
#'             use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
#'             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup'),
#'             main='Expression for Top drivers',scale='row',
#'             cluster_rows=TRUE,cluster_columns=TRUE,
#'             clustering_distance_rows='pearson',
#'             clustering_distance_columns='pearson',
#'             row_names_gp = gpar(fontsize = 12),
#'             pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
#' draw.heatmap(mat=ac_mat,use_genes=ms_tab[rownames(sig_driver),'originalID_label'],
#'              use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
#'              use_samples=colnames(ac_mat),
#'              use_sample_label=phe_info[colnames(ac_mat),'geo_accession'],
#'              phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup'),
#'              main='Activity for Top drivers',scale='row',
#'              cluster_rows=TRUE,cluster_columns=TRUE,
#'              clustering_distance_rows='pearson',
#'              clustering_distance_columns='pearson',
#'              row_names_gp = gpar(fontsize = 6),
#'              pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
#'
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1/driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' exp_mat <- Biobase::exprs(analysis.par$cal.eset) ## rownames matches originalID
#' ac_mat <- Biobase::exprs(analysis.par$merge.ac.eset) ## rownames matches originalID_label
#' phe_info <- Biobase::pData(analysis.par$cal.eset) ## phenotype information
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' draw.heatmap(mat=exp_mat,use_genes=ms_tab[rownames(sig_driver),'originalID'],
#'             use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
#'             use_samples=colnames(exp_mat),
#'             use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
#'             phenotype_info=phe_info,
#'             use_phe=c('gender','pathology','subgroup'),
#'             main='Expression for Top drivers',scale='row',
#'             cluster_rows=TRUE,cluster_columns=TRUE,
#'             clustering_distance_rows='pearson',
#'             clustering_distance_columns='pearson',
#'             row_names_gp = gpar(fontsize = 12),
#'             pdf_file=sprintf('%s/heatmap_demo1.pdf',
#'             analysis.par$out.dir.PLOT),
#'             pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
#' draw.heatmap(mat=ac_mat,use_genes=ms_tab[rownames(sig_driver),'originalID_label'],
#'              use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
#'              use_samples=colnames(ac_mat),
#'              use_sample_label=phe_info[colnames(ac_mat),'geo_accession'],
#'              phenotype_info=phe_info,
#'              use_phe=c('gender','pathology','subgroup'),
#'              main='Activity for Top drivers',scale='row',
#'              cluster_rows=TRUE,cluster_columns=TRUE,
#'              clustering_distance_rows='pearson',
#'              clustering_distance_columns='pearson',
#'              row_names_gp = gpar(fontsize = 6),
#'              pdf_file=sprintf('%s/heatmap_demo2.pdf',
#'              analysis.par$out.dir.PLOT),
#'              pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
#'}
#' @export
draw.heatmap <- function(mat=NULL,use_genes=rownames(mat),use_gene_label=use_genes,use_samples=colnames(mat),use_sample_label=use_samples,
                         phenotype_info=NULL,use_phe=NULL,main="",scale='none',pdf_file=NULL,
                         cluster_rows=TRUE,cluster_columns=TRUE,
                         show_row_names=TRUE,show_column_names=TRUE,
                         clustering_distance_rows='pearson',clustering_distance_columns='pearson',
                         use_color=NULL,pre_define=NULL,
                         ...){
  #
  all_input_para <- c('mat','use_genes','use_gene_label','use_samples','use_sample_label','main','scale',
                      'cluster_rows','cluster_columns','show_row_names','show_column_names','clustering_distance_rows','clustering_distance_columns')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('cluster_rows',c(TRUE,FALSE),envir=environment()),
                 check_option('cluster_columns',c(TRUE,FALSE),envir=environment()),
                 check_option('show_row_names',c(TRUE,FALSE),envir=environment()),
                 check_option('show_column_names',c(TRUE,FALSE),envir=environment()),
                 check_option('scale',c("none", "row",'column'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  names(use_gene_label) <- use_genes
  names(use_sample_label) <- use_samples
  if(is.null(rownames(phenotype_info))==FALSE){
    ori_phenotype_info <- phenotype_info
    phenotype_info <- as.data.frame(phenotype_info[colnames(mat),],stringsAsFactors=FALSE)
    colnames(phenotype_info) <- colnames(ori_phenotype_info)
  }
  for(i in colnames(phenotype_info)){
    phenotype_info[,i] <- clean_charVector(phenotype_info[,i])
  }
  if(exists('row_names_gp')==FALSE) row_names_gp <- gpar(fontsize = 12)
  if(exists('column_names_gp')==FALSE) column_names_gp <- gpar(fontsize = 12)
  use_genes <- base::intersect(use_genes,rownames(mat))
  use_samples <- base::intersect(use_samples,colnames(mat))
  use_mat <- mat[use_genes,use_samples]
  rownames(use_mat) <- use_gene_label[rownames(use_mat)]
  colnames(use_mat) <- use_sample_label[colnames(use_mat)]
  row_names_max_width <- base::max(strwidthMod(rownames(use_mat),'inches',cex=row_names_gp[[1]]/7))
  row_names_max_width <- unit(row_names_max_width,'inches')
  column_names_max_height <- base::max(strwidthMod(colnames(use_mat),'inches',cex=column_names_gp[[1]]/7))
  column_names_max_height <- unit(column_names_max_height,'inches')
  if(scale=='row'){use_mat <- t(apply(use_mat,1,do.std))}
  if(scale=='column'){use_mat <- apply(use_mat,2,do.std)}
  if(base::length(use_phe)==0){
    if(scale=='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main,name='Raw value',
                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                     show_row_names=show_row_names,show_column_names=show_column_names,
                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
    if(scale!='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main, name='Z value',
                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                     show_row_names=show_row_names,show_column_names=show_column_names,
                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
  }else{
    if(base::length(use_phe)==1){
      use_phe_info <- as.data.frame(phenotype_info[,use_phe],stringsAsFactors=FALSE)
      rownames(use_phe_info) <- rownames(phenotype_info)
      colnames(use_phe_info) <- gsub(' ','.',use_phe)
    }else{
      use_phe_info <- phenotype_info[,use_phe]
      colnames(use_phe_info) <- gsub(' ','.',use_phe)
    }
    use_phe <- colnames(use_phe_info)
    l2c <- get.class.color(base::unique(as.character(as.matrix(use_phe_info))),use_color=use_color,pre_define=pre_define)
    use_col <- lapply(use_phe,function(x)l2c[base::unique(use_phe_info[,x])])
    names(use_col) <- use_phe
    ha_column <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(use_phe_info),col = use_col)
    if(scale=='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main, top_annotation = ha_column,name='Raw value',
                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                     show_row_names=show_row_names,show_column_names=show_column_names,
                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
    if(scale!='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main, top_annotation = ha_column,name='Z value',
                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                     show_row_names=show_row_names,show_column_names=show_column_names,
                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
  }
  ht_list <- ht1
  if(is.null(pdf_file)==FALSE){
    ww <- 1.25*column_names_gp[[1]]/72*ncol(use_mat)+base::max(strwidthMod(rownames(use_mat),'inches',ori=TRUE))+5
    hh <- 1.25*row_names_gp[[1]]/72*nrow(use_mat)+base::max(strwidthMod(colnames(use_mat),'inches',ori=TRUE))+3
    pdf(pdf_file,width=ww,height=hh)
  }
  ComplexHeatmap::draw(ht_list,heatmap_legend_side='left',annotation_legend_side='right')
  if(is.null(pdf_file)==FALSE) {while (!is.null(dev.list()))  dev.off();}
  return(TRUE)
}

################################ Function enrichment related functions

##
#' Merge Selected Major GeneSets from MsigDB
#'
#' \code{merge_gs} combines selected major gene set collections (e.g. "H", "C1" ) together, and return a list object with sub-collections as elements.
#' Each element contains a vector of genes belong to that sub-collection gene set.
#'
#' @param all_gs2gene list, the list returned by \code{gs.preload()}.
#' @param use_gs a vector of characters, names of major gene set collections. Users can call \code{all_gs2gene_info} to see all the available collections.
#' Default is c("H", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG").
#' @return Return a list object with sub-collection gene sets as elements. Each element contains a vector of genes.
#' @examples
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,
#'                        use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
#'
#' @export
merge_gs <- function(all_gs2gene=all_gs2gene,use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5')){
  #
  all_input_para <- c('all_gs2gene')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(use_gs)==TRUE | 'all' %in% use_gs){ ## all gene sets
    use_gs <- base::unique(all_gs2gene_info$Category)
  }
  #
  use_gs <- base::unique(use_gs)
  use_gs <- base::intersect(use_gs,names(all_gs2gene))
  if(base::length(use_gs)==0){
    message('Wrong setting for use_gs, no intersection with names(all_gs2gene), please check and re-try!');return(FALSE);
  }
  message(sprintf('%s will be merged !',paste(use_gs,collapse=';')))
  nn <- unlist(lapply(all_gs2gene[use_gs],names))
  use_gs2gene <- unlist(all_gs2gene[use_gs],recursive = FALSE)
  names(use_gs2gene)<-nn
  if(base::length(base::unique(nn))<base::length(nn)){
    message('duplicate names observed, will merge by name!')
    nn1 <- base::unique(nn)
    mod_use_gs2gene <- list()
    for(i in nn1){
      mod_use_gs2gene[[i]] <- base::unique(unlist(use_gs2gene[which(nn==i)]))
    }
    names(mod_use_gs2gene) <- nn1
    use_gs2gene <- mod_use_gs2gene
  }
  use_gs2gene
}

# simple functions
list2mat <- function(input_list){
  all_x <- base::unique(unlist(input_list))
  all_y <- base::unique(names(input_list))
  mat1 <- matrix(0,nrow=base::length(all_x),ncol=base::length(all_y))
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
    input_v <- clean_charVector(input_v); names(input_v) <- input_vn
    for(i in 1:base::length(input_v)){
      if(input_v[i] %in% names(tmp2)){
        tmp2[[input_v[i]]] <- c(tmp2[[input_v[i]]],names(input_v)[i])
      }else{
        tmp2[[input_v[i]]] <- names(input_v)[i]
      }
    }
  }else{
    tmp1 <- stats::aggregate(names(input_v),list(input_v),function(x)base::paste(x,collapse=sep))
    tmp2 <- tmp1$x; names(tmp2) <- tmp1$Group.1
  }
  tmp2
}

#' Gene Set Enrichment Analysis by Fisher's Exact Test
#'
#' \code{funcEnrich.Fisher} performs gene set enrichment analysis to the input gene list, by using the Fisher's Exact Test.
#' Background gene list is accepeted.
#'
#' @param input_list a vector of characters, a vector of gene symbols. If gene symbols are not available, users can call \code{get_IDtransfer}
#' and \code{get_name_transfertab} for ID conversion.
#' @param bg_list a vector of characters, a vector of background gene symbols. If NULL, genes in \code{gs2gene} will be used as background.
#' Default is NULL.
#' @param gs2gene list, a list contains elements of gene sets.
#' The name of the element is gene set, each element contains a vector of genes in that gene set.
#' If NULL, will use \code{all_gs2gene}, which is created by function \code{gs.preload}. Default is NULL.
#' @param use_gs a vector of characters, the names of gene sets.
#' If \code{gs2gene} is NULL, \code{all_gs2gene} will be used. The \code{use_gs} must be the subset of \code{names(all_gs2gene)}.
#' If "all", all the gene sets in \code{gs2gene} will be used.
#' If user input his own \code{gs2gene} list, \code{use_gs} will be set to "all" as default.
#' Default is c("H", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG").
#' @param min_gs_size numeric, the minimum size of gene set to analysis. Default is 5.
#' @param max_gs_size numeric, the maximum size of gene set to analysis. Default is 500.
#' @param Pv_adj character, method to adjust P-value. Default is "fdr".
#' For details, please check \code{p.adjust.methods}.
#' @param Pv_thre numeric, threshold for the adjusted P-values. Default is 0.1.
#'
#' @return Return a data.frame, contains gene sets with significant enrichment statistics. Column details are as follows,
#'
#' \item{#Name}{Name of the enriched gene set}
#' \item{Total_item}{Background size}
#' \item{Num_item}{Number of genes in the gene set (filtered by the background list)}
#' \item{Num_list}{Number of input genes for testing (filtered by the background list)}
#' \item{Num_list_item}{Number of input genes annotated by the gene set (filtered by the background list)}
#' \item{Ori_P}{Original P-value from Fisher's Exact Test}
#' \item{Adj_P}{Adjusted P-value}
#' \item{Odds_Ratio}{Odds ratio from the 2*2 matrix used for Fisher's Exact Test}
#' \item{Intersected_items}{A vector of the intersected genes, collapsed by ';'. Number is equal to Num_list_item}
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
#'                                bg_list=ms_tab[,'geneSymbol'],
#'                                use_gs=c('H','C5'),
#'                                Pv_thre=0.1,Pv_adj = 'none')
#' \dontrun{
#' }
#' @export
funcEnrich.Fisher <- function(input_list=NULL,bg_list=NULL,
                              use_gs=NULL,
                              gs2gene=NULL,
                              min_gs_size=5,max_gs_size=500,Pv_adj='fdr',Pv_thre=0.1){
  #
  all_input_para <- c('input_list','min_gs_size','max_gs_size','Pv_adj','Pv_thre')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('Pv_adj',c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(gs2gene)==TRUE){ ## use inner gs2gene
    if(is.null(use_gs)==TRUE){
      use_gs <- c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG')
    }else{
      if(use_gs[1] == 'all'){
        use_gs <- c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)
      }
    }
    if(base::length(base::setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)))>0){
      message(sprintf('Input %s not in all_gs2gene, please check all_gs2gene_info (items in Category or Sub-Category) and re-try!',
                      base::paste(base::setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)),collapse=';')));
      return(FALSE)
    }
    if(base::length(use_gs)>1){
      gs2gene <- merge_gs(all_gs2gene,use_gs = use_gs)
    }else{
      gs2gene <- all_gs2gene[[use_gs]]
    }
  }else{
    if(is.null(use_gs)==TRUE){
      use_gs <- 'all'
    }
    if(base::length(use_gs)>1){
      if(class(gs2gene[[1]])=='list') gs2gene <- merge_gs(gs2gene,use_gs = use_gs)
    }else{
      if(use_gs == 'all'){
        if(class(gs2gene[[1]])=='list') gs2gene <- merge_gs(gs2gene,use_gs = names(gs2gene))
      }else{
        if(class(gs2gene[[1]])=='list') gs2gene <- gs2gene[[use_gs]]
      }
    }
  }
  all_gs <- names(gs2gene)
  input_list <- base::unique(input_list)
  bg_list <- base::unique(bg_list)
  if(!is.null(bg_list)){
    use_gs2gene <- lapply(gs2gene,function(x){base::intersect(x,bg_list)})
    names(use_gs2gene) <- names(gs2gene)
  }else{
    use_gs2gene <- gs2gene
  }
  bg_list <- base::unique(unlist(use_gs2gene))
  ## size selection
  s1 <- unlist(lapply(use_gs2gene,length))
  w1 <- which(s1>=min_gs_size & s1<=max_gs_size)
  use_gs2gene <- use_gs2gene[w1]
  all_gs <- names(use_gs2gene) ## all tested gene set number
  ## input filter
  input_list <- base::intersect(input_list,bg_list)
  bg_or <- base::length(input_list)/base::length(bg_list)
  s1 <- unlist(lapply(use_gs2gene,function(x){
    base::length(base::intersect(input_list,x))/base::length(x)
  }))
  w1 <- which(s1>bg_or)
  use_gs2gene <- use_gs2gene[w1]
  empty_vec <- as.data.frame(matrix(NA,ncol=9));colnames(empty_vec) <- c('#Name','Total_item','Num_item','Num_list','Num_list_item','Ori_P','Adj_P','Odds_Ratio','Intersected_items')
  if(base::length(w1)==0) return(empty_vec)
  ## fisher~
  pv <- lapply(use_gs2gene,function(x){
    n11 <- base::length(base::intersect(input_list,x))
    n12 <- base::length(base::intersect(input_list,base::setdiff(bg_list,x)))
    n21 <- base::length(base::setdiff(x,input_list))
    n22 <- base::length(base::setdiff(bg_list,base::unique(c(input_list,x))))
    ft <- fisher.test(base::cbind(c(n11,n12),c(n21,n22)))$p.value
    or <- n11/n12/(n21/n22)
    c(base::length(bg_list),base::length(x),base::length(input_list),n11,ft,or,base::paste(base::intersect(input_list,x),collapse=';'))
  })
  pv <- do.call(rbind,pv)
  pv <- as.data.frame(pv,stringsAsFactors=FALSE)
  colnames(pv) <- c('Total_item','Num_item','Num_list','Num_list_item','Ori_P','Odds_Ratio','Intersected_items')
  pv[1:6] <- lapply(pv[1:6],as.numeric)
  pv$Adj_P <- p.adjust(pv$Ori_P,method=Pv_adj,n=base::length(all_gs))
  pv$`#Name` <- rownames(pv)
  pv <- pv[,c(9,1:5,8,6:7)]
  pv <- pv[order(pv$Ori_P),]
  use_pv <- pv[which(pv$Adj_P<=Pv_thre),]
  return(use_pv)
}

#' Bar Plot for Gene Set Enrichment Analysis Result
#'
#' \code{draw.funcEnrich.bar} draws a horizontal bar plot to visualize the gene set enrichment analysis.
#' Users can choose to display P-values and the top intersected genes from each gene set.
#'
#' @param funcEnrich_res data.frame, containing the result of functional enrichment analysis.
#' It is highly suggested to use \code{funcEnrich.Fisher} to create this data frame.
#' If users decided to prepare the data.frame on their own, please make sure the column names match the following parameters.
#' @param top_number numeric, the number of top enriched gene sets to be displayed. Default is 30.
#' @param Pv_col character, the name of the column in \code{funcEnrich_res} which contains P-value. Default is "Ori_P".
#' @param name_col character, the name of the column in \code{funcEnrich_res} which contains gene set name. Default is "#Name".
#' @param item_col character, the name of the column in \code{funcEnrich_res} which contains intersected genes collapsed with ";".
#' Default is "Intersected_items".
#' @param Pv_thre numeric, threshold of P-values. Genes or drivers with P-values lower than the threshold will be kept. Default is 0.1.
#' @param display_genes logical, if TRUE, the intersected genes will be displayed. Default is FALSE.
#' @param gs_cex numeric, giving the amount by which the text of gene sets names should be magnified relative to the default. Default is 0.5.
#' @param gene_cex numeric, giving the amount by which the text of gene symbols should be magnified relative to the default. Default is 0.5.
#' @param main character, an overall title for the plot.
#' @param bar_col character, the color code used to plot the bars. Default is brewer.pal(8,'RdBu')[7].
#' @param eg_num numeric, the number of intersected gene symbols to display on the right side of the bar. Default is 5.
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
#'                           bg_list=ms_tab[,'geneSymbol'],
#'                           use_gs=c('H','C5'),Pv_thre=0.1,
#'                           Pv_adj = 'none')
#' draw.funcEnrich.bar(funcEnrich_res=res1,top_number=5,
#'                    main='Function Enrichment for Top drivers',
#'                    gs_cex=0.4,gene_cex=0.5)
#' draw.funcEnrich.bar(funcEnrich_res=res1,top_number=3,
#'                    main='Function Enrichment for Top drivers',
#'                    display_genes = TRUE,eg_num=3,
#'                    gs_cex=0.3,gene_cex=0.3)
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
#'                           bg_list=ms_tab[,'geneSymbol'],
#'                           use_gs=c('H','C5'),Pv_thre=0.1,Pv_adj = 'none')
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' draw.funcEnrich.bar(funcEnrich_res=res1,top_number=30,
#'                     main='Function Enrichment for Top drivers',
#'                     pdf_file=sprintf('%s/funcEnrich_bar_nogene.pdf',
#'                     analysis.par$out.dir.PLOT))
#' draw.funcEnrich.bar(funcEnrich_res=res1,top_number=30,
#'                     main='Function Enrichment for Top drivers',
#'                     display_genes = TRUE,gs_cex=0.6,
#'                     pdf_file=sprintf('%s/funcEnrich_bar_withgene.pdf',
#'                     analysis.par$out.dir.PLOT))
#' }
#' @export
draw.funcEnrich.bar <- function(funcEnrich_res=NULL,top_number=30,
                                Pv_col='Ori_P',item_col='Intersected_items',
                                Pv_thre=0.1,display_genes=FALSE,name_col='#Name',
                                gs_cex=0.5,gene_cex=0.5,main="",bar_col=brewer.pal(8,'RdBu')[7],eg_num=5,
                                pdf_file=NULL){
  #
  all_input_para <- c('funcEnrich_res','Pv_col','item_col','Pv_thre','display_genes','name_col',
                      'gs_cex','gene_cex','main','bar_col','eg_num')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('display_genes',c(TRUE,FALSE),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(top_number)==TRUE) top_number <- nrow(funcEnrich_res)
  funcEnrich_res <- funcEnrich_res[which(funcEnrich_res[,Pv_col]<=Pv_thre),]
  if(nrow(funcEnrich_res)>top_number) funcEnrich_res <- funcEnrich_res[1:top_number,]
  pv_val <- funcEnrich_res[,Pv_col]
  s1 <- funcEnrich_res[[item_col]]
  s1 <- sapply(s1,function(x){
    x1 <- unlist(strsplit(x,';'))
    if(base::length(x1)>eg_num){x2 <- sprintf('Total %d items, e.g: %s',base::length(x1),base::paste(x1[1:5],collapse=';'));x<-x2;}
    return(x)
  })
  ## plot design (inch)
  # W: |textWidth|main(4)|textWidth1|
  # H: |1|textHeight|1
  plot_part <- function(ori=FALSE,before_off=FALSE){
    #print(strwidth('W',units = 'inch'))
    textWidth   <- base::max(strwidthMod(funcEnrich_res[,name_col],units='inch',cex=gs_cex,ori=ori))+par.char2inch()[1]
    textWidth1  <- base::max(strwidthMod(s1,units='inch',cex=gene_cex,ori=ori))+par.char2inch()[1]
    each_textHeight <- base::max(strheightMod(funcEnrich_res[,name_col],units='inch',cex=gs_cex))
    if(display_genes==TRUE) each_textHeight <- base::max(each_textHeight,base::max(strheightMod(s1,units='inch',cex=gene_cex)))
    textHeight <- nrow(funcEnrich_res)*1.5*each_textHeight*1.04 # default space=0.2, set 0.5 here
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE){
      if(display_genes==TRUE){
        ww <- textWidth+4+textWidth1; hh <- 2+textHeight
      }else{
        ww <- textWidth+4+1; hh <- 2+textHeight
      }
      #print(sprintf('pdf_file with width: %s and height %s',signif(ww,2),signif(hh,2)))
      pdf(pdf_file,width=ww,height=hh)
    }
    if(display_genes==TRUE) par(mai=c(1,textWidth,1,textWidth1)) else par(mai=c(1,textWidth,1,1))
    a<-graphics::barplot(rev(-log10(pv_val)),horiz=TRUE,border = NA,col=bar_col,main=main,xaxt='n',
               xlim=c(0,max(-log10(pv_val))),space=0.5,width=each_textHeight*par.inch2pos()[2])
    graphics::mtext(side=1,'P-value',line=0.5/par.lineHeight2inch(),xpd=TRUE) ## middle position
    graphics::axis(side=1,at=0:round(par()$usr[2]),labels=10^-(0:round(par()$usr[2])),xpd=TRUE)
    graphics::text(par()$usr[1]-0.5*par.char2pos()[1],a,rev(funcEnrich_res[,name_col]),adj=1,xpd=TRUE,cex=gs_cex);
    if(display_genes==TRUE) graphics::text(rev(-log10(pv_val))+0.5*par.char2pos()[1],a,rev(s1),adj=0,xpd=TRUE,cex=gene_cex)
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}
  return(TRUE)
}

#' Cluster Plot for Gene Set Enrichment Analysis Result
#'
#' \code{draw.funcEnrich.cluster} draws a cluster plot based on binary matrix, to visualize the existence of genes in the enriched gene sets.
#' The P-value of enrichment is also displayed in the plot.
#'
#' @param funcEnrich_res data.frame,  containing the result of functional enrichment analysis.
#' It is highly suggested to use \code{funcEnrich.Fisher} to create this data frame.
#' If users decided to prepare the data.frame on their own, please make sure the column names match the following parameters.
#' @param top_number numeric, the number of top enriched gene sets to be displayed. Default is 30.
#' @param Pv_col character, the name of the column in \code{funcEnrich_res} which contains P-value. Default is "Ori_P".
#' @param name_col character, the name of the column in \code{funcEnrich_res} which contains gene set name. Default is "#Name".
#' @param item_col character, the name of the column in \code{funcEnrich_res} which contains intersected genes collapsed with ";".
#' Default is "Intersected_items".
#' @param Pv_thre numeric, threshold of P-values. Genes or drivers with P-values lower than the threshold will be kept. Default is 0.1.
#' @param gs_cex numeric, giving the amount by which the text of gene sets names should be magnified relative to the default. Default is 0.5.
#' @param gene_cex numeric, giving the amount by which the text of gene symbols should be magnified relative to the default. Default is 0.8.
#' @param pv_cex numeric, giving the amount by which the text of P-values should be magnified relative to the default. Default is 0.7.
#' @param main character, an overall title for the plot.
#' @param h numeric, the height where the cluster tree should be cut. The same parameter as \code{cutree}. Default is 0.95.
#' @param inner_color character, the color code for the indexing box inside the main plot region.
#' Could be one character or a character vector.
#' If want to set the color by gene, could input the color code character with names set as genes.
#' If want to set the color by Z-score, could use `z2col` to generate the color code.
#' Default is brewer.pal(9,'Reds')[3].
#' @param cluster_gs logical, if TRUE, gene sets will be clustered. Default is TRUE.
#' @param cluster_gene logical, if TRUE, genes will be clustered. Default is TRUE.
#' @param use_genes a vector of characters, a vector of gene symbols to display.
#' If NULL, all the genes in the top enriched gene sets will be displayed. Default is NULL.
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#' @param return_mat logical, if TRUE, return a binary matrix. Rows are gene sets, columns are genes. Default if FALSE.
#'
#' @return If \code{return_mat==FALSE}, return a logical value. If TRUE, plot has been created successfully.
#' If \code{return_mat == TRUE}, return a binary matrix of the cluster. Rows are gene sets, columns are genes.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
#'                           bg_list=ms_tab[,'geneSymbol'],
#'                           use_gs=c('H','C5'),Pv_thre=0.1,Pv_adj = 'none')
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.5,
#'                        gene_cex=0.9,pv_cex=0.8)
#' DA_Z <- z2col(ms_tab[rownames(sig_driver),'Z.G4.Vs.others_DA'],
#'         blue_col=brewer.pal(9,'Blues')[3],
#'         red_col=brewer.pal(9,'Reds')[3],col_max_thre=6)
#' names(DA_Z) <- ms_tab[rownames(sig_driver),'geneSymbol']
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.5,
#'                        gene_cex=0.9,pv_cex=0.8,inner_color=DA_Z)
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=10,gs_cex = 0.6,
#'                        gene_cex=1,pv_cex=1,
#'                        cluster_gs=TRUE,cluster_gene = TRUE)
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=15,gs_cex = 0.8,
#'                        gene_cex=0.9,pv_cex=0.8,
#'                        cluster_gs=TRUE,cluster_gene = FALSE)
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=20,gs_cex = 0.8,
#'                         gene_cex=0.9,pv_cex=0.8,
#'                         cluster_gs=FALSE,cluster_gene = TRUE)
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=20,gs_cex = 1,
#'                         gene_cex=1,pv_cex=0.8,
#'                         cluster_gs=FALSE,cluster_gene = FALSE)
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
#'                           bg_list=ms_tab[,'geneSymbol'],
#'                           use_gs=c('H','C5'),Pv_thre=0.1,Pv_adj = 'none')
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.8,
#'                         gene_cex=0.9,pv_cex=0.8,
#'                         pdf_file = sprintf('%s/funcEnrich_cluster.pdf',
#'                         analysis.par$out.dir.PLOT))
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 1.4,
#'                         gene_cex=1.5,pv_cex=1.2,
#'                         pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',
#'                         analysis.par$out.dir.PLOT),
#'                         cluster_gs=TRUE,cluster_gene = TRUE)
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.8,
#'                         gene_cex=0.9,pv_cex=0.8,
#'                         pdf_file = sprintf('%s/funcEnrich_clusterGS.pdf',
#'                         analysis.par$out.dir.PLOT),
#'                         cluster_gs=TRUE,cluster_gene = FALSE)
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.8,
#'                         gene_cex=0.9,pv_cex=0.8,
#'                         pdf_file = sprintf('%s/funcEnrich_clusterGENE.pdf',
#'                         analysis.par$out.dir.PLOT),
#'                         cluster_gs=FALSE,cluster_gene = TRUE)
#' draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 1.5,
#'                         gene_cex=1.4,pv_cex=1.2,
#'                         pdf_file = sprintf('%s/funcEnrich_clusterNO.pdf',
#'                         analysis.par$out.dir.PLOT),
#'                         cluster_gs=FALSE,cluster_gene = FALSE)
#' }
#' @export
draw.funcEnrich.cluster <- function(funcEnrich_res=NULL,top_number=30,
                                    Pv_col='Ori_P',name_col='#Name',item_col='Intersected_items',Pv_thre=0.1,
                                    gs_cex=0.7,gene_cex=0.8,pv_cex=0.7,
                                    main="",h=0.95,inner_color=brewer.pal(9,'Reds')[3],
                                    cluster_gs=TRUE,cluster_gene=TRUE,
                                    pdf_file=NULL,use_genes=NULL,return_mat=FALSE){
  #
  all_input_para <- c('funcEnrich_res','Pv_col','item_col','Pv_thre','name_col',
                      'gs_cex','gene_cex','pv_cex','main','h','inner_color',
                      'cluster_gs','cluster_gene','return_mat')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('cluster_gs',c(TRUE,FALSE),envir=environment()),
                 check_option('cluster_gene',c(TRUE,FALSE),envir=environment()),
                 check_option('return_mat',c(TRUE,FALSE),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(top_number)==TRUE) top_number <- nrow(funcEnrich_res)
  funcEnrich_res <- funcEnrich_res[which(funcEnrich_res[,Pv_col]<=Pv_thre),]
  if(nrow(funcEnrich_res)>top_number) funcEnrich_res <- funcEnrich_res[1:top_number,]
  pv_val <- funcEnrich_res[,Pv_col]; names(pv_val) <- rownames(funcEnrich_res)
  all_g2s <- lapply(funcEnrich_res[,item_col],function(x1)unlist(strsplit(x1,';')))
  names(all_g2s) <- funcEnrich_res[,name_col]
  mat1 <- t(list2mat(all_g2s))
  mat1 <- mat1[rev(funcEnrich_res[,name_col]),]
  if(is.null(use_genes)==FALSE) mat1 <- mat1[,base::intersect(colnames(mat1),use_genes)]
  if(ncol(mat1)==0){message('No genes left, please check and re-try!');return(FALSE)}
  #mat1 <- mat1[,order(colnames(mat1))]
  h_gs <- hclust(dist(mat1,method='binary'))
  h_gene <- hclust(dist(t(mat1),method='binary'))
  gs_cluster <- cutree(h_gs,h=h)
  gene_cluster <- cutree(h_gene,h=h)
  if(cluster_gs==FALSE){gs_cluster <- rep(1,length.out=nrow(mat1));names(gs_cluster)<-rownames(mat1);}
  if(cluster_gene==FALSE){gene_cluster <- rep(1,length.out=ncol(mat1));names(gene_cluster)<-colnames(mat1);}
  cc1 <- grDevices::colorRampPalette(brewer.pal(8,'Dark2'))(base::length(base::unique(gs_cluster)))
  cc2 <- grDevices::colorRampPalette(brewer.pal(9,'Pastel1'))(base::length(base::unique(gene_cluster)))
  cc3 <- grDevices::colorRampPalette(brewer.pal(9,'Reds')[3:9])(100)
  #inner_color <- cc3[1]
  # get gs order
  if(cluster_gs==TRUE) gs_cluster <- gs_cluster[h_gs$order]
  tmp2 <- vec2list(gs_cluster,sep=NULL)
  tmp2 <- tmp2[rev(order(unlist(lapply(tmp2,function(x)base::min(funcEnrich_res[x,Pv_col])))))]
  mat1 <- mat1[unlist(tmp2),]
  if(cluster_gene==TRUE) mat1 <- mat1[,h_gene$order]
  gs_cluster <- gs_cluster[rownames(mat1)]
  gene_cluster <- gene_cluster[colnames(mat1)]
  #a <- heatmap(mat1,scale='none',col=c('white','red'),ColSideColors = cc2[gene_cluster[colnames(mat1)]],
  # RowSideColors = cc1[gs_cluster[rownames(mat1)]],distfun = function(x){dist(x,method='binary')},margins=c(5,20),Rowv=NA,Colv=NA)
  pv <- pv_val[rownames(mat1)]
  pv1 <- format(pv,scientific = TRUE,digits = 3)
  plot_part <- function(ori=TRUE,before_off=FALSE){
    gsWidth <- base::max(strwidthMod(rownames(mat1),units='inch',cex=gs_cex,ori=ori))+par.char2inch()[1]
    gsHeight <- base::max(strheightMod(rownames(mat1),units='inch',cex=gs_cex))*nrow(mat1)*1.75
    geneWidth <- base::max(strheightMod(colnames(mat1),units='inch',cex=gene_cex))*ncol(mat1)*1.5
    geneHeight <- base::max(strwidthMod(colnames(mat1),units='inch',cex=gene_cex,ori=ori))*1.05+par.char2inch()[2]*1.1 ## add bar height
    pvWidth   <- base::max(strwidthMod(pv1,units='inch',cex=pv_cex,ori=ori))+par.char2inch()[1]
    pvHeight   <- base::max(strheightMod(pv1,units='inch',cex=pv_cex,ori=ori))*nrow(mat1)*1.75
    gsHeight <- base::max(gsHeight,pvHeight)
    ##
    ##
    mr <- 1/pvWidth
    geneWidth1 <- ceiling((geneWidth+0.5)*mr)
    pvWidth1 <- 1
    gsWidth1 <- ceiling((gsWidth+0.5)*mr)
    if(geneWidth1+pvWidth1+gsWidth1>200){ ## avoid too many graphics::layout values
      mr <- 180/(geneWidth1+pvWidth1+gsWidth1)
      geneWidth1 <- round(geneWidth1*mr)
      pvWidth1 <- round(pvWidth1*mr)
      if(pvWidth1<1){
        pvWidth1 <- 1;
      }
      gsWidth1 <- ceiling(gsWidth1*mr)
    }
    #####
    ww <- (gsWidth+pvWidth+geneWidth)+1
    hh <- geneHeight+gsHeight+0.5
    #print(c(ww,hh))
    ## pdf
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE) pdf(file=pdf_file,width=ww,height=hh)
    ## graphics::layout: 1|2|3
    graphics::layout(t(matrix(c(rep(1,geneWidth1),rep(2,pvWidth1),rep(3,gsWidth1)),byrow=TRUE)))
    #print(t(matrix(c(rep(1,geneWidth1),rep(2,pvWidth1),rep(3,gsWidth1)),byrow=TRUE)))
    #print(ww);print(hh)
    par(mai=c(0.5,0.5,geneHeight,0));
    if(length(inner_color)==1 | is.null(names(inner_color)==TRUE) | length(intersect(colnames(mat1),names(inner_color)))<1){ ## length 1 or without names
      graphics::image(t(mat1),col=c('white',inner_color[1]),xaxt='n',yaxt='n',bty='n')
    }else{
      gene_order <- colnames(mat1)
      mat2 <- mat1;
      for(i in 1:ncol(mat2)){
        mat2[which(mat2[,i]==1),i] <- i
      }
      w1 <- setdiff(names(inner_color),gene_order)
      inner_color[w1] <- 'grey'
      inner_color_mod <- inner_color[gene_order]
      graphics::image(t(mat2),col=c('white',inner_color_mod),xaxt='n',yaxt='n',bty='n')
    }
    pp <- par()$usr;
    gs_cs <- cumsum(base::table(gs_cluster)[base::unique(gs_cluster)])
    gene_cs <- cumsum(base::table(gene_cluster)[base::unique(gene_cluster)])
    xx <- (pp[2]-pp[1])/base::length(gene_cluster);
    yy <- (pp[4]-pp[3])/base::length(gs_cluster)
    graphics::abline(h=gs_cs*yy+pp[3],col='black',lwd=0.25)
    graphics::abline(v=gene_cs*xx+pp[1],col='black',lwd=0.25)
    graphics::abline(v=pp[1:2],col='black',lwd=0.5);
    graphics::abline(h=pp[3:4],col='black',lwd=0.5)
    ## draw gene name
    yy <- par.char2pos()[2]*0.9
    if(cluster_gene==TRUE){
      graphics::text(c(1:base::length(gene_cluster))*xx+pp[1]-xx/2,pp[4]+yy,colnames(mat1),xpd=TRUE,adj=0,cex=gene_cex,srt=90)
      graphics::rect(xleft=c(1:base::length(gene_cluster))*xx+pp[1]-xx,xright=c(1:base::length(gene_cluster))*xx+pp[1],ybottom=pp[4],ytop=pp[4]+yy*0.8,
           col=cc2[gene_cluster[colnames(mat1)]],xpd=TRUE,border=NA)
    }else{
      graphics::text(c(1:base::length(gene_cluster))*xx+pp[1]-xx/2,pp[4]+0.5*yy,colnames(mat1),xpd=TRUE,adj=0,cex=gene_cex,srt=90)
    }
    # draw p-value
    #pp <- par()$usr;
    par(mai=c(0.5,0,geneHeight,0));
    graphics::plot(1,xaxt='n',yaxt='n',bty='n',xlim=c(pp[1],pp[2]),ylim=c(pp[3],pp[4]),col='white',xlab="",ylab="")
    pp <- par()$usr;
    yy <- (pp[4]-pp[3])/base::length(gs_cluster)
    pv_c <- z2col(-qnorm(pv))
    graphics::rect(xleft=pp[1],xright=pp[2],ybottom=c(1:base::length(gs_cluster))*yy+pp[3]-yy,
         ytop=c(1:base::length(gs_cluster))*yy+pp[3],
         col=pv_c,border = NA)
    graphics::text(0.5,c(1:base::length(gs_cluster))*yy+pp[3]-yy/2,pv1,xpd=TRUE,adj=0.5,cex=pv_cex)
    graphics::abline(v=pp[1],col='black',lwd=2);
    # draw gs name
    par(mai=c(0.5,0,geneHeight,0.5));
    zz <- base::min(c(xx,yy))
    graphics::plot(1,xaxt='n',yaxt='n',bty='n',xlim=c(pp[1],pp[2]),ylim=c(pp[3],pp[4]),col='white',xlab="",ylab="")
    pp <- par()$usr;
    yy <- (pp[4]-pp[3])/base::length(gs_cluster)
    graphics::text(pp[1]+zz*0.2,c(1:base::length(gs_cluster))*yy+pp[3]-yy/2,rownames(mat1),xpd=TRUE,adj=0,cex=gs_cex)
    # get region for p-value
    #graphics::abline(v=pp[1:2],col='black',lwd=0.5);
    graphics::abline(h=pp[3:4],col='black',lwd=0.5);
    graphics::abline(v=pp[1],col='black',lwd=0.5);
    graphics::abline(h=gs_cs*yy+pp[3],col='black',lwd=0.25)
    ##
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}
  if(return_mat==TRUE){
    return(list(mat=mat1,gene_cluster=gene_cluster,gs_cluster=gs_cluster))
  }else{
    return(TRUE)
  }
}
#' Heat Bubble Matrix Plot for Top Drivers in NetBID2 Analysis
#'
#' \code{draw.bubblePlot} combines the matrix bubble chart and the heat map, using bubble color to compare P-values (performed by Fisher's Exact Test) and bubble size to compare the intersected size for target genes.
#' Rows are enriched gene set, columns are top drivers. Users can also check number of protein-coding genes targetted by each driver.
#'
#'
#' @param driver_list a vector of characters, the names of top drivers.
#' @param show_label a vector of characters, the names of top drivers to be displayed in the plot.
#' If NULL, the names in \code{driver_list} will be displayed. Default is NULL.
#' @param Z_val a vector of numerics, the Z statistics of the \code{driver_list}.
#' It is highly suggested to assign names to this vector. If the vector is nameless, the function will use the names of \code{driver_list} by default.
#' @param driver_type a vector of characters, the biotype or other characteristics of the driver.
#' In the demo, we use \code{"gene_biotype"} column in the master table as input.
#' It is highly suggested to assign names to this vector. If the vector is nameless, the function will use the names of \code{driver_list} by default.
#' Default is NULL.
#' @param target_list list, the driver-to-target list object. The names of the list elements are drivers.
#' Each element is a data frame, usually contains at least three columns. "target", target gene names; "MI", mutual information; "spearman", spearman correlation coefficient.
#' Users can call \code{get_net2target_list} to create this list and follow the suggested pipeline.
#' @param transfer2symbol2type data.frame, the ID-conversion table for converting the original ID into gene symbol and gene biotype (at gene level),
#' or into transcript symbol and transcript biotype (at transcript level).
#' It is highly suggested to use \code{get_IDtransfer2symbol2type} to create this ID-conversion table.
#' @param gs2gene list, a list contains elements of gene sets. The name of the element is gene set, each element contains a vector of genes in that gene set.
#' If NULL, will use \code{all_gs2gene}, which is created by function \code{gs.preload}. Default is NULL.
#' @param use_gs a vector of characters, the names of gene sets. If \code{gs2gene} is NULL, \code{all_gs2gene} will be used. And the \code{use_gs} must be the subset of names(all_gs2gene).
#' Please check \code{all_gs2gene_info} for detailed cateogory description.
#' Default is c("H", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG").
#' @param display_gs_list a vector of characters, the names of gene sets to be displayed in the plot.
#' If NULL, all the gene sets will be displayed in descending order of their significance. Default is NULL.
#' @param bg_list a vector of characters, a vector of background gene symbols. If NULL, genes in \code{gs2gene} will be used as background. Default is NULL.
#' @param min_gs_size numeric, the minimum size of gene set to analysis. Default is 5.
#' @param max_gs_size numeric, the maximum size of gene set to analysis, Default is 500.
#' @param Pv_adj character, method to adjust P-value. Default is "none". For details, please check \code{p.adjust.methods}.
#' @param Pv_thre numeric, threshold for the adjusted P-values. Default is 0.1.
#' @param top_geneset_number integer, the number of top enriched gene sets to be displayed in the plot. Default is 30.
#' @param top_driver_number integer, the number of top significant drivers to be displayed in the plot. Default is 30.
#' @param main character, an overall title for the plot.
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be save. Default is NULL.
#' @param mark_gene a vector of characters, a vector of gene symbols to be highlighted red in the plot. Default is NULL.
#' @param driver_cex numeric, giving the amount by which the text of driver symbols should be magnified relative to the default. Default is 1.
#' @param gs_cex numeric, giving the amount by which the text of gene set names should be magnified relative to the default. Default is 1.
#' @param only_return_mat logicial, if TRUE, the function will only return the gene set Vs. driver matrix with value representing the Z-statistics of the significance test;
#' and the plot will not be generated. Default is FALSE.
#'
#' @return Return a logical value if only_return_mat=FALSE. If TRUE, the plot has been created successfully.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' db.preload(use_level='gene',use_spe='human',update=FALSE)
#' use_genes <- base::unique(analysis.par$merge.network$network_dat$target.symbol)
#' transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',
#'                                            use_genes=use_genes,
#'                                            dataset='hsapiens_gene_ensembl')
#' ## get transfer table !!!
#' draw.bubblePlot(driver_list=rownames(sig_driver),
#'                show_label=ms_tab[rownames(sig_driver),'gene_label'],
#'                Z_val=ms_tab[rownames(sig_driver),'Z.G4.Vs.others_DA'],
#'                driver_type=ms_tab[rownames(sig_driver),'gene_biotype'],
#'                target_list=analysis.par$merge.network$target_list,
#'                transfer2symbol2type=transfer_tab,
#'                min_gs_size=5,
#'                max_gs_size=500,use_gs=c('H'),
#'                top_geneset_number=5,top_driver_number=5,
#'                main='Bubbleplot for top driver targets',
#'                gs_cex = 0.4,driver_cex = 0.5)
#'  ## the cex is set just in case of figure margin too large,
#'  ## in real case, user could set cex larger or input pdf file name
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' use_genes <- base::unique(analysis.par$merge.network$network_dat$target.symbol)
#' transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',
#'                                            use_genes=use_genes,
#'                                            dataset='hsapiens_gene_ensembl')
#' ## get transfer table !!!
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' mark_gene <- c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D') ## marker for Group4
#' draw.bubblePlot(driver_list=rownames(sig_driver),
#'                show_label=ms_tab[rownames(sig_driver),'gene_label'],
#'                Z_val=ms_tab[rownames(sig_driver),'Z.G4.Vs.others_DA'],
#'                driver_type=ms_tab[rownames(sig_driver),'gene_biotype'],
#'                target_list=analysis.par$merge.network$target_list,
#'                transfer2symbol2type=transfer_tab,
#'                min_gs_size=5,max_gs_size=500,
#'                use_gs=use_gs=c('CP:KEGG','CP:BIOCARTA','H'),
#'                top_geneset_number=30,top_driver_number=50,
#'                pdf_file = sprintf('%s/bubbledraw.pdf',
#'                analysis.par$out.dir.PLOT),
#'                main='Bubbleplot for top driver targets',
#'                mark_gene=ms_tab[which(ms_tab$geneSymbol %in% mark_gene),
#'                'originalID_label'])
#' }
#' @export
draw.bubblePlot <- function(driver_list=NULL,show_label=driver_list,Z_val=NULL,driver_type=NULL,
                            target_list=NULL,transfer2symbol2type=NULL,
                            bg_list=NULL,min_gs_size=5,max_gs_size=500,
                            gs2gene=NULL,use_gs=NULL,
                            display_gs_list=NULL,
                            Pv_adj='none',Pv_thre=0.1,
                            top_geneset_number=30,top_driver_number=30,
                            pdf_file=NULL,main="",mark_gene=NULL,driver_cex=1,gs_cex=1,only_return_mat=FALSE){
  #
  all_input_para <- c('driver_list','show_label','Z_val','target_list','transfer2symbol2type',
                      'min_gs_size','max_gs_size','Pv_adj','Pv_thre','top_geneset_number','top_driver_number',
                      'driver_cex','gs_cex','only_return_mat')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('only_return_mat',c(TRUE,FALSE),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(names(show_label))==TRUE){names(show_label) <- driver_list}
  if(is.null(names(Z_val))==TRUE){names(Z_val) <- driver_list}
  if(is.null(driver_type)==FALSE){
    if(is.null(names(driver_type))==TRUE){names(driver_type) <- driver_list}
  }
  driver_list <- driver_list[order(abs(Z_val),decreasing = TRUE)]
  if(base::length(driver_list)>top_driver_number){
    driver_list <- driver_list[1:top_driver_number]
  }
  driver_list <- driver_list[order(Z_val[driver_list],decreasing=TRUE)]
  ## get target gene for driver_list
  transfer_tab <- transfer2symbol2type
  if(base::length(base::intersect(c('gene_biotype'),colnames(transfer_tab)))==0 & base::length(base::intersect(c('transcript_biotype'),colnames(transfer_tab)))==0){
    message('Input transfer table must contain the biotype information, please use get_IDtransfer2symbol2type() function to generate it!');return(FALSE);
  }
  #rownames(transfer_tab) <- transfer_tab[,1]
  target_gene <- lapply(driver_list,function(x){
    x1 <- target_list[[x]]$target
    w1 <- which.max(unlist(lapply(transfer_tab,function(x)length(intersect(x,x1)))))
    w2 <- which(colnames(transfer_tab) %in% c('gene_biotype','transcript_biotype'))[1]
    x1 <- x1[which(x1 %in% transfer_tab[,w1])]
    x2 <- transfer_tab[which(transfer_tab[,w1] %in% x1),]
    x3 <- x2[which(x2[,w2]=='protein_coding'),]
    target <- base::unique(x2[,'external_gene_name'])
    target_pc <- base::unique(x3[,'external_gene_name'])
    return(list(target,target_pc))
  })
  names(target_gene) <- driver_list
  ##
  f_res <- lapply(target_gene,function(x){
    funcEnrich.Fisher(input_list=x[[1]],bg_list=bg_list,gs2gene=gs2gene,use_gs=use_gs,
                      min_gs_size=min_gs_size,max_gs_size=max_gs_size,
                      Pv_adj='none',Pv_thre=Pv_thre)
  })
  names(f_res) <- names(target_gene)
  ## get display matrix
  all_path <- base::unique(unlist(lapply(f_res,function(x){x[[1]]})))
  all_path <- all_path[which(is.na(all_path)==FALSE)] ## get all sig path
  if(is.null(display_gs_list)==FALSE){
    all_path <- base::intersect(display_gs_list,all_path)
    if(base::length(all_path)<3){
      message('The number for passed gene sets is smaller than 3, please check the display_gs_list and re-try!')
      return(FALSE)
    }
  }
  f_mat <- lapply(f_res,function(x){
    as.data.frame(x)[all_path,5:6]
  })
  f_mat2 <- do.call(rbind,lapply(f_mat,function(x)-qnorm(x[[2]])))
  #  print(do.call(rbind,lapply(f_mat,function(x)x[[2]])))
  f_mat2[which(is.na(f_mat2)==TRUE | f_mat2==-Inf)] <- 0
  f_mat1 <- do.call(rbind,lapply(f_mat,function(x)x[,1]))
  colnames(f_mat1) <- all_path
  colnames(f_mat2) <- all_path
  f_mat3 <- t(apply(f_mat2,1,function(x){
    o1 <- order(abs(x),decreasing = TRUE)[1:3];
    x1<-x;x1[base::setdiff(1:base::length(x1),o1)] <- 0;x1
  }))
  ## use top
  min_path_num <- 5
  max_path_num <- top_geneset_number
  all_path_order <- apply(f_mat3,2,max)
  all_path_order <- sort(all_path_order,decreasing = TRUE)
  if(base::length(all_path_order) > max_path_num){
    w1 <- 1:base::length(all_path_order)
    if(base::length(w1)>=min_path_num & base::length(w1)<=max_path_num){
      all_path <- names(sort(all_path_order[w1],decreasing=TRUE))
    }else{
      if(base::length(w1)<min_path_num){
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
  if(only_return_mat==TRUE) return(f_mat2) ####
  plot_part <- function(ori=FALSE,before_off=FALSE){
    gsWidth  <- base::max(strwidthMod(colnames(f_mat1),'inches',cex=gs_cex,ori=ori))+par.char2inch()[1]
    gsHeight <- base::max(strheightMod(colnames(f_mat1),'inches',cex=gs_cex)*nrow(f_mat1))
    driverWidth  <- base::max(strwidthMod(show_label[rownames(f_mat1)],'inches',cex=driver_cex,ori=ori))+par.char2inch()[2]
    driverHeight <- base::max(strheightMod(show_label[rownames(f_mat1)],'inches',cex=driver_cex)*ncol(f_mat1))
    if(is.null(driver_type)==FALSE){
      rw <- base::max(strwidthMod(driver_type[driver_list],'inches',cex=gs_cex,ori=ori))+6*par.char2inch()[1]
    }else{
      rw <- 15*par.char2inch()[1]
    }
    gsWidth <- base::max(gsWidth,+par.char2inch()[1]*10+strwidthMod('target_size\n(protein_coding)','inches',cex=0.8,ori=ori))
    ## output to pdf
    # width:gsWidth|nc*0.5|rw
    # height:driverWidth|2*0.5|(nr)*0.5|0.5|1
    ww <- gsWidth + nc*0.5 +rw
    hh <- driverWidth + 2*0.5 + nr*0.5 +0.5*1.5+1
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=ww,height=hh)
    graphics::layout(1);par(mai=c(driverWidth,gsWidth,1,rw))
    graphics::plot(1,bty='n',col='white',xlim=c(0,nc),ylim=c(-2,nr+1.5),xaxt='n',yaxt='n',xlab="",ylab="",main=main,xaxs='i',yaxs='i')
    graphics::segments(x0=0,x1=nc,y0=0:nr,y1=0:nr,col='dark grey',xpd=TRUE)
    graphics::segments(x0=0:nc,x1=0:nc,y0=0,y1=nr,col='dark grey',xpd=TRUE)
    graphics::segments(x0=0:nc,x1=0:nc,y0=0,y1=-2,col='grey',xpd=TRUE)
    #
    pp <- par()$usr
    graphics::text(pp[1]-par.char2pos()[1],1:nr-0.5,colnames(f_mat1),xpd=TRUE,srt=0,adj=1,cex=gs_cex) ## sig pathways
    if(is.null(mark_gene)==TRUE){
      graphics::text(1:nc-0.5,-2-par.char2pos()[2],show_label[rownames(f_mat1)],xpd=TRUE,srt=90,adj=1,cex=driver_cex) ## sig regulators
    }else{
      bc <- rep('black',length.out=base::length(rownames(f_mat1)))
      bc[which(rownames(f_mat1) %in% mark_gene)] <- 'red'
      print(base::table(bc))
      graphics::text(1:nc-0.5,-2-par.char2pos()[2],show_label[rownames(f_mat1)],xpd=TRUE,srt=90,adj=1,col=bc,cex=driver_cex) ## sig regulators
    }
    ## draw circle
    max_size <- base::max(f_mat1,na.rm=TRUE)
    f_mat1 <- f_mat1/max_size
    cc_r <- matrix(z2col(f_mat2,n_len=30,sig_thre=qnorm(1-0.1)),ncol=ncol(f_mat2),byrow = FALSE)
    for(i in 1:nrow(f_mat1)){
      for(j in 1:ncol(f_mat1)){
        plotrix::draw.circle(i-0.5,j-0.5,radius=f_mat1[i,j]/2,col=cc_r[i,j])
      }
    }
    ## draw circle legend
    legend_size <- base::unique(round(base::seq(1,max_size,length.out=base::min(5,-1+nrow(f_mat1)))))
    for(i in 1:base::length(legend_size)){
      draw.circle(base::length(legend_size)-i+1.5,nr+0.5,radius=0.5*legend_size[i]/max_size)
      graphics::text(base::length(legend_size)-i+1.5,nr+1,legend_size[i],xpd=TRUE,pos=3)
    }
    graphics::text(0.5,nr+0.5,'Size ')
    ## draw p-value legend
    pp <- par()$usr
    p_label <- c(1,0.1,0.05,0.01,0.001,1e-4,1e-10)
    p_thre <- qnorm(1-c(1,0.1,0.05,0.01,0.001,0.0001,1e-10))
    p_col  <- z2col(p_thre,n_len=30,sig_thre=qnorm(1-0.1))
    p_col_m  <- z2col(-p_thre,n_len=30,sig_thre=qnorm(1-0.1))
    dd <- par.char2pos()[2]*1.5
    ybottom <- base::seq(nr-6*dd,nr-dd,length.out=1+base::length(p_thre))[1:base::length(p_thre)]
    ytop <- base::seq(nr-6*dd,nr-dd,length.out=1+base::length(p_thre))[2:(base::length(p_thre)+1)]
    for(i in 1:base::length(p_thre)){
      graphics::rect(pp[2]+par.char2pos()[1]*1,ybottom[i],pp[2]+par.char2pos()[1]*3,ytop[i],col=p_col[i],xpd=TRUE)
      graphics::text(pp[2]+par.char2pos()[1]*3.5,(ybottom[i]+ytop[i])/2,pos=4,p_label[i],xpd=TRUE)
    }
    ybottom <- base::seq(nr-11*dd,nr-6*dd,length.out=1+base::length(p_thre))[1:base::length(p_thre)]
    ytop <- base::seq(nr-11*dd,nr-6*dd,length.out=1+base::length(p_thre))[2:(base::length(p_thre)+1)]
    for(i in 2:base::length(p_thre)){
      graphics::rect(pp[2]+par.char2pos()[1]*1,ybottom[i],pp[2]+par.char2pos()[1]*3,ytop[i],col=rev(p_col_m)[i-1],xpd=TRUE)
      graphics::text(pp[2]+par.char2pos()[1]*3.5,(ybottom[i]+ytop[i])/2,pos=4,rev(p_label)[i-1],xpd=TRUE)
    }
    graphics::text(pp[2]+par.char2pos()[1]*2,nr-0.5*dd,'P-Value',xpd=TRUE,adj=0)
    # target size
    max_target_size <- base::max(unlist(lapply(target_gene,function(x)base::max(unlist(lapply(x,length))))))
    ori_size <- unlist(lapply(target_gene,function(x)base::length(x[[1]])))
    pro_size <- unlist(lapply(target_gene,function(x)base::length(x[[2]])))
    graphics::rect(0:(nc-1)+0.15,-2,0:(nc-1)+0.45,1.5*ori_size/max_target_size-2,col='blue')
    graphics::rect(0:(nc-1)+0.55,-2,0:(nc-1)+0.85,1.5*pro_size/max_target_size-2,col='green')
    graphics::segments(y0=-2,y1=-0.5,x0=0,x1=0,xpd=TRUE)
    graphics::segments(y0=seq(-2,-0.5,length.out=3),y1=seq(-2,-0.5,length.out=3),x0=-0.25,x1=0,xpd=TRUE)
    text(-0.3,seq(-2,-0.5,length.out=3),round(base::seq(0,max_target_size,length.out=3)),adj=1,xpd=TRUE)
    legend(-5,-0.5,c('target_size','target_size\n(protein_coding)'),fill=c('blue','green'),border=NA,bty='n',cex=0.8,xpd=TRUE)
    # add sig color
    sig_col <- z2col(Z_val[driver_list],n_len=30)
    graphics::rect(0:(nc-1)+0.35,-0.4,0:(nc-1)+0.65,-0.1,col=sig_col)
    # add driver type !!!
    if(is.null(driver_type)==FALSE){
      cc_tmp <- get.class.color(base::unique(driver_type[driver_list]))
      graphics::points(x=0:(nc-1)+0.5,y=rep(-2-0.75*par.char2pos()[1],length.out=nc),col=cc_tmp[driver_type[driver_list]],pch=16,cex=2,xpd=TRUE)
      #graphics::legend(pp[2],0.5,names(cc_tmp),fill=cc_tmp,border=NA,bty='n',cex=0.8,xpd=TRUE)
      #print(base::seq(from=0.5,by= -par.char2pos()[2]*1.5,length.out=base::length(cc_tmp)))
      graphics::points(x=pp[2]+par.char2pos()[1],y=base::seq(from=-2,by= par.char2pos()[2]*1.5,length.out=base::length(cc_tmp)),col=cc_tmp,cex=1.5,xpd=TRUE,pch=16)
      graphics::text(x=pp[2]+par.char2pos()[1]*2,y=base::seq(from=-2,by= par.char2pos()[2]*1.5,length.out=base::length(cc_tmp)),names(cc_tmp),pos=4,xpd=TRUE)
    }
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}
  return(TRUE)
}


#' Gene Set Enrichment Analysis by GSEA
#'
#' \code{funcEnrich.GSEA} performs gene set enrichment analysis to the input gene list, by using the GSEA Test.
#'
#' @param rank_profile a named vector of numerics, the differential values (DE or DA) calculated from a sample comparison (e.g. "G4 vs. Others").
#' Names of the vector must be gene names.
#' For the DA, user could use 'processDriverProfile()' to convert the DA profile into gene-name based profile.
#' The differential values can be "logFC" or "t-statistics".
#' @param gs2gene list, a list contains elements of gene sets.
#' The name of the element is gene set, each element contains a vector of genes in that gene set.
#' If NULL, will use \code{all_gs2gene}, which is created by function \code{gs.preload}. Default is NULL.
#' @param use_gs a vector of characters, the names of gene sets.
#' If \code{gs2gene} is NULL, \code{all_gs2gene} will be used. The \code{use_gs} must be the subset of \code{names(all_gs2gene)}.
#' If "all", all the gene sets in \code{gs2gene} will be used.
#' If user input his own \code{gs2gene} list, \code{use_gs} will be set to "all" as default.
#' Default is c("H", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG").
#' @param min_gs_size numeric, the minimum size of gene set to analysis. Default is 5.
#' @param max_gs_size numeric, the maximum size of gene set to analysis. Default is 500.
#' @param Pv_adj character, method to adjust P-value. Default is "fdr".
#' For details, please check \code{p.adjust.methods}.
#' @param Pv_thre numeric, threshold for the adjusted P-values. Default is 0.1.
#' @param test_strategy choose from "KS" and "GSEA". Default is "GSEA".
#' If "KS", will perform a Kolmogorov-Smirnov test to get the significance value.
#' @param nperm numeric, number of random permutations. Default is 1000.
#' This function only do gene-label based permutation reshuffling.
#' @param use_seed integer, the random seed. Default is 999.
#'
#' @return Return a data.frame, contains gene sets with significant enrichment statistics.
#' Column details are as follows (test_strategy=GSEA),
#'
#' \item{#Name}{Name of the enriched gene set}
#' \item{Total_item}{Size in the profile}
#' \item{Num_item}{Number of genes in the gene set (filtered by the profile list)}
#' \item{Ori_P}{Original P-value from GSEA Test}
#' \item{Adj_P}{Adjusted P-value}
#' \item{ES}{Enrichment Score}
#' \item{NES}{normalized Enrichment Score}
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' ## get significant gene set by driver's DA profile
#' DA_profile <- processDriverProfile(Driver_name=ms_tab$gene_label,
#'                                     Driver_profile=ms_tab$logFC.G4.Vs.others_DA,
#'                                     choose_strategy='absmax',
#'                                     return_type ='gene_statistics')
#' res1 <- funcEnrich.GSEA(rank_profile=DA_profile,
#'                          use_gs=c('H'),
#'                          Pv_thre=0.1,Pv_adj = 'none')
#' \dontrun{
#' }
#' @export
funcEnrich.GSEA <- function(rank_profile=NULL,
                            use_gs=NULL,
                            gs2gene=NULL,
                            min_gs_size=5,max_gs_size=500,
                            Pv_adj='fdr',Pv_thre=0.1,
                            test_strategy='GSEA',
                            nperm=1000,use_seed=999){
  #
  all_input_para <- c('rank_profile','min_gs_size','max_gs_size','Pv_adj','Pv_thre',
                      'test_strategy','nperm','use_seed')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('Pv_adj',c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('test_strategy',c("GSEA","KS"),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}

  #
  if(is.null(gs2gene)==TRUE){ ## use inner gs2gene
    if(is.null(use_gs)==TRUE){
      use_gs <- c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG')
    }else{
      if(use_gs[1] == 'all'){
        use_gs <- c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)
      }
    }
    if(base::length(base::setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)))>0){
      message(sprintf('Input %s not in all_gs2gene, please check all_gs2gene_info (items in Category or Sub-Category) and re-try!',
                      base::paste(base::setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)),collapse=';')));
      return(FALSE)
    }
    if(base::length(use_gs)>1){
      gs2gene <- merge_gs(all_gs2gene,use_gs = use_gs)
    }else{
      gs2gene <- all_gs2gene[[use_gs]]
    }
  }else{
    if(is.null(use_gs)==TRUE){
      use_gs <- 'all'
    }
    if(base::length(use_gs)>1){
      if(class(gs2gene[[1]])=='list') gs2gene <- merge_gs(gs2gene,use_gs = use_gs)
    }else{
      if(use_gs == 'all'){
        if(class(gs2gene[[1]])=='list') gs2gene <- merge_gs(gs2gene,use_gs = names(gs2gene))
      }else{
        if(class(gs2gene[[1]])=='list') gs2gene <- gs2gene[[use_gs]]
      }
    }
  }
  all_gs <- names(gs2gene)
  bg_list <- names(rank_profile)
  if(!is.null(bg_list)){
    use_gs2gene <- lapply(gs2gene,function(x){base::intersect(x,bg_list)})
    names(use_gs2gene) <- names(gs2gene)
  }else{
    use_gs2gene <- gs2gene
  }
  bg_list <- base::unique(unlist(use_gs2gene))
  ## size selection
  s1 <- unlist(lapply(use_gs2gene,length))
  w1 <- which(s1>=min_gs_size & s1<=max_gs_size)
  use_gs2gene <- use_gs2gene[w1]
  all_gs <- names(use_gs2gene) ## all tested gene set number

  ## input filter
  empty_vec <- as.data.frame(matrix(NA,ncol=7));
  colnames(empty_vec) <- c('#Name','Total_item','Num_item','Ori_P','Adj_P','ES','NES')
  if(length(bg_list)==0) return(empty_vec)
  ## GSEA
  pf1 <- sort(rank_profile,decreasing = T)
  if(test_strategy=='GSEA'){
    set.seed(seed = use_seed, kind = NULL)
    pv <- lapply(names(use_gs2gene),function(x0){
      message(sprintf('Calculate permutation for %s',x0))
      x <- use_gs2gene[[x0]]
      es <- get_ES(pf1,x)$ES ## ES
      es.p <- unlist(lapply(1:nperm,function(x1){
        get_ES(pf1,sample(names(pf1),size = length(x)))$ES
      }))
      if(es>0){
        es.p.pos <- es.p[which(es.p>0)]
        nes <- es/mean(es.p.pos)
        pv  <- length(which(es.p.pos>=es))/length(es.p.pos)
      }else{
        es.p.neg <- es.p[which(es.p<0)]
        nes <- es/abs(mean(es.p.neg))
        pv  <- length(which(es.p.neg<=es))/length(es.p.neg)
      }
      return(c(length(x),pv,es,nes))
    })
    pv <- do.call(rbind,pv)
    pv <- as.data.frame(pv,stringsAsFactors=FALSE)
    colnames(pv) <-  c('Num_item','Ori_P','ES','NES')
    rownames(pv) <- names(use_gs2gene)
    pv$Adj_P <- p.adjust(pv$Ori_P,method=Pv_adj,n=base::length(all_gs))
    pv$`#Name` <- rownames(pv)
    pv$Total_item <- length(pf1)
    pv <- pv[,c('#Name','Total_item','Num_item','Ori_P','Adj_P','ES','NES')]
    pv <- pv[order(pv$Ori_P),]
    use_pv <- pv[which(pv$Adj_P<=Pv_thre),]
  }
  if(test_strategy=='KS'){
    pv <- lapply(names(use_gs2gene),function(x0){
      message(sprintf('Calculate KS for %s',x0))
      x <- use_gs2gene[[x0]]
      res <- ks.test(pf1[x],pf1)
      return(c(res$statistic,res$p.value,length(x)))
    })
    pv <- do.call(rbind,pv)
    pv <- as.data.frame(pv,stringsAsFactors=FALSE)
    rownames(pv) <- names(use_gs2gene)
    colnames(pv) <-  c('D-statistics','Ori_P','Num_item')
    pv$Adj_P <- p.adjust(pv$Ori_P,method=Pv_adj,n=base::length(all_gs))
    pv$`#Name` <- rownames(pv)
    pv$Total_item <- length(pf1)
    pv <- pv[,c('#Name','Total_item','Num_item','Ori_P','Adj_P','D-statistics')]
    pv <- pv[order(pv$Ori_P),]
    use_pv <- pv[which(pv$Adj_P<=Pv_thre),]
  }
  return(use_pv)
}

#' GSEA (Gene Set Enrichment Analysis) Plot for one Gene Set or one Driver
#'
#' \code{draw.GSEA} draws a GSEA plot to analyze one gene set (with gene list annotated) or one driver (with list of target genes).
#'
#'
#' @param rank_profile a named vector of numerics, the differential values (DE or DA) calculated from a sample comparison (e.g. "G4 vs. Others").
#' Names of the vector must be gene names.
#' For the DA, user could use `processDriverProfile()` to convert the DA profile into gene-name based profile.
#' The differential values can be "logFC" or "t-statistics".
#' @param use_genes a vector of characters, a vector of genes to display. The genes can either be annotated genes in gene set or the targe genes from a specific driver.
#' The gene names must be a subset of \code{names(rank_profile)}.
#' @param use_direction a vector of numeric 1s and -1s, 1 is positive regulation from driver, -1 is negative regulation from driver.
#' Users can get this vector by converting the signs of "spearman". If NULL, no regulation direction will be displayed. Default is NULL.
#' @param main character, an overall title for the plot. Default is "".
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#' @param annotation character, the annotation set by users for easier reference.
#' Normally the annotation is the P-value or other statistics to show the significance of the interested gene set or driver.
#' If NULL, will perform a Kolmogorov-Smirnov test to get the significance value.
#' If want to get the statistics from GSEA test, could use `funcEnrich.GSEA()` to get the statistics first.
#' Default is NULL.
#' @param annotation_cex numeric, giving the amount by which the text of annotation should be magnified relative to the default. Default is 1.2.
#' @param left_annotation character, annotation displayed on the left of the figure, representing left condition of the \code{rank_profile}. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure, representing right condition of the \code{rank_profile}. Default is "".
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#'
#' @examples
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#'
#' ## draw for the most significant gene set
#' # by driver's DA profile
#' DA_profile <- processDriverProfile(Driver_name=ms_tab$gene_label,
#'                                     Driver_profile=ms_tab$logFC.G4.Vs.others_DA,
#'                                     choose_strategy='absmax',
#'                                     return_type ='gene_statistics')
#' res1 <- funcEnrich.GSEA(rank_profile=DA_profile,
#'                          use_gs=c('H'),
#'                          Pv_thre=0.1,Pv_adj = 'none')
#' top_gs <- res1[1,'#Name'] ## draw for the top 1
#' annot <- sprintf('NES: %s \nAdjusted P-value: %s',
#'           signif(res1[1,'NES'],2),
#'           signif(res1[1,'Adj_P'],2))
#' draw.GSEA(rank_profile=DA_profile,
#'           use_genes=all_gs2gene$H[[top_gs]],
#'           main=sprintf('GSEA plot for gene set %s',
#'           top_gs),
#'           annotation=annot,annotation_cex=1.2,
#'           left_annotation='high in G4',
#'           right_annotation='high in others')
#'
#' ## draw for the most significant driver
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' DE_profile <- analysis.par$DE[[1]]$`Z-statistics`;
#' names(DE_profile) <- rownames(analysis.par$DE[[1]])
#' use_driver <- driver_list[1]
#' use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
#' use_target_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman) ## 1/-1
#' annot <- sprintf('P-value: %s',signif(ms_tab[use_driver,'P.Value.G4.Vs.others_DA'],2))
#'
#' ## draw for the driver
#' draw.GSEA(rank_profile=DE_profile,
#'           use_genes=use_target_genes,
#'           use_direction=use_target_direction,
#'           main=sprintf('GSEA plot for driver %s',
#'           ms_tab[use_driver,'gene_label']),
#'           annotation=annot,annotation_cex=1.2,
#'           left_annotation='high in G4',
#'           right_annotation='high in others')
#' draw.GSEA(rank_profile=DE_profile,
#'           use_genes=use_target_genes,
#'           use_direction=NULL,
#'           main=sprintf('GSEA plot for driver %s',
#'           ms_tab[use_driver,'gene_label']),
#'           annotation=NULL,annotation_cex=1.2,
#'           left_annotation='high in G4',
#'           right_annotation='high in others')
#'
#' ## draw for the gene set
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' use_target_genes <- all_gs2gene[[1]][[1]]
#' draw.GSEA(rank_profile=DE_profile,
#'          use_genes=use_target_genes,
#'          main=sprintf('GSEA plot for %s',names(all_gs2gene[[1]][1])),
#'          left_annotation='high in G4',
#'          right_annotation='high in others')
#'
#' \dontrun{
#' #' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' DE_profile <- analysis.par$DE[[1]]$`Z-statistics`;
#' names(DE_profile) <- rownames(analysis.par$DE[[1]])
#' use_driver <- driver_list[1]
#' use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
#' use_target_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman) ## 1/-1
#' annot <- sprintf('P-value: %s',signif(ms_tab[use_driver,'P.Value.G4.Vs.others_DA'],2))
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,
#'           use_direction=use_target_direction,
#'           main=sprintf('GSEA plot for driver %s',ms_tab[use_driver,'gene_label']),
#'           pdf_file = sprintf('%s/GSEA_driver.pdf',
#'           analysis.par$out.dir.PLOT),
#'           annotation=annot,annotation_cex=1.2,
#'           left_annotation='high in G4',
#'           right_annotation='high in others')
#'
#' ## draw for the gene set
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' use_target_genes <- all_gs2gene[[1]][[1]]
#' draw.GSEA(rank_profile=DE_profile,
#'          use_genes=use_target_genes,
#'          main=sprintf('GSEA plot for %s',names(all_gs2gene[[1]][1])),
#'          pdf_file = sprintf('%s/GSEA_GS_each.pdf',analysis.par$out.dir.PLOT),
#'          left_annotation='high in G4',
#'          right_annotation='high in others')
#'}
#' @export
draw.GSEA <- function(rank_profile=NULL,use_genes=NULL,use_direction=NULL,main="",pdf_file=NULL,
                      annotation=NULL,annotation_cex=1.2,left_annotation=NULL,right_annotation=NULL){
  #
  all_input_para <- c('rank_profile','use_genes','main','annotation_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  #### start plot
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  if(is.null(use_direction)==FALSE){
    w1 <- which(use_genes %in% names(rank_profile))
    w2 <- unique(use_direction[w1])
    if(length(w2)==1){
      if(w2==1) use_direction <- NULL;
    }
  }
  if(is.null(use_direction)==FALSE){
    new_rank_profile <- c(rank_profile,-rank_profile)
    names(new_rank_profile) <- c(paste0('POS_',names(rank_profile)),paste0('NEG_',names(rank_profile)))
    new_use_genes <- use_genes
    pos_genes <- new_use_genes[which(use_direction==1)]
    neg_genes <- new_use_genes[which(use_direction==-1)]
    new_use_genes[which(use_direction==1)] <- paste0('POS_',new_use_genes[which(use_direction==1)])
    new_use_genes[which(use_direction==-1)] <- paste0('NEG_',new_use_genes[which(use_direction==-1)])
  }
  rank_profile <- sort(rank_profile,decreasing = TRUE)
  r_len <- base::length(rank_profile)
  use_pos <- which(names(rank_profile) %in% use_genes)
  plot_part <- function(ori=FALSE,before_off=FALSE){
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE){
      pdf(pdf_file,width=10,height=10)
    }
    # height: 4|3|2|1
    graphics::layout(matrix(c(rep(4,10),3,2,rep(1,10)),ncol=1))
    # plot
    ## rank for all, 1
    par(mar=c(10,6,0.5,2))
    mm <- base::max(abs(rank_profile))
    y1 <- seq(-mm,mm,length.out=7); y1 <- round(y1,1)
    unit <- r_len/10; unit <- round(unit/100)*100
    x1 <- base::seq(0,r_len,by=unit);x1 <- base::unique(x1); x1 <- c(x1,base::max(x1)+unit)
    par(usr=c(0,base::max(x1),-mm,mm))
    graphics::plot(rank_profile,col='grey',pch=16,xaxt='n',xlab="",ylab="",bty='l',type='n',ylim=c(-mm,mm),xaxs='i',las=1,yaxs='i')
    graphics::polygon(x=c(0,1:r_len,r_len),y=c(0,rank_profile,0),col='grey',border=NA)
    if(is.null(left_annotation)==FALSE) graphics::text(0,mm-par.char2pos()[2],pos=4,left_annotation,col='red',xpd=TRUE,cex=annotation_cex)
    if(is.null(right_annotation)==FALSE) graphics::text(r_len,-mm+par.char2pos()[2],pos=2,right_annotation,col='blue',xpd=TRUE,cex=annotation_cex)
    pp <- par()$usr
    graphics::mtext(side=2,line = 3,'Ranked list metric (PreRanked)',cex=annotation_cex)
    graphics::mtext(side=1,line = 3.5,'Rank in Ordered Dataset',cex=annotation_cex)
    x1 <- x1[which(x1<base::length(rank_profile))]
    x1[base::length(x1)] <- base::length(rank_profile)
    graphics::segments(x0=x1,x1=x1,y0=pp[3]-par.char2pos()[2]/5,y1=pp[3],xpd=TRUE)
    graphics::text(x1,pp[3]-par.char2pos()[2]/2,get_label_manual(x1),adj=1,xpd=TRUE)
    # get zero cross
    w1 <- base::which.min(abs(rank_profile))
    graphics::abline(v=w1,lty=2,col='grey')
    graphics::text(w1,-mm/4,sprintf('Zero cross at %d',w1),adj=0.5)
    if(is.null(use_direction)==FALSE){
      graphics::legend(pp[1]/2+pp[2]/2,pp[3]-(pp[4]-pp[3])/4,lty=1,lwd=2,c('Enrichment profile/Hits (positive)','Enrichment profile/Hits (negative)','Ranking metric scores'),
             col=c(pos_col,neg_col,'grey'),xpd=TRUE,horiz=TRUE,xjust=0.5,cex=annotation_cex,bty='n')
    }else{
      graphics::legend(pp[1]/2+pp[2]/2,pp[3]-(pp[4]-pp[3])/4,lty=1,lwd=2,c('Enrichment profile','Hits','Ranking metric scores'),
             col=c('green','black','grey'),xpd=TRUE,horiz=TRUE,xjust=0.5,cex=annotation_cex,bty='n')
    }
    pm <- par()$usr
    ## get image bar, 2
    par(mar=c(0,6,0,2))
    use_col <- z2col(rank_profile,sig_thre = 0,n_len = 100,blue_col='blue',red_col='red')
    graphics::image(x=as.matrix(1:r_len),col=use_col,bty='n',xaxt='n',yaxt='n',xlim=c(pm[1],pm[2])/r_len)
    graphics::abline(v=use_pos/r_len,col='grey')
    pp <- par()$usr
    graphics::text(0,pp[3]/2+pp[4]/2,pos=2,sprintf('Size:%s',base::length(use_pos)),xpd=TRUE)
    ## mark gene position; 3
    par(mar=c(0,6,0,2))
    graphics::plot(1,col='white',xlab="",ylab="",bty='n',xlim=c(1,r_len),ylim=c(0,1),xaxt='n',yaxt='n',xaxs='i')
    if(is.null(use_direction)==FALSE){
      use_pos_P <- which(names(rank_profile) %in% pos_genes)
      use_pos_N <- which(names(rank_profile) %in% neg_genes)
      if(length(use_pos_P)>0) graphics::segments(y0=1,y1=0.5,x0=use_pos_P,x1=use_pos_P,col=pos_col)
      if(length(use_pos_N)>0) graphics::segments(y0=0,y1=0.5,x0=use_pos_N,x1=use_pos_N,col=neg_col)
      graphics::abline(h=0.5,col='light grey')
      graphics::text(0,0.75,pos=2,sprintf('Pos_Size:%s',base::length(use_pos_P)),xpd=TRUE)
      graphics::text(0,0.25,pos=2,sprintf('Neg_Size:%s',base::length(use_pos_N)),xpd=TRUE)
    }else{
      graphics::abline(v=use_pos)
    }
    ## GSEA ES, 4
    par(mar=c(0,6,5,2))
    # get ES score
    es <- ''
    if(is.null(use_direction)==FALSE){
      es_res_pos <- get_ES(rank_profile,pos_genes)
      es_res_neg <- get_ES(rank_profile,neg_genes)
      es_res <- es_res_pos
      #print(es_res_pos$RES);print(es_res_neg$RES)
      y2 <- base::seq(base::min(c(es_res_pos$RES,es_res_neg$RES),na.rm=T),base::max(c(es_res_pos$RES,es_res_neg$RES),na.rm=T),length.out=7); y2 <- round(y2,1)
      graphics::plot(es_res_pos$RES,col=pos_col,xaxt='n',xlab="",ylab="",bty='n',
           xlim=c(1,r_len),type='l',lwd=3,ylim=c(base::min(c(es_res_pos$RES,es_res_neg$RES),na.rm=T),base::max(y2)),main=main,xpd=TRUE,xaxs='i',las=1,cex.main=annotation_cex)
      lines(es_res_neg$RES,col=neg_col,lwd=3,xpd=TRUE)
      w1 <- base::which.max(abs(es_res_pos$RES));
      if(length(w1)==1) graphics::segments(x0=w1,x1=w1,y0=0,y1=es_res_pos$RES[w1],lty=2,col='grey')
      w1 <- base::which.max(abs(es_res_neg$RES));
      if(length(w1)==1) graphics::segments(x0=w1,x1=w1,y0=0,y1=es_res_neg$RES[w1],lty=2,col='grey')
    }else{
      es_res <- get_ES(rank_profile,use_genes)
      y2 <- base::seq(base::min(es_res$RES),base::max(es_res$RES),length.out=7); y2 <- round(y2,1)
      graphics::plot(es_res$RES,col='green',xaxt='n',xlab="",ylab="",bty='n',
           xlim=c(1,r_len),type='l',lwd=3,ylim=c(base::min(es_res$RES),base::max(y2)),main=main,xpd=TRUE,xaxs='i',las=1)
      w1 <- base::which.max(abs(es_res$RES));
      graphics::segments(x0=w1,x1=w1,y0=0,y1=es_res$RES[w1],lty=2,col='grey')
      es <- sprintf('ES: %s',format(es_res$RES[w1],digits=3))
    }
    graphics::abline(h=0,lty=2,col='dark grey')
    pp <- par()$usr
    graphics::mtext(side=2,line = 3,'Enrichment score (ES)',cex=annotation_cex)
    # add annotation
    if(is.null(annotation)==TRUE){
      if(is.null(use_direction)==FALSE){
        pv <- ks.test(new_rank_profile,new_rank_profile[new_use_genes])$p.value
      }else{
        pv <- ks.test(rank_profile,rank_profile[use_genes])$p.value
      }
      #print(t.test(rank_profile,rank_profile[use_genes]))
      if(pv==0){
        pv <- '<2.2e-16'
      }else{
        if(pv<0.01) pv <- format(pv,digits = 3,scientific = TRUE) else pv <- signif(pv,3)
      }
      annotation <- sprintf("%s\nKS test p-value:%s",es,pv)
    }
    pp <- par()$usr
    if(es_res$RES[base::which.max(abs(es_res$RES))]>0)
      graphics::text(r_len,pp[4]-2*par.char2pos()[2],annotation,pos=2,cex=annotation_cex,xpd=TRUE)
    else
      graphics::text(0,pp[3]+2*par.char2pos()[2],annotation,pos=4,cex=annotation_cex,xpd=TRUE)
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);dev.off()} else {plot_part()}
  return(TRUE)
}

## get enrichment score
get_ES <- function(rank_profile=NULL,use_genes=NULL,weighted.score.type=1){
  gene.list <- names(rank_profile)
  correl.vector <- rank_profile
  tag.indicator <- sign(match(gene.list, use_genes, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- base::length(gene.list)
  Nh <- base::length(use_genes)
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
  max.ES <- base::max(RES)
  min.ES <- base::min(RES)
  if(is.na(max.ES)==TRUE) max.ES <- 0
  if(is.na(min.ES)==TRUE) min.ES <- 0
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- base::which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- base::which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}
##
get_z2p_each <- function(x,use_star=FALSE,digit_num=2,twosided=T){
  x <- abs(x)
  if(twosided==T) use_pv <- pnorm(x,lower.tail = F)*2 ## two-tail
  if(twosided==F) use_pv <- pnorm(x,lower.tail = F) ## one-tail
  x_star <- ''
  if(as.character(use_pv)!='0'){
    use_p <- format(use_pv,digits=digit_num,scientific = TRUE)
  }else{
    low_p <- .Machine$double.xmin
    low_z <- sapply(10^(-(1:(1+-log10(low_p)))),function(xx)combinePvalVector(xx,twosided = twosided))
    #low_z <- sapply(10^(-(1:(1+-log10(low_p)))),function(xx)ifelse(xx == 0, 0, combinePvalVector(xx,twosided = twosided)))
    use_pv <- low_z[2,which(low_z[1,]>=x)[1]]
    use_p <- format(use_pv, digits=3,scientific = TRUE)
    use_p[which(use_p=='NA')] <- '<1e-308'
    use_p <- as.character(use_p)
    x_star <- '***'
  }
  x_star[which(use_pv<0.05)] <-'*'
  x_star[which(use_pv<0.01)] <-'**'
  x_star[which(use_pv<0.001)] <-'***'
  if(use_star==TRUE) use_p<-paste0(use_p,x_star)
  return(use_p)
}
get_z2p <- function(x,use_star=FALSE,digit_num=2,twosided=T){
  x[which(is.na(x)==TRUE)] <- 0
  x <- abs(x)
  x[which(is.na(x)==TRUE)] <- 0 ##
  use_p <- unlist(lapply(x,function(xx)get_z2p_each(xx,use_star=use_star,digit_num=digit_num,twosided=twosided)))
  return(use_p)
}

#' Draw GSEA (gene set enrichment analysis) Plot with NetBID Analysis of Drivers
#'
#' \code{draw.GSEA.NetBID} creates a GSEA plot for drivers with more NetBID analysis information. Such as number of target genes, ranking of target genes in
#' differential expressed file, differential expression (DE) and differential activity (DA) values.
#'
#' @param DE data.frame, a data.frame created either by function \code{getDE.limma.2G} or \code{getDE.BID.2G}. Row names are gene/driver names,
#' columns must include gene/driver name and calculated differencial values (e.g. "ID", "logFC", "AveExpr", "P.Value" etc.).
#' @param name_col character, the name of the column in \code{DE} contains gene names. If NULL, will use the row names of \code{DE}.
#' Default is NULL.
#' @param profile_col character, the name of the column in \code{DE} contains calculated differencial value (e.g. "logFC" or "P.Value").
#' If \code{DE} is created by \code{getDE.limma.2G} or \code{getDE.BID.2G}, this parameter should be set to "logFC" or "t".
#' @param profile_trend character, users can choose between "pos2neg" and "neg2pos". "pos2neg" means high \code{profile_col} in target group will be shown on the left.
#' "neg2pos" means high \code{profile_col} in control group will be shown on the left. Default is "pos2neg".
#' For details, please check online tutorial.
#' @param driver_list a vector of characters, the names of top drivers.
#' @param show_label a vector of characters, the names of top drivers.
#' If NULL, will display the names in \code{driver_list}. Default is NULL.
#' @param driver_DA_Z a vector of numerics, the Z statistics of differential activity (DA) value of the \code{driver_list}.
#' It is highly suggested to give names to the vector, otherwise the names of \code{driver_list} will be used.
#' @param driver_DE_Z a vector of numerics, the Z statistics of differential expressed (DE) value of the \code{driver_list}.
#' It is highly suggested to give names to the vector, otherwise the names of \code{driver_list} will be used.
#' @param target_list list, the driver-to-target list object. The names of the list elements are drivers.
#' Each element is a data frame, usually contains at least three columns. "target", target gene names; "MI", mutual information; "spearman", spearman correlation coef- ficient.
#' Users can call \code{get_net2target_list} to create this list.
#' @param top_driver_number numeric, number for the top significant drivers to be displayed in the plot. Default is 30.
#' @param target_nrow numeric, users can choose between 1 and 2. Number of panels to mark the ranking of target genes.
#' If 1, the ranking of target genes will be marked in one panel.
#' If 2, the ranking of target genes will be marked in two panels. Upper panel for positively-regulated, lower panel for negatively-regulated.
#' Default is 2. For details, please check online tutorial.
#' @param target_col character, name of the color palette used for display marker line in the panel. Users can choose between "black" and "RdBu".
#' If "black", the marker line in the panel is black.
#' If "RdBu", the marker line in the panel is Red to Blue.
#' If \code{target_col_type} is set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If \code{target_col_type} is set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' with significant high set for red and low for blue. The significant threshold is set by \code{profile_sig_thre}.
#' Default is 'RdBu'.
#' @param target_col_type character, name of the color palette used for display target genes. This parameter works only when \code{target_col} is set as "RdBu".
#' Users can choose between "PN" and "DE".
#' If "PN", positively-regulated genes will be colored red and negatively-regulated genes will be colored blue.
#' If "DE", the color shades is decided by its differentiated value.
#' Default is "PN".
#' @param left_annotation character, annotation on the left of profile curve, indicating high in control group or target group.
#' Default is "".
#' @param right_annotation character, annotation on the right of profile curve, indicating high in the opposite group of \code{left_annotation}.
#' Default is "".
#' @param main character, an overall title for the plot. Default is "".
#' @param profile_sig_thre numeric, threshold value for target genes. This parameter works only when \code{target_col_type} is set as "DE" and \code{target_col} is set as "RdBu".
#' Non-significant target genes will be colored grey. Default is 0.
#' @param Z_sig_thre numeric, threshold value of Z statistics from \code{driver_DA_Z} and \code{driver_DE_Z}. Significant values will have background color.
#' Default is 1.64.
#' @param pdf_file character, the file path to save figure as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.

#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' comp <- 'G4.Vs.others'
#' DE <- analysis.par$DE[[comp]]
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' draw.GSEA.NetBID(DE=DE,profile_col='t',
#'                  name_col='ID',
#'                  profile_trend='neg2pos',
#'                  driver_list = driver_list,
#'                  show_label=ms_tab[driver_list,'gene_label'],
#'                  driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
#'                  driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
#'                  target_list=analysis.par$merge.network$target_list,
#'                  top_driver_number=5,
#'                  target_nrow=2,target_col='RdBu',
#'                  left_annotation = 'test_left',
#'                  right_annotation = 'test_right',
#'                  main='test',target_col_type='DE',
#'                  Z_sig_thre=1.64,
#'                  profile_sig_thre = 1.64)
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' comp <- 'G4.Vs.others'
#' DE <- analysis.par$DE[[comp]]
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',
#'                  driver_list = driver_list,
#'                  show_label=ms_tab[driver_list,'gene_label'],
#'                  driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
#'                  driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
#'                  target_list=analysis.par$merge.network$target_list,
#'                  top_driver_number=30,
#'                  target_nrow=2,
#'                  target_col='RdBu',
#'                  left_annotation = 'test_left',
#'                  right_annotation = 'test_right',
#'                  main='test',
#'                  target_col_type='DE',
#'                  Z_sig_thre=1.64,
#'                  profile_sig_thre = 1.64,
#'                  pdf_file=sprintf('%s/NetBID_GSEA_demo1.pdf',
#'                  analysis.par$out.dir.PLOT))
#'draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='neg2pos',
#'                  driver_list = driver_list,
#'                  show_label=ms_tab[driver_list,'gene_label'],
#'                  driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
#'                  driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
#'                  target_list=analysis.par$merge.network$target_list,
#'                  top_driver_number=30,
#'                  target_nrow=1,
#'                  target_col='RdBu',
#'                  left_annotation = 'test_left',
#'                  right_annotation = 'test_right',
#'                  main='test',target_col_type='PN',
#'                  Z_sig_thre=1.64,profile_sig_thre = 1.64,
#'                  pdf_file=sprintf('%s/NetBID_GSEA_demo2.pdf',
#'                  analysis.par$out.dir.PLOT))
#'}
#' @export
draw.GSEA.NetBID <- function(DE=NULL,name_col=NULL,profile_col=NULL,profile_trend='pos2neg',
                             driver_list=NULL,show_label=driver_list,driver_DA_Z=NULL,driver_DE_Z=NULL,target_list=NULL,
                             top_driver_number=30,target_nrow=2,target_col='RdBu',target_col_type='PN',
                             left_annotation="",right_annotation="",main="",
                             profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  #
  all_input_para <- c('DE','profile_col','profile_trend','driver_list','show_label',
                      'driver_DA_Z','driver_DE_Z','target_list','top_driver_number','target_nrow',
                      'target_col','target_col_type','left_annotation','right_annotation','main',
                      'profile_sig_thre','Z_sig_thre')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('profile_trend',c('pos2neg','neg2pos'),envir=environment()),
                 check_option('target_nrow',c(1,2),envir=environment()),
                 check_option('target_col',c('RdBu','black'),envir=environment()),
                 check_option('target_col_type',c('PN','DE'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  if(!profile_col %in% colnames(DE)){
    message(sprintf('%s not in colnames of DE, please check and re-try!',profile_col))
    return(FALSE)
  }
  if(is.null(names(driver_DA_Z))) names(driver_DA_Z) <- driver_list
  if(is.null(names(driver_DE_Z))) names(driver_DE_Z) <- driver_list
  if(is.null(names(show_label))) names(show_label) <- driver_list
  if(base::length(driver_list)>top_driver_number){
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
  ## calculate graphics::layout
  if(is.null(name_col)==TRUE){
    DE <- base::cbind(DE[,base::setdiff(colnames(DE),'name')],name=rownames(DE),stringsAsFactors=FALSE)
    name_col <- 'name'
  }
  w1 <- which(is.na(DE[,profile_col])==FALSE)
  DE <- DE[w1,]
  if(profile_trend=='pos2neg') DE <- DE[order(DE[,profile_col],decreasing = TRUE),] else DE <- DE[order(DE[,profile_col],decreasing = FALSE),]
  DE_profile <- DE[,profile_col]
  DE_profile_name <- DE[,name_col]
  n_gene <- base::length(DE_profile)

  plot_part <- function(ori=FALSE,before_off=FALSE){
    if(target_nrow==2){
      n_driver <- base::length(driver_list)*2
      ratio1 <- ceiling(n_driver/15) ## profile height to rows
      ratio2 <- 4 ## width of profile to DA/DE
      rr <- 1
    } else {
      n_driver <- base::length(driver_list)
      ratio1 <- ceiling(1.2*n_driver/15) ## profile height to rows
      ratio2 <- 4 ## width of profile to DA/DE
      rr <- 1
    }
    #
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE){
      pdf(pdf_file,width=(rr*2+ratio2)*1.5,height=(ratio1+rr)*1.5)
    }
    # get graphics::layout
    graphics::layout(matrix(c(rep(0,length.out=rr),rep(1,length.out=ratio2),rep(0,length.out=rr*1),
                    rep(c(rep(4,length.out=rr),rep(2,length.out=ratio2),rep(3,length.out=rr*1)),
                        length.out=ratio1*(ratio2+rr*2))),
                  ncol=c(ratio2+rr*2),byrow=TRUE))
    ## graphics::layout
    # |0|1|0
    # |4|2|3
    ## plot 1
    par(mar=c(1.5,1.5,4,0))
    mm <- quantile(DE_profile,probs=c(0.0001,0.9999));
    mm <- base::max(abs(mm)); mm <- c(-mm,mm)
    y1 <- base::seq(mm[1],mm[2],length.out=5); y1 <- round(y1,1)
    unit <- n_gene/10; unit <- round(unit/100)*100
    x1 <- base::seq(0,n_gene,by=unit);x1 <- base::unique(x1); x1 <- c(x1[1:(base::length(x1)-1)],n_gene)
    par(usr=c(0,base::length(DE_profile),mm[1],mm[2]))
    graphics::plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
    pp <- par()$usr; rr <- (pp[2]-pp[1])/n_gene
    graphics::polygon(x=c(pp[1],c(1:n_gene)*rr+pp[1],pp[2]),y=c(0,DE_profile,0),col='grey',border='grey',xpd=TRUE,lwd=0.3)
    if(profile_trend=='pos2neg'){
      if(is.null(left_annotation)==FALSE) graphics::text(pp[1]+(pp[2]-pp[1])/100,mm[2]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
      if(is.null(right_annotation)==FALSE) graphics::text(pp[2]-(pp[2]-pp[1])/100,mm[1]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
    }else{
      if(is.null(left_annotation)==FALSE) graphics::text(pp[1]+(pp[2]-pp[1])/100,mm[1]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
      if(is.null(right_annotation)==FALSE) graphics::text(pp[2]-(pp[2]-pp[1])/100,mm[2]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    }
    graphics::axis(side=2,at=y1,labels=y1)
    graphics::mtext(side=2,line = 2.5,profile_col,cex=1)
    graphics::segments(pp[1],mm[1],pp[2],mm[1],xpd=TRUE)
    graphics::segments(x1*rr+pp[1],mm[1]-(mm[2]-mm[1])/30,x1*rr+pp[1],mm[1],xpd=TRUE)
    graphics::text(x1*rr+pp[1],mm[1]-(mm[2]-mm[1])/10,get_label_manual(x1),adj=1,xpd=TRUE)
    ## plot2
    par(mar=c(2,1.5,2,0))
    graphics::plot(1,col='white',xlab="",ylab="",xlim=c(0,n_gene),xaxt='n',yaxt='n')
    pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
    yy1 <- base::seq(from=pp[3],to=pp[4],length.out=n_driver+1)
    yy2 <- base::seq(from=pp[3],to=pp[4],length.out=base::length(driver_list)+1)
    #graphics::segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
    #graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
    #graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey')
    # add columns
    use_target_list <- target_list[driver_list]
    if(target_col_type=='DE'){
      cc <- z2col(DE_profile,sig_thre=profile_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[5:9],blue_col=brewer.pal(9,'Blues')[5:9],
                  col_max_thre=base::max(abs(DE_profile)))
      #names(cc) <- DE_profile_name
      cc[which(cc=='white')] <- 'light grey'
    }
    if(target_nrow==1){
      for(i in 1:base::length(driver_list)){
        t1 <- use_target_list[[driver_list[[i]]]]
        w0 <- which(DE_profile_name %in% t1$target)
        w1 <- w0*rr+pp[1]
        if(target_col=='black'){
          graphics::segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            graphics::segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],lwd=1.5,
                     col=cc[w0])
          }else{
            graphics::segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],lwd=1.5,
                     col=z2col(t1$spearman,sig_thre=0,col_max_thre=1,col_min_thre=0.01,
                               red_col = pos_col,blue_col=neg_col))
          }
        }
      }
    }
    if(target_nrow==2){
      for(i in 1:base::length(driver_list)){
        t1 <- use_target_list[[driver_list[[i]]]]
        t11 <- t1[which(t1$spearman>=0),]$target
        t12 <- t1[which(t1$spearman<0),]$target
        w0 <- which(DE_profile_name %in% t11)
        w1 <- w0*rr+pp[1]
        if(base::length(w1)>0){
          if(target_col=='black'){
            graphics::segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col='black',lwd=1)
          }else{
            if(target_col_type=='DE'){
              graphics::segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col=cc[w0],lwd=1.5)
            }else{
              graphics::segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col=pos_col,lwd=1.5)
            }
          }
        }
        w0 <- which(DE_profile_name %in% t12)
        w1 <- w0*rr+pp[1]
        if(base::length(w1)>0){
          if(target_col=='black'){
            graphics::segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col='black',lwd=1)
          }else{
            if(target_col_type=='DE'){
              graphics::segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col=cc[w0],lwd=1.5)
            }else{
              graphics::segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col=neg_col,lwd=1.5)
            }
          }
        }
      }
    }
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey')
    ## plot 3
    par(mar=c(2,0.5,2,2))
    graphics::plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
    pp <- par()$usr
    graphics::rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
    yy2 <- base::seq(from=pp[3],to=pp[4],length.out=base::length(driver_list)+1)
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white',xpd=TRUE)
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey',xpd=TRUE)
    graphics::abline(v=c(pp[1],(pp[1]+pp[2])/2,pp[2]))
    ## add text
    mm_min <- base::min(base::min(abs(driver_DA_Z[driver_list]),na.rm=TRUE)*0.9,base::min(abs(driver_DE_Z[driver_list]),na.rm=TRUE)*0.9)
    mm_min <- base::max(mm_min,Z_sig_thre)
    mm_max <- base::max(base::max(abs(driver_DA_Z[driver_list]),na.rm=TRUE)*1.1,base::max(abs(driver_DE_Z[driver_list]),na.rm=TRUE)*1.1)
    c1 <- z2col(driver_DA_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
                col_min_thre=mm_min,col_max_thre=mm_max)
    c2 <- z2col(driver_DE_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
                col_min_thre=mm_min,col_max_thre=mm_max)
    for(i in 1:base::length(driver_list)){
      z1 <- driver_DA_Z[driver_list[i]]
      z2 <- driver_DE_Z[driver_list[i]]
      p1 <- get_z2p(z1)
      p2 <- get_z2p(z2)
      graphics::rect(xleft=pp[1],xright=(pp[1]+pp[2])/2,ybottom=yy2[i],ytop=yy2[i+1],col=c1[i],border='dark grey',xpd=TRUE)
      graphics::rect(xright=pp[2],xleft=(pp[1]+pp[2])/2,ybottom=yy2[i],ytop=yy2[i+1],col=c2[i],border='dark grey',xpd=TRUE)
      graphics::text(x=(pp[1]+(pp[1]+pp[2])/2)/2,y=(yy2[i]+yy2[i+1])/2,p1,adj=0.5)
      graphics::text(x=(pp[2]+(pp[1]+pp[2])/2)/2,y=(yy2[i]+yy2[i+1])/2,p2,adj=0.5)
    }
    text((pp[1]+(pp[1]+pp[2])/2)/2,pp[4],'DA',xpd=TRUE,cex=1.5,pos=3)
    text((pp[2]+(pp[1]+pp[2])/2)/2,pp[4],'DE',xpd=TRUE,cex=1.5,pos=3)
    ## plot 4
    par(mar=c(2,6,2,0.2))
    graphics::plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
    pp <- par()$usr
    #yy1 <- base::seq(from=pp[3],to=pp[4],length.out=n_driver+1)
    #yy11 <- (yy1[1:(base::length(yy1)-1)]+yy1[2:base::length(yy1)])/2
    yy2 <- base::seq(from=pp[3],to=pp[4],length.out=base::length(driver_list)+1)
    yy22 <- (yy2[1:(base::length(yy2)-1)]+yy2[2:base::length(yy2)])/2
    dyy <- yy2[2]-yy2[1]
    graphics::text(show_label,x=(pp[1]+pp[2])/2,y=yy22,xpd=TRUE,adj=1)
    # add target size
    target_size <- do.call(rbind,lapply(use_target_list,function(x){
      x1 <- base::length(which(x$spearman>=0))
      x2 <- base::length(which(x$spearman<0))
      c(x1,x2)
    }))
    if(target_nrow==2){
      mm <- base::max(target_size)
      tt <- pp[2]-(pp[1]+pp[2])*0.55
      for(i in 1:base::length(driver_list)){
        graphics::rect(xleft=(pp[1]+pp[2])*0.55,xright=(pp[1]+pp[2])*0.55+target_size[i,1]/mm*tt,
             ybottom=yy22[i],ytop=yy22[i]+dyy*0.35,col=pos_col,border=NA)
        graphics::rect(xleft=(pp[1]+pp[2])*0.55,xright=(pp[1]+pp[2])*0.55+target_size[i,2]/mm*tt,
             ytop=yy22[i],ybottom=yy22[i]-dyy*0.35,col=neg_col,border=NA)
      }
      graphics::segments(x0=(pp[1]+pp[2])*0.55,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
      sst <- round(base::seq(0,mm,length.out=3))
      ss <- sst*tt/mm+(pp[1]+pp[2])*0.55
      graphics::segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
      graphics::text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
      text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
      graphics::legend(2,pp[3],c('Positively-regulated','Negatively-regulated'),fill=c(pos_col,neg_col),border = NA,bty='n',cex=0.6,xpd=TRUE,xjust=1,yjust=1)
    }
    #
    if(target_nrow==1){
      target_size <- base::rowSums(target_size)
      mm <- base::max(target_size)
      tt <- pp[2]-(pp[1]+pp[2])*0.55
      for(i in 1:base::length(driver_list)){
        graphics::rect(xleft=(pp[1]+pp[2])*0.55,xright=(pp[1]+pp[2])*0.55+target_size[i]/mm*tt,
             ybottom=yy22[i]-dyy*0.2,ytop=yy22[i]+dyy*0.2,col='dark grey',border=NA)
      }
      graphics::segments(x0=(pp[1]+pp[2])*0.55,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
      sst <- round(base::seq(0,mm,length.out=3))
      ss <- sst*tt/mm+(pp[1]+pp[2])*0.55
      graphics::segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
      graphics::text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
      text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
    }
  }
  ##
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);dev.off()} else {plot_part()}
  return(TRUE)
}
###
#' Draw GSEA (gene set enrichment analysis) Plot with NetBID Analysis of Gene Sets
#'
#' \code{draw.GSEA.NetBID.GS} creates a GSEA plot for gene sets with more NetBID analysis information.
#' Such as, number of genes in each gene set, marking the rank of annotated genes in the differential expression profile and differential activity (DA) values.
#'
#' @param DE data.frame, a data.frame created either by function \code{getDE.limma.2G} or \code{getDE.BID.2G}.
#' Row names are gene names, columns must include the calculated differencial values (e.g. "ID", "logFC", "AveExpr", "P.Value" etc.).
#' @param name_col character, the name of the column in \code{DE} contains gene names. If NULL, will use the row names of \code{DE}. Default is NULL.
#' @param profile_col character, the name of the column in \code{DE} contains calculated differencial value (e.g. "logFC" or "P.Value").
#' If \code{DE} is created by \code{getDE.limma.2G} or \code{getDE.BID.2G}, this parameter should be set to "logFC" or "t".
#' @param profile_trend character, users can choose between "pos2neg" and "neg2pos". "pos2neg" means high \code{profile_col} in target group will be shown on the left.
#' "neg2pos" means high \code{profile_col} in control group will be shown on the left. Default is "pos2neg".
#' @param use_gs2gene list, contains elements of gene sets. Element name is gene set name, each element contains a vector of genes belong to that gene set.
#' This list can be created by calling one element from \code{all_gs2gene}, or merge several gene sets into one by using \code{merge_gs}.
#' @param sig_gs_list a vector of characters, the names of top gene sets.
#' @param gs_DA_Z a vector of numerics, the Z-statistics of differentail activity (DA) values for the \code{sig_gs_list}.
#' It is highly suggested to assign name to the vector, otherwise will use name of \code{sig_gs_list}.
#' @param top_gs_number integer, the number of top significant gene sets to be displayed. Default is 30.
#' @param target_col character, name of the color palette used for display marker line in the panel.
#' Users can choose between "black" and "RdBu". If "black", the marker line in the panel is black. If "RdBu", the marker line in the panel is Red to Blue.
#' The color shade of the marker line is decided by each gene's significance of differentiation. High in red, low in blue.
#' Default is "RdBu".
#' @param left_annotation character, annotation on the left of profile curve, indicating high in control group or target group. Default is "".
#' @param right_annotation character, annotation on the right of profile curve, indicating high in the opposite group of \code{left_annotation}. Default is "".
#' @param main character, an overall title for the plot. Default is "".
#' @param profile_sig_thre numeric, threshold value for target genes. This parameter works only when \code{target_col_type} is set as "DE" and \code{target_col} is set as "RdBu".
#' Non-significant target genes will be colored grey. Default is 0.
#' @param Z_sig_thre numeric, threshold value of Z-statistics from \code{driver_DA_Z} and \code{driver_DE_Z}. Significant values will have background color. Default is 1.64.
#' @param pdf_file character, the file path to save figure as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#' @examples
#' \dontrun{
#' db.preload(use_level='transcript',use_spe='human',update=FALSE)
#'
#' ## get all_gs2gene
#'
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#'
#' ms_tab <- analysis.par$final_ms_tab
#' comp <- 'G4.Vs.others'
#' DE <- analysis.par$DE[[comp]]
#' analysis.par$out.dir.PLOT <- getwd()
#'
#' ## directory for saving the pdf files
#' exp_mat_gene <- Biobase::exprs(analysis.par$cal.eset)
#'
#' ## calculate activity for all genesets
#' use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,
#'                        use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
#' ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,cal_mat = exp_mat_gene)
#'
#' ## get DA for the gene set
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')]
#' # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')]
#' # get sample list for G1
#' DA_gs <- getDE.limma.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,
#'                         G1_name='G4',G0_name='others')
#' ## or use: DA_gs <- getDE.BID.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,
#'                         G1_name='G4',G0_name='others')
#' ## draw vocalno plot for top sig-GS
#' sig_gs <- draw.volcanoPlot(dat=DA_gs,
#'                            label_col='ID',
#'                            logFC_col='logFC',
#'                            Pv_col='P.Value',
#'                            logFC_thre=0.25,
#'                            Pv_thre=1e-4,
#'                            main='Volcano Plot for gene sets',
#'                            show_label=TRUE,
#'                            label_type = 'distribute',
#'                            label_cex = 0.5,
#'                            pdf_file=sprintf('%s/vocalno_GS_DA.pdf',
#'                            analysis.par$out.dir.PLOT))
#' ## GSEA plot for the significant gene sets
#' draw.GSEA.NetBID.GS(DE=DE,name_col='ID',
#'                     profile_col='t',profile_trend='pos2neg',
#'                     sig_gs_list = sig_gs$ID,
#'                     gs_DA_Z=DA_gs[sig_gs$ID,'Z-statistics'],
#'                     use_gs2gene = use_gs2gene,
#'                     top_gs_number=5,target_col='RdBu',
#'                     left_annotation = 'test_left',
#'                     right_annotation = 'test_right',
#'                     main='test',Z_sig_thre=1.64,profile_sig_thre = 0,
#'                     pdf_file=sprintf('%s/NetBID_GSEA_GS_demo1.pdf',
#'                     analysis.par$out.dir.PLOT))
#'}
#' @export
draw.GSEA.NetBID.GS <- function(DE=NULL,name_col=NULL,profile_col=NULL,profile_trend='pos2neg',
                                sig_gs_list=NULL,gs_DA_Z=NULL,use_gs2gene=NULL,
                                top_gs_number=30,target_col='RdBu',
                                left_annotation="",right_annotation="",main="",
                                profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  #
  all_input_para <- c('DE','name_col','profile_col','profile_trend','sig_gs_list',
                      'gs_DA_Z','top_gs_number',
                      'target_col','left_annotation','right_annotation','main',
                      'profile_sig_thre','Z_sig_thre')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('profile_trend',c('pos2neg','neg2pos'),envir=environment()),
                 check_option('target_col',c('RdBu','black'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
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
    DE <- base::cbind(DE[,base::setdiff(colnames(DE),'name')],name=rownames(DE),stringsAsFactors=FALSE)
    name_col <- 'name'
  }
  if(is.null(names(gs_DA_Z))) names(gs_DA_Z) <- sig_gs_list
  if(base::length(sig_gs_list)>top_gs_number){
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
  ## calculate graphics::layout
  w1 <- which(is.na(DE[,profile_col])==FALSE)
  DE <- DE[w1,]
  if(profile_trend=='pos2neg') DE <- DE[order(DE[,profile_col],decreasing = TRUE),] else DE <- DE[order(DE[,profile_col],decreasing = FALSE),]
  DE_profile <- DE[,profile_col]
  DE_profile_name <- DE[,name_col]
  use_gs2gene <- lapply(use_gs2gene,function(x)base::intersect(x,DE_profile_name))
  use_gs2gene <- use_gs2gene[sig_gs_list]
  #names(DE_profile) <- rownames(DE)
  n_gene <- base::length(DE_profile)
  n_driver <- base::length(sig_gs_list)
  plot_part <- function(ori=FALSE,before_off=FALSE){
    gswidth <- base::max(strwidthMod(sig_gs_list,units='inches',cex=1))
    gsheight <- base::max(strheightMod(sig_gs_list,units='inches',cex=1))*length(sig_gs_list)*2
    gswidth <- ceiling(gswidth)
    ratio1 <- ceiling(gsheight/1.5) ## profile height to rows
    ratio2 <- 6 ## width of main profile
    rr1 <- ceiling(gswidth/0.5)+1
    rr2 <- 2
    # width: 0.2|rr1|0.2|ratio2|0.1|rr2|0.2
    # 4: gswidth,0.5
    # height: 0.2|ratio1|0.3|1|0.2
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE){
      pdf(pdf_file,width=0.7+(gswidth+0.5)/rr1*(rr1+rr2+ratio2),height=0.5+(gsheight/ratio1)*(ratio1+1))
    }
    # get graphics::layout
    graphics::layout(matrix(c(rep(0,length.out=rr1),rep(1,length.out=ratio2),rep(0,length.out=rr2*1),
                    rep(c(rep(4,length.out=rr1),rep(2,length.out=ratio2),rep(3,length.out=rr2*1)),
                        length.out=ratio1*(ratio2+rr1+rr2))),
                  ncol=c(ratio2+rr1+rr2),byrow=TRUE))
    #print(matrix(c(rep(0,length.out=rr1),rep(1,length.out=ratio2),rep(0,length.out=rr2*1),
    #               rep(c(rep(4,length.out=rr1),rep(2,length.out=ratio2),rep(3,length.out=rr2*1)),
    #                   length.out=ratio1*(ratio2+rr1+rr2))),
    #             ncol=c(ratio2+rr1+rr2),byrow=TRUE))
    # 0|1|0
    # 4|2|3
    ## plot 1
    par(mai=c(0.1,0.1,0.2,0.1))
    mm <- quantile(DE_profile,probs=c(0.0001,0.9999));
    mm <- base::max(abs(mm)); mm <- c(-mm,mm)
    y1 <- base::seq(mm[1],mm[2],length.out=5); y1 <- round(y1,1)
    unit <- n_gene/10; unit <- round(unit/100)*100
    x1 <- base::seq(0,n_gene,by=unit);x1 <- base::unique(x1); x1 <- c(x1[1:(base::length(x1)-1)],n_gene)
    par(usr=c(0,base::length(DE_profile),mm[1],mm[2]))
    graphics::plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
    pp <- par()$usr; rr <- (pp[2]-pp[1])/n_gene
    graphics::polygon(x=c(pp[1],c(1:n_gene)*rr+pp[1],pp[2]),y=c(0,DE_profile,0),col='grey',border='grey',xpd=TRUE,lwd=0.3)
    if(profile_trend=='pos2neg'){
      if(is.null(left_annotation)==FALSE) graphics::text(pp[1]+(pp[2]-pp[1])/100,mm[2]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
      if(is.null(right_annotation)==FALSE) graphics::text(pp[2]-(pp[2]-pp[1])/100,mm[1]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
    }else{
      if(is.null(left_annotation)==FALSE) graphics::text(pp[1]+(pp[2]-pp[1])/100,mm[1]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
      if(is.null(right_annotation)==FALSE) graphics::text(pp[2]-(pp[2]-pp[1])/100,mm[2]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    }
    graphics::axis(side=2,at=y1,labels=y1)
    graphics::mtext(side=2,line = 2.5,profile_col,cex=1)
    graphics::segments(pp[1],mm[1],pp[2],mm[1],xpd=TRUE)
    graphics::segments(x1*rr+pp[1],mm[1]-(mm[2]-mm[1])/30,x1*rr+pp[1],mm[1],xpd=TRUE)
    graphics::text(x1*rr+pp[1],mm[1]-(mm[2]-mm[1])/10,get_label_manual(x1),adj=1,xpd=TRUE,cex=0.6)
    ## plot2
    par(mai=c(0.2,0.1,0.2,0.1))
    graphics::plot(1,col='white',xlab="",ylab="",xlim=c(0,n_gene),xaxt='n',yaxt='n')
    pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
    yy1 <- base::seq(from=pp[3],to=pp[4],length.out=n_driver+1)
    yy2 <- base::seq(from=pp[3],to=pp[4],length.out=base::length(sig_gs_list)+1)
    # add columns
    use_target_list <- use_gs2gene[sig_gs_list]
    cc <- z2col(DE_profile,sig_thre=profile_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[5:9],blue_col=brewer.pal(9,'Blues')[5:9],
                col_max_thre=base::max(abs(DE_profile)))
    cc[which(cc=='white')] <- 'light grey'
    for(i in 1:base::length(sig_gs_list)){
      t1 <- use_target_list[[sig_gs_list[i]]]
      w0 <- which(DE_profile_name %in% t1)
      w1 <- w0*rr+pp[1]
      if(target_col=='black'){
        graphics::segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],col='black',lwd=1)
      }else{
        graphics::segments(x0=w1,x1=w1,y0=yy1[i],y1=yy1[i+1],lwd=1.5,col=cc[w0])
      }
    }
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey')
    ## plot 3
    par(mai=c(0.2,0,0.2,0.2))
    graphics::plot(1,col='white',xlab="",ylab="",xlim=c(0,1),xaxt='n',yaxt='n',bty='n')
    pp <- par()$usr
    graphics::rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
    yy2 <- base::seq(from=pp[3],to=pp[4],length.out=base::length(sig_gs_list)+1)
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white',xpd=TRUE)
    graphics::segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey',xpd=TRUE)
    ## add text
    mm_min <- base::min(abs(gs_DA_Z[sig_gs_list]),na.rm=TRUE)*0.9
    mm_min <- base::max(mm_min,Z_sig_thre)
    mm_max <- base::max(abs(gs_DA_Z[sig_gs_list]),na.rm=TRUE)*1.1
    c1 <- z2col(gs_DA_Z[sig_gs_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
                col_min_thre=mm_min,col_max_thre=mm_max)
    for(i in 1:base::length(sig_gs_list)){
      z1 <- gs_DA_Z[sig_gs_list[i]]
      p1 <- get_z2p(z1)
      graphics::rect(xleft=pp[1],xright=pp[2],ybottom=yy2[i],ytop=yy2[i+1],col=c1[i],border='dark grey',xpd=TRUE)
      graphics::text(x=(pp[1]+pp[2])/2,y=(yy2[i]+yy2[i+1])/2,p1,adj=0.5)
    }
    text((pp[1]+pp[2])/2,pp[4],'DA',xpd=TRUE,cex=1.5,pos=3)
    ## plot 4
    par(mai=c(0.2,0.2,0.2,0.1))
    graphics::plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
    pp <- par()$usr
    yy2 <- base::seq(from=pp[3],to=pp[4],length.out=base::length(sig_gs_list)+1)
    yy22 <- (yy2[1:(base::length(yy2)-1)]+yy2[2:base::length(yy2)])/2
    dyy <- yy22[2]-yy22[1]
    # add target size
    target_size <- unlist(lapply(use_gs2gene,length))
    mm <- base::max(target_size)
    rr <- ceiling(gswidth)*1.5
    tt_left <- pp[2]-(pp[2]-pp[1])/(1+rr)
    tt <- (pp[2]-pp[1])/(1+rr)
    graphics::text(show_label,x=tt_left-tt/25,y=yy22,xpd=TRUE,adj=1)
    for(i in 1:base::length(sig_gs_list)){
      graphics::rect(xleft=tt_left,xright=tt_left+target_size[i]/mm*tt,
           ybottom=yy22[i]-dyy*0.2,ytop=yy22[i]+dyy*0.2,col='dark grey',border=NA)
    }
    graphics::segments(x0=tt_left,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(base::seq(0,mm,length.out=3))
    ss <- sst*tt/mm+tt_left
    graphics::segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+0.3*(pp[4]-pp[3])/length(sig_gs_list),xpd=TRUE)
    graphics::text(x=ss,y=pp[4]+0.4*(pp[4]-pp[3])/length(sig_gs_list),srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Size',x=tt_left-tt/10,y=pp[4]+(pp[4]-pp[3])/length(sig_gs_list),adj=1,xpd=TRUE,cex=0.8)
    ##
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off();} else {plot_part()}
  graphics::layout(1);
  return(TRUE)
}



#' Merge Target Gene List for Two Drivers
#'
#' \code{merge_target_list} merges target gene list for two drivers together. Shared target genes with high "MI (mutual information)" statistics will be kept in the final target list.
#'
#' @param driver1 character, the name of the first driver.
#' @param driver2 character, the name of the second driver.
#' @param target_list list, the driver-to-target list object. The names of the list elements are drivers (e.g. driver1 and driver2).
#' Each element is a data frame, usually contains at least three columns. "target", target gene names; "MI", mutual information; "spearman", spearman correlation coefficient.
#' Users can call \code{get_net2target_list} to create this list.
#' @return Return a data.frame with rows of target genes, column of "target", "MI", "spearman".
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' driver1 <- ms_tab[1,'originalID_label']
#' driver2 <- ms_tab[2,'originalID_label']
#' m1 <- merge_target_list(driver1=driver1,driver2=driver2,
#'                           target_list=analysis.par$merge.network$target_list)
#' @export
merge_target_list <- function(driver1=NULL,driver2=NULL,target_list=NULL){
  #
  all_input_para <- c('driver1','driver2','target_list')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  t1 <- target_list[driver1][[1]]
  t2 <- target_list[driver2][[1]]
  ov <- base::intersect(t1$target,t2$target)
  rownames(t1) <- t1$target
  rownames(t2) <- t2$target
  if(base::length(ov)==0){
    t_out <- base::rbind(t1,t2)
  }else{
    t_out1 <- base::rbind(t1[base::setdiff(rownames(t1),ov),],t2[base::setdiff(rownames(t2),ov),])
    t_out2 <- list()
    for(each_ov in ov){
      x1 <- t1[each_ov,2:3]
      x2 <- t2[each_ov,2:3]
      if(x1[1]>x2[1]) t_out2[[each_ov]] <- t1[each_ov,] else t_out2[[each_ov]] <- t2[each_ov,]
    }
    t_out2 <- do.call(rbind,t_out2)
    t_out <- base::rbind(t_out1,t_out2)
  }
  return(t_out)
}
#
#' Scatter Box Plot of Driver's Expression Values and Activity Values across Subgroup of Samples
#'
#' \code{draw.categoryValue} draws a scatter box plot to visualize one selected driver's expression value and activity value across different phenotype subgroups of samples.
#' Two side-by-side scatter box plots will be created. The left plot shows driver's activity values in different phenotype subgroups, each point is a sample.
#' The right plot shows driver's expression value in different phenotype subgroups, each point is a sample.
#'
#' @param ac_val a vector of numerics, the activity values of the selected driver across all samples.
#' @param exp_val a vector of numerics, the expression values of the selected driver across all samples.
#' @param use_obs_class a vector of characters, the category of sample. The order of samples here must match the order in \code{ac_val} and \code{exp_val}.
#' Users can call \code{get_obs_label} to create this vector.
#' @param class_order a vector of characters, the order of catefory (subgroup).
#' If NULL, will use the alphabetical order of the category (subgroup). Default is NULL.
#' @param category_color a vector of characters, a vector of color codes for each category in \code{class_order}.
#' If NULL, will call \code{get.class.color} to create the vector. Default is NULL.
#' @param stripchart_color character, the color of the scatter of points. Default is "black" with transparency alpha equals 0.7.
#' @param strip_cex numeric, giving the amount by which the size of scattered points should be magnified relative to the default. Default is 1.
#' @param class_srt numeric, rotation angle of the column labels (subgroup labels) at the bottom of the box plot. Default is 90.
#' @param class_cex numeric, giving the amount by which the text of category (subgroup) labels should be magnified relative to the default. Default is 1.
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be save. Default is NULL.
#' @param main_ac character, the main title of the plot to show activity values. Default is "".
#' @param main_exp character, the main title of the plot to show expression values. Default is "".
#' @param main_cex numeric, giving the amount by which the text of the main title should be magnified relative to the default. Default is 1.
#' @param use_color a vector of color codes, colors to be assigned to each member of display label. Default is brewer.pal(9, 'Set1').
#' @param pre_define a vector of characters, pre-defined color codes for a certain input (e.g. c("blue", "red") with names c("A", "B")). Default is NULL.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' use_driver <- driver_list[3]
#' exp_mat <- Biobase::exprs(analysis.par$cal.eset)
#' ## expression,the rownames could match originalID
#' ac_mat  <- Biobase::exprs(analysis.par$merge.ac.eset)
#' ## activity,the rownames could match originalID_label
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
#' draw.categoryValue(ac_val=ac_mat[use_driver,],
#'                    exp_val=exp_mat[ms_tab[use_driver,'originalID'],],
#'                    use_obs_class=use_obs_class,
#'                    class_order=c('WNT','SHH','G4'),
#'                    class_srt=30,
#'                    main_ac = ms_tab[use_driver,'gene_label'],
#'                    main_exp=ms_tab[use_driver,'geneSymbol'],
#'                    pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' use_driver <- driver_list[3]
#' exp_mat <- Biobase::exprs(analysis.par$cal.eset)
#' ## rownames could match originalID
#' ac_mat  <- Biobase::exprs(analysis.par$merge.ac.eset)
#' ## rownames could match originalID_label
#' phe_info <- Biobase::pData(analysis.par$cal.eset)
#' use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' draw.categoryValue(ac_val=ac_mat[use_driver,],
#'                    exp_val=exp_mat[ms_tab[use_driver,'originalID'],],
#'                    use_obs_class=use_obs_class,
#'                    class_order=c('WNT','SHH','G4'),
#'                    class_srt=30,
#'                    main_ac = ms_tab[use_driver,'gene_label'],
#'                    main_exp=ms_tab[use_driver,'geneSymbol'],
#'                    pdf_file=sprintf('%s/categoryValue_demo1.pdf',
#'                    analysis.par$out.dir.PLOT))
#'}
#' @export
draw.categoryValue <- function(ac_val=NULL,exp_val=NULL,use_obs_class=NULL,category_color=NULL,
                               stripchart_color=get_transparent('black',0.7),strip_cex=1,class_order=NULL,class_srt=90,class_cex=1,pdf_file=NULL,
                               main_ac="",main_exp="",main_cex=1,
                               use_color=NULL,pre_define=NULL){
  #
  all_input_para <- c('ac_val','use_obs_class','strip_cex',
                      'class_srt','class_cex','main_ac','main_exp','main_cex')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  use_obs_class <- clean_charVector(use_obs_class)
  #
  if(is.null(names(use_obs_class))==FALSE){
    if(is.null(ac_val)==FALSE) use_obs_class <- use_obs_class[names(ac_val)]
    if(is.null(exp_val)==FALSE) use_obs_class <- use_obs_class[names(exp_val)]
  }
  if(is.null(class_order)){
    class_order <- sort(base::unique(use_obs_class))
  }
  if(is.null(category_color)==TRUE){
    class_col  <- get.class.color(class_order,use_color=use_color,pre_define=pre_define)
    class_col1 <- get_transparent(class_col,0.5)
  }else{
    class_col1 <- category_color
  }
  c1 <- 0
  if(is.null(ac_val)==FALSE){c1 <- c1+1}
  if(is.null(exp_val)==FALSE){c1 <- c1+1}
  labelWidth <- base::max(strwidthMod(class_order,'inches',cex=class_cex)*sin(class_srt*pi/180))
  if(is.null(pdf_file)==FALSE){
    hh <- 5+1.5+labelWidth
    if(c1==1) pdf(pdf_file,width=1.5+3,height=hh)
    if(c1==2) pdf(pdf_file,width=1.5+3*2,height=hh)
  }
  if(c1>1) graphics::layout(t(matrix(1:c1)))
  par(mai=c(labelWidth+0.5,1,1,0.5))
  if(is.null(ac_val)==FALSE){
    ddf <- data.frame(data=ac_val,class=factor(use_obs_class,levels=class_order))
    a <- boxplot(data~class,data=ddf,ylab='Activity Value',col=class_col1,outline=FALSE,border='dark grey',cex.lab=1.2,names=NA,bty='n',
                 ylim=c(base::min(ddf$data),base::max(ddf$data)),main=main_ac,cex.main=main_cex,xlab='')
    graphics::text(1:base::length(class_order),par()$usr[3]-(par()$usr[4]-par()$usr[3])/20,adj=0.5+class_srt/180,class_order,srt=class_srt,xpd=TRUE,cex=class_cex)
    graphics::stripchart(data~class,data=ddf,add=TRUE,pch=16,method='jitter',vertical=TRUE,col=stripchart_color,cex=strip_cex)
  }
  if(is.null(exp_val)==FALSE){
    ddf <- data.frame(data=exp_val,class=factor(use_obs_class,levels=class_order))
    a <- boxplot(data~class,data=ddf,col=class_col1,ylab='Expression Value',outline=FALSE,border='dark grey',cex.lab=1.2,names=NA,bty='n',
                 ylim=c(base::min(ddf$data),base::max(ddf$data)),main=main_exp,cex.main=main_cex,xlab='')
    graphics::text(1:base::length(class_order),par()$usr[3]-(par()$usr[4]-par()$usr[3])/20,adj=0.5+class_srt/180,class_order,srt=class_srt,xpd=TRUE,cex=class_cex)
    graphics::stripchart(data~class,data=ddf,add=TRUE,pch=16,method='jitter',vertical=TRUE,col=stripchart_color,cex=strip_cex)
  }
  graphics::layout(1)
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}

### simple functions
get_transparent <- function(x,alpha=0.1){
  grDevices::rgb(t(grDevices::col2rgb(x)/255),alpha=alpha)
}

get_label_manual <- function(x){
  x1 <- sapply(x,function(x2){
    x3 <- unlist(strsplit(as.character(x2),""))
    x4 <- base::length(x3)%/%3 ## add number
    if(x4>0){
      pp <- base::length(x3)-base::seq(1,x4)*3; x3[pp] <- paste0(x3[pp],','); base::paste(x3,collapse="")
    }else{
      x2
    }
  })
  unlist(x1)
}

#' Target Network Structure Plot for One Driver
#'
#' \code{draw.targetNet} draws a network structure to display the target genes of one selected driver. Edges of positively-regulated target genes are orange,
#' edges of negatively-regulated target genes are green. The width of the edges shows the strength of regulation.
#'
#' @param source_label character, the label of selected one driver.
#' @param source_z numeric, the Z-statistic of the selected driver. The color shade of driver's node in the network is decided by this Z-statistic.
#' If NULL, the driver node will be colored grey. Default is NULL.
#' @param edge_score a named vector of numerics, indicating the correlation between the driver and its target genes. The range of the numeric value is from -1 to 1.
#' Positive value means it is positively-regulated by driver and vice versa. The names of the vector are gene names.
#' @param label_cex numeric, giving the amount by which the text of target gene names should be magnified relative to the default. Default is 0.7.
#' @param source_cex numeric, giving the amount by which the text of driver name should be magnified relative to the default. Default is 1.
#' @param arrow_direction character, users can choose between "in" and "out". Default is "out".
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#' @param n_layer integer, number of circle layers to display. Default is 1.
#' @param alphabetical_order logical, if TRUE, the targe gene names will be sorted alphabetically. If FALSE, will be sorted by statistics. Default is FALSE.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#'
#' @examples
#' source_label <- 'test1'
#' source_z <- 1.96
#' edge_score <- (sample(1:200,size=100,replace=TRUE)-100)/100
#' names(edge_score) <- paste0('G',1:100)
#' draw.targetNet(source_label=source_label,source_z=source_z,
#'                edge_score=edge_score)
#' draw.targetNet(source_label=source_label,source_z=source_z,
#'                edge_score=edge_score,n_layer=2)
#' draw.targetNet(source_label=source_label,source_z=source_z,
#'                edge_score=edge_score,
#'                arrow_direction='in',
#'                source_cex=2)
#' \dontrun{
#' source_label <- 'test1'
#' source_z <- 1.96
#' edge_score <- (sample(1:200,size=100,replace=TRUE)-100)/100
#' names(edge_score) <- paste0('G',1:100)
#' analysis.par <- list()
#' analysis.par$out.dir.PLOT <- getwd()
#' draw.targetNet(source_label=source_label,source_z=source_z,
#'                edge_score=edge_score,
#'                pdf_file=sprintf('%s/targetNet.pdf',
#'                analysis.par$out.dir.PLOT))
#'}
#' @export
draw.targetNet <- function(source_label="",source_z=NULL,edge_score=NULL,
                           label_cex=0.7,source_cex=1,
                           pdf_file=NULL,arrow_direction='out',n_layer=1,alphabetical_order=FALSE){
  #
  all_input_para <- c('source_label','edge_score','label_cex','source_cex',
                      'arrow_direction','n_layer','alphabetical_order')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('alphabetical_order',c(TRUE,FALSE),envir=environment()),
                 check_option('arrow_direction',c('in','out'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  edge_score<- sort(edge_score)
  tmp1 <- sapply(base::unique(names(edge_score)),function(x){
    x1 <- edge_score[which(names(edge_score)==x)]
    x1[base::which.max(abs(x1))]
  })
  names(tmp1) <- base::unique(names(edge_score))
  edge_score <- tmp1
  edge_score<- edge_score[order(edge_score,decreasing = TRUE)]
  g1 <- names(edge_score)
  ec <- z2col(edge_score*100,sig_thre=0,n_len=base::length(edge_score),red_col=pos_col,blue_col=neg_col);names(ec) <- names(edge_score)
  ec <- get_transparent(ec,alpha=0.8); names(ec) <- names(edge_score)
  ew <- 2*label_cex*(abs(edge_score)-base::min(abs(edge_score)))/(base::max(abs(edge_score))-base::min(abs(edge_score)))+label_cex/2; names(ew) <- names(edge_score)
  if(base::max(abs(edge_score))-base::min(abs(edge_score))==0){ew <- rep(label_cex/2,length.out=length(edge_score)); names(ew) <- names(edge_score)}
  t2xy <- function(tt,radius=1) {
    t2p <- pi*2 * tt
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  if(alphabetical_order==TRUE)  g1 <- sort(g1)
  plot_part <- function(ori=FALSE,before_off=FALSE){
    geneWidth <- base::max(strwidthMod(g1,'inches',cex=label_cex))
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=6+2*geneWidth*n_layer,height=6+2*geneWidth*n_layer)
    par(mai=c(1,1,1,1))
    graphics::plot(1,xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white',xlab="",ylab="",xaxt='n',yaxt='n')
    pp <- par()$usr
    ## add target label
    #tt <- seq(-0.5,0.5,length.out=base::length(g1)+1)[-1]; ## g1
    rad_v <- base::seq(from=0.5,to=1,length.out=n_layer)
    if(n_layer==1) rad_v=0.8
    #tmp1 <- ceiling(base::length(g1)/sum(rad_v)*rad_v)
    uu <- ceiling(base::length(g1)/base::length(rad_v))
    tmp1 <- rep(uu,length.out=base::length(rad_v))
    if(n_layer>1) tmp1[base::length(tmp1)] <- base::length(g1)-sum(tmp1[1:(base::length(tmp1)-1)]) else tmp1 <- base::length(g1)
    tmp1<-cumsum(tmp1)
    all_g1 <- list()
    for(i in 1:n_layer){
      if(i==1) all_g1[[i]] <- g1[1:tmp1[i]] else all_g1[[i]]  <- g1[(tmp1[i-1]+1):tmp1[i]]
    }
    if(alphabetical_order==FALSE) all_g1 <- lapply(all_g1,function(x)x[order(edge_score[x])])
    #all_tt <- lapply(1:n_layer,function(i)seq(-0.5,0.5,length.out=base::length(all_g1[[i]])+1)[-1])
    all_tt <- lapply(1:n_layer,function(i){
      if(i==1) return(seq(-0.5,0.5,length.out=uu+1)[-1][1:base::length(all_g1[[i]])])
      if(i>1) return(seq(-0.5-(i-1)/(n_layer*uu),0.5-(i-1)/(n_layer*uu),length.out=uu+1)[-1][1:base::length(all_g1[[i]])])
    })
    all_p <- lapply(1:n_layer,function(i)t2xy(all_tt[[i]],radius=rad_v[i]))
    # add line
    for(i in 1:n_layer){
      each_v <- rad_v[i]
      p1 <- all_p[[i]]
      g1_use <- all_g1[[i]]
      tt <- all_tt[[i]]
      geneWidth <- strwidthMod(source_label,'user',cex=source_cex)
      if(arrow_direction=='out'){
        p2<-t2xy(tt,radius=each_v-label_cex/36);
        p3<-t2xy(tt,radius=each_v-label_cex/48);
        graphics::arrows(x0=0,y0=0,x1=p2$x,y1=p2$y,col=ec[g1_use],lwd=ew[g1_use],angle=10,length=0.1*label_cex,xpd=TRUE);
      }else{
        p2<-t2xy(tt,radius=each_v-label_cex/36);
        p3<-t2xy(tt,radius=each_v-label_cex/36);
        p4<-t2xy(tt,radius=geneWidth/2);
        graphics::arrows(x0=p2$x,y0=p2$y,x1=p4$x,y1=p4$y,col=ec[g1_use],lwd=ew[g1_use],angle=5,length=0.1*label_cex,xpd=TRUE);
      }
      graphics::points(p3$x,p3$y,pch=16,col='dark grey',cex=label_cex)
    }
    # add label
    for(j in 1:n_layer){
      p1 <- all_p[[j]]
      g1_use <- all_g1[[j]]
      for(i in 1:base::length(p1$x)) graphics::text(p1$x[i],p1$y[i],g1_use[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    }
    ## add source label
    if(is.null(source_z)==TRUE){
      draw.ellipse(0,0,a=1.05*geneWidth/2,b=1.05*geneWidth/2,col='light grey',border=NA)
    }else{
      draw.ellipse(0,0,a=1.05*geneWidth/2,b=1.05*geneWidth/2,col=z2col(source_z),border=NA)
    }
    graphics::text(0,0,source_label,adj=0.5,xpd=TRUE,cex=source_cex)
    graphics::legend(x=pp[1],y=pp[3],fill=c(pos_col,neg_col),c('Positively-regulated','Negatively-regulated'),bty='n',xpd=T,border=NA,cex=label_cex,horiz = FALSE)
  }
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}
  return(TRUE)
}


#' Target Ntwork Structure Plot for Two Drivers
#'
#' \code{draw.targetNet.TWO} draws a network structure to display the target genes of two selected drivers. Edges of positively-regulated target genes are orange,
#' edges of negatively-regulated target genes are green. The width of the edges shows the strength of regulation.
#' It will also print out the number of shared and unique targe genes for each driver, with P-value and odds ratio.
#'
#' @param source1_label character, the label of the first selected driver (to be displayed on the left).
#' @param source2_label character, the label of the second selected driver (to be displayed on the right).
#' @param source1_z numeric, the Z-statistic of the first driver. The color shade of driver’s node in the network is decided by this Z-statistic.
#' If NULL, the driver will be colored grey. Default is NULL.
#' @param source2_z numeric, the Z-statistic of the second driver.The color shade of driver’s node in the network is decided by this Z-statistic.
#' If NULL, the driver will be colored grey. Default is NULL.
#' @param edge_score1 a named vector of numerics, indicating the correlation between the first driver and its target genes.
#' The range of the numeric value is from -1 to 1. Positive value means it is positively-regulated by driver and vice versa. The names of the vector are gene names.
#' @param edge_score2 a named vector of numerics, indicating the correlation between the seconde driver and its target genes.
#' The range of the numeric value is from -1 to 1. Positive value means it is positively-regulated by driver and vice versa. The names of the vector are gene names.
#' @param arrow_direction1 character, the arrow direction for first driver. Users can choose between "in" and "out". Default is "out".
#' @param arrow_direction2 character, the arrow direction for second driver. Users can choose between "in" and "out". Default is "out".
#' @param label_cex numeric, giving the amount by which the text of target gene names should be magnified relative to the default. Default is 0.7.
#' @param source_cex numeric, giving the amount by which the text of driver name should be magnified relative to the default. Default is 1.
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#' @param total_possible_target numeric or a vector of characters. If input is numeric, it is the total number of possible target genes.
#' If input is a vector of characters, it is the background list of all possible target genes.
#' This parameter will be passed to function \code{test.targetNet.overlap} to test whether the target genes of the two drivers are significantly intersected.
#' If NULL, will do not perform this test. Default is NULL.
#' @param show_test logical, if TRUE, the test result will be printed and returned. Default is FALSE.
#' @param n_layer integer, number of circle layers to display. Default is 1.
#' @param alphabetical_order logical, if TRUE, the targe gene names will be sorted alphabetically. If FALSE, will be sorted by statistics. Default is FALSE.
#' @return If \code{show_test}==FALSE, will return a logical value indicating whether the plot has been successfully generated,
#' otherwise will return the statistics of testing when total_possible_target is not NULL.
#'
#' @examples
#' source1_label <- 'test1'
#' source1_z <- 1.96
#' edge_score1 <- (sample(1:160,size=80,replace=TRUE)-80)/80
#' names(edge_score1) <- sample(paste0('G',1:1000),size=80)
#' source2_label <- 'test2'
#' source2_z <- -2.36
#' edge_score2 <- (sample(1:240,size=120,replace=TRUE)-120)/120
#' names(edge_score2) <- sample(paste0('G',1:1000),size=120)
#' draw.targetNet.TWO(source1_label=source1_label,
#'                source2_label=source2_label,
#'                source1_z=source1_z,source2_z=source2_z,
#'                edge_score1=edge_score1,edge_score2=edge_score2,
#'                total_possible_target=paste0('G',1:1000),
#'                show_test=TRUE,label_cex=0.6)
#' draw.targetNet.TWO(source1_label=source1_label,
#'                source2_label=source2_label,
#'                source1_z=source1_z,source2_z=source2_z,
#'                edge_score1=edge_score1,edge_score2=edge_score2,
#'                total_possible_target=paste0('G',1:1000),
#'                show_test=TRUE,label_cex=0.6,n_layer=2)
#'
#' \dontrun{
#' source1_label <- 'test1'
#' source1_z <- 1.96
#' edge_score1 <- (sample(1:160,size=100,replace=TRUE)-80)/80
#' names(edge_score1) <- sample(paste0('G',1:1000),size=100)
#' source2_label <- 'test2'
#' source2_z <- -2.36
#' edge_score2 <- (sample(1:240,size=100,replace=TRUE)-120)/120
#' names(edge_score2) <- sample(paste0('G',1:1000),size=100)
#' analysis.par <- list()
#' analysis.par$out.dir.PLOT <- getwd()
#' draw.targetNet.TWO(source1_label=source1_label,
#'                source2_label=source2_label,
#'                source1_z=source1_z,source2_z=source2_z,
#'                edge_score1=edge_score1,edge_score2=edge_score2,
#'                total_possible_target=paste0('G',1:1000),show_test=TRUE,
#'                pdf_file=sprintf('%s/targetNetTWO.pdf',
#'                analysis.par$out.dir.PLOT))
#' }
#' @export
draw.targetNet.TWO <- function(source1_label="",source2_label="",
                               source1_z=NULL,source2_z=NULL,
                               edge_score1=NULL,edge_score2=NULL,
                               arrow_direction1='out',arrow_direction2='out',
                               label_cex=0.7,source_cex=1,pdf_file=NULL,
                               total_possible_target=NULL,show_test=FALSE,n_layer=1,alphabetical_order=FALSE){
  #
  all_input_para <- c('source1_label','source2_label','source1_z','source2_z',
                      'edge_score1','edge_score2','label_cex','source_cex',
                      'arrow_direction1','arrow_direction2',
                      'n_layer','alphabetical_order','show_test')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('alphabetical_order',c(TRUE,FALSE),envir=environment()),
                 check_option('show_test',c(TRUE,FALSE),envir=environment()),
                 check_option('arrow_direction1',c('in','out'),envir=environment()),
                 check_option('arrow_direction2',c('in','out'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  tmp1 <- sapply(base::unique(names(edge_score1)),function(x){
    x1 <- edge_score1[which(names(edge_score1)==x)];x1[base::which.max(abs(x1))]
  })
  names(tmp1) <- base::unique(names(edge_score1));edge_score1 <- tmp1
  tmp1 <- sapply(base::unique(names(edge_score2)),function(x){
    x1 <- edge_score2[which(names(edge_score2)==x)];x1[base::which.max(abs(x1))]
  })
  names(tmp1) <- base::unique(names(edge_score2));edge_score2 <- tmp1
  edge_score1<- sort(edge_score1,decreasing = FALSE)
  edge_score2<- sort(edge_score2,decreasing = TRUE)
  g12 <- base::intersect(names(edge_score1),names(edge_score2))
  g1  <- base::setdiff(names(edge_score1),names(edge_score2))
  g2  <- base::setdiff(names(edge_score2),names(edge_score1))
  ec1 <- z2col(edge_score1*100,sig_thre=0,n_len=base::length(edge_score1),red_col=pos_col,blue_col=neg_col);names(ec1) <- names(edge_score1)
  ec2 <- z2col(edge_score2*100,sig_thre=0,n_len=base::length(edge_score2),red_col=pos_col,blue_col=neg_col);names(ec2) <- names(edge_score2)
  if(base::max(abs(edge_score1)) - base::min(abs(edge_score1))==0){
	  ew1 <- rep(2 * label_cex+ label_cex/2,length.out=length(edge_score1))
  }else{
	  ew1 <- 2 * label_cex * (abs(edge_score1) - base::min(abs(edge_score1)))/(base::max(abs(edge_score1)) -
										   base::min(abs(edge_score1))) + label_cex/2
  }
  names(ew1) <- names(edge_score1)
  if(base::max(abs(edge_score2)) - base::min(abs(edge_score2))==0){
	  ew2 <- rep(2 * label_cex+ label_cex/2,length.out=length(edge_score2))
  }else{
          ew2 <- 2 * label_cex * (abs(edge_score2) - base::min(abs(edge_score2)))/(base::max(abs(edge_score2)) -                                                              base::min(abs(edge_score2))) + label_cex/2
  }
  names(ew2) <- names(edge_score2)
  #ew1 <- 2*label_cex*(abs(edge_score1)-base::min(abs(edge_score1)))/(base::max(abs(edge_score1))-base::min(abs(edge_score1)))+label_cex/2; names(ew1) <- names(edge_score1)
  #ew2 <- 2*label_cex*(abs(edge_score2)-base::min(abs(edge_score2)))/(base::max(abs(edge_score2))-base::min(abs(edge_score2)))+label_cex/2; names(ew2) <- names(edge_score2)
  t2xy <- function(tt,radius=1,init.angle=0) {
    t2p <- pi*2 * tt + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  #g1|g12|g2
  if(alphabetical_order==TRUE){
    g1 <- sort(g1);
    g2 <- sort(g2);
    g12 <- sort(g12);
  }
  plot_part <- function(ori=FALSE,before_off=FALSE){
    geneWidth <- base::max(strwidthMod(g1,'inches',cex=label_cex))
    if(before_off==TRUE) dev.off()
    if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=10+4*geneWidth,height=8+2*geneWidth)
    graphics::plot(1,xlim=c(-1.4,1.4),ylim=c(-1,1),bty='n',col='white',xlab="",ylab="",xaxt='n',yaxt='n')
    par(mai=c(1,1,1,1))

    geneWidth1 <- strwidthMod(source1_label,'user',cex=source_cex)
    geneWidth2 <- strwidthMod(source2_label,'user',cex=source_cex)
    geneWidth <- base::max(geneWidth1,geneWidth2)

    ag <- 0.245-0.01*n_layer
    lp <- 0.4
    rad_v <- base::seq(from=0.6,to=1,length.out=n_layer)
    if(n_layer==1) rad_v=0.8

    if(base::length(g1)>0){
      # get info
      uu <- ceiling(base::length(g1)/base::length(rad_v))
      tmp1 <- rep(uu,length.out=base::length(rad_v))
      if(n_layer>1) tmp1[base::length(tmp1)] <- base::length(g1)-sum(tmp1[1:(base::length(tmp1)-1)]) else tmp1 <- base::length(g1)
      tmp1<-cumsum(tmp1)
      all_g1 <- list()
      for(i in 1:n_layer){
        if(i==1) all_g1[[i]] <- g1[1:tmp1[i]] else all_g1[[i]]  <- g1[(tmp1[i-1]+1):tmp1[i]]
      }
      #all_g1 <- lapply(all_g1,function(x)x[order(edge_score[x])])
      all_tt <- lapply(1:n_layer,function(i){
        if(i==1) return(seq(-ag,ag,length.out=uu)[1:base::length(all_g1[[i]])])
        if(i>1) return(seq(-ag-2*ag*(i-1)/(n_layer*uu),ag-2*ag*(i-1)/(n_layer*uu),length.out=uu)[1:base::length(all_g1[[i]])])
      })
      all_p <- lapply(1:n_layer,function(i)t2xy(all_tt[[i]],radius=rad_v[i],init.angle= -180))

      # add line
      for(i in 1:n_layer){
        each_v <- rad_v[i]
        p1 <- all_p[[i]]
        g1_use <- all_g1[[i]]
        tt <- all_tt[[i]]
        if(arrow_direction1=='out'){
          p2<-t2xy(tt,radius=each_v-label_cex/36,init.angle= -180);
          p3<-t2xy(tt,radius=each_v-label_cex/48,init.angle= -180);
          graphics::arrows(x0=-lp,y0=0,x1=p2$x-lp,y1=p2$y,col=ec1[g1_use],lwd=ew1[g1_use],angle=10,length=0.1*label_cex,xpd=TRUE);
        }else{
          p2<-t2xy(tt,radius=each_v-label_cex/36,init.angle= -180);
          p3<-t2xy(tt,radius=each_v-label_cex/36,init.angle= -180);
          p4<-t2xy(tt,radius=geneWidth/2,init.angle= -180);
          graphics::arrows(x0=p2$x-lp,y0=p2$y,x1=p4$x-lp,y1=p4$y,col=ec1[g1_use],lwd=ew1[g1_use],angle=5,length=0.1*label_cex,xpd=TRUE);
        }
        graphics::points(p3$x-lp,p3$y,pch=16,col='dark grey',cex=label_cex)
      }
      # add label
      for(j in 1:n_layer){
        p1 <- all_p[[j]]
        g1_use <- all_g1[[j]]
        for(i in 1:base::length(p1$x)) graphics::text(p1$x[i]-lp,p1$y[i],g1_use[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
      }
    }
    ##
    if(base::length(g2)>0){
      # get info
      uu <- ceiling(base::length(g2)/base::length(rad_v))
      tmp1 <- rep(uu,length.out=base::length(rad_v))
      if(n_layer>1) tmp1[base::length(tmp1)] <- base::length(g2)-sum(tmp1[1:(base::length(tmp1)-1)]) else tmp1 <- base::length(g2)
      tmp1<-cumsum(tmp1)
      all_g1 <- list()
      for(i in 1:n_layer){
        if(i==1) all_g1[[i]] <- g2[1:tmp1[i]] else all_g1[[i]]  <- g2[(tmp1[i-1]+1):tmp1[i]]
      }
      #all_g1 <- lapply(all_g1,function(x)x[order(edge_score[x])])
      all_tt <- lapply(1:n_layer,function(i){
        if(i==1) return(seq(-ag,ag,length.out=uu)[1:base::length(all_g1[[i]])])
        if(i>1) return(seq(-ag-2*ag*(i-1)/(n_layer*uu),ag-2*ag*(i-1)/(n_layer*uu),length.out=uu)[1:base::length(all_g1[[i]])])
      })
      all_p <- lapply(1:n_layer,function(i)t2xy(all_tt[[i]],radius=rad_v[i],init.angle=0))

      # add line
      for(i in 1:n_layer){
        each_v <- rad_v[i]
        p1 <- all_p[[i]]
        g1_use <- all_g1[[i]]
        tt <- all_tt[[i]]
        if(arrow_direction2=='out'){
          p2<-t2xy(tt,radius=each_v-label_cex/36);
          p3<-t2xy(tt,radius=each_v-label_cex/48);
          graphics::arrows(x0=lp,y0=0,x1=p2$x+lp,y1=p2$y,col=ec2[g1_use],lwd=ew2[g1_use],angle=10,length=0.1*label_cex,xpd=TRUE);
        }else{
          p2<-t2xy(tt,radius=each_v-label_cex/36);
          p3<-t2xy(tt,radius=each_v-label_cex/36);
          p4<-t2xy(tt,radius=geneWidth/2);
          graphics::arrows(x0=p2$x+lp,y0=p2$y,x1=p4$x+lp,y1=p4$y,col=ec2[g1_use],lwd=ew2[g1_use],angle=5,length=0.1*label_cex,xpd=TRUE);
        }
        graphics::points(p3$x+lp,p3$y,pch=16,col='dark grey',cex=label_cex)
      }
      # add label
      for(j in 1:n_layer){
        p1 <- all_p[[j]]
        g1_use <- all_g1[[j]]
        for(i in 1:base::length(p1$x)) graphics::text(p1$x[i]+0.4,p1$y[i],g1_use[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
      }
    }
    ##
    if(base::length(g12)>0){
      rm <- base::min(0.1*base::length(g12),1)
      dd <- par.char2pos()[2]/2; nr <-ceiling(base::length(g12)/(rm*2/dd));each_col_n<-ceiling(base::length(g12)/nr)
      tt <- base::seq(rm,-rm,length.out=each_col_n);
      xx <- seq(-lp+geneWidth,lp-geneWidth,length.out=nr)
      if(nr==1) xx<-0
      tt <- unlist(lapply(tt,function(x)rep(x,length.out=nr)))[1:base::length(g12)]
      xx <- rep(xx,length.out=base::length(xx)*each_col_n)[1:base::length(g12)]

      if(arrow_direction1=='out'){
        graphics::arrows(x0=-lp,y0=0,x1=xx,y1=tt,col=ec1[g12],lwd=ew1[g12],angle=10,length=0.1*label_cex,xpd=TRUE)
      }else{
        p4<-t2xy(tt,radius=geneWidth/2);
        graphics::arrows(x0=xx,y0=tt,x1=-lp,y1=0,col=ec1[g12],lwd=ew1[g12],angle=5,length=0.1*label_cex,xpd=TRUE);
      }
      if(arrow_direction2=='out'){
        graphics::arrows(x0=lp,y0=0,x1=xx,y1=tt,col=ec2[g12],lwd=ew2[g12],angle=10,length=0.1*label_cex,xpd=TRUE)
      }else{
        p4<-t2xy(tt,radius=geneWidth/2);
        graphics::arrows(x0=xx,y0=tt,x1=lp,y1=0,col=ec1[g12],lwd=ew1[g12],angle=5,length=0.1*label_cex,xpd=TRUE);
      }
      boxtext(xx,tt,labels=g12,col.bg=get_transparent('light grey',0.3),cex=label_cex)
    }
    if(is.null(source2_z)==TRUE)
      draw.ellipse(lp,0,a=geneWidth/2,b=geneWidth/2,col='light grey',border=NA)
    else
      draw.ellipse(lp,0,a=geneWidth/2,b=geneWidth/2,col=z2col(source2_z),border=NA)

    graphics::text(lp,0,source2_label,adj=0.5,cex=source_cex)

    if(is.null(source1_z)==TRUE)
      draw.ellipse(-lp,0,a=geneWidth/2,b=geneWidth/2,col='light grey',border=NA)
    else
      draw.ellipse(-lp,0,a=geneWidth/2,b=geneWidth/2,col=z2col(source1_z),border=NA)
    text(-lp,0,source1_label,adj=0.5,cex=source_cex)
    pp <- par()$usr
    graphics::legend(x=pp[1],y=pp[3],fill=c(pos_col,neg_col),c('Positively-regulated','Negatively-regulated'),bty='n',xpd=T,border=NA,cex=label_cex,horiz = TRUE)
  }
  # fisher test for target
  if(is.null(pdf_file)==FALSE){plot_part(ori=TRUE);plot_part(ori=TRUE,before_off=TRUE);while (!is.null(dev.list()))  dev.off()} else {plot_part()}

  if(is.null(total_possible_target)==FALSE & show_test==TRUE){
    res <- test.targetNet.overlap(source1_label,source2_label,names(edge_score1),names(edge_score2),total_possible_target)
    return(res)
  }
  return(TRUE)
}

#' Test for Intersection of Target Genes between Two Drivers
#'
#' \code{test.targetNet.overlap} performs Fisher's exact test to see whether the target genes from two drivers are significantly intersected.
#'
#'
#' @param source1_label character, the label of the first selected driver.
#' @param source2_label character, the label of the second selected driver.
#' @param target1 a vector of characters, the list of target genes from the first driver.
#' @param target2 a vector of characters, the list of target genes from the second driver.
#' @param total_possible_target numeric or a vector of characters. If input is numeric, it is the total number of possible target genes.
#' If input is a vector of characters, it is the background list of all possible target genes.
#'
#' @return Return statistics of the testing, including the \code{P.Value}, \code{Odds_Ratio} and \code{Intersected_Number}.
#'
#' @examples
#' source1_label <- 'test1'
#' target1 <- sample(paste0('G',1:1000),size=80)
#' source2_label <- 'test2'
#' target2 <- sample(paste0('G',1:1000),size=120)
#' test.targetNet.overlap(source1_label=source1_label,source2_label=source2_label,
#'                target1=target1,target2=target2,
#'                total_possible_target=paste0('G',1:1000))
#' \dontrun{
#' }
#' @export
test.targetNet.overlap <- function(source1_label=NULL,source2_label=NULL,
                                   target1=NULL,target2=NULL,
                                   total_possible_target=NULL){
  #
  all_input_para <- c('source1_label','source2_label','target1','target2')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  t1  <- base::unique(target1)
  t2  <- base::unique(target2)
  print(sprintf('%s has %d unique targets !',source1_label,base::length(t1)))
  print(sprintf('%s has %d unique targets !',source2_label,base::length(t2)))
  n11 <- base::length(base::intersect(t1,t2))
  n12 <- base::length(base::setdiff(t1,t2))
  n21 <- base::length(base::setdiff(t2,t1))
  if(class(total_possible_target) %in% c('integer','numeric')){
    n22 <- total_possible_target-base::length(union(t1,t2))
  }else{
    n22 <- base::length(base::setdiff(total_possible_target,c(t1,t2)))
  }
  mm  <- base::cbind(c(n11,n21),c(n12,n22))
  ft  <- fisher.test(mm)$p.value
  or  <- n11/n12/(n21/n22)
  rownames(mm) <- c(sprintf('In %s target',source1_label),sprintf('Not in %s target',source1_label))
  colnames(mm) <- c(sprintf('In %s target',source2_label),sprintf('Not in %s target',source2_label))
  print(mm)
  res <- c('P.Value'=ft,'Odds_Ratio'=or,'Intersected_Number'=mm[1,1])
  return(res)
}

#' Convert Pairwise Network Data Frame to Driver-to-Target List
#'
#' \code{get_net2target_list} is a helper function in the \code{get.SJAracne.network}.
#' But if users have their own pairwise gene network files, they can convert it to driver-to-target list object.
#'
#' @param net_dat data.frame, must contain two columns with column names "source" (driver) and "target" (target genes).
#' "MI" (mutual information) and "spearman" (spearman correlation coefficient) columns are optional, but strongly suggested to use.
#' If "MI" and "spearman" columns are missing, errors may occur in some following steps (e.g. es.method='weightedmean' in \code{cal.Activity}).
#'
#' @return Return a list. The names of the list elements are drivers.
#' Each element is a data frame, contains three columns. "target", target gene names;
#' "MI", mutual information; "spearman", spearman correlation coefficient.
#'
#' @examples
#' tf.network.file <- sprintf('%s/demo1/network/SJAR/project_2019-02-14/%s/%s',
#'                    system.file(package = "NetBID2"),
#'                    'output_tf_sjaracne_project_2019-02-14_out_.final',
#'                    'consensus_network_ncol_.txt')
#' net_dat      <- read.delim(file=tf.network.file,stringsAsFactors = FALSE)
#' target_list  <- get_net2target_list(net_dat)
#' @export
get_net2target_list <- function(net_dat=NULL) {
  all_input_para <- c('net_dat')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  all_source <- base::unique(net_dat$source)
  all_target <- lapply(all_source, function(x) {
    n1 <- net_dat[which(net_dat$source == x), base::intersect(c('target', 'MI', 'spearman'),colnames(net_dat))]
    if(class(n1)=='character') n1 <- data.frame('target'=n1,'MI'=1,'spearman'=1,stringsAsFactors=F)
    n1 <- unique(n1)
    if(length(unique(n1$target))!=length(n1$target)){ ## multiple
      t1 <- table(n1$target)
      w1 <- names(which(t1==1)); w2 <- names(which(t1>1))
      w21 <- n1[which(n1$target %in% w1),]
      w22 <- do.call(rbind,lapply(w2,function(x){
        x1 <- n1[which(n1$target==x),]
        x1 <- x1[which.max(x1$MI),]
      }))
      n1 <- rbind(w21,w22)
    }
    rownames(n1) <- n1$target
    return(n1)
  })
  names(all_target) <- all_source
  return(all_target)
}

#' Read SJARACNe Network Result and Return it as List Object
#'
#' \code{get.SJAracne.network} reads SJARACNe network construction result and returns a list object
#' with network data frame, driver-to-target list and igraph object wrapped inside.
#'
#' In the demo, "consensus_network_ncol_.txt" file will be read and convert into a list object.
#' This list contains three elements, \code{network_data}, \code{target_list} and \code{igraph_obj}.
#' \code{network_dat} is a data.frame, contains all the information of the network SJARACNe constructed.
#' \code{target_list} is a driver-to-target list object. Please check details in \code{get_net2target_list}.
#' \code{igraph_obj} is an igraph object used to save this directed and weighted network.
#' Each edge of the network has two attributes, \code{weight} and \code{sign}.
#' \code{weight} is the "MI (mutual information)" value and \code{sign} is the sign of the spearman
#' correlation coefficient (1, positive regulation; -1, negative regulation).
#'
#' @param network_file character, the path for storing network file. For the output of SJAracne, the name of the network file will be "consensus_network_ncol_.txt" under the output directory.
#'
#' @return Return a list containing three elements, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#'
#' @examples
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             project_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
#'
#' \dontrun{
#' }
#' @export
get.SJAracne.network <- function(network_file=NULL){
  all_input_para <- c('network_file')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  net_dat      <- read.delim(file=network_file,stringsAsFactors = FALSE)
  target_list  <- get_net2target_list(net_dat)
  igraph_obj   <- graph_from_data_frame(net_dat[,c('source','target')],directed=TRUE) ## add edge weight ???
  if('MI' %in% colnames(net_dat)) igraph_obj   <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  if('spearman' %in% colnames(net_dat)) igraph_obj   <- set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}

#' Update Network List Object Using Constraints
#'
#' \code{update_SJAracne.network} updates the network object created by \code{get.SJAracne.network}, using constraints like statistical thresholds and interested gene list.
#'
#' @param network_list list, the network list object created by \code{get.SJAracne.network}.
#' The list contains three elements, \code{network_dat}, \code{target_list} and \code{igraph_obj}. For details, please check \code{get.SJAracne.network}.
#' @param all_possible_drivers a vector of characters, all possible drivers used to filter the network.
#' If NULL, will use drivers from \code{network_list}. Default is NULL.
#' @param all_possible_targets a vector of characters, all possible target genes used to filter the network.
#' If NULL, will use targets from \code{network_list}. Default is NULL.
#' @param force_all_drivers logical, if TRUE, will include all drivers from \code{all_possible_drivers} into the final network.
#' For \code{network_dat} and \code{target_list} in the network list object, all genes in \code{all_possible_drivers} will not be filtered using the following statistical thresholds.
#' For \code{igraph_obj} in the network list object, all genes in \code{all_possible_drivers} that don't exist in the original network, will be kept as vertices.
#' Default is TRUE.
#' @param force_all_targets logical, if TRUE, will include all genes from \code{all_possible_targets} into the final network.
#' For \code{network_dat} and \code{target_list} in the network list object, all genes in \code{all_possible_drivers} will not be filtered using the following statistical thresholds.
#' For \code{igraph_obj} in the network list object, all genes in \code{all_possible_drivers} that don't exist in the original network, will be kept as vertices.
#' Default is TRUE.
#' @param min_MI numeric, minimum threshold for "MI (mutual information)". Default is 0.
#' @param max_p.value numeric, maximum threshold for P-value. Default is 1.
#' @param min_spearman_value numeric, minimum threshold for spearman absolute value. Default is 0.
#' @param min_pearson_value numeric, minimum threshold for pearson absolute value. Default is 0.
#' @param spearman_sign_use a vector of numerics, users can choose from 1, -1 and c(1, -1). 1 means only positve spearman values will be used.
#' -1 means only negative spearman values will be used. Default is c(1,-1).
#' @param pearson_sign_use a vector of numerics, users can choose from 1, -1 and c(1, -1). 1 means only positve pearson values will be used.
#' -1 means only negative pearson values will be used. Default is c(1,-1).
#' @param directed logical, if TRUE, the network is a directed graph. Default is TRUE.
#' @param weighted logical, if TRUE, the network is weighted. Default is TRUE.
#'
#' @return Return a list containing three elements, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#'
#' @examples
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             project_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
#' all_possible_drivers <- c(names(analysis.par$tf.network$target_list)[1:1000],
#'                            'addition_driver_1','addition_driver_2')
#' tf.network.update <- update_SJAracne.network(network_list=analysis.par$tf.network,
#'                                              all_possible_drivers=all_possible_drivers,
#'                                              force_all_drivers=TRUE,
#'                                              force_all_targets=FALSE,
#'                                              pearson_sign_use=1)
#' print(base::intersect(c('addition_driver_1','addition_driver_2'),
#'       names(V(tf.network.update$igraph_obj)))) ## check
#' \dontrun{
#' }
#' @export update_SJAracne.network
update_SJAracne.network <- function(network_list=NULL,
                                    all_possible_drivers=NULL,
                                    all_possible_targets=NULL,
                                    force_all_drivers=TRUE,
                                    force_all_targets=TRUE,
                                    min_MI=0,max_p.value=1,
                                    min_spearman_value=0,
                                    min_pearson_value=0,
                                    spearman_sign_use=c(1,-1),
                                    pearson_sign_use=c(1,-1),
                                    directed=TRUE,weighted=TRUE
){
  #
  all_input_para <- c('network_list','force_all_drivers','force_all_targets',
                      'min_MI','max_p.value','min_spearman_value','min_pearson_value',
                      'spearman_sign_use','pearson_sign_use','directed','weighted')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('force_all_drivers',c(TRUE,FALSE),envir=environment()),
                 check_option('force_all_targets',c(TRUE,FALSE),envir=environment()),
                 check_option('directed',c(TRUE,FALSE),envir=environment()),
                 check_option('weighted',c(TRUE,FALSE),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  n1 <- names(network_list)
  n2 <- c('network_dat','target_list','igraph_obj')
  if(base::length(base::setdiff(n2,n1))>0){
    message(sprintf('%s not included in the network_list, please chech and re-try !',base::paste(base::setdiff(n2,n1),collapse=';')));
    return(FALSE)
  }
  ori_net_dat <- network_list$network_dat
  net_dat <- ori_net_dat
  if(is.null(all_possible_drivers)==TRUE) all_possible_drivers <- base::unique(net_dat$source)
  if(is.null(all_possible_targets)==TRUE) all_possible_targets <- base::unique(net_dat$target)
  # basic statistics filter
  w1 <- which(net_dat$MI>=min_MI & net_dat$p.value<=max_p.value
              & abs(net_dat$pearson)>=min_pearson_value &
                abs(net_dat$spearman)>=min_spearman_value)
  all_w1 <- w1 ## use rows
  # sign choose
  if(!1 %in% spearman_sign_use){ ## do not use the positive ones
    w1 <- which(net_dat$spearman<=0)
    all_w1 <- base::setdiff(all_w1,w1)
  }
  if(!-1 %in% spearman_sign_use){ ## do not use the negative ones
    w1 <- which(net_dat$spearman>=0)
    all_w1 <- base::setdiff(all_w1,w1)
  }
  if(!1 %in% pearson_sign_use){ ## do not use the positive ones
    w1 <- which(net_dat$pearson<=0)
    all_w1 <- base::setdiff(all_w1,w1)
  }
  if(!-1 %in% pearson_sign_use){ ## do not use the negative ones
    w1 <- which(net_dat$pearson>=0)
    all_w1 <- base::setdiff(all_w1,w1)
  }
  # nodes filter
  all_possible_nodes <- base::unique(c(base::unique(net_dat$source),base::unique(net_dat$target))) ## original all possible
  if(force_all_drivers==TRUE){
    w1 <- which(net_dat$source %in% all_possible_drivers)
    all_w1 <- base::unique(c(all_w1,w1))
    all_possible_nodes <- base::unique(c(all_possible_nodes,all_possible_drivers))
  }
  if(force_all_targets==TRUE){
    w1 <- which(net_dat$target %in% all_possible_targets)
    all_w1 <- base::unique(c(all_w1,w1))
    all_possible_nodes <- base::unique(c(all_possible_nodes,all_possible_targets))
  }
  message(sprintf('%d from %d edges are kept in the network !',base::length(all_w1),nrow(ori_net_dat)))
  message(sprintf('%d nodes will be used to generate the igraph!',base::length(all_possible_nodes)))
  # keep all genes in all* to be in the igraph
  net_dat <- net_dat[which(net_dat$source %in% all_possible_drivers & net_dat$target %in% all_possible_targets),] ## filter by all
  igraph_obj   <- graph_from_data_frame(net_dat[,c('source','target')],directed=directed,vertices = all_possible_nodes) ##
  if(weighted==TRUE & 'MI' %in% colnames(net_dat)) igraph_obj <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  if(directed==TRUE & 'spearman' %in% colnames(net_dat)) igraph_obj   <- set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
  target_list  <- get_net2target_list(net_dat)
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}


#' QC Tables and Plots for Network Object
#'
#' \code{draw.network.QC} creates tables and plots for showing some basic statistics,
#' driver statistics and scale-free checking of the target network.
#'
#' @param igraph_obj igraph object, created by \code{get.SJAracne.network}.
#' @param outdir character, the output directory for saving QC tables and plots.
#' @param prefix character, the prefix of output QC figures names. Default is "".
#' @param directed logical, if TRUE, this network will be treated as a directed network. Default is TRUE.
#' @param weighted logical, if TRUE, this network will be treated as a weighted network. Default is FALSE.
#' @param generate_html logical, if TRUE, a html file will be created by R Markdown.
#' If FALSE, plots will be save as separated PDFs.
#' Default is TRUE.
#' @param html_info_limit logical, if TRUE, the statistics for network QC html will be limited. Default is TRUE.
#' @return Return a logical value. If TRUE, success in creating QC tables and plots.
#'
#' @examples
#' \dontrun{
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             project_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' draw.network.QC(analysis.par$tf.network$igraph_obj,
#'                 outdir=analysis.par$out.dir.QC,prefix='TF_net_')
#' }
#' @export
draw.network.QC <- function(igraph_obj,outdir=NULL,prefix="",directed=TRUE,weighted=FALSE,generate_html=TRUE,html_info_limit=TRUE){
  #
  all_input_para <- c('igraph_obj','prefix','directed','weighted','generate_html','html_info_limit')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('directed',c(TRUE,FALSE),envir=environment()),
                 check_option('weighted',c(TRUE,FALSE),envir=environment()),
                 check_option('generate_html',c(TRUE,FALSE),envir=environment()),
                 check_option('html_info_limit',c(TRUE,FALSE),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message(paste0("The output directory: \"", outdir, "\" is created!"))
  }else
    message(paste0("The output will overwrite the files in directory: \"",outdir,"\""))
  if(class(igraph_obj)!='igraph'){
    message('Should input igraph object ! ');return(FALSE);
  }
  net <- igraph_obj
  if(generate_html==TRUE){
    if(rmarkdown::pandoc_available()==FALSE){
      stop('pandoc not available, please set Sys.setenv(RSTUDIO_PANDOC=$pandoc_installed_path), or set generate_html=FALSE')
    }
    directed <- directed
    weighted <- weighted
    output_rmd_file <- sprintf('%s/%snetQC.Rmd',outdir,prefix)
    file.copy(from=system.file('Rmd/net_QC.Rmd',package = "NetBID2"),to=output_rmd_file)
    rmarkdown::render(output_rmd_file, rmarkdown::html_document(toc = TRUE))
    return(TRUE)
  }
  ###
  deg <- igraph::degree(net,mode='out')
  source_list <- names(deg)[which(deg>0)]
  #
  res_file <- sprintf('%s/%snetwork_info.pdf',outdir,prefix)
  pdf(res_file,width=8,height=8);
  #
  par(mar=c(6,6,6,8))
  d_out <- igraph::degree(net,mode = 'all')
  a <- graphics::hist(d_out,breaks = 20,xlab='Degree',cex.lab=1.2,cex.axis=1.2,cex.main=1.2,
            main=sprintf('Density plot of degree distribution for all nodes \n (network node:%d, network edge:%d)',base::length(V(net)),base::length(E(net))));
  d1 <- stats::density(d_out)
  mm <- base::max(a$counts)/base::max(d1$y);mm1 <- base::seq(0,base::max(d1$y),length.out = 5);mm2 <- format(mm1,scientific=TRUE,digits=3);mm2[1]<-'0';
  lines(x=d1$x,y=d1$y*mm,col=get_transparent('red',0.5),lwd=1.5)
  graphics::axis(side=4,at=mm*mm1,labels=mm2,las=2);
  graphics::mtext(side=4,line = 6,'Density',cex=1.2)
  d_out <- igraph::degree(net,mode = 'out')[source_list]
  a <- graphics::hist(d_out,breaks = 20,xlab='Degree',cex.lab=1.2,cex.axis=1.2,cex.main=1.2,
            main=sprintf('Density plot of target size for all %s drivers \n (Size from %s to %s; mean Size: %s, median Size: %s )',base::length(source_list),base::min(d_out),base::max(d_out),format(base::mean(d_out),digits=5),stats::median(d_out)));
  d1 <- stats::density(d_out)
  mm <- base::max(a$counts)/base::max(d1$y);mm1 <- base::seq(0,base::max(d1$y),length.out = 5);mm2 <- format(mm1,scientific=TRUE,digits=3);mm2[1]<-'0';
  lines(x=d1$x,y=d1$y*mm,col=get_transparent('red',0.5),lwd=1.5)
  graphics::axis(side=4,at=mm*mm1,labels=mm2,las=2);
  graphics::mtext(side=4,line = 6,'Density',cex=1.2)
  #
  res1 <- check_scalefree(net)
  dev.off()
  return(TRUE)
}

## functions to check the scale free feature of the network
check_scalefree <- function(igraph_obj) {
  gr1 <- igraph_obj
  fp1 <- igraph::degree_distribution(gr1)
  dd <- as.data.frame(base::cbind(k = 1:base::max(igraph::degree(gr1)), pk = fp1[-1]))
  dd$pk <- dd$pk + 1 / base::length(V(gr1))
  r2 <-
    stats::lm(log10(dd$pk) ~ log10(dd$k))
  r3 <- summary(r2)$adj.r.squared
  if(base::length(dd$k)>100) graphics::plot(pk ~ k,data = dd,log = 'xy',main = sprintf('R2:%s', format(r3,digits=4)),pch=16,col=get_transparent('dark grey',0.8),cex.lab=1.4,cex.axis=1.2)
  if(base::length(dd$k)<=100) graphics::plot(pk ~ k,data = dd,log = 'xy',main = sprintf('R2:%s', format(r3,digits=4)),pch=16,col=get_transparent('black',0.8),cex.lab=1.4,cex.axis=1.2)
  graphics::abline(a=r2$coefficients[1],b=r2$coefficients[2],col=get_transparent('red',0.5),lwd=2)
  return(r3)
}


#' Prepare Data Files for Running SJARACNe
#'
#' \code{SJAracne.prepare} prepares data files for running SJAracne. SJARACNe is a scalable software tool for gene network reverse engineering from big data.
#' Detailed description and how to run SJARACNe can be found in its GitHub repository.
#' The usage of SJARACNe may be updated and the bash file generated by this function may not fit the version in use (not suit for SJAracne 0.2.0).
#' Please check \url{https://github.com/jyyulab/SJARACNe/} for details.
#'
#' @param eset an ExpressionSet class object, which contains the expression matrix.
#' @param use.samples a vector of characters, the list of sample used to run SJARACNe.
#' @param TF_list a vector of characters, the TF list.
#' @param SIG_list a vector of characters, the SIG list.
#' @param SJAR.main_dir character, the path to save the results generated by SJARACNe.
#' @param SJAR.project_name character, the project name used to label the output directory.
#' @param IQR.thre numeric, the IQR filter threshold to filter all non-driver genes.
#' @param IQR.loose_thre numeric, the IQR filter threshold to filter for all driver(TF/SIG) genes.
#' @param geneSymbol_column character, the column name in fdata(eset) which contains gene symbol.
#' If NULL, will use the main ID in exprs(eset) to fulfill the "geneSymbol" column. Default is NULL.
#' @param add_options additional option for running SJARACNe.
#' @examples
#' \dontrun{
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' db.preload(use_level='gene',use_spe='human',update=FALSE)
#' use_gene_type <- 'external_gene_name' ## this should user-defined !!!
#' use_genes <- rownames(Biobase::fData(network.par$net.eset))
#' use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
#' #select sample for analysis
#' phe <- Biobase::pData(network.par$net.eset)
#' use.samples <- rownames(phe) ## use all samples, or choose to use some samples
#' prj.name <- network.par$project.name # can use other names, if need to run different use samples
#' network.par$out.dir.SJAR <- 'test' ## set the directory
#' SJAracne.prepare(eset=network.par$net.eset,
#'                  use.samples=use.samples,
#'                  TF_list=use_list$tf,
#'                  SIG_list=use_list$sig,
#'                  IQR.thre = 0.5,IQR.loose_thre = 0.1,
#'                  SJAR.project_name=prj.name,
#'                  SJAR.main_dir=network.par$out.dir.SJAR)
#' }
#' @export
SJAracne.prepare <-
  function(eset,use.samples = rownames(Biobase::pData(eset)),
           TF_list=NULL,SIG_list=NULL,
           SJAR.main_dir='',
           SJAR.project_name = "",
           IQR.thre=0.5,IQR.loose_thre=0.1,add_options='',geneSymbol_column=NULL) {
    #
    all_input_para <- c('eset','use.samples','TF_list','SIG_list','SJAR.main_dir','SJAR.project_name','IQR.thre','IQR.loose_thre')
    check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
    if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
    #
    SJAR.outdir <- file.path(SJAR.main_dir, SJAR.project_name)
    if (!file.exists(SJAR.outdir)) {
      dir.create(SJAR.outdir, recursive = TRUE)
    }
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
    d <- Biobase::exprs(eset)[, use.samples]
    # filter genes with count=0
    d <- d[!apply(d == 0, 1, all), ]
    # filter genes with IQR
    choose1 <- IQR.filter(d, rownames(d),thre = IQR.thre,loose_thre=IQR.loose_thre,loose_gene=base::unique(c(TF_list,SIG_list)))
    d <- d[choose1, ]
    use.genes <- rownames(d)
    use.genes <- use.genes[which(is.na(use.genes)==FALSE)]
    w1 <- which(rownames(d) %in% use.genes)
    d <- d[w1,]
    use.genes <- rownames(d)
    # write exp data to exp format
    use.genes.symbol <- use.genes
    if(is.null(geneSymbol_column)==FALSE){
      if(geneSymbol_column %in% colnames(Biobase::fData(eset))){use.genes.symbol <- Biobase::fData(eset)[use.genes,geneSymbol_column]}
    }
    use.genes <- clean_charVector(use.genes)
    use.genes.symbol <- clean_charVector(use.genes.symbol)
    expdata <- data.frame(isoformId = use.genes, geneSymbol = use.genes.symbol, d, stringsAsFactors=FALSE)
    #
    write.table(
      expdata,
      file = SJAR.expression_matrix,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    ##
    cat(base::intersect(use.genes, TF_list),file = SJAR.hub_genes.tf,sep = '\n')
    cat(base::intersect(use.genes, SIG_list),file = SJAR.hub_genes.sig,sep = '\n')
    # write scripts to bash file for tf
    network_project_name <- SJAR.project_name
    tf_out_path <- SJAR.main_dir
    tf_pj <- sprintf('%s_TF',network_project_name)
    sig_out_path <- SJAR.main_dir
    sig_pj <- sprintf('%s_SIG',network_project_name)
    cmd_tf <- sprintf('sjaracne %s %s %s %s %s',tf_pj,SJAR.expression_matrix,SJAR.hub_genes.tf,tf_out_path,add_options)
    cmd_sig <- sprintf('sjaracne %s %s %s %s %s',sig_pj,SJAR.expression_matrix,SJAR.hub_genes.sig,sig_out_path,add_options)
    cat(cmd_tf,file = SJAR.bash_file.tf,sep = '\n')
    cat(cmd_sig,file = SJAR.bash_file.sig,sep = '\n')
    message(sprintf('The running command is: %s for TF, and the bash file is generated in %s.',cmd_tf,SJAR.bash_file.tf))
    message(sprintf('The running command is: %s for SIG, and the bash file is generated in %s.',cmd_sig,SJAR.bash_file.sig))
    message('The bash files only suit for SJAracne 0.1.0, if you are using other versions (e.g 0.2.0), please check the usage online and prepare the bash scripts by yourself !')
    return(TRUE)
  }

####
#' Calculate Differential Expression (DE) or Differential Activity (DA) by Using Bayesian Inference
#'
#' \code{bid} calculates the differential expression (DE) / differential activity (DA) by using Bayesian Inference method.
#' Users can choose different regression models and pooling strategies.
#'
#' It is a core function inside \code{getDE.BID.2G}.
#' This function allows users to have access to more options when calculating the statistics using Bayesian Inference method.
#' In some cases, the input expression matrix could be at probe/transcript level, but DE/DA calculated at gene level is expected.
#' By setting pooling strategy, users can successfully solve the special cases.
#' The P-value is estimated by the posterior distribution of the coefficient.
#'
#' @param mat matrix, the expression/activity matrix of IDs (gene/transcript/probe) from one gene. Rows are IDs, columns are samples.
#' It is strongly suggested to contain rownames of IDs and column names of samples. Example, geneA has two probes A1 and A2 across all 6 samples (Case-rep1, Case-rep2, Case-rep3, Control-rep1, Control-rep2 and Control-rep3).
#' The \code{mat} of geneA is a 2*6 numeric matrix. Likewise, if geneA has only one probe, the \code{mat} is a one-row matrix.
#' @param use_obs_class a vector of characters, the category of sample.
#' If the vector names are not available, the order of samples in \code{use_obs_class} must be the same as in \code{mat}.
#' Users can call \code{get_obs_label} to create this vector.
#' @param class_order a vector of characters, the order of the sample's category.
#' The first class in this vector will be considered as the control group by default.
#' If NULL, the order will be assigned using alphabetical order. Default is NULL.
#' @param class_ordered logical, if TRUE, the \code{class_order} will be ordered. And the order must be consistent with the phenotypic trend,
#' such as "low", "medium", "high". Default is TRUE.
#' @param method character, users can choose between "MLE" and "Bayesian".
#' "MLE", the maximum likelihood estimation, will call generalized linear model(glm/glmer) to perform data regression.
#' "Bayesian", will call Bayesian generalized linear model (bayesglm) or multivariate generalized linear mixed model (MCMCglmm) to perform data regression.
#' Default is "Bayesian".
#' @param family character or family function or the result of a call to a family function.
#' This parameter is used to define the model's error distribution. See \code{?family} for details.
#' Currently, options are gaussian, poisson, binomial(for two-group sample classes)/category(for multi-group sample classes)/ordinal(for multi-group sample classes with class_ordered=TRUE).
#' If set with gaussian or poission, the response variable in the regression model will be the expression level, and the independent variable will be the sample's phenotype.
#' If set with binomial, the response variable in the regression model will be the sample phenotype, and the independent variable will be the expression level.
#' For binomial, category and ordinal input, the family will be automatically reset, based on the sample's class level and the setting of \code{class_ordered}.
#' Default is gaussian.
#' @param pooling character, users can choose from "full","no" and "partial".
#' "full", use probes as independent observations.
#' "no", use probes as independent variables in the regression model.
#' "partial", use probes as random effect in the regression model.
#' Default is "full".
#' @param prior.V.scale numeric, the V in the parameter "prior" used in \code{MCMCglmm}.
#' It is meaningful to set when one choose "Bayesian" as method and "partial" as pooling.
#' Default is 0.02.
#' @param prior.R.nu numeric, the R-structure in the parameter "prior" used in \code{MCMCglmm}.
#' It is meaningful to set when one choose "Bayesian" as method and "partial" as pooling.
#' Default is 1.
#' @param prior.G.nu numeric, the G-structure in the parameter "prior" used in \code{MCMCglmm}.
#' It is meaningful to set when one choose "Bayesian" as method and "partial" as pooling.
#' Default is 2.
#' @param nitt numeric, the parameter "nitt" used in \code{MCMCglmm}.
#' It is meaningful to set when one choose "Bayesian" as method and "partial" as pooling.
#' Default is 13000.
#' @param burnin numeric, the parameter "burnin" used in \code{MCMCglmm}.
#' It is meaningful to set when one choose "Bayesian" as method and "partial" as pooling.
#' Default is 3000.
#' @param thin numeric, the parameter "thin" used in \code{MCMCglmm}.
#' It is meaningful to set when one choose "Bayesian" as method and "partial" as pooling.
#' Default is 10.
#' @param std logical, if TRUE, the expression matrix will be normalized by column. Default is TRUE.
#' @param logTransformed logical, if TRUE, log transformation has been performed. Default is TRUE.
#' @param log.base numeric, the base of log transformation when \code{do.logtransform} is set to TRUE. Default is 2.
#' @param average.method character, the method applied to calculate FC (fold change). Users can choose between "geometric" and "arithmetic".
#' Default is "geometric".
#' @param pseudoCount integer, the integer added to avoid "-Inf" showing up during log transformation in the FC (fold change) calculation.
#' @param return_model logical, if TRUE, the regression model will be returned; Otherwise, just return basic statistics from the model. Default is FALSE.
#' @param use_seed integer, the random seed. Default is 999.
#' @param verbose logical, if TRUE, print out additional information during calculation. Default is FALSE.
#'
#' @return Return a one-row data frame with calculated statistics for one gene/gene set if \code{return_model} is FALSE.
#' Otherwise, the regression model will be returned.
#'
#' @examples
#' mat <- matrix(c(0.50099,1.2108,1.0524,-0.34881,-0.13441,-0.87112,
#'                 1.84579,2.0356,2.6025,1.62954,1.88281,1.29604),
#'                 nrow=2,byrow=TRUE)
#' rownames(mat) <- c('A1','A2')
#' colnames(mat) <- c('Case-rep1','Case-rep2','Case-rep3',
#'                   'Control-rep1','Control-rep2','Control-rep3')
#' res1 <- bid(mat=mat,
#'            use_obs_class = c(rep('Case',3),rep('Control',3)),
#'            class_order = c('Control','Case'))
#' \dontrun{
#' }
#' @export
bid <- function(mat=NULL,use_obs_class=NULL,class_order=NULL,class_ordered=TRUE,
                method='Bayesian',family=gaussian,pooling='full',
                prior.V.scale=0.02,prior.R.nu=1,prior.G.nu=2,nitt = 13000,burnin =3000,thin=10,
                std=TRUE,logTransformed=TRUE,log.base=2,
                average.method='geometric',pseudoCount=0,return_model=FALSE,use_seed=999,verbose=FALSE){
  #check input
  #
  all_input_para <- c('mat','use_obs_class','class_order','class_ordered',
                      'method','family','pooling',
                      'prior.V.scale','prior.R.nu','prior.G.nu',
                      'nitt','burnin','thin','std','log.base',
                      'logTransformed','average.method','pseudoCount',
                      'return_model','use_seed','verbose')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('class_ordered',c(TRUE,FALSE),envir=environment()),
                 check_option('std',c(TRUE,FALSE),envir=environment()),
                 check_option('logTransformed',c(TRUE,FALSE),envir=environment()),
                 check_option('return_model',c(TRUE,FALSE),envir=environment()),
                 check_option('verbose',c(TRUE,FALSE),envir=environment()),
                 check_option('method',c('Bayesian','MLE'),envir=environment()),
                 check_option('pooling',c('full','no','partial'),envir=environment()),
                 check_option('average.method',c('arithmetic','geometric'),envir=environment())
                 )
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  use_obs_class <- clean_charVector(use_obs_class)
  #
  if (is.character(family)){
    if(family %in% c('category','ordinal')) family <- 'binomial' ## use automatically judging
    if (!family %in% c('gaussian','binomial', 'poisson')) {
      message("Only gaussian,poisson, and binomial/category/ordinal are supported, please check and re-try!");return(FALSE);
    }
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family))
    family <- family()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (!family$family %in% c('gaussian','binomial', 'poisson')) {
    message("Only gaussian,poisson, and binomial/category/ordinal are supported, please check and re-try!");return(FALSE);
  }
  method<-match.arg(method)
  #sample name
  if(is.null(colnames(mat))==TRUE & is.null(names(use_obs_class))==TRUE){
    colnames(mat) <- paste0('Sample',1:ncol(mat))
    names(use_obs_class) <- paste0('Sample',1:ncol(mat))
  }
  if(is.null(colnames(mat))==TRUE & is.null(names(use_obs_class))==FALSE){
    colnames(mat) <- names(use_obs_class)
  }
  if(is.null(colnames(mat))==FALSE & is.null(names(use_obs_class))==TRUE){
    names(use_obs_class) <- colnames(mat)
  }
  use_obs_class <- use_obs_class[colnames(mat)]
  ##check sample class
  class_order  <- base::intersect(class_order,use_obs_class)
  if(base::length(class_order)<=1){message('Sample class order is smaller than two classes, please check and re-try!');return(FALSE);}
  if(verbose==TRUE) message(sprintf('%d sample classes will be used in calculation and %s will be treated as control',base::length(class_order),class_order[1]))
  ##generate comp
  comp <- factor(use_obs_class,levels=class_order)
  ##check input matrix
  #remove NAs
  d <- mat
  nna<-apply(!is.na(d),2,all)
  if(dim(d)[1]==1){ d <- t(as.matrix(d[,nna])); rownames(d) <- rownames(mat) }else{ d <- d[,nna]}
  ##generate data frame
  d   <- reshape::melt(t(d))
  dat <- data.frame(response=d$value,treatment=rep(comp,base::length(base::unique(d$X2))),probe=d$X2)
  ##calculate FC
  n.levels<-base::length(base::unique(dat$probe)) ##
  n.treatments<-base::length(base::unique(dat$treatment)) ## 2grp or mgrp
  AveExpr<-base::mean(dat$response)
  if(n.treatments==2){
    FC.val<-FC(dat$response,dat$treatment,
               logTransformed=logTransformed,log.base=log.base,
               average.method=average.method,pseudoCount=pseudoCount)
    res<-c(FC=FC.val,AveExpr=AveExpr,n.levels=as.integer(n.levels))
  }else{
    res<-c(AveExpr=AveExpr,n.levels=as.integer(n.levels))
  }
  ##re-strandarize the input data
  if(std==TRUE & family$family=='Poisson'){
    message('negative values not allowed for the Poisson family, please set std=FALSE and re-try!');return(FALSE)
  }
  if(std==TRUE & sd(dat$response)>0)
    dat$response<-0.5*(dat$response-base::mean(dat$response))/sd(dat$response)

  if(class_ordered==TRUE) dat$treatment <- as.ordered(dat$treatment) ## for multi-groups
  if(n.treatments==2) class_ordered=FALSE
  ## main part !!!
  # inner functions
  get_info_model<-function(M,dat,method=NULL,df_sta=NULL){
    sum.tmp<-summary(M)
    if('summary.glm' %in% class(sum.tmp) | 'summary.clm'  %in% class(sum.tmp)){
      if('summary.glm' %in% class(sum.tmp) ) w1 <- grep('treatment|response',rownames(sum.tmp$coef))
      if('summary.clm' %in% class(sum.tmp) ) w1 <- grep('\\||response',rownames(sum.tmp$coef))
      if(is.null(df_sta)==TRUE){
        if(method=='MLE') df_sta<-sum.tmp$df.residual else df_sta <- summary(M)$df.residual-summary(M)$df[1] ## ???
      }
      if('z value' %in% colnames(sum.tmp$coef)){ ###
        z_sta <- sum.tmp$coef[w1,'z value']
        p_sta <- sum.tmp$coef[w1,4]
        t_sta <- NA
      } else{
        t_sta <- sum.tmp$coef[w1,'t value']
        if(base::length(grep('^Pr',colnames(sum.tmp$coef)))>0){
          p_sta <- sum.tmp$coef[w1,grep('^Pr',colnames(sum.tmp$coef))]
        }else{
          p_sta<-2*pt(abs(t_sta),lower.tail=FALSE,df=df_sta)
        }
        z_sta<-sign(t_sta)*abs(qnorm(p_sta/2))
      }
      if(verbose==TRUE) print(sum.tmp$coef)
      rs<-c(t=t_sta,'P.Value'=p_sta,'Z-statistics'=z_sta)
    }
    if('summary.merMod' %in% class(sum.tmp)){
      if(verbose==TRUE) print(sum.tmp$coef)
      w1 <- grep('treatment|response',rownames(sum.tmp$coef))
      t_sta<-sum.tmp$coef[w1,3]
      if(is.null(df_sta)==TRUE) df_sta<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
      p_sta<-2*pt(abs(t_sta),lower.tail=FALSE,df=df_sta)
      z_sta<-sign(t_sta)*abs(qnorm(p_sta/2))
      rs<-c(t=t_sta,'P.Value'=p_sta,'Z-statistics'=z_sta)
    }
    if('summary.MCMCglmm' %in% class(sum.tmp)){
      if(verbose==TRUE) print(sum.tmp$sol)
      t_sta<-sum.tmp$sol[2,1]/sd(M$Sol[,2]) ## t
      p_sta<-2*pt(-abs(t_sta),df=df_sta)
      if(is.na(p_sta)) p_sta<-sum.tmp$sol[2,5]
      z_sta<-sign(t_sta)*abs(qnorm(p_sta/2))
      rs<-c(t=t_sta,'P.Value'=p_sta,'Z-statistics'=z_sta)
    }
    return(rs)
  }
  ## check sd
  rs_NA <- c(t=0,'P.Value'=1,'Z-statistics'=0)
  if(sd(dat$response)==0){rs <- rs_NA; return(rs);}
  df_sta <- NULL
  ################### in total: 2(MLE/Bayesian)*3(family)*3(pooling)*2(n.levels)=36 *2(random_effect)=72
  ################### MLE, 6 conditions
  #### gaussian/poisson
  if(method=='MLE' & family$family %in% c('gaussian','poisson') & (n.levels==1 | pooling=='full')){ ## n.levels==1/full
    if(class_ordered==TRUE){message('If need to get p-value between each group, try to set family to "ordinial"!');}
    M <- stats::glm(response ~ treatment, data=dat, family=family)
  }
  if(method=='MLE' & family$family %in% c('gaussian','poisson') & n.levels>1 & pooling=='no'){
    if(class_ordered==TRUE){message('If need to get p-value between each group, try to set family to "ordinial"!');}
    M <- stats::glm(response ~ treatment + probe, data=dat,family=family)
  }
  if(method=='MLE' & family$family %in% c('gaussian','poisson') & n.levels>1 & pooling=='partial'){
    if(class_ordered==TRUE){message('If need to get p-value between each group, try to set family to "ordinial"!');}
    M <- lme4::lmer(response ~ treatment + (treatment + 1 | probe), data=dat)
  }
  #### binomial
  if(method=='MLE' & family$family=='binomial' & (n.levels==1 | pooling=='full')){ ## n.levels==1
    if(class_ordered==FALSE) M <- stats::glm(treatment ~ response, data=dat, family=family)
    if(class_ordered==TRUE) M <- ordinal::clm(treatment ~ response, data=dat)
  }
  if(method=='MLE' & family$family=='binomial' & n.levels>1 & pooling=='no'){
    if(class_ordered==FALSE) M <- stats::glm(treatment ~ response + probe, data=dat,family=family)
    if(class_ordered==TRUE) M <- ordinal::clm(treatment ~ response + probe, data=dat)
  }
  if(method=='MLE' & family$family=='binomial' & n.levels>1 & pooling=='partial'){
    if(class_ordered==FALSE) M <- lme4::lmer(treatment ~ response + (response + 1 | probe), data=dat,family=family)
    if(class_ordered==TRUE){
      if(base::length(levels(dat$probe))<3){message('Random-effect terms has less than three levels, treat as fix effect!');M <- ordinal::clm(treatment ~ response + probe, data=dat)}
      if(base::length(levels(dat$probe))>=3) M <- ordinal::clmm(treatment ~ response + (response+1|probe), data=dat) ## clmm must have grouping factor larger than 3
    }
  }
  ################### Bayesian, 8 conditions
  set.seed(use_seed)
  if(family$link=='logit'){prior.scale<-2.5;glmm.family<-'categorial';}
  if(family$link=='probit'){prior.scale<-2.5*1.6;glmm.family<-'ordinal';}
  #### gaussian/poisson
  if(method=='Bayesian' & family$family %in% c('gaussian','poisson') & (pooling=='full' | pooling=='no' & n.levels==1)){
    if(class_ordered==TRUE){message('If need to get p-value between each group, try to set family to "ordinial"!');}
    M <- arm::bayesglm(response ~ treatment, data=dat, family=family)
  }
  if(method=='Bayesian' & family$family %in% c('gaussian','poisson') & n.levels>1 & pooling=='no'){
    if(class_ordered==TRUE){message('If need to get p-value between each group, try to set family to "ordinial"!');}
    M <- arm::bayesglm(response ~ treatment + probe, data=dat,family=family)
  }
  #
  prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu)) ## default prior
  if(method=='Bayesian' & family$family %in% c('gaussian','poisson') & n.levels==1 & pooling=='partial'){
    if(class_ordered==TRUE){message('If need to get p-value between each group, try to set family to "ordinial"!');}
    M <- MCMCglmm::MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family=family$family)
    df_sta<-nrow(dat)-n.treatments
  }
  if(method=='Bayesian' & family$family %in% c('gaussian','poisson') & n.levels>1 & pooling=='partial'){
    prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(n.treatments)*prior.V.scale, nu = prior.G.nu)))
    if(class_ordered==TRUE){message('If need to get p-value between each group, try to set family to "ordinial"!');}
    M <- MCMCglmm::MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family=family$family)
    df_sta <-nrow(dat)-(n.levels+1)*n.treatments
  }
  #### binomial
  if(method=='Bayesian' & family$family=='binomial' & (pooling=='full' | pooling=='no' & n.levels==1)){ ###
    if(class_ordered==FALSE) M <- arm::bayesglm(treatment ~ response, data=dat, family=family, prior.scale = prior.scale)
    if(class_ordered==TRUE){
      M <- MCMCglmm::MCMCglmm(treatment ~ response, data=dat, family='ordinal',prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
      df_sta<-nrow(dat)-2
    }
  }
  if(method=='Bayesian' & family$family=='binomial' & n.levels>1 & pooling=='no'){ ###
    if(class_ordered==FALSE) M <- arm::bayesglm(treatment ~ response + probe, data=dat, family=family,prior.scale = prior.scale)
    if(class_ordered==TRUE){
      M <- MCMCglmm::MCMCglmm(treatment ~ response + probe, data=dat, family='ordinal',prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
      df_sta <-nrow(dat)-(n.levels+1)*2
    }
  }
  # for partial
  if(method=='Bayesian' & family$family=='binomial' & n.levels==1 & pooling=='partial'){
    if(family$link=='logit'){prior<-list(R = list(V = prior.V.scale, nu=n.levels))}
    if(family$link=='probit'){prior<-list(R = list(V = prior.V.scale, nu=n.levels+1))}
    if(!family$link %in% c('probit','logit')){message('For partial pooling with Binomial family and Bayeisan method, only logit and probit model are supported!');return(FALSE)}
    if(class_ordered==FALSE) M <- MCMCglmm::MCMCglmm(treatment~response,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
    if(class_ordered==TRUE) M <- MCMCglmm::MCMCglmm(treatment~response,family='ordinal', prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
    df_sta<-nrow(dat)-2
  }
  if(method=='Bayesian' & family$family=='binomial' & n.levels>1 & pooling=='partial'){
    if(family$link=='logit'){prior<-list(R = list(V = prior.V.scale, nu=n.levels), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))}
    if(family$link=='probit'){prior<-list(R = list(V = prior.V.scale, nu=n.levels+1), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))}
    if(!family$link %in% c('probit','logit')){message('For partial pooling with Binomial family and Bayeisan method, only logit and probit model are supported!');return(FALSE)}
    if(class_ordered==FALSE) M <- MCMCglmm::MCMCglmm(treatment~response, random=~idh(1+response):probe,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
    if(class_ordered==TRUE) M <- MCMCglmm::MCMCglmm(treatment~response, random=~idh(1+response):probe,family='ordinal', prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
    df_sta <-nrow(dat)-(n.levels+1)*2
  }
  ##
  if(return_model==TRUE) return(M)
  rs <- get_info_model(M=M,dat=dat,method=method,df_sta=df_sta); rs <- c(res,rs)
  return(rs)
}

##
#fold change function, inner function
#first class is 0.
#positive: class1/class0
#negative: class0/class1
FC <- function(x,cl,logTransformed = TRUE,
               log.base = 2,average.method = c('geometric', 'arithmetic'),
               pseudoCount = 0) {
  x.class0 <- x[(cl == levels(cl)[1])] + pseudoCount
  x.class1 <- x[(cl == levels(cl)[2])] + pseudoCount
  if (missing(average.method))
    average.method <- 'geometric'
  if (logTransformed) {
    if (is.na(log.base) | log.base < 0)
      stop('Please specify log.base !\n')
    logFC <- base::mean(x.class1) - base::mean(x.class0)
    FC.val <- sign(logFC) * log.base ^ abs(logFC)
  } else{
    logFC <-
      ifelse(average.method == 'arithmetic',
             log(base::mean(x.class1)) - log(base::mean(x.class0)),
             base::mean(log(x.class1) - base::mean(log(x.class0))))
    FC.val <- sign(logFC) * exp(abs(logFC))
  }
  FC.val[FC.val == 0 | is.na(FC.val)] <- 1
  FC.val
}
################## internal function for graphic
par.pos2inch <- function(){
  user.range <- par("usr")[c(2,4)] - par("usr")[c(1,3)]
  region.pin <- par("pin")
  return(region.pin/user.range)
}
par.inch2pos <- function(){return(1/par.pos2inch())}
par.char2inch <- function(){return(par()$cin)} ## letter W
par.devSize2xypos <- function(){
  pp <- par()$usr
  rr <- par("din")/par.pos2inch()
  mm <- c(pp[1]/2+pp[2]/2,pp[3]/2+pp[4]/2)
  xyp <- c(mm-(rr/2),mm+(rr/2))
  xyp <- xyp[c(1,3,2,4)]
  return(xyp)
}
par.lineHeight2inch <- function(){
  lheight <- par()$lheight
  y1 <- par.char2inch()[2]*lheight ## line height in inches
  y1
}
par.char2pos <- function(){par()$cxy}
strheightMod <- function(s, units = "inch", cex = 1,ori=TRUE,mod=FALSE){
  s <- s[which(is.na(s)==F)]
  if(ori==TRUE) return(strheight(s=s,units=units,cex=cex))
  if(units=='user') return(par.char2pos()[2]*cex)
  if(units=='inch' | units=='inches') return(par.char2inch()[2]*cex)
}
strwidthMod <- function(s, units = "inch", cex = 1,ori=TRUE,mod=FALSE){
  s <- s[which(is.na(s)==F)]
  if(ori==TRUE) return(strwidth(s=s,units=units,cex=cex))
  if(mod==TRUE){
    plot.new()
    rt <- strwidth(s,units=units)/strwidth('W',units=units); rt <- ceiling(rt)
    if(units=='user') r1 <- par.char2pos()[1]*cex*rt
    if(units=='inch') r1 <- par.char2inch()[1]*cex*rt
    dev.off(); return(r1)
  }else{
    if(units=='user') return(par.char2pos()[1]*cex*nchar(s))
    if(units=='inch'| units=='inches') return(par.char2inch()[1]*cex*nchar(s))
  }
}

##
#' Lazy mode for NetBID2 result visualization
#'
#' \code{NetBID.lazyMode.DriverVisualization} is an integrated function to draw visualization plots for top drivers.
#'
#' User need to strictly follow the NetBID2 pipeline to get the complicated list object analysis.par.
#' "intgroup", "use_comp" should be specified;
#' "transfer_tab" could be set to NULL with "main_id_type" specified, but it is suggested to input by hand if available.
#'
#' @param analysis.par list, stores all related datasets from driver analysis step.
#' @param intgroup character, one interested phenotype group from the \code{analysis.par$cal.eset}.
#' @param use_comp character, the name of the comparison of interest, should be included in the colnames of \code{analysis.par$DA} and \code{analysis.par$DE}.
#' @param main_id_type character, the type of driver's ID. It comes from the attribute name in biomaRt package.
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' For details, user can call \code{biomaRt::listAttributes()} to display all available attributes in the selected dataset.
#' @param transfer_tab data.frame, the ID conversion table. Users can call \code{get_IDtransfer} to get this table.
#' @param use_gs a vector of characters, names of major gene set collections. Users can call \code{all_gs2gene_info} to see all the available collections.
#' Default is c("H", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG").
#' @param min_Size numeric, minimum size for the target genes. Default is 30.
#' @param max_Size numeric, maximum size for the target genes. Default is 1000.
#' @param top_number number for the top significant genes/drivers in the combine results to be displayed on the plot.
#' Default is 30.
#' @param top_strategy character, choose from "Both", "Up", "Down".
#' If set to "Both", top drivers with highest absolute Z statistics will be displayed.
#' If set to "Up", only top up-regulated drivers will be displayed.
#' If set to "Down", only top down-regulated drivers will be displayed.
#' Default is "Both".
#' @param logFC_thre numeric, the threshold of logFC. Genes or drivers with absolute logFC value higher than the threshold will be kept.
#' Default is 0.05.
#' @param Pv_thre numeric, the threshold of P-values. Genes or drivers with P-values lower than the threshold will be kept.
#' Default is 0.05.
#'
#' @examples
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' analysis.par$out.dir.PLOT <- 'test/'
#' NetBID.lazyMode.DriverVisualization(analysis.par=analysis.par,
#'                                     intgroup='subgroup',use_comp='G4.Vs.others',
#'                                     transfer_tab=analysis.par$transfer_tab,
#'                                     logFC_thre=0.2,Pv_thre=1e-4)
#' }
#' @export
NetBID.lazyMode.DriverVisualization <- function(analysis.par=NULL,intgroup=NULL,use_comp=NULL,
                                                main_id_type='external_gene_name',
                                                transfer_tab=NULL,use_gs=c("H", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG"),
                                                min_Size=30,max_Size=1000,top_number=30,top_strategy='Both',
                                                logFC_thre=0.05,Pv_thre=0.05){
  all_input_para <- c('analysis.par','intgroup','use_comp','main_id_type',
                      'min_Size','max_Size','top_number','top_strategy','logFC_thre','Pv_thre')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  if('final_ms_tab' %in% names(analysis.par)){
    ms_tab <- analysis.par$final_ms_tab
  }else{
    message('final_ms_tab not included in the analysis.par list, please check and re-try!');return(FALSE)
  }
  if('cal.eset' %in% names(analysis.par)){
    exp_mat <- Biobase::exprs(analysis.par$cal.eset)
  }else{
    message('cal.eset not included in the analysis.par list, please check and re-try!');return(FALSE)
  }
  if('merge.ac.eset' %in% names(analysis.par)){
    ac_mat <- Biobase::exprs(analysis.par$merge.ac.eset)
  }else{
    message('merge.ac.eset not included in the analysis.par list, please check and re-try!');return(FALSE)
  }
  if(!'DE' %in% names(analysis.par)){
    message('DE not included in the analysis.par list, please check and re-try!');return(FALSE)
  }
  if(!'DA' %in% names(analysis.par)){
    message('DA not included in the analysis.par list, please check and re-try!');return(FALSE)
  }
  if(!use_comp %in% names(analysis.par$DA)){
    message(sprintf('%s not included in the analysis.par$DA list, please check and re-try!',use_comp));return(FALSE)
  }
  if(!use_comp %in% names(analysis.par$DE)){
    message(sprintf('%s not included in the analysis.par$DE list, please check and re-try!',use_comp));return(FALSE)
  }
  if(!intgroup %in% colnames(Biobase::pData(analysis.par$cal.eset))){
    message(sprintf('%s not included in the Biobase::pData(analysis.par$cal.eset), please check and re-try!',intgroup));return(FALSE)
  }
  if(exists('db_info')==FALSE){
    message('Please run db.preload() first to load in db_info !');return(FALSE)
  }
  if(is.null(transfer_tab)==TRUE & 'transfer_tab' %in% names(analysis.par)){
    transfer_tab <- analysis.par$transfer_tab
  }
  use_genes = base::unique(c(analysis.par$final_ms_tab$originalID,
                             analysis.par$merge.network$network_dat$target))
  if(is.null(transfer_tab)==TRUE){
    message('Begin get transfer table at gene level (for gene function enrichment analysis )')
    transfer_tab <- get_IDtransfer2symbol2type(from_type=main_id_type,use_genes = use_genes,use_level = 'gene',ignore_version = TRUE)
  }else{
    if(!'external_gene_name' %in% colnames(transfer_tab) & !'external_transcript_name' %in% colnames(transfer_tab)){
      message('Begin get transfer table at gene level (for gene function enrichment analysis )')
      transfer_tab <- get_IDtransfer2symbol2type(from_type=main_id_type,use_genes = use_genes,use_level = 'gene',ignore_version = TRUE)
    }
    if(!'external_gene_name' %in% colnames(transfer_tab) & 'external_transcript_name' %in% colnames(transfer_tab)){
      transfer_tab$external_gene_name <- gsub('(.*)-.*','\\1',transfer_tab$external_transcript_name)
    }
  }
  print(str(transfer_tab))
  if(exists('all_gs2gene')==FALSE){
    message('Please run gs.preload() first to load in all_gs2gene or prepare the same format of this object!');return(FALSE)
  }
  phe <- Biobase::pData(analysis.par$cal.eset)
  DE <- analysis.par$DE[[use_comp]];DA <- analysis.par$DA[[use_comp]]
  G1 <- gsub('Ave.(.*)','\\1',colnames(DE)[grep('^Ave\\.',colnames(DE))])[2]
  G0 <- gsub('Ave.(.*)','\\1',colnames(DE)[grep('^Ave\\.',colnames(DE))])[1]
  #
  w1 <- which(ms_tab$Size>=min_Size & ms_tab$Size<=max_Size)
  message(sprintf('%s out of %s drivers passed the size filteration!',base::length(w1),nrow(ms_tab)))
  if(length(w1)==0){
    message('No drivers passed, will not use size filteration !');
  }else{
    ms_tab <- ms_tab[w1,]
  }
  print('Begin Volcano plot by draw.volcanoPlot() ')
  res1 <- draw.volcanoPlot(dat=ms_tab,label_col = 'gene_label',
                           logFC_col = sprintf('logFC.%s_DA',use_comp),Pv_col = sprintf('P.Value.%s_DA',use_comp),
                           logFC_thre = logFC_thre,Pv_thre = Pv_thre,
                           pdf_file=sprintf('%s/%s_Volcano.pdf',analysis.par$out.dir.PLOT,use_comp))
  print('Finish Volcano plot ')
  if(top_strategy=='Up') res1 <- res1[which(res1[,sprintf('logFC.%s_DA',use_comp)]>0),,drop=FALSE]
  if(top_strategy=='Down') res1 <- res1[which(res1[,sprintf('logFC.%s_DA',use_comp)]<0),,drop=FALSE]
  driver_list <- rownames(res1)
  driver_DE_Z <- ms_tab[driver_list,sprintf('Z.%s_DE',use_comp)]
  driver_DA_Z <- ms_tab[driver_list,sprintf('Z.%s_DA',use_comp)]
  if(base::length(driver_list)>=top_number){
    driver_list <- driver_list[order(abs(driver_DA_Z),decreasing = TRUE)[1:top_number]]
  }
  driver_DE_Z <- ms_tab[driver_list,sprintf('Z.%s_DE',use_comp)];names(driver_DE_Z) <- driver_list
  driver_DA_Z <- ms_tab[driver_list,sprintf('Z.%s_DA',use_comp)];names(driver_DA_Z) <- driver_list
  if(nrow(res1)<5){
    message(sprintf('%s drivers remained, too few for draw top drivers, will skip plots for top drivers !',nrow(res1)));
  }else{
    message(sprintf('%s drivers (%s) will be displayed!',base::length(driver_list),base::paste(driver_list,collapse=';')))
    print('Begin TopDA_GSEA plot by draw.GSEA.NetBID() ')
    draw.GSEA.NetBID(DE=DE,profile_col = 'logFC',name_col = 'ID',
                     driver_list = driver_list,show_label=ms_tab[driver_list,'gene_label'],
                     driver_DE_Z = driver_DE_Z,
                     driver_DA_Z = driver_DA_Z,
                     target_list = analysis.par$merge.network$target_list,
                     pdf_file=sprintf('%s/%s_TopDA_GSEA.pdf',analysis.par$out.dir.PLOT,use_comp),
                     top_driver_number = top_number,main=use_comp)
    print(sprintf('Finish TopDA_GSEA plot by draw.GSEA.NetBID(), please check %s',sprintf('%s/%s_TopDA_GSEA.pdf',analysis.par$out.dir.PLOT,use_comp)))
    print('Begin TopDA_Heatmap_ac plot by draw.heatmap() ')
    draw.heatmap(mat=ac_mat,use_genes = ms_tab[driver_list,'originalID_label'],use_gene_label = ms_tab[driver_list,'gene_label'],
                 phenotype_info  = Biobase::pData(analysis.par$merge.ac.eset),use_phe=intgroup,
                 pdf_file=sprintf('%s/%s_TopDA_Heatmap_ac.pdf',analysis.par$out.dir.PLOT,use_comp),scale='row')
    print(sprintf('Finish TopDA_Heatmap_ac plot by draw.heatmap(), please check %s',sprintf('%s/%s_TopDA_Heatmap_ac.pdf',analysis.par$out.dir.PLOT,use_comp)))
    print('Begin TopDA_Heatmap_exp plot by draw.heatmap() ')
    draw.heatmap(mat=exp_mat,use_genes = ms_tab[driver_list,'originalID'],use_gene_label = ms_tab[driver_list,'geneSymbol'],
                 phenotype_info  = Biobase::pData(analysis.par$merge.ac.eset),use_phe=intgroup,
                 pdf_file=sprintf('%s/%s_TopDA_Heatmap_exp.pdf',analysis.par$out.dir.PLOT,use_comp),scale='row')
    print(sprintf('Finish TopDA_Heatmap_exp plot by draw.heatmap(), please check %s',sprintf('%s/%s_TopDA_Heatmap_exp.pdf',analysis.par$out.dir.PLOT,use_comp)))
    print('Begin TopDA_FuncEnrich plot by funcEnrich.Fisher() and draw.funcEnrich.cluster()')
    res1 <- funcEnrich.Fisher(input_list = ms_tab[driver_list,'geneSymbol'],bg_list = ms_tab$geneSymbol,
                              Pv_thre = 0.1,use_gs=use_gs,Pv_adj = 'none')
    out2excel(res1,out.xlsx = sprintf('%s/%s_TopDA_FuncEnrich.xlsx',analysis.par$out.dir.PLOT,use_comp))
    print(sprintf('Finish TopDA_FuncEnrich plot by funcEnrich.Fisher(), please check %s',
                  sprintf('%s/%s_TopDA_FuncEnrich.xlsx',analysis.par$out.dir.PLOT,use_comp)))
    if(nrow(res1)>=3){
      draw.funcEnrich.cluster(res1,pdf_file=sprintf('%s/%s_TopDA_FuncEnrich.pdf',analysis.par$out.dir.PLOT,use_comp),top_number = top_number,
                              Pv_thre = 0.1)
      print(sprintf('Finish TopDA_FuncEnrich plot by draw.funcEnrich.cluster(), please check %s',
                    sprintf('%s/%s_TopDA_FuncEnrich.pdf',analysis.par$out.dir.PLOT,use_comp)))
    }else{
      message('Too few results for Function enrichment analysis for top drivers, will pass the FuncEnrich Plot!')
    }
    print('Begin TopDA_BubblePlot plot by draw.bubblePlot() ')
    draw.bubblePlot(driver_list = driver_list,show_label = ms_tab[driver_list,'gene_label'],
                    transfer2symbol2type = transfer_tab,
                    target_list=analysis.par$merge.network$target_list,
                    Z_val=driver_DA_Z,Pv_thre=0.1,top_geneset_number = top_number,
                    top_driver_number = top_number,use_gs=use_gs,
                    pdf_file=sprintf('%s/%s_TopDA_BubblePlot.pdf',analysis.par$out.dir.PLOT,use_comp))
    print(sprintf('Finish TopDA_BubblePlot plot by draw.bubblePlot(), please check %s',sprintf('%s/%s_TopDA_BubblePlot.pdf',analysis.par$out.dir.PLOT,use_comp)))
  }
  ## detailed for top
  pf <- DE$logFC; names(pf)<-DE$ID
  print('Begin Each TopDA_GSEA plot by draw.GSEA() ')
  pdf(sprintf('%s/%s_EachTopDA_GSEA.pdf',analysis.par$out.dir.PLOT,use_comp),width=8,height=8)
  for(each_driver in driver_list){
    print(each_driver)
    use_direction <- base::sign(analysis.par$merge.network$target_list[[each_driver]]$spearman)
    if(length(unique(use_direction))==1){
      if(unique(use_direction)==1){
        use_direction <- NULL
      }
    }
    draw.GSEA(rank_profile = pf,use_genes = analysis.par$merge.network$target_list[[each_driver]]$target,
              use_direction = use_direction,
              annotation=sprintf('P.Value:%s',get_z2p(driver_DA_Z[each_driver])),
              left_annotation = sprintf('High in %s',G1),
              right_annotation = sprintf('High in %s',G0),main=ms_tab[each_driver,'gene_label'])
  }
  dev.off()
  print(sprintf('Finish EachTopDA_GSEA plot by draw.GSEA(), please check %s',sprintf('%s/%s_EachTopDA_GSEA.pdf',analysis.par$out.dir.PLOT,use_comp)))
  #
  print('Begin EachTopDA_CateBox plot by draw.categoryValue() ')
  pdf(sprintf('%s/%s_EachTopDA_CateBox.pdf',analysis.par$out.dir.PLOT,use_comp),width=8,height=8)
  for(each_driver in driver_list){
    if(ms_tab[each_driver,'originalID'] %in% rownames(exp_mat)){
      draw.categoryValue(ac_val   = ac_mat[ms_tab[each_driver,'originalID_label'],],
                         exp_val  = exp_mat[ms_tab[each_driver,'originalID'],],
                         main_ac  = sprintf("%s\n(P.Value:%s)",ms_tab[each_driver,'gene_label'],get_z2p(driver_DA_Z[each_driver])),
                         main_exp = sprintf("%s\n(P.Value:%s)",ms_tab[each_driver,'geneSymbol'],get_z2p(driver_DE_Z[each_driver])),
                         use_obs_class=get_obs_label(phe,intgroup))
    }else{
      draw.categoryValue(ac_val   = ac_mat[ms_tab[each_driver,'originalID_label'],],
                         main_ac  = sprintf("%s\n(P.Value:%s)",ms_tab[each_driver,'gene_label'],get_z2p(driver_DA_Z[each_driver])),
                         use_obs_class=get_obs_label(phe,intgroup))
    }
  }
  dev.off()
  print(sprintf('Finish EachTopDA_CateBox plot by draw.categoryValue(), please check %s',sprintf('%s/%s_EachTopDA_CateBox.pdf',analysis.par$out.dir.PLOT,use_comp)))
  #
  print('Begin EachTopDA_TargetNet plot by draw.targetNet() ')
  pdf(sprintf('%s/%s_EachTopDA_TargetNet.pdf',analysis.par$out.dir.PLOT,use_comp),width=8,height=8)
  for(each_driver in driver_list){
    if('spearman' %in% colnames(analysis.par$merge.network$target_list[[each_driver]])){
      es <- analysis.par$merge.network$target_list[[each_driver]]$MI*sign(analysis.par$merge.network$target_list[[each_driver]]$spearman);
    }else{
      es <- rep(1,length.out=length(analysis.par$merge.network$target_list[[each_driver]]$target))
    }
    names(es) <- analysis.par$merge.network$target_list[[each_driver]]$target
    if(main_id_type!='external_gene_name'){
      es <- base::cbind(es,es)
      colnames(es) <- c('Sample1','Sample2')
      tmp_eset <- generate.eset(es)
      if('external_transcript_name' %in% names(transfer_tab)){
        tmp_eset <- update_eset.feature(tmp_eset,use_feature_info = transfer_tab,
                                        from_feature = main_id_type,to_feature = 'external_transcript_name')
      }else{
        tmp_eset <- update_eset.feature(tmp_eset,use_feature_info = transfer_tab,
                                        from_feature = main_id_type,to_feature = 'external_gene_name')
      }
      es <- Biobase::exprs(tmp_eset)[,1];names(es) <- rownames(Biobase::exprs(tmp_eset))
    }
    if(base::length(es)<30){n_layer=1;label_cex=1;}
    if(base::length(es)<50 & base::length(es)>=30){n_layer=1;label_cex=0.8;}
    if(base::length(es)>=50){n_layer <- 1+floor(base::length(es)/50);label_cex=0.7;}
    if(n_layer>=4){n_layer <- n_layer/2; label_cex<-label_cex-0.1;}
    if(n_layer>=4){n_layer <- n_layer/2; label_cex<-label_cex-0.1;}
    draw.targetNet(source_label = ms_tab[each_driver,'gene_label'],source_z = driver_DA_Z[each_driver],
                   edge_score=es,n_layer = n_layer,label_cex=label_cex)
  }
  dev.off()
  print(sprintf('Finish EachTopDA_TargetNet plot by draw.targetNet(), please check %s',sprintf('%s/%s_EachTopDA_TargetNet.pdf',analysis.par$out.dir.PLOT,use_comp)))
  message(sprintf('Finish All !!! Check %s',analysis.par$out.dir.PLOT))
  return(TRUE)
}

##
#' Lazy mode for NetBID2 driver estimation
#'
#' \code{NetBID.lazyMode.DriverEstimation} is an integrated function for NetBID2 driver estimation.
#'
#' The function will return the complicated list object analysis.par if set return_analysis.par=TRUE.
#' Meanwhile, the master table and the RData containing analysis.par will be automatically saved.
#'
#' @param project_main_dir character, name or absolute path of the main working directory for driver analysis.
#' @param project_name character, name of the project folder.
#' @param tf.network.file character, the path of the TF network file (e.g. "XXX/consensus_network_ncol_.txt").
#' @param sig.network.file character, the path of the SIG network file (e.g. "XXX/consensus_network_ncol_.txt").
#' @param cal.eset ExpressionSet class, the ExpressionSet for analysis.
#' @param main_id_type character, the type of driver's ID. It comes from the attribute name in biomaRt package.
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' For details, user can call \code{biomaRt::listAttributes()} to display all available attributes in the selected dataset.
#' @param cal.eset_main_id_type character, the type of cal.eset's ID. It comes from the attribute name in biomaRt package.
#' Such as "ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version" or "refseq_mrna".
#' For details, user can call \code{biomaRt::listAttributes()} to display all available attributes in the selected dataset.
#' @param use_level character, users can choose "transcript" or "gene". Default is "gene".
#' @param transfer_tab data.frame, the ID conversion table. Users can call \code{get_IDtransfer} to get this table.
#' Only useful when "cal.eset_main_id_type" does not equal to "main_id_type".
#' This is mainly for converting ID for cal.eset, must include "cal.eset_main_id_type" and "main_id_type".
#' If NULL, will automatically generate it.
#' @param intgroup character, one interested phenotype group from the \code{cal.eset}.
#' @param G1_name character, the name of experimental group (e.g. "Male"), must be the character in \code{intgroup}.
#' @param G0_name character, the name of control group (e.g. "Female"), must be the character in \code{intgroup}.
#' @param comp_name character, the name of the comparison of interest.
#' @param do.QC logical, if TRUE, will perform network QC and activity eSet QC plots. Default is TRUE.
#' @param DE_strategy character, use limma or bid to calculate differentiated expression/activity. Default is 'bid'.
#' @param return_analysis.par logical, if TRUE, will return the complicated list object analysis.par.
#'
#' @examples
#' \dontrun{
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' tf.network.file <- sprintf('%s/SJAR/%s/output_tf_sjaracne_%s_out_.final/%s',
#'                 network.dir,network.project.name,network.project.name,
#'                 'consensus_network_ncol_.txt')
#' sig.network.file <- sprintf('%s/SJAR/%s/output_sig_sjaracne_%s_out_.final/%s',
#'                 network.dir,network.project.name,network.project.name,
#'                 'consensus_network_ncol_.txt')
#' load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir))
#' cal.eset <- network.par$net.eset
#' db.preload(use_level='gene')
#' analysis.par <- NetBID.lazyMode.DriverEstimation(project_main_dir=project_main_dir,
#'                 project_name=project_name,
#'                 tf.network.file=tf.network.file,
#'                 sig.network.file=sig.network.file,
#'                 cal.eset=cal.eset,
#'                 main_id_type='external_gene_name',
#'                 cal.eset_main_id_type='external_gene_name',
#'                 intgroup='subgroup',
#'                 G1_name='G4',G0_name='SHH',
#'                 comp_name='G4.Vs.SHH',
#'                 do.QC=FALSE,return_analysis.par=TRUE)
#' }
#' @export
NetBID.lazyMode.DriverEstimation <- function(project_main_dir=NULL,project_name=NULL,
                                             tf.network.file=NULL,sig.network.file=NULL,
                                             cal.eset=NULL,
                                             main_id_type=NULL,cal.eset_main_id_type=NULL,use_level='gene',
                                             transfer_tab=NULL,
                                             intgroup=NULL,G1_name=NULL,G0_name=NULL,comp_name=NULL,
                                             do.QC=TRUE,DE_strategy='bid',return_analysis.par=TRUE){
  #
  if(exists('analysis.par')==TRUE){
    stop('analysis.par is occupied in the current session,please manually run: rm(analysis.par) and re-try, otherwise will not change !');
  }
  all_input_para <- c('project_main_dir','project_name','tf.network.file','sig.network.file','cal.eset',
                      'main_id_type','cal.eset_main_id_type','intgroup','G1_name','G0_name','DE_strategy','use_level')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}

  check_res <- c(check_option('DE_strategy',c('limma','bid'),envir=environment()),
                 check_option('use_level',c('gene','transcript'),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}

  if(exists('tf_sigs')==FALSE){
    message('Please run db.preload() first to load in db_info !');return(FALSE)
  }
  if('ensembl_transcript_id' %in% names(tf_sigs$tf)){
    message('You setting is at transcript level! If this setting is not TRUE, please run db.preload() again !')
  }else{
    message('You setting is at gene level! If this setting is not TRUE, please run db.preload() again !')
  }
  print('Current db info:')
  print(db_info)
  if(!intgroup %in% colnames(Biobase::pData(cal.eset))){
    message(sprintf('%s not included in the Biobase::pData(cal.eset), please check and re-try!',intgroup));return(FALSE)
  }
  phe <- Biobase::pData(cal.eset)
  G1 <- rownames(phe)[which(phe[,intgroup]==G1_name)];G0 <- rownames(phe)[which(phe[,intgroup]==G0_name)];
  if(base::length(G1)==0){message(sprintf('NO Sample annotated by %s in %s, please check and re-try!',G1_name,intgroup));return(FALSE)}
  if(base::length(G0)==0){message(sprintf('NO Sample annotated by %s in %s, please check and re-try!',G0_name,intgroup));return(FALSE)}
  if(is.null(comp_name)==TRUE) comp_name <- sprintf('%s.Vs.%s',G1_name,G0_name)
  if(file.exists(tf.network.file)==FALSE){message('%s not exists, please check and re-try!');return(FALSE)}
  if(file.exists(sig.network.file)==FALSE){message('%s not exists, please check and re-try!');return(FALSE)}
  # create workspace
  print('Begin create workspace by NetBID.analysis.dir.create() ')
  analysis.par <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
                                             project_name=project_name,
                                             tf.network.file=tf.network.file,
                                             sig.network.file=sig.network.file)
  print('Finish create workspace')
  # read in network
  print('Begin read in network from file by get.SJAracne.network() ')
  analysis.par$tf.network <- get.SJAracne.network(analysis.par$tf.network.file)
  analysis.par$sig.network <- get.SJAracne.network(analysis.par$sig.network.file)
  print('Finish read in network from file')
  if(do.QC==TRUE & nrow(analysis.par$tf.network$network_dat)>0) draw.network.QC(analysis.par$tf.network$igraph_obj,outdir = analysis.par$out.dir.QC,prefix = 'TF_')
  if(do.QC==TRUE & nrow(analysis.par$sig.network$network_dat)>0) draw.network.QC(analysis.par$sig.network$igraph_obj,outdir = analysis.par$out.dir.QC,prefix = 'SIG_')
  analysis.par$merge.network <- merge_TF_SIG.network(analysis.par$tf.network,analysis.par$sig.network)
  if(cal.eset_main_id_type!=main_id_type){
    message(sprintf('The ID type for the network is %s, and the ID type for the cal.eset is %s, need to transfer the ID type of cal.eset from %s to %s',
                    main_id_type,cal.eset_main_id_type,cal.eset_main_id_type,main_id_type))
    if(is.null(transfer_tab)==FALSE){
      w1 <- base::setdiff(c(cal.eset_main_id_type,main_id_type),colnames(transfer_tab))
      if(base::length(w1)>0){
        message(sprintf('Wrong ID type for the input transfer_tab, missing %s, try to prepare the correct transfer_tab or set it to NULL',base::paste(w1,collapse=';')));
        return(FALSE)
      }
    }
    if(is.null(transfer_tab)==FALSE){
      transfer_tab1 <- transfer_tab
    }else{
      transfer_tab1 <- get_IDtransfer(from_type = cal.eset_main_id_type,to_type=main_id_type,use_genes=rownames(Biobase::exprs(cal.eset)),
                                      ignore_version=TRUE)
    }
    cal.eset <- update_eset.feature(cal.eset,use_feature_info = transfer_tab1,
                                    from_feature = cal.eset_main_id_type,to_feature = main_id_type)
    message('Finish transforming the cal.eset ID ! ')
  }
  es.method <- 'mean'
  if('MI' %in% colnames(analysis.par$merge.network$network_dat) & 'spearman' %in% colnames(analysis.par$merge.network$network_dat)) es.method <- 'weightedmean'
  print(sprintf('Begin calculate activity by cal.Activity(), the setting for es.method is %s ',es.method))
  ac_mat <- cal.Activity(igraph_obj = analysis.par$merge.network$igraph_obj,
                         cal_mat = Biobase::exprs(cal.eset),es.method=es.method)
  analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info = phe)
  if(do.QC==TRUE) draw.eset.QC(analysis.par$merge.ac.eset,outdir = analysis.par$out.dir.QC,prefix = 'AC_')
  print('Finish calculate activity ')
  print(sprintf('Begin prepare ID transfer table to external_%s_name',use_level))
  use_genes = base::unique(c(analysis.par$final_ms_tab$originalID,
                             analysis.par$merge.network$network_dat$target))
  transfer_tab <- get_IDtransfer2symbol2type(from_type=main_id_type,
                                             use_genes = use_genes,
                                             use_level=use_level,ignore_version=TRUE)
  print(str(transfer_tab))
  print('Finish prepare ID transfer table ')
  # DE/DA
  print(sprintf('Begin get DE/DA for %s by getDE.limma.2G() ',comp_name))
  if(DE_strategy=='limma'){
    de <- getDE.limma.2G(cal.eset,G1 = G1,G0=G0,G1_name=G1_name,G0_name=G0_name)
    da <- getDE.limma.2G(analysis.par$merge.ac.eset,G1 = G1,G0=G0,G1_name=G1_name,G0_name=G0_name)
  }else{
    de <- getDE.BID.2G(cal.eset,G1 = G1,G0=G0,G1_name=G1_name,G0_name=G0_name)
    da <- getDE.BID.2G(analysis.par$merge.ac.eset,G1 = G1,G0=G0,G1_name=G1_name,G0_name=G0_name)
  }
  DE <- list(de); DA <- list(da); names(DE) <- names(DA) <- comp_name
  print(sprintf('Finish get DE/DA for %s',comp_name))
  #
  print('Begin generate master table by generate.masterTable() ')
  ms_tab <- generate.masterTable(use_comp=comp_name,DE=DE,DA=DA,target_list = analysis.par$merge.network$target_list,
                                 main_id_type = main_id_type,transfer_tab=transfer_tab,tf_sigs = tf_sigs)
  print('Finish generate master table')
  out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
  out2excel(ms_tab,out.xlsx = out_file)
  message(sprintf('Finish output master table to excel file, please check %s',out_file))
  analysis.par$DE <- DE;analysis.par$DA <- DA;
  analysis.par$final_ms_tab <- ms_tab;analysis.par$cal.eset <- cal.eset
  analysis.par$transfer_tab <- transfer_tab
  NetBID.saveRData(analysis.par = analysis.par,step = 'ms-tab')
  message(sprintf('Finish save RData to file, please check %s/analysis.par.Step.ms-tab.RData',analysis.par$out.dir.DATA))
  if(return_analysis.par==TRUE) return(analysis.par) else return(TRUE)
}

#' Draw Oncoprint Plot to display the mutation information
#'
#' \code{draw.oncoprint} plots the heatmap to display the mutation information for samples.
#'
#' @param phenotype_info data.frame, phenotype of samples. Users can call \code{Biobase::pData(eset)} to create.
#' @param Missense_column character, column names from \code{phenotype_info}.
#' @param Missense_label, character, gene label for the Missense_column.
#' @param Amplification_column character, column names from \code{phenotype_info}.
#' @param Amplification_label, character, gene label for the Amplification_column.
#' @param Deletion_column character, column names from \code{phenotype_info}.
#' @param Deletion_label, character, gene label for the mDeletion_column.
#' @param Sample_column character, column names from \code{phenotype_info}. If not NULL, will show sample names in the figure.
#' @param main character, an overall title for the plot. Default is "".
#' @param pdf_file character, the file path to save plot as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#' @param ..., for more options, please check \code{?oncoPrint} for more details.
#'
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#'
#' @examples
#' all_sample <- sprintf('Sample%s',1:30)
#' group1_sample <- sample(all_sample,18) ## demo sample for KRAS missense mutation
#' group2_sample <- sample(all_sample,12) ## demo sample for MYC amplification
#' group3_sample <- sample(all_sample,10) ## demo sample for MYC missense mutation
#' group4_sample <- sample(all_sample,1) ## demo sample for MYC deletion
#' phenotype_info_demo <-
#'   data.frame(sample =sprintf('Sample%s',1:30),
#'              KRAS_MIS=ifelse(all_sample %in% group1_sample,1,0),
#'              MYC_AMP=ifelse(all_sample %in% group2_sample,1,0),
#'              MYC_MIS=ifelse(all_sample %in% group3_sample,1,0),
#'              MYC_DEL=ifelse(all_sample %in% group4_sample,1,0))
#' draw.oncoprint(phenotype_info=phenotype_info_demo,
#'                Missense_column=c('KRAS_MIS','MYC_MIS'),Missense_label=c('KRAS','MYC'),
#'                Amplification_column=c('MYC_AMP'),Amplification_label=c('MYC'),
#'                Deletion_column=c('MYC_DEL'),Deletion_label=c('MYC'),
#'                Sample_column='sample',
#'                main="OncoPrint for the demo dataset")
#' \dontrun{
#'}
#' @export
draw.oncoprint <- function(phenotype_info=NULL,
                           Missense_column=NULL,Missense_label=NULL,
                           Amplification_column=NULL,Amplification_label=NULL,
                           Deletion_column=NULL,Deletion_label=NULL,
                           Sample_column=NULL,
                           main="",pdf_file=NULL,
                           ...){
  #
  all_input_para <- c('phenotype_info')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(Missense_column)==TRUE & is.null(Amplification_column)==TRUE
     & is.null(Deletion_column)==TRUE){
    stop("At least one of the Missense/Amplication/Deletion columns is required")
  }
  if(length(Missense_column)!=length(Missense_label)){
    stop('the length of Missense_column and Missense_label must be the same')
  }
  if(length(Amplification_column)!=length(Amplification_label)){
    stop('the length of Amplification_column and Amplification_label must be the same')
  }
  if(length(Deletion_column)!=length(Deletion_label)){
    stop('the length of Deletion_column and Deletion_label must be the same')
  }
  sample_name <- NULL;
  if(is.null(Sample_column)==FALSE){
    if(Sample_column %in% colnames(phenotype_info)){
      sample_name <- phenotype_info[,Sample_column[1]]
    }
  }
  ##
  use_col <- rep(0,length.out=3);names(use_col) <- c('AMP','MIS','DEL')
  all_col <- c(Missense_label,Amplification_label,Deletion_label)
  uni_col <- unique(all_col)
  ## each row is one label(gene), each column is one sample
  mat <- matrix("  ",nrow=length(uni_col),ncol=nrow(phenotype_info));
  colnames(mat) <- rownames(phenotype_info);
  rownames(mat) <- uni_col;
  if(is.null(Missense_column)==FALSE){
    use_col['MIS'] <- 1;
    for(i in 1:length(Missense_column)){
      w1 <- which(phenotype_info[,Missense_column[i]]!=0)
      for(w2 in w1){
        if(mat[Missense_label[i],w2]==''){
          mat[Missense_label[i],w2] <- 'MIS;'
        }else{
          mat[Missense_label[i],w2] <- sprintf('%s%s',mat[Missense_label[i],w2],'MIS;')
        }
      }
    }
  }
  if(is.null(Amplification_column)==FALSE){
    use_col['AMP'] <- 1;
    for(i in 1:length(Amplification_column)){
      w1 <- which(phenotype_info[,Amplification_column[i]]!=0)
      for(w2 in w1){
        if(mat[Amplification_label[i],w2]==''){
          mat[Amplification_label[i],w2] <- 'AMP;'
        }else{
          mat[Amplification_label[i],w2] <- sprintf('%s%s',mat[Amplification_label[i],w2],'AMP;')
        }
      }
    }
  }
  if(is.null(Deletion_column)==FALSE){
    use_col['DEL'] <- 1;
    for(i in 1:length(Deletion_column)){
      w1 <- which(phenotype_info[,Deletion_column[i]]!=0)
      for(w2 in w1){
        if(mat[Deletion_label[i],w2]==''){
          mat[Deletion_label[i],w2] <- 'DEL;'
        }else{
          mat[Deletion_label[i],w2] <- sprintf('%s%s',mat[Deletion_label[i],w2],'DEL;')
        }
      }
    }
  }
  use_colname <- names(use_col)[which(use_col==1)]
  use_colname_full <- c('Amplification','Missense','Deletion')[which(use_col==1)]
  #
  col = c(AMP=brewer.pal(8,'Set1')[1], ## red
          MIS=brewer.pal(8,'Set1')[2], ## blue
          DEL=brewer.pal(8,'Set1')[4], ## purple
          OTHER='light grey')
  #########################################################
  if(is.null(pdf_file)==FALSE){
    ww <- 0.15*ncol(mat)+4
    hh <- nrow(mat)+2
    pdf(pdf_file,width=ww,height=hh)
  }
  #
  dxx <- 0.3
  xx <- seq(1,ncol(mat));
  if(is.null(sample_name)==FALSE){
    plot(1,col='white',xlim=c(1,ncol(mat)*1.2),ylim=c(0,1+nrow(mat)),bty='n',
         xaxt='n',yaxt='n',xlab='',ylab='',main=main)
    text(x=xx,y=1-0.05,
         sample_name,srt=60,xpd=T,adj=1)
  }else{
    plot(1,col='white',xlim=c(1,ncol(mat)*1.2),ylim=c(1,1+nrow(mat)),bty='n',
         xaxt='n',yaxt='n',xlab='',ylab='',main=main)
  }
  for(i in 1:nrow(mat)){
    rect(xleft=xx-dxx,xright=xx+dxx,
         ybottom=i,ytop=i+0.8,xpd=TRUE,border=NA,col=col['OTHER'])
    text(1,i+0.4,rownames(mat)[i],pos=2,xpd=T)
    mm <- mat[i,]
    w1 <- grep('AMP',mm)
    if(length(w1)>0) rect(xleft=xx[w1]-dxx,xright=xx[w1]+dxx,
                          ybottom=i,ytop=i+0.8,xpd=TRUE,border=NA,col=col['AMP'])
    w1 <- grep('DEL',mm)
    if(length(w1)>0) rect(xleft=xx[w1]-dxx,xright=xx[w1]+dxx,
                          ybottom=i,ytop=i+0.8,xpd=TRUE,border=NA,col=col['DEL'])
    w1 <- grep('MIS',mm)
    if(length(w1)>0) rect(xleft=xx[w1]-dxx,xright=xx[w1]+dxx,
                          ybottom=i+0.3,ytop=i+0.5,xpd=TRUE,border=NA,col=col['MIS'])
  }
  #
  pp <- par()$usr
  legend(y=nrow(mat)/2+0.5,x=1+ncol(mat),
         fill = col[use_colname],use_colname_full,xpd=T,border=NA,bty='n',
         yjust=0.5,xjust=0)
  #
  if(is.null(pdf_file)==FALSE) {while (!is.null(dev.list()))  dev.off();}
  return(TRUE)
}
