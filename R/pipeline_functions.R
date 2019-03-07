
#' @import Biobase limma tximport igraph biomaRt openxlsx msigdbr ConsensusClusterPlus
#' @importFrom GEOquery getGEO
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plot3D scatter3D
#' @importFrom vsn meanSdPlot
#' @importFrom plotrix draw.ellipse
#' @importFrom impute impute.knn
#' @importFrom umap umap
#' @importFrom rhdf5 H5Fopen H5Fclose
#' @importFrom plyr ddply
#' @importFrom DESeq2 DESeqDataSetFromTximport DESeq
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom reshape melt
#' @importFrom graphics plot
#' @importFrom aricode clustComp

##   ‘MCMCglmm’ ‘arm’ ‘reshape’
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
#library(gplots) ## for some plot like venn
library(igraph) ## for network related functions
library(plotrix) ## for draw.ellipse

library(biomaRt) ## for gene id conversion
library(openxlsx) ## for output into excel
library(impute) ## for impute

library(msigdbr) ## for msigDB gene sets
library(ComplexHeatmap) ## for complex heatmap
library(umap) ## for umap visualization
library(rhdf5) ## for read in MICA results
#source('pipeline_functions_bid.R')

####
#' Load database files used for NetBID2 into R workspace.
#'
#' \code{db.preload} returns the TF (transcription factors) and Sig (signaling factors) list (tf_sigs)
#' and biomart database information (db_info) by input interested species name and make choice from
#' gene or transcript level.
#'
#' This is a pre-processing function for NetBID2, user could input the species name (e.g human, mouse),
#' analysis level (transcript or gene level) and optionally input TF list or SIG list
#' (otherwise will use list from package data). The function could automatically download information
#' from biomart and save into RData under the db/ directory with specified species name and analysis level.
#'
#' @param use_level character, either 'transcript' or 'gene', default is 'gene'
#' @param use_spe character, input the species name (e.g 'human', 'mouse', 'rat'), default is 'human'
#' @param update logical,whether to update if previous RData has been generated, default FALSE
#' @param TF_list a character vector,input the list of TF names, if NULL, will use pre-defined list in the package, default NULL
#' @param SIG_list a character vector,input the list of SIG names, if NULL, will use pre-defined list in the package, default NULL
#' @param input_attr_type character, input the type for the list of TF_list, SIG_list. If no input for TF_list, SIG_list, just leave it to NULL, default NULL.
#' See biomaRt \url{https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html} for more details.
#' @param main.dir character, main file path for NetBID2,
#' if NULL, will set to \code{system.file(package = "NetBID2")}. Default is NULL.
#' @param db.dir character, file path for saving the RData, default is \code{db} directory under \code{main.dir} when setting for \code{main.dir}.
#'
#' @return Reture TRUE if success and FALSE if not. Will load two variables into R workspace, tf_sigs and db_info.
#'
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
                       db.dir=sprintf("%s/db/",main.dir)){
  ## load annotation info, including: TF/Sig list, gene info
  if(is.null(main.dir)==TRUE){
    main.dir <- system.file(package = "NetBID2")
    message(sprintf('main.dir not set, will use package directory: %s',main.dir))
  }
  if(is.null(db.dir)==TRUE){
    db.dir <- sprintf("%s/db/",main.dir)
  }
  message(sprintf('Will use directory %s as the db.dir',db.dir))
  use_spe <- toupper(use_spe)
  output.db.dir <- sprintf('%s/%s',db.dir,use_spe)
  if(!file.exists(output.db.dir)){
    dir.create(output.db.dir)
  }
  RData.file <- sprintf('%s/%s_%s.RData', output.db.dir,use_spe,use_level)
  if (update == TRUE | !file.exists(RData.file)) {
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
                                 'gene_biotype','gene_biotype',
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
      tf_sigs$tf[[each_id_type]] <- setdiff(unique(TF_info[[each_id_type]]),"")
      tf_sigs$sig[[each_id_type]] <- setdiff(unique(SIG_info[[each_id_type]]),"")
    }
    db_info <- all_ds[w1,]
    save(tf_sigs,db_info=db_info,file = RData.file)
  }
  load(RData.file,.GlobalEnv)
  return(TRUE)
}

#' Get transcription factor (TF) and signaling factor (SIG) list for the input gene/transcript type
#'
#' \code{get.TF_SIG.list} is a gene ID conversion function to get the TF/SIG list
#' for the input gene list with selected gene/transcript type.
#'
#' @param use_genes a vector of characters, all possible genes used in network generation.
#' If NULL, will not filter the TF/SIG list by this gene list. Default is NULL.
#' @param use_gene_type character, attribute name from the biomaRt package,
#' such as 'ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version', 'refseq_mrna'.
#' Full list could be obtained by
#' \code{listAttributes(mart)$name}, where \code{mart} is the output of \code{useMart} function.
#' (e.g mart <- useMart('ensembl',db_info[1]))
#' The type must be the gene type for the input \code{use_genes}. Default is 'ensembl_gene_id'.
#' @param dataset character, name for the dataset used for ID conversion,
#' such as 'hsapiens_gene_ensembl'.
#' If NULL, will use \code{db_info[1]} if run \code{db.preload} brefore. Default is NULL.
#'
#'
#' @return This function will return a list containing two parts,
#' for $tf saving the TF list for the input gene type and $sig saving the SIG list.
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
                            use_gene_type='ensembl_gene_id',
                            dataset=NULL){
  if(is.null(dataset)==TRUE){
    dataset <- db_info[1]
  }
  if(is.null(tf_sigs)==TRUE){
    message('tf_sigs not loaded yet, please run db.preload() before processing !');return(FALSE);
  }
  n1 <- names(tf_sigs$tf)[-1]
  if(use_gene_type %in% n1){
    if(is.null(use_genes)==TRUE){
      TF_list <- unique(tf_sigs$tf[[use_gene_type]])
      SIG_list <- unique(tf_sigs$sig[[use_gene_type]])
    }else{
      TF_list <- unique(intersect(use_genes,tf_sigs$tf[[use_gene_type]]))
      SIG_list <- unique(intersect(use_genes,tf_sigs$sig[[use_gene_type]]))
    }
  }else{
    mart <- useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
    #filters <- listFilters(mart)
    attributes <- listAttributes(mart)
    if(!use_gene_type %in% attributes$name){
      message(sprintf('%s not in the attributes for %s, please check and re-try !',use_gene_type,dataset));return(FALSE)
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

#' Gene ID conversion related functions.
#'
#' \code{get_IDtransfer} will generate a transfer table for ID conversion by input the from-gene type and to-gene type.
#'
#' @param from_type character, attribute name from the biomaRt package,
#' such as 'ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version', 'refseq_mrna'.
#' Full list could be obtained by
#' \code{listAttributes(mart)$name}, where \code{mart} is the output of \code{useMart} function.
#' The type must be the gene type for the input \code{use_genes}.
#' @param to_type character, character, attribute name from the biomaRt package,
#' the gene type will be convereted to.
#' @param add_type character, character, attribute name from the biomaRt package,
#' additional type in the final table
#' @param use_genes a vector of characters, gene list used for ID conversion.
#' If NULL, will extract all possible genes.
#' @param dataset character, name for the dataset used for ID conversion, such as 'hsapiens_gene_ensembl'.
#' If NULL, will use \code{db_info[1]} if run \code{db.preload} brefore. Default is NULL.
#'
#' @return
#' \code{get_IDtransfer} will return a data.frame for the transfer table.
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
get_IDtransfer <- function(from_type=NULL,to_type=NULL,add_type=NULL,use_genes=NULL,dataset=NULL){
  if(is.null(dataset)==TRUE){
    if(exists('db_info')==FALSE){
      message('db_info not found, please run db.preload() first!');return(FALSE);
    }
    dataset <- db_info[1]
  }
  mart <- useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
  attributes <- listAttributes(mart)

  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,dataset));return(FALSE)
  }
  if(!to_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',to_type,dataset));return(FALSE)
  }
  if(is.null(use_genes)==TRUE | length(use_genes)>100){
    tmp1 <- getBM(attributes=c(from_type,to_type,add_type),values=1,mart=mart,filters='strand')
    tmp2 <- getBM(attributes=c(from_type,to_type,add_type),values=-1,mart=mart,filters='strand')
    tmp1 <- rbind(tmp1,tmp2)
    if(is.null(use_genes)==FALSE){
      tmp1 <- tmp1[which(tmp1[,1] %in% use_genes),]
    }
  }else{
    tmp1 <- getBM(attributes=c(from_type,to_type,add_type),values=use_genes,mart=mart,filters=from_type)
  }
  w1 <- apply(tmp1,1,function(x)length(which(is.na(x)==TRUE | x=="")))
  transfer_tab <- tmp1[which(w1==0),]
  return(transfer_tab)
}

#' Gene ID conversion related functions.
#'
#' \code{get_IDtransfer_betweenSpecies} will generate a transfer table for ID conversion between species.
#'
#' @param from_spe character, input the species name (e.g 'human', 'mouse', 'rat') that the \code{use_genes} belong to, default is 'human'
#' @param to_spe character, input the species name (e.g 'human', 'mouse', 'rat') that need to transfered to, default is 'mouse'
#' @param from_type character, attribute name from the biomaRt package,
#' such as 'ensembl_gene_id',
#' 'ensembl_gene_id_version'
#' 'ensembl_transcript_id',
#' 'ensembl_transcript_id_version',
#' 'refseq_mrna'.
#' Full list could be obtained by
#' \code{listAttributes(mart)$name}, where \code{mart} is the output of \code{useMart} function.
#' The type must be the gene type for the input \code{use_genes}.
#' @param to_type character, character, attribute name from the biomaRt package, the gene type will be convereted to.
#' @param use_genes a vector of characters, gene list used for ID conversion. Must be the genes with \code{from_type} in \code{from_spe}.
#' If NULL, will output all possible genes in transfer table. Default is NULL.
#'
#' @return
#' \code{get_IDtransfer_betweenSpecies} will return a data.frame for the transfer table.
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
                                          use_genes=NULL){
  from_spe <- toupper(from_spe)
  to_spe <- toupper(to_spe)
  ensembl <- useMart("ensembl")
  all_ds  <- listDatasets(ensembl)
  w1 <- grep(sprintf("^%s GENES",from_spe),toupper(all_ds$description))
  if(length(w1)==1){
    from_spe_ds <- all_ds[w1,1]
    mart1 <- useMart(biomart="ensembl", dataset=from_spe_ds) ## get id for input spe
  }
  if(length(w1)==0){
    message(sprintf('Check input from_spe parameter: %s, not included in the ensembl database',from_spe))
    return(FALSE)
  }
  if(length(w1)>1){
    w2 <- paste(all_ds[w1,2],collapse=';')
    message(sprintf('Check input from_spe parameter: %s, more than one species match in ensembl database : %s,
                    please check and re-try',from_spe,w2))
    return(FALSE)
  }
  w1 <- grep(sprintf("^%s GENES",to_spe),toupper(all_ds$description))
  if(length(w1)==1){
    to_spe_ds <- all_ds[w1,1]
    mart2 <- useMart(biomart="ensembl", dataset=to_spe_ds) ## get id for input spe
  }
  if(length(w1)==0){
    message(sprintf('Check input to_spe parameter: %s, not included in the ensembl database',to_spe))
    return(FALSE)
  }
  if(length(w1)>1){
    w2 <- paste(all_ds[w1,2],collapse=';')
    message(sprintf('Check input to_spe parameter: %s, more than one species match in ensembl database : %s,
                    please check and re-try',to_spe,w2))
    return(FALSE)
  }
  #### mart1 mart2
  attributes <- listAttributes(mart1)
  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,from_spe));return(FALSE)
  }
  attributes <- listAttributes(mart2)
  if(!to_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',to_type,to_spe));return(FALSE)
  }
  ## get homolog between from_spe to to_spe
  cn1 <- gsub('(.*)_gene_ensembl','\\1',from_spe_ds)
  cn2 <- attributes$name ## attribute names in mart2
  cn3 <- cn2[grep(sprintf('%s_homolog_associated_gene_name',cn1),cn2)]
  if(length(cn3)!=1){
    message('No homolog info found in Biomart, sorry !');return(FALSE)
  }
  tmp1 <- get_IDtransfer(from_type=from_type,to_type='external_gene_name',use_genes=use_genes,dataset=from_spe_ds)
  tmp2 <- getBM(attributes=c(cn3,'external_gene_name'),values=TRUE,
                mart=mart2,filters=sprintf('with_%s_homolog',cn1))
  colnames(tmp1) <- sprintf('%s_%s',colnames(tmp1),from_spe)
  colnames(tmp2) <- sprintf('%s_%s',colnames(tmp2),to_spe)
  tmp3 <- merge(tmp1,tmp2,by.x=sprintf('external_gene_name_%s',from_spe),by.y=sprintf('%s_%s',cn3,to_spe))
  transfer_tab <- tmp3[,c(2,3,1)]
  if(to_type != 'external_gene_name'){
    tmp4 <- get_IDtransfer(from_type='external_gene_name',to_type=to_type,use_genes=tmp3[,3],dataset=to_spe_ds)
    colnames(tmp4) <- sprintf('%s_%s',colnames(tmp4),to_spe)
    tmp5 <- merge(tmp3,tmp4)
    transfer_tab <- tmp5[,c(3,4,2,1)]
  }
  return(transfer_tab)
}


#' Gene ID conversion related functions.
#' \code{get_IDtransfer2symbol2type} will generate the transfer table for the original ID to the gene symbol and gene biotype (at gene level)
#' or transcript symbol and transcript biotype (at transcript level).
#'
#' @param from_type character, attribute name from the biomaRt package,
#' such as 'ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version', 'refseq_mrna'.
#' Full list could be obtained by
#' \code{listAttributes(mart)$name}, where \code{mart} is the output of \code{useMart} function.
#' The type must be the gene type for the input \code{use_genes}.
#' @param use_genes a vector of characters, gene list used for ID conversion.
#' If NULL, will output all possible genes in transfer table. Default is NULL.
#' @param dataset character, name for the dataset used for ID conversion, such as 'hsapiens_gene_ensembl'.
#' If NULL, will use \code{db_info[1]} if run \code{db.preload} brefore. Default is NULL.
#' @param use_level character, either 'transcript' or 'gene', default is 'gene'
#'
#' @return
#' \code{get_IDtransfer2symbol2type} will return a data.frame for the transfer table with gene/transcript biotype.
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
get_IDtransfer2symbol2type <- function(from_type=NULL,use_genes=NULL,dataset=NULL,use_level='gene'){
  if(is.null(dataset)==TRUE){
    dataset <- db_info[1]
  }
  message(sprintf('Your setting is at %s level',use_level))
  mart <- useMart(biomart="ensembl", dataset=dataset) ## get mart for id conversion !!!! db_info is saved in db RData
  attributes <- listAttributes(mart)
  if(!from_type %in% attributes$name){
    message(sprintf('%s not in the attributes for %s, please check and re-try !',from_type,dataset));return(FALSE)
  }
  if(use_level=='gene') tmp1 <- get_IDtransfer(from_type=from_type,to_type='external_gene_name',add_type='gene_biotype',use_genes=use_genes,dataset=dataset)
  if(use_level=='transcript')   tmp1 <- get_IDtransfer(from_type=from_type,to_type='external_gene_name',add_type='transcript_biotype',use_genes=use_genes,dataset=dataset)
  transfer_tab <- tmp1
  return(transfer_tab)
}

#' Gene ID conversion related functions.
#'
#' \code{get_name_transfertab} will get the transfered ID by input the original ID and transfer table.
#'
#' @param use_genes a vector of characters, gene list used for ID conversion.
#' @param transfer_tab data.frame, the transfer table for ID conversion, could be obtained by \code{get_IDtransfer}.
#' @param from_type character, attribute name from the biomaRt package,
#' such as 'ensembl_gene_id',
#' 'ensembl_gene_id_version'
#' 'ensembl_transcript_id',
#' 'ensembl_transcript_id_version',
#' 'refseq_mrna'.
#' Full list could be obtained by
#' \code{listAttributes(mart)$name}, where \code{mart} is the output of \code{useMart} function.
#' The type must be the gene type for the input \code{use_genes}.
#'  If NULL, will use the first column in the transfer table.
#' @param to_type character, character, attribute name from the biomaRt package,
#'  the gene type will be convereted to.
#'  If NULL, will use the second column in the transfer table.
#'
#' @return
#' \code{get_name_transfertab} will return the list for the converted IDs.
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
get_name_transfertab <- function(use_genes=NULL,transfer_tab=NULL,from_type=NULL,to_type=NULL){
  if(is.null(from_type)==TRUE){from_type=colnames(transfer_tab)[1];}
  if(is.null(to_type)==TRUE){to_type=colnames(transfer_tab)[2];}
  transfer_tab <- unique(transfer_tab[,c(from_type,to_type)])
  x <- use_genes;
  t1 <- unique(transfer_tab[which(transfer_tab[,from_type] %in% x),])
  c1 <- unique(t1[,from_type])
  if(length(c1)<nrow(t1)){
    message('Gene ID in from type contain multiple items!');return(FALSE)
  }
  rownames(t1) <- t1[,from_type]
  x1 <- t1[x,to_type]
  w1 <- which(is.na(x1)==TRUE)
  x1[w1] <- x[w1]
  x1
}

#' Create directory for the network generation part. Suggested but not required.
#'
#' \code{NetBID.network.dir.create} will generate a working directory structure for the network generation part in NetBID2.
#' This function aims to assist researchers to organize the working directory. It is suggested but not required.
#'
#' This function need to input the main directory for the project with the project name.
#' It will generate three sub-directories, QC/ for the QC-related plots, DATA/ for saving the RData, and SJAR/ for running SJAracne.
#' Such organization is suggested but not required to use all functions in NetBID2.
#' This function will return a variable called \code{network.par}.
#' This variable is strongly suggested to use in the network generation part of NetBID2.
#' It will be used to store all related datasets during calculation.
#'
#' @param project_main_dir character, main directory for the project.
#' @param prject_name character, project name.
#'
#' @return a list called network.par, including main.dir,project.name, out.dir, out.dir.QC, out.dir.DATA, out.dir.SJAR

#' @examples
#'
#' \dontrun{
#' NetBID.network.dir.create(project_main_dir='demo1/',project_name='network_test')
#' }
#' @export
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

#' Create directory for the driver analysis part. Suggested but not required.
#'
#' \code{NetBID.analysis.dir.create} will generate a working directory structure for the driver analysis part in NetBID2.
#' This function aims to assist researchers to organize the working directory. It is suggested but not required.
#'
#' This function need to input the main directory for the analysis project with the project name,
#' and need to input the project directory for the network generation with the network_project_name represents the name of the network for use (SJAR.project_name in SJ.SJAracne.prepare)
#' It will generate three sub-directories, QC/ for the QC-related plots, DATA/ for saving the RData, and PLOT/ for saving visulization plots.
#' Such organization is suggested but not required to use all functions in NetBID2.
#' This function will return a variable called \code{analysis.par}.
#' This variable is strongly suggested to use in the driver analysis part of NetBID2.
#' It will be used to store all related datasets during calculation.
#'
#' @param project_main_dir character, main directory for the project in the driver analysis part.
#' @param prject_name character, project name.
#' @param network_dir character, main directory for the project in the network generation part.
#' @param network_project_name character, the project name of the network for use (SJAR.project_name in SJ.SJAracne.prepare or SJAracne.prepare);
#' This parameter is optional. If previously has not follow the NetBID2 suggested pipeline, could leave this to NULL
#' but set the real path to tf.network.file and sig.network.file if want to follow the driver analysis part of NetBID2 suggested pipeline.
#' @param tf.network.file character, file path of the TF network (XXX/consensus_network_ncol_.txt). Optional, if do not set network_project_name.
#' @param sig.network.file character, file path of the SIG network (XXX/consensus_network_ncol_.txt). Optional, if do not set network_project_name.
#'
#' @return a list called analysis.par, including main.dir,project.name, out.dir, out.dir.QC, out.dir.DATA, out.dir.PLOT
#'
#' @examples
#'
#' \dontrun{
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' #
#' project_main_dir <- 'demo1/'
#' project_name <- 'driver_test'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             prject_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' }
#' @export
NetBID.analysis.dir.create <- function(project_main_dir=NULL,prject_name=NULL,
                                       network_dir=NULL,
                                       network_project_name=NULL,
                                       tf.network.file=NULL,
                                       sig.network.file=NULL){
  if(exists('analysis.par')==TRUE){message('analysis.par is occupied in the current session,please manually run: rm(analysis.par) and re-try, otherwise will not change !');
    return(analysis.par)}
  if(is.null(project_main_dir)==TRUE){message('project_main_dir required, please input and re-try!');return(FALSE)}
  if(is.null(prject_name)==TRUE){message('prject_name required, please input and re-try!');return(FALSE)}
  if(is.null(network_dir)==TRUE){message('network_dir required, please input and re-try!');return(FALSE)}
  #if(is.null(network_project_name)==TRUE){message('network_project_name required, please input and re-try!');return(FALSE)}
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
  analysis.par$tf.network.file <- NULL
  analysis.par$sig.network.file <- NULL
  if(is.null(tf.network.file)==FALSE){
    analysis.par$tf.network.file  <- tf.network.file
  }else{
    tf_net1 <- sprintf('%s/SJAR/%s/output_tf_sjaracne_%s_out_.final/consensus_network_ncol_.txt',
                      network_dir,network_project_name,network_project_name) ## old version of sjaracne
    if(file.exists(tf_net1)) analysis.par$tf.network.file <- tf_net1
    tf_net2 <- sprintf('%s/SJAR/SJARACNE_%s_TF/SJARACNE_out.final/consensus_network_ncol_.txt',
                      network_dir,network_project_name) ## new version of sjaracne
    if(file.exists(tf_net2)) analysis.par$tf.network.file <- tf_net2
  }
  if(is.null(tf.network.file)==FALSE){
    analysis.par$sig.network.file <- sig.network.file
  }else{
    sig_net1 <- sprintf('%s/SJAR/%s/output_sig_sjaracne_%s_out_.final/consensus_network_ncol_.txt',
                       network_dir,network_project_name,network_project_name) ## old version of sjaracne
    if(file.exists(sig_net1)) analysis.par$sig.network.file <- sig_net1
    sig_net2 <- sprintf('%s/SJAR/SJARACNE_%s_SIG/SJARACNE_out.final/consensus_network_ncol_.txt',
                       network_dir,network_project_name) ## new version of sjaracne
    if(file.exists(sig_net2)) analysis.par$sig.network.file <- sig_net2
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

#' Automatically save RData for NetBID2. Suggested in the NetBID2 pipeline analysis but not required.
#'
#' \code{NetBID.saveRData} is a function strongly suggested to control the pipeline analysis in NetBID2.
#'
#' Users could save two complicate list object, network.par in the network generation part and analysis.par in the driver analysis part,
#' into the data directory (network.par$out.dir.DATA or analysis.par$out.dir.DATA), with the name of the RData marked by \code{step} name.
#' The two lists could make user to save the whole related dataset in each step \code{NetBID.saveRData} and easy to get them back by using \code{NetBID.loadRData}.
#' The RData saved from each step could be used to run the following analysis without repeating the former steps.
#'
#' @param network.par list, store all related datasets during network generation part.
#' @param analysis.par list, store all related datasets during driver analysis part.
#' @param step character, name for the step.
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

#' Automatically load RData for NetBID2. Suggested in the NetBID2 pipeline analysis but not required.
#'
#' \code{NetBID.loadRData} is a function strongly suggested to control the pipeline analysis in NetBID2.
#'
#' Users could save two complicate list object, network.par in the network generation part and analysis.par in the driver analysis part,
#' into the data directory (network.par$out.dir.DATA or analysis.par$out.dir.DATA), with the name of the RData marked by \code{step} name.
#' The two lists could make user to save the whole related dataset in each step \code{NetBID.saveRData} and easy to get them back by using \code{NetBID.loadRData}.
#' The RData saved from each step could be used to run the following analysis without repeating the former steps.
#'
#' @param network.par list, store all related datasets during network generation part.
#' @param analysis.par list, store all related datasets during driver analysis part.
#' @param step character, name for the step.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#'
#' @export
NetBID.loadRData <- function(network.par=NULL,analysis.par=NULL,step='exp-load'){
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

#' Load gene expression set from GEO database
#'
#' \code{load.exp.GEO} is a function to get GEO dataset, save into ExpressionSet class object and save to RData.
#' If the RData exists, user could choose to directly load from it.
#'
#' This function aims to provide simple way to get expression set from GEO database and manage the related RData.
#'
#' @param out.dir character, the file path to save the RData or if the RData has already generated, load from it.
#' Default is network.par$out.dir.DATA if follow the NetBID2 suggested pipeline.
#' @param GSE character, the GSE ID to download.
#' @param GPL character, the GPL ID to download.
#' @param getGPL logical, whether to download the GPL file or not. Default is TRUE.
#' @param update logical, whether to update the RData if the data is already in the out.dir.Default is FALSE
#'
#' @return ExpressionSet class object if success.
#' @examples
#'
#' \dontrun{
#' net_eset <- load.exp.GEO(out.dir='test/',
#'                          GSE='GSE116028',
#'                          GPL='GPL6480',
#'                          getGPL=TRUE,
#'                          update=FALSE)
#' }
#' @export
load.exp.GEO <- function(out.dir=network.par$out.dir.DATA,GSE = NULL,GPL = NULL,getGPL=TRUE,update = FALSE){
  if(is.null(out.dir)==TRUE){
    message('out.dir required, please re-try')
    return(FALSE)
  }
  if (is.null(GSE) | is.null(GPL)){
    message('GSE and GPL required, please re-try')
    return(FALSE)
  }
  expRData_dir <- sprintf('%s/%s_%s.RData', out.dir, GSE,GPL)
  if (file.exists(expRData_dir) & update == FALSE) {
    message(sprintf('RData exist in %s and update==TRUE, will directly load from RData .',expRData_dir))
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

#' Load gene expression set from Salmon output (demo version).
#'
#' \code{load.exp.RNASeq.demoSalmon} is a function to read in salmon results and convert to eSet/DESeqDataSet class object.
#'
#' This function assist to read in RNASeq results from salmon.
#' Due to the complicated condition (e.g reference sequence) in running salmon, this function is just a demo version that may not suit for all conditions.
#'
#' @param salmon_dir character, the main file directory to save the salmon results.
#' @param tx2gene data.frame or NULL, this parameter will be passed to \code{tximport}, check for detail.
#' If NULL, will read the transcript name in one of the salmon output
#' (this only works when using e.g gencode.v29.transcripts.fa from GENCODE as reference).
#' @param use_phenotype_info data.frame, phenotype information dataframe, must contain the columns \code{use_sample_col} and \code{use_design_col}.
#' @param use_sample_col character, the column name to indicate which column in \code{use_phenotype_info} should be used as the sample name.
#' @param use_design_col character, the column name to indicate which column in \code{use_phenotype_info} should be used as the design feature of the samples.
#' @param return_type character, the class of the return object, choose from 'eset','dds'. 'eset' is the ExpressionSet class object,
#' 'dds' is the DESeqDataSet class object.
#' Default is 'eset'
#' @param merge_level character, choose from 'gene' and 'transcript',
#' if choose 'gene' and original salmon results is mapped to the transcriptome,
#' expression matrix will be merged to gene level.
#' (this only works when using e.g gencode.v29.transcripts.fa from GENCODE as reference).
#' @export
load.exp.RNASeq.demoSalmon <- function(salmon_dir = "",tx2gene=NULL,
                                       use_phenotype_info = NULL,
                                       use_sample_col=NULL,
                                       use_design_col=NULL,
                                       return_type='eset',
                                       merge_level='gene') {
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
  eset <- load.exp.RNASeq.demo(files,type='salmon',
                               tx2gene=tx2gene,
                               use_phenotype_info=use_phenotype_info,
                               use_sample_col=use_sample_col,
                               use_design_col=use_design_col,
                               return_type=return_type,
                               merge_level=merge_level)
  return(eset)
}

#' Load gene expression set from RNASeq results (demo version).
#'
#' \code{load.exp.RNASeq.demo} is a function to read in RNASeq results and convert to eSet/DESeqDataSet class object.
#'
#' This function assist to read in RNASeq results from different resources.
#' Due to the complicated condition (e.g reference sequence) in running RNASeq,
#' this function is just a demo version that may not suit for all conditions.
#'
#' @param files a vector of characters,filenames for the transcript-level abundances, will be passed to \code{tximport}.
#' Check for detail.
#' @param type character, he type of software used to generate the abundances, will be passed to \code{tximport}.
#' Check for detail.
#' @param tx2gene data.frame or NULL, this parameter will be passed to \code{tximport}, check for detail.
#' @param use_phenotype_info data.frame, phenotype information dataframe, must contain the columns \code{use_sample_col} and \code{use_design_col}.
#' @param use_sample_col character, the column name to indicate which column in \code{use_phenotype_info} should be used as the sample name.
#' @param use_design_col character, the column name to indicate which column in \code{use_phenotype_info} should be used as the design feature of the samples.
#' @param return_type character, the class of the return object, choose from 'eset','dds'. 'eset' is the ExpressionSet class object,
#' 'dds' is the DESeqDataSet class object.
#' Default is 'eset'
#' @param merge_level character, choose from 'gene' and 'transcript',
#' if choose 'gene' and original salmon results is mapped to the transcriptome,
#' expression matrix will be merged to gene level.
#' (this only works when using e.g gencode.v29.transcripts.fa from GENCODE as reference).
#' @export
load.exp.RNASeq.demo <- function(files,type='salmon',
                                 tx2gene=NULL,
                                 use_phenotype_info = NULL,
                                 use_sample_col=NULL,
                                 use_design_col=NULL,
                                 return_type='eset',
                                 merge_level='gene') {
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

#' Generate expression Set object
#'
#' \code{generate.eset} is a function to generate eSet object by input expression matrix (required),
#' phenotype and feature information (optional).
#'
#' This function aims to assist the generation of eSet object,
#' especially for sometimes only expression matrix is available.
#'
#' @param exp_mat data matrix, the expression data matrix with each row a gene/transcript and each column a sample.
#' @param phenotype_info data.frame, the phenotype information for the samples in \code{exp_mat},
#' the rownames must match the colnames of \code{exp_mat}. If NULL, will generate a single column dataframe.
#' Default is NULL.
#' @param feature_info data.frame,the feature information for the genes/transcripts/probes in \code{exp_mat},
#' the rownames must match the rownames of \code{exp_mat}. If NULL, will generate a single column dataframe.
#' Default is NULL.
#' @param annotation_info character, the annotation for the eSet. Default is "".
#'
#' @return an ExressionSet object.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000),nrow=1000,ncol=10)
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' eset <- generate.eset(exp_mat=mat1)
#' @export
generate.eset <- function(exp_mat=NULL, phenotype_info=NULL, feature_info=NULL, annotation_info="") {
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
  if(class(phenotype_info)=='character' | is.null(dim(phenotype_info))==TRUE){
    phenotype_info <- data.frame(group = phenotype_info, stringsAsFactors = FALSE)
    rownames(phenotype_info) <- colnames(exp_mat)
  }
  if(class(feature_info)=='character' | is.null(dim(feature_info))==TRUE){
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

#' Merge two ExpressionSet class object.
#'
#' \code{merge_eset} is a function to merge two ExpressionSet class objects.
#' If the genes in two objects are same, the expression matrix will be directly merged, otherwise,
#' Z-transformation is performed before merging.
#'
#' @param eset1 ExpressionSet class, the first dataset to merge.
#' @param eset2 ExpressionSet class, the second dataset to merge.
#' @param group1 character, name for the first group.
#' @param group2 character, name for the second group.
#' @param use_col a vector of characters, the column names in the phenotype information to be kept in the merged ExpressionSet.
#' If NULL, will use the intersected column names between the two datasets. Default is NULL.
#' @param group_col_name, character, name for the column indicate the originate of the samples in the phenotype information dataframe in the merged ExpressionSet.
#' Default is 'original_group'.
#' @param remove_batch logical, indicate whether or not to remove batch effect between two sample set.
#' Default is FALSE.
#'
#' @return an ExressionSet object.
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
                       remove_batch = FALSE) {
  mat1 <- exprs(eset1)
  mat2 <- exprs(eset2)
  w1 <- intersect(rownames(mat1), rownames(mat2))
  if(length(w1)==0){
    message('No overlap genes between two eSet, please check and re-try!');return(FALSE);
  }
  if(length(w1)<nrow(mat1) & length(w1)<nrow(mat2)){
    ## z-transformation
    mat1 <- std(mat1)
    mat2 <- std(mat2)
  }
  rmat <- cbind(mat1[w1, ], mat2[w1, ])
  #choose1 <- apply(rmat <= quantile(rmat, probs = 0.05), 1, sum) <= ncol(rmat) * 0.90 ## low expressed genes
  #rmat <- rmat[choose1, ]
  phe1 <- pData(eset1)
  phe2 <- pData(eset2)
  if(length(use_col)==0){
    use_col <- intersect(colnames(phe1),colnames(phe2))
  }
  if(length(use_col)>1)
    rphe <- rbind(phe1[colnames(mat1), use_col], phe2[colnames(mat2), use_col])
  if(length(use_col)==1){
    rphe <- c(phe1[colnames(mat1), use_col], phe2[colnames(mat2), use_col])
    rphe <- data.frame(rphe,stringsAsFactors=FALSE); colnames(rphe) <- use_col;
    rownames(rphe) <- colnames(rmat)
  }
  if(is.null(group1)==TRUE) group1 <- 'group1'
  if(is.null(group2)==TRUE) group2 <- 'group2'
  rphe[[group_col_name]]<- c(rep(group1, ncol(mat1)), rep(group2, ncol(mat2)))
  if (remove_batch == TRUE) {
    rmat <- removeBatchEffect(rmat,batch=rphe[[group_col_name]])
  }
  reset <- generate.eset(rmat,phenotype_info = rphe, annotation_info = 'combine')
  return(reset)
}

#' Update ExpressionSet feature information, mainly for gene ID conversion
#'
#' \code{update_eset.feature} is a function to update the feature information in the ExpressionSet object.
#'
#' This function is designed mainly for gene ID conversion.
#' User could input the transfer table for ID conversion, which could be obtained from the original feature table
#' (if use \code{load.exp.GEO} and set getGPL==TRUE) or by running the function \code{get_IDtransfer}.
#' The relationship betweeen the original ID and the target ID could be as follows:
#' 1) one->one, just change the ID label.
#' 2) multiple->one, the expression value for the target ID will be the merge of the original ID, with user-defined choice of merge_method.
#' 3) one->mulitple, the expression value for the original ID will be distributed to the matched target IDs, with user-defined choice of distribute_method.
#' warning message will be presented in this condition.
#' 4) multiple-->multiple, the function will do distribute step 3) first and merge 2).
#'
#' @param use_eset ExpressionSet class object, the original ExpressionSet to update.
#' @param use_feature_info data.frame, the transfer table for ID conversion.
#' @param from_feature character,the column name in \code{use_feature_info},must be the same type of the rownames for the expression matrix in \code{use_eset}.
#' @param to_feature character, the column in \code{use_feature_info}, the target ID will be convereted to.
#' @param merge_method character, startegy to merge the expression value, choose from 'median','mean','max','min'. Default is "median".
#' @param distribute_method character, strategy to distribute the expression value, choose from 'mean', 'equal'. Default is "equal".
#'
#' @return an ExressionSet object.
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
#' print(exprs(eset)[test_transfer_table$Gene,])
#' print(exprs(new_eset))
#'
#' @export update_eset.feature
update_eset.feature <- function(use_eset=NULL,use_feature_info=NULL,from_feature=NULL,to_feature=NULL,
                                merge_method='median',distribute_method='equal'){
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
  w1 <- which(use_feature_info[,1]!="" & use_feature_info[,2]!="" & is.na(use_feature_info[,1])==FALSE & is.na(use_feature_info[,2])==FALSE)
  use_feature_info <- use_feature_info[w1,]
  g1 <- rownames(mat) ## rownames for the expmat
  f1 <- as.character(use_feature_info[,from_feature]) ## from feature info
  t1 <- as.character(use_feature_info[,to_feature]) ## to feature info
  w1 <- which(f1 %in% g1); f1 <- f1[w1]; t1 <- t1[w1]; ## only consider features in the rownames of expmat
  if(length(w1)==0){
    message(sprintf('Rownames of the expression matrix was not included in the %s column, please check and re-try !',from_feature))
    return(use_eset)
  }
  message(sprintf('%d transfer pairs related with %d rows from original expression matrix will be keeped !',length(w1),length(g1)))
  fc1 <- table(f1); tc1 <- table(t1); fw1 <- which(fc1>1); tw1 <- which(tc1>1); ## check duplicate records
  if(length(fw1)>0){
    message(sprintf('Original feature %s has %d items with duplicate records, will distribute the original values equal to all related items !
                    if do not want this, please check and retry !',from_feature,length(fw1)))
    #return(use_eset)
    w2 <- which(f1 %in% names(fw1)) ## need to distribute
    w0 <- setdiff(1:length(f1),w2) ## do not need to distribute
    if(distribute_method=='equal'){
      v1 <- mat[f1[w2],]; ## distribute equal
    }
    if(distribute_method=='mean'){
      v1 <- mat[f1[w2],]; ## distribute mean
      tt <- as.numeric(table(f1[w2])[f1[w2]])
      v1 <- v1/tt;
    }
    rownames(v1) <- paste0(f1[w2],'-',t1[w2]);
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

#' Update ExpressionSet phenotype information, mainly for changing sample name and extract useful phenotype information for sample clustering.
#'
#' \code{update_eset.phenotype} is a function to update the phenotype information in the ExpressionSet object.
#'
#' This function is designed mainly for extracting useful phenotype information for samples.
#' Especially for phenotype information directly get from GEO database.
#'
#' @param use_eset ExpressionSet class object, the original ExpressionSet to update.
#' @param use_phenotype_info data.frame, the phenotype information dataframe.
#' @param use_sample_col character,the column name in \code{use_phenotype_info}, which contains the sample name.
#' @param use_col character, the columns in \code{use_phenotype_info} to be kept.
#' If set to 'auto', wil extract columns with unique sample feature ranges from 2 to sample size-1.
#' If set to 'GEO-auto', will extract columns: 'geo_accession','title','source_name_ch1',and columns end with ':ch1'.
#' Default is "auto".
#' @return an ExressionSet object.
#' @examples
#' \dontrun{
#' net_eset <- load.exp.GEO(out.dir='test/',
#'                          GSE='GSE116028',
#'                          GPL='GPL6480',
#'                          getGPL=TRUE,
#'                          update=FALSE)
#' net_eset <- update_eset.phenotype(use_eset=net_eset,
#'                                   use_phenotype_info=pData(net_eset),
#'                                   use_sample_col='geo_accession',
#'                                   use_col='GEO-auto')
#' }
#' @export update_eset.phenotype
update_eset.phenotype <- function(use_eset=NULL,use_phenotype_info=NULL,use_sample_col=NULL,use_col='auto'){
  if(is.null(use_eset)){
    message('use_eset required, please re-try !');
    return(use_eset)
  }
  if(is.null(use_phenotype_info)) use_phenotype_info <- pData(use_eset)
  if(is.null(use_sample_col)==FALSE){
    if(!use_sample_col %in% colnames(use_phenotype_info)){
      message(sprintf('%s not in the colnames of use_phenotype_info, please re-try!',use_sample_col));return(use_eset)
    }
  }
  if(is.null(use_col)) use_col <- colnames(use_phenotype_info)
  mat <- exprs(use_eset)
  s1 <- colnames(mat) ## all samples
  if(is.null(use_sample_col)==TRUE){
    p1 <- rownames(use_phenotype_info)
  }else{
    p1 <- use_phenotype_info[,use_sample_col] ## sample in the phenotype info
  }
  w1 <- which(p1 %in% s1); p1 <- p1[w1]; ## only consider samples in the colnames of expmat
  if(length(w1)==0){
    if(is.null(use_sample_col)==TRUE){
      message('Colnames of the expression matrix was not included in the rownames of use_phenotype_info, please check and re-try !')
    }else{
      message(sprintf('Colnames of the expression matrix was not included in the %s column, please check and re-try !',use_sample_col))
    }
    return(use_eset)
  }
  message(sprintf('%d out of %d samples from the expression matrix will be keeped !',length(w1),length(s1)))
  mat_new <- mat[,p1]
  use_phenotype_info <- use_phenotype_info[w1,]
  n1 <- colnames(use_phenotype_info)
  if(use_col[1] == 'GEO-auto'){
    w1 <- c('geo_accession','title','source_name_ch1',n1[grep(':ch1',n1)])
    p1 <- use_phenotype_info[,w1]
    colnames(p1)[4:ncol(p1)] <- gsub('(.*):ch1','\\1',colnames(p1)[4:ncol(p1)])
    colnames(p1)[3] <- gsub('(.*)_ch1','\\1',colnames(p1)[3])
    rownames(p1) <- use_phenotype_info[,use_sample_col]
    p1 <- as.data.frame(apply(p1,2,function(x){if(class(x)=='factor'){as.character(x)}else{x}}),stringsAsFactors=FALSE)
    new_phenotype_info <- p1
  }else{
    if(use_col[1] == 'auto'){
      u1 <- apply(use_phenotype_info,2,function(x)length(unique(x)))
      w1 <- which(u1>=2 & u1<=nrow(use_phenotype_info)-1)
      if(length(w1)==0){
        message('No column could match the auto criteria, please check and re-try!');return(FALSE)
      }
      new_phenotype_info <- use_phenotype_info[,w1]
    }else{
      if(length(setdiff(use_col,n1))>0){
        message(sprintf('%s not in use_phenotype_info, please re-try!',paste(setdiff(use_col,n1),collapse=';')));return(FALSE)
      }
      p1 <- use_phenotype_info[,use_col]
      p1 <- as.data.frame(apply(p1,2,function(x){if(class(x)=='factor'){as.character(x)}else{x}}),stringsAsFactors=FALSE)
      new_phenotype_info <- p1
    }
  }
  #print(new_phenotype_info)
  new_eset <- generate.eset(exp_mat=mat_new, phenotype_info=new_phenotype_info, feature_info=fData(use_eset), annotation_info=annotation(use_eset))
  return(new_eset)
}

#' IQR (interquartile range) filter for the genes in the expression matrix
#'
#' \code{IQR.filter} is a function to extract genes by their IQR value.
#'
#' This function aims to extract out most variable genes (defined by the IQR value).
#' This step will be used to perform sample cluster and to prepare the input for SJAracne.
#'
#' @param exp_mat matrix, the gene expression matrix, with each row a gene/transcript/probe and each column a sample.
#' @param use_genes a vector of characters, the gene list to report. Default is the rownames of \code{exp_mat}.
#' @param thre numeric, the quantile threshold of IQR. Default is 0.5.
#' @param loose_gene a vector of characters, the gene list that only need to pass the \code{loose_thre}.
#' This parameter is designed for inputing possible drivers used in SJAracne. Default is NULL.
#' @param loose_thre numeric,the quantile threshold of IQR for the genes in \code{loose_gene}. Default is 0.1.
#' @return a vector with logical values indicate which genes should be kept.
#' @examples
#' mat1 <- matrix(rnorm(15000),nrow=1500,ncol=10)
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' choose1 <- IQR.filter(mat1,thre=0.5,
#'                      loose_gene=paste0('Gene',1:100))
#' @export
IQR.filter <- function(exp_mat,use_genes=rownames(exp_mat),thre = 0.5,loose_gene=NULL,loose_thre=0.1) {
  use_genes <- intersect(use_genes,rownames(exp_mat))
  use_genes <- setdiff(use_genes,"")
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

#' Simple function to normalize RNASeq Read Count data.
#'
#' \code{RNASeqCount.normalize.scale} is a simple version function to normalize the RNASeq read count data.
#'
#' Strongly suggest to follow the DESeq2 pipeline to process the RNASeq datasets.
#' \code{load.exp.RNASeq.demo} and \code{load.exp.RNASeq.demoSalmon} is included in NetBID2 but may not suit for all conditions.
#'
#' @param mat matrix, the original input matrix for the read data, each row is a gene/transcript with each column a sample.
#' @param total integer, total read counts, if NULL will use the mean of the column sum. Default is NULL.
#' @param pseudoCount integer, pseudo count to add for all read counts to avoid -Inf in following log transformation.
#' Default is 1.
#'
#' @return a matrix containing the normalized RNASeq count.
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
  d <- mat
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

#' Calculate activity value for all possible drivers
#'
#' \code{cal.Activity} is a function to calculate activity for all possible drivers by
#' input the target list for drivers and the expression matrix for target genes.
#'
#' @param target_list a list for the target gene information for the drivers.
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list} could be automatically generated
#' by \code{get_net2target_list} by running
#' \code{get.SJAracne.network}.
#' @param cal_mat numeric matrix,the expression matrix for all genes/transcripts for calculation.
#' @param es.method character, strategy to calculate the activity for the driver,
#' choose from mean, weighted mean, maxmean, absmean;
#' in which the weighted in the weighted mean is the MI (mutual information) value * sign of correlation (use the spearman correlation sign).
#' Default is 'mean'.
#'
#' @return the activity matrix with each row a driver and each column a sample (same sample order with cal_mat)
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ac_mat <- cal.Activity(target_list=analysis.par$tf.network$target_list,
#'                        cal_mat=exprs(analysis.par$cal.eset),
#'                        es.method='weightedmean')
#' @export
cal.Activity <- function(target_list=NULL, cal_mat=NULL, es.method = 'mean') {
  ## mean, absmean, maxmean, weightedmean
  use_genes <- row.names(cal_mat)
  all_target <- target_list
  #all_target <- all_target[intersect(use_genes, names(all_target))] ## if the driver is not included in cal_mat but its target genes are included, will also calculate activity
  ac.mat <-
    matrix(NA, ncol = ncol(cal_mat), nrow = length(all_target)) ## generate activity matrix, each col for sample, each row for source target
  #z-normalize each sample
  cal_mat <- apply(cal_mat, 2, std)
  for (i in 1:length(all_target)) {
    x <- names(all_target)[i]
    x1 <- all_target[[x]]
    x2 <- unique(intersect(rownames(x1), use_genes)) ## filter target by cal genes
    x1 <- x1[x2, ] ## target info
    target_num <- length(x2)
    if (target_num == 0)
      next
    if (target_num == 1){
      if (es.method == 'weightedmean') ac.mat[i, ] <- cal_mat[x2,]
      if (es.method != 'weightedmean') ac.mat[i, ] <- cal_mat[x2,]*x1$MI * sign(x1$spearman)
      next
    }
    if (es.method != 'weightedmean')
      ac.mat[i, ] <- apply(cal_mat[x2,], 2, es, es.method)
    if (es.method == 'weightedmean') {
      weight <- x1$MI * sign(x1$spearman)
      ac.mat[i, ] <- apply(cal_mat[x2,] * weight, 2, es, 'mean')
    }
  }
  rownames(ac.mat) <- names(all_target)
  colnames(ac.mat) <- colnames(cal_mat)
  return(ac.mat)
}

#' Calculate activity value for gene sets.
#'
#' \code{cal.Activity.GS} is a function to calculate activity for all gene sets.
#' @param use_gs2gene a list for geneset to genes, the name for the list is the gene set name and the content in each list is the vector for genes belong to that gene set.
#' Default is all_gs2gene[c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG')] with all_gs2gene loaded by using \code{gs.preload}.
#' @param cal_mat numeric matrix,the expression matrix for all genes/transcripts for calculation.
#' @param es.method character, strategy to calculate the activity for the driver,
#' choose from mean, absmean, maxmean;
#' Default is 'mean'.
#' @return the activity matrix with each row a gene set and each column a sample (same sample order with cal_mat)
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,
#'                         use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
#' exp_mat_gene <- exprs(analysis.par$cal.eset)
#' ## each row is a gene symbol, if not, must convert ID first
#' ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,
#'                         cal_mat = exp_mat_gene)
#' @export
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

#' Get differential expression (DE)/differential activity (DA) between case and control sample groups by Bayesian Inference.
#'
#' \code{getDE.BID.2G} is a function aims to get DE/DA genes/drivers with detailed statistical information between case (G1) and control (G0) groups.
#'
#' @param eset ExpressionSet class object, the input gene expression or driver activity.
#' @param G1 a vector of characters, the sample list used as the case.
#' @param G0 a vecotr of characters, the sample list used as the control.
#' @param G1_name character, the group name for the samples in G1, default is "G1".
#' @param G0_name character, the group name for the samples in G0, default is "G0".
#' @param verbose logical, whether or not to print verbose information during calculation.
#' Default is TRUE.
#'
#' @return
#' A dataframe for all genes with the columns as the output of \code{topTable} in limma.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- pData(analysis.par$cal.eset)
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
#' @export
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
  de<- ddply(d,'id','combRowEvid.2grps',comp=comp,family=gaussian,method='Bayesian',n.iter=5000,nitt=25000,burnin=5000,thin=1,pooling=c('full'),logTransformed=TRUE,restand=FALSE,average.method=c('geometric'))
  names(de)<-gsub('.full',"",names(de))
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
  tT <- tT[rownames(new_mat),]
  exp_G1 <- rowMeans(new_mat[,G1])
  exp_G0 <- rowMeans(new_mat[,G0])
  tT <- cbind(tT,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  tT <- tT[order(tT$P.Value),]
  return(tT)
}

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

#' Get differential expression (DE)/differential activity (DA) between case and control sample groups by limma related functions.
#'
#' \code{getDE.limma.2G} is a function aims to get DE/DA genes/drivers with detailed statistical information between case (G1) and control (G0) groups.
#'
#' @param eset ExpressionSet class object, the input gene expression or driver activity.
#' @param G1 a vector of characters, the sample list used as the case.
#' @param G0 a vecotr of characters, the sample list used as the control.
#' @param G1_name character, the group name for the samples in G1, default is "G1".
#' @param G0_name character, the group name for the samples in G0, default is "G0".
#' @param verbose logical, whether or not to print verbose information during calculation.Default is TRUE.
#' @param random_effect a vector of characters, vector or factor specifying a blocking variable.
#' Default is NULL indicating no random effect will be considered.
#'
#' @return
#' A dataframe for all genes with the columns as the output of \code{topTable} in limma.
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- pData(analysis.par$cal.eset)
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
  tT <- topTable(fit2, adjust.method = "fdr", number = Inf,coef=1)
  tT <- tT[order(tT$P.Value, decreasing = FALSE), ]
  tT <- cbind(ID=rownames(tT),tT,stringsAsFactors=FALSE)
  tT <- tT[rownames(new_mat),]
  exp_G1 <- rowMeans(new_mat[,G1])
  exp_G0 <- rowMeans(new_mat[,G0])
  w1 <- which(tT$P.Value<=0);
  if(length(w1)>0) tT$P.Value[w1] <- .Machine$double.xmin;
  z_val <- sapply(tT$P.Value*sign(tT$logFC),function(x)combinePvalVector(x,twosided = TRUE)[1])
  if(is.null(random_effect)==TRUE){
    tT <- cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  }else{
    tT <- cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1,
                'Ave.G0_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G0'],
                'Ave.G1_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G1'])
  }
  if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  return(tT)
}

#' Combine Pvalues by Fisher or Stouffer's method
#'
#' \code{combinePvalVector} is a function to combine P-values by Stouffer of Fisher method.
#'
#' @param pvals, a vector of numeric values, the P-value values need to combine.
#' The sign of the P-value will be added to indicate the direction of testing if signed is set to TRUE.
#' @param method character, choose from Stouffer or Fisher, default is "Stouffer".
#' @param signed logical, whether the sign of the pvals will be considered in calculation.
#' Default is TRUE.
#' @param twosided logical, whether the pvalues are from two-sided test or not.
#' If not, pvalues must between 0 and 0.5.
#' Default is TRUE.
#' @return a vector contains the 'Z-statistics' and 'P.Value'
#' @examples
#' combinePvalVector(c(0.1,1e-3,1e-5))
#' combinePvalVector(c(0.1,1e-3,-1e-5))
#' @export
combinePvalVector <-
  function(pvals,
           method = c('Stouffer', 'Fisher'),
           signed = TRUE,
           twosided = TRUE) {
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

#' Merge activity value ExpressionSet for TF (transcription factors) and Sig (signaling factors)
#'
#' \code{merge_TF_SIG.AC} is a function to merge the activity value ExpressionSet for TF (transcription factors) and Sig (signaling factors).
#'
#' This function aims to merge the driver activity ExpressionSet, by automatically add "_TF"/"_SIG" suffix for drivers.
#'
#' @param TF_AC ExpressionSet object containing the activity value for all TFs.
#' @param SIG_AC ExpressionSet object containing the activity value for all SIGs.
#'
#' @return an ExpressionSet object
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ac_mat_TF <- cal.Activity(target_list=analysis.par$tf.network$target_list,
#'                        cal_mat=exprs(analysis.par$cal.eset),
#'                        es.method='weightedmean')
#' ac_mat_SIG <- cal.Activity(target_list=analysis.par$tf.network$target_list,
#'                        cal_mat=exprs(analysis.par$cal.eset),
#'                        es.method='weightedmean')
#' analysis.par$ac.tf.eset  <- generate.eset(exp_mat=ac_mat_TF,
#'                                           phenotype_info=pData(analysis.par$cal.eset))
#' analysis.par$ac.sig.eset <- generate.eset(exp_mat=ac_mat_SIG,
#'                                           phenotype_info=pData(analysis.par$cal.eset))
#' analysis.par$merge.ac.eset <- merge_TF_SIG.AC(TF_AC=analysis.par$ac.tf.eset,
#'                                           SIG_AC=analysis.par$ac.sig.eset)
#' @export
merge_TF_SIG.AC <- function(TF_AC=NULL,SIG_AC=NULL){
  mat_TF <- exprs(TF_AC)
  mat_SIG <- exprs(SIG_AC)
  funcType <- c(rep('TF',nrow(mat_TF)),rep('SIG',nrow(mat_SIG)))
  rn <- c(rownames(mat_TF),rownames(mat_SIG))
  rn_label <- paste(rn,funcType,sep='_')
  mat_combine <- rbind(mat_TF,mat_SIG[,colnames(mat_TF)])
  rownames(mat_combine) <- rn_label
  eset_combine <- generate.eset(exp_mat=mat_combine,phenotype_info=pData(TF_AC)[colnames(mat_combine),],
                                feature_info=NULL,annotation_info='activity in dataset')
  return(eset_combine)
}

#' Merge target network for TF (transcription factors) and Sig (signaling factors).
#'
#' \code{merge_TF_SIG.network} is a function to merge the target network from TF and Sig.
#'
#' @param TF_network TF network obtained by \code{get.SJAracne.network}.
#' @param SIG_network SIG network obtained by \code{get.SJAracne.network}.
#' @return
#' This function will return the same structure object as TF_network/SIG_network,
#' which is a list containing three items, \code{network_dat},
#' \code{target_list} and \code{igraph_obj}.
#' @examples
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             prject_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
#' analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,
#'                                                    SIG_network=analysis.par$sig.network)
#' @export
merge_TF_SIG.network <- function(TF_network=NULL,SIG_network=NULL){
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

#' Generate final master table for drivers.
#'
#' \code{generate.masterTable.TF_SIG} is a function to automatically generate master tables
#' by input TF (transcription factor)/Sig (signaling factor) results separately.
#'
#' This function is designed for automatically generate master table, with :
#' DE (differentiated expression), DA_TF (differentiated activity for TF), DA_SIG (differentiated activity for SIG),
#' TF_network (a list for the target gene information for the TF), SIG_network (a list for the target gene information for the SIG),
#' and necessary additional information (main_id_type, tf_sigs, z_col and display_col).
#' If results from TF/SIG have merged before, please use \code{generate.masterTable}.
#'
#' @param use_comp a vector of characters, the comparison name used to display in the master table.
#' @param DE a list of DE results, the list name must contain the items in \code{use_comp}.
#' @param DA_TF a list of TF DA results, the list name must contain the items in \code{use_comp}.
#' @param DA_SIG a list of SIG DA results, the list name must contain the items in \code{use_comp}.
#' @param TF_network a list for the target gene information for the TFs.
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{TF_network} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param SIG_network a list for the target gene information for the SIGs.
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{TF_network} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param main_id_type character, the main gene id type. The attribute name from the biomaRt package,
#' such as 'ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version', 'refseq_mrna'. Full list could be obtained by
#' \code{listAttributes(mart)$name}, where \code{mart} is the output of \code{useMart} function.
#' @param tf_sigs list, which contain the detailed information for the TF and SIGs, this can be obtained by running \code{db.preload}.
#' @param z_col character, column name in \code{DE}, \code{DA_TF}, \code{DA_SIG} that contains the Z statistics. Default is 'Z-statistics'.
#' @param display_col character,column name in \code{DE}, \code{DA_TF}, \code{DA_SIG} to be kept in the master table. Default is c('logFC','P.Value').
#' @return a data.frame that contains the information for all tested drivers
#' The column "originalID" and "originalID_label" contain the ID same with the original dataset.
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' #analysis.par$final_ms_tab ## this is master table generated before
#' ac_mat <- cal.Activity(target_list=analysis.par$tf.network$target_list,
#'                        cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
#' analysis.par$ac.tf.eset  <- generate.eset(exp_mat=ac_mat,
#'                                           phenotype_info=pData(analysis.par$cal.eset))
#' ac_mat <- cal.Activity(target_list=analysis.par$sig.network$target_list,
#'                        cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
#' analysis.par$ac.sig.eset  <- generate.eset(exp_mat=ac_mat,
#'                                           phenotype_info=pData(analysis.par$cal.eset))
#' phe_info <- pData(analysis.par$cal.eset)
#' all_subgroup <- unique(phe_info$subgroup) ##
#' for(each_subtype in all_subgroup){
#'  comp_name <- sprintf('%s.Vs.others',each_subtype) ## each comparison must give a name !!!
#'  G0  <- rownames(phe_info)[which(phe_info$`subgroup`!=each_subtype)] # get sample list for G0
#'  G1  <- rownames(phe_info)[which(phe_info$`subgroup`==each_subtype)] # get sample list for G1
#'  DE_gene_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,
#'                                  G1_name=each_subtype,G0_name='other')
#'  analysis.par$DE[[comp_name]] <- DE_gene_limma
#'  DA_driver_limma <- getDE.limma.2G(eset=analysis.par$ac.tf.eset,G1=G1,G0=G0,
#'                                    G1_name=each_subtype,G0_name='other')
#'  analysis.par$DA_TF[[comp_name]] <- DA_driver_limma
#'  DA_driver_limma <- getDE.limma.2G(eset=analysis.par$ac.sig.eset,G1=G1,G0=G0,
#'                                    G1_name=each_subtype,G0_name='other')
#'  analysis.par$DA_SIG[[comp_name]] <- DA_driver_limma
#' }
#' all_comp <- names(analysis.par$DE) ## get all comparison name for output
#' db.preload(use_level='gene',use_spe='human',update=FALSE);
#' test_ms_tab <- generate.masterTable.TF_SIG(use_comp=all_comp,
#'                                            DE=analysis.par$DE,
#'                                            DA_TF=analysis.par$DA_TF,
#'                                            DA_SIG=analysis.par$DA_SIG,
#'                                            TF_network=analysis.par$tf.network$target_list,
#'                                            SIG_network=analysis.par$sig.network$target_list,
#'                                            tf_sigs=tf_sigs,
#'                                            z_col='Z-statistics',
#'                                            display_col=c('logFC','P.Value'),
#'                                            main_id_type='external_gene_name')
#' @export
generate.masterTable.TF_SIG <- function(use_comp=NULL,DE=NULL,DA_TF=NULL,DA_SIG=NULL,
                                        TF_network=NULL,SIG_network=NULL,
                                        main_id_type=NULL,
                                        tf_sigs=tf_sigs,
                                        z_col='Z-statistics',display_col=c('logFC','P.Value')){
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
    x1 <- x[which(x!="")]
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

#' Generate final master table for drivers.
#'
#' \code{generate.masterTable} is a function to automatically generate master tables
#' by input the merged TF/SIG results.
#'
#' This function is designed for automatically generate master table, with :
#' DE (differentiated expression), DA(differentiated activity for drivers),network (a list for the target gene information for the drivers),
#' and necessary additional information (main_id_type, tf_sigs, z_col and display_col).
#' If results from TF/SIG do not have merged before, please use \code{generate.masterTable.TF_SIG}.
#'
#' @param use_comp a vector of characters, the comparison name used to display in the master table.
#' @param DE a list of DE results, the list name must contain the items in \code{use_comp}.
#' @param DA a list of DA results, the list name must contain the items in \code{use_comp}.
#' @param network a list for the target gene information for the drivers.
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{TF_network} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param main_id_type character, the main gene id type. The attribute name from the biomaRt package,
#' such as 'ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version', 'refseq_mrna'. Full list could be obtained by
#' \code{listAttributes(mart)$name}, where \code{mart} is the output of \code{useMart} function.
#' @param tf_sigs list, which contain the detailed information for the TF and SIGs, this can be obtained by running \code{db.preload}.
#' @param z_col character, column name in \code{DE}, \code{DA} that contains the Z statistics. Default is 'Z-statistics'.
#' @param display_col character,column name in \code{DE}, \code{DA} to be kept in the master table. Default is c('logFC','P.Value').
#' @return a data.frame that contains the information for all tested drivers.
#' The column "originalID" and "originalID_label" contain the ID same with the original dataset.
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' #analysis.par$final_ms_tab ## this is master table generated before
#' ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,
#'                        cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
#' analysis.par$ac.merge.eset  <- generate.eset(exp_mat=ac_mat,
#'                                              phenotype_info=pData(analysis.par$cal.eset))
#' phe_info <- pData(analysis.par$cal.eset)
#' all_subgroup <- unique(phe_info$subgroup) ##
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
#'                                            network=analysis.par$merge.network$target_list,
#'                                            tf_sigs=tf_sigs,
#'                                            z_col='Z-statistics',
#'                                            display_col=c('logFC','P.Value'),
#'                                            main_id_type='external_gene_name')
#' @export
generate.masterTable <- function(use_comp=NULL,DE=NULL,DA=NULL,
                                 network=NULL,main_id_type=NULL,
                                 tf_sigs=tf_sigs,
                                 z_col='Z-statistics',display_col=c('logFC','P.Value')){
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
    x1 <- x[which(x!="")]
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

#' Output master table to excel files.
#'
#' \code{out2excel} is a function to output dataframe into an excel file, mainly for output the master table generated by \code{generate.masterTable}.
#'
#' @param all_ms_tab list or dataframe, if dataframe, it is the master table generated by \code{generate.masterTable}.
#' If it is a list, each item in list could be a dataframe containing the master table.
#' The name of the list will be the sheet name in the excel file.
#' @param out.xlsx character, file name for the output excel.
#' @param mark_gene list, list of marker genes, additional info to add in the master table.The name of the list is the name for the mark group.
#' @param mark_col character, the color used to label the marker genes. If NULL, will use \code{get.class.color} to get the colors.
#' @param mark_strategy character, choose from 'color','add_column'. 'Color' means the mark_gene will be displayed by its background color;
#' 'add_column' means the mark_gene will be displayed in separate columns with content TRUE/FALSE indicating whether the genes belong to each mark group.
#' @param workbook_name character, workbook name for the output excel. Default is 'ms_tab'.
#' @param only_z_sheet logical, if TRUE will generate a separate sheet only contain Z related columns in DA/DE. Default is FALSE.
#' @param z_column character, the column name that contain the Z-statistics. If NULL, will find columns start with "Z.".
#' Default is NULL.
#' @param sig_thre numeric, threshold for the Z-statistics, values passed the threshold will be marked by the color automatically generated by \code{z2col}.
#' Default is 1.64.
#' @return logical value indicating whether the file has been sucessfully generated or not.
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
#' mark_col <- get.class.color(names(mark_gene))
#' outfile <- 'test_out.xlsx'
#' out2excel(ms_tab,out.xlsx = out_file,mark_gene,mark_col)
#' }
#' @export
out2excel <- function(all_ms_tab,out.xlsx,
                      mark_gene=NULL,
                      mark_col=NULL,
                      mark_strategy='color',
                      workbook_name='ms_tab',
                      only_z_sheet=FALSE,
                      z_column=NULL,sig_thre=1.64){
  wb <- createWorkbook(workbook_name)
  if(!mark_strategy %in% c('color','add_column')){
    message('mark_strategy must be color or add_column, please check and re-try!');return(FALSE);
  }
  if(!'list' %in% class(all_ms_tab)){
    all_ms_tab <- list('Sheet1'=as.data.frame(all_ms_tab))
  }
  if(only_z_sheet==TRUE){
    nn <- names(all_ms_tab)
    all_ms_tab <- lapply(all_ms_tab,function(x){
      w1 <- grep('.*_[DE|DA]',colnames(x))
      w2 <- setdiff(colnames(x)[w1],colnames(x)[w1][grep('^Z.*',colnames(x)[w1])])
      w3 <- setdiff(colnames(x),w2)
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
        r1 <- do.call(cbind,lapply(mark_gene,function(x1){
          ifelse(g1 %in% x1,'TRUE','FALSE')
        }))
        new_x <- cbind(x,r1)
        colnames(new_x) <- c(colnames(x),new_col_name)
        new_x
      })
    }
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
  saveWorkbook(wb, out.xlsx, overwrite = TRUE)
  return(TRUE)
  ##
}

#' Load MSigDB for NetBID2 into R workspace.
#'
#' \code{gs.preload} will load two variables into R workspace, the list for gene set to genes (all_gs2gene)
#'  and a dataframe (all_gs2gene_info) for detailed description for gene sets.
#'
#' This is a pre-processing function for NetBID2 advanced analysis, user only need to input the species name (e.g Homo sapiens, Mus musculus).
#' The function could automatically download information from MSigDB by the functions in \code{msigdbr} and save into RData under the db/ directory
#' with specified species name.
#'
#' @param use_spe character, input the species name (e.g 'Homo sapiens', 'Mus musculus'). Full list of available species name could be found by \code{msigdbr_show_species()}.
#' Default is 'Homo sapiens'
#' @param update logical,whether to update if previous RData has been generated, default FALSE
#' @param main.dir character, main file path for NetBID2, if NULL, will set to system.file(package = "NetBID2"). Default is NULL.
#' @param db.dir character, file path for saving the RData, default is \code{db} directory under \code{main.dir} when setting for \code{main.dir}.
#'
#' @return Reture TRUE if success and FALSE if not. Will load two variables into R workspace, all_gs2gene and all_gs2gene_info
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
  ## only support geneSymbol (because pipeline-generated master table will contain geneSymbol column)
  if(is.null(main.dir)==TRUE){
    main.dir <- system.file(package = "NetBID2")
    message(sprintf('main.dir not set, will use package directory: %s',main.dir))
  }
  if(is.null(db.dir)==TRUE){
    db.dir <- sprintf("%s/db/",main.dir)
  }
  message(sprintf('Will use directory %s as the db.dir',db.dir))
  all_spe <- msigdbr_show_species()
  if(!use_spe %in% all_spe){
    message(sprintf('%s not in %s, please check and re-try !',use_spe,paste(all_spe,collapse=';')))
    return(FALSE)
  }
  use_spe1 <- gsub(' ','_',use_spe)
  out_file <- sprintf('%s/%s_gs2gene.RData',db.dir,use_spe1)
  if(file.exists(out_file)==FALSE | update==TRUE){
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
    all_gs_subcat <- setdiff(unique(all_gs_info$gs_subcat),"")
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
    all_gs2gene_info[,c(1,2,4,5)] <- as.data.frame(apply(all_gs2gene_info[,c(1,2,4,5)],2,function(x){x[which(is.na(x)==TRUE)] <- "";x}),stringsAsFactors=FALSE)
    save(all_gs2gene,all_gs2gene_info,file=out_file)
  }
  load(out_file,.GlobalEnv)
  message('all_gs2gene loaded, you could see all_gs2gene_info to check the details !')
  return(TRUE)
}

########################### visualization functions
## simple functions to get info
#' Generate a vector for sample category.
#'
#' \code{get_obs_label} will generate a vector for sample categories with names to the vector representing the sample name.
#'
#' This is a simple function to generate sample category vector from a data.frame.
#' Mainly used for input preparation in the visualization plots.
#'
#' @param phe_info data.frame, phenotype dataframe for the samples with sample names in rownames, e.g from \code{pData(eset)}.
#' @param use_col a vector of numeric or character, the column index or column name for extraction to get the sample category vector.
#' @param collapse character, character string to separate the results when the length use_col is more than 1. default is "|".
#'
#' @return
#' Will return a vector for sample categories with names to the vector representing the sample name.

#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' phe_info <- pData(analysis.par$cal.eset)
#' use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
#' print(use_obs_class)
#' \dontrun{
#'}
#' @export
get_obs_label <- function(phe_info,use_col,collapse='|'){
  obs_label<-phe_info[,use_col];
  if(length(use_col)>1){
    obs_label<-apply(obs_label,1,function(x)paste(x,collapse=collapse))
  }
  names(obs_label) <- rownames(phe_info);
  obs_label
}

#' Get interested groups from the pData of an ExpressionSet object.
#'
#' \code{get_int_group} is a simple function to extract columns with unique sample feature ranges from 2 to sample size-1.
#'
#' @param eset, an ExpressionSet object, the input object for analysis.
#' @return a vector of characters, the column names which could be used for sample cluster analysis
#' @examples
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' intgroups <- get_int_group(network.par$net.eset)
#' @export
get_int_group <- function(eset){
  phe <- pData(eset)
  feature_len <- apply(phe,2,function(x)length(unique(x)))
  intgroup <- colnames(phe)[which(feature_len>1 & feature_len<nrow(phe))]
  return(intgroup)
}

#' Get Score between predicted label and observed label.
#'
#' \code{get_clustComp} is a function to compare the predicted label and observed label and return the score.
#'
#' @param pred_label a vector of characters, the predicted label.
#' @param obs_label a vector of characters, the observed label
#' @param strategy character, the strategy to compare with labels,
#' choose from 'ARI (adjusted rand index)', 'NMI (normalized mutual information)', 'Jaccard'. Default is 'ARI'.
#' @return score for the comparison
#' @examples
#' obs_label <- c('A','A','A','B','B','C','D')
#' pred_label  <- c(1,1,1,1,2,2,2)
#' get_clustComp(pred_label,obs_label)
#' @export
get_clustComp <- function(pred_label, obs_label,strategy='ARI') {
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
    names(obs_label) <- as.character(1:length(obs_label))
    names(pred_label) <- names(obs_label)
  }
  if(!strategy %in% c('ARI','NMI','Jaccard')){
    message('Only support RI,ARI,NMI and Jaccard!');return(FALSE);
  }
  if(strategy=='Jaccard') res1 <- get_jac(pred_label, obs_label) else res1 <- clustComp(pred_label, obs_label)[[strategy]]
  return(res1)
}

# get jaccard accuracy
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

#' Draw the cluster comparison between predicted label and observed label.
#'
#' \code{draw.clustComp} is a function to draw the comparison between the predicted label and observed label.
#'
#' @param pred_label a vector of characters, the predicted label.
#' @param obs_label a vector of characters, the observed label
#' @param strategy character, the strategy to compare with labels,
#' choose from 'ARI (adjusted rand index)', 'NMI (normalized mutual information)', 'Jaccard'. Default is 'ARI'.
#' @param use_col logical, whether or not to use color in the plot. Default is TRUE.
#' @param low_K integer, the lowest number to display on the figures. Number smaller than this will directly display the sample name.
#' Default is 5
#' @param highlight_clust a vector of characters, the cluster need to be highlighted on the plot.
#' @param main character, the title for the plot.
#' @param clust_cex numeric, cex for the cluster label. Default is 1
#' @param outlier_cex numeric, cex for the sample names. Default is 0.3
#' @return a matrix for the number in the plot
#' @examples
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' mat <- exprs(network.par$net.eset)
#' phe <- pData(network.par$net.eset)
#' intgroup <- 'subgroup'
#' pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,
#'                              obs_label=get_obs_label(phe,intgroup),
#'                              kmeans_strategy='consensus')
#' draw.clustComp(pred_label,get_obs_label(phe,intgroup),outlier_cex=1,low_K=2,use_col=TRUE)
#' draw.clustComp(pred_label,get_obs_label(phe,intgroup),outlier_cex=1,low_K=2,use_col=FALSE)
#' @export
draw.clustComp <- function(pred_label, obs_label,strategy='ARI',
                           use_col=TRUE,low_K=5,
                           highlight_clust=NULL,
                           main=NULL,clust_cex=1,outlier_cex=0.3) {
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
    names(obs_label) <- as.character(1:length(obs_label))
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
  t1 <- table(list(pred_label[nn],obs_label[nn]))
  layout(1)
  par(mar=c(2,12,3,4))
  if(use_col==TRUE){
    image(t1,col=c('white',colorRampPalette(brewer.pal(8,'Reds'))(length(unique(as.numeric(t1)))-1)),bty='n',xaxt='n',yaxt='n',
          main=mm)
  }else{
    image(t1,bty='n',xaxt='n',yaxt='n',
          main=mm,col='white')
  }

  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])
  xx <- seq(pp[1],pp[2],length.out = nrow(t1)+1)
  yy <- seq(pp[3],pp[4],length.out = ncol(t1)+1)
  xxx <- (xx[1:(length(xx)-1)]+xx[2:length(xx)])/2
  yyy <- (yy[1:(length(yy)-1)]+yy[2:length(yy)])/2

  abline(h=yy);abline(v=xx)
  text(pp[1],yyy,colnames(t1),adj=1,xpd=TRUE,col=ifelse(colnames(t1) %in% highlight_clust,2,1),cex=clust_cex)
  text(xxx,pp[4]+max(strheight(rownames(t1),units='inches',cex=1))/12,rownames(t1),adj=0.5,srt=0,xpd=TRUE,
       col=ifelse(rownames(t1) %in% highlight_clust,2,1),cex=clust_cex)

  for(i in 1:nrow(t1)){
    for(j in 1:ncol(t1)){
      v1 <- t1[i,j]
      if(v1==0) next
      if(v1>low_K){
        text(xxx[i],yyy[j],v1,cex=clust_cex)
      }else{
        v2 <- names(obs_label)[which(pred_label==rownames(t1)[i] & obs_label==colnames(t1)[j])]
        v2 <- paste(v2,collapse='\n')
        text(xxx[i],yyy[j],v2,cex=outlier_cex)
      }
    }
  }
  return(t1)
}


#' Get the color for the input Z statistics.
#'
#' \code{z2col} is a function to transfer the input Z statistics to a color bar.
#'
#' @param x a vector of numeric values. The input Z statistics.
#' @param n_len integer, number of unique colors. Default is 60.
#' @param sig_thre numeric, the threshold for significance (absolute Z statistics), values do not pass the threshold will be colored in 'white'.
#' @param col_min_thre numeric, the threshold for the lowest values used to generate the color bar. Default is 0.01.
#' @param col_max_thre numeric, the threshold for the maximum values used to generate the color bar. Default is 3.
#' @param blue_col a vector of characters, the blue colors used for the negative values in Z statistics. Default is brewer.pal(9,'Set1')[2].
#' @param red_col a vector of characters, the red lors used for the negative values in Z statistics. Default is brewer.pal(9,'Set1')[1].
#' @return a vector of color characters
#' @examples
#' t1 <- sort(rnorm(mean=0,sd=2,n=100))
#' image(as.matrix(t1),col=z2col(t1))
#' @export
z2col <- function(x,n_len=60,sig_thre=0.01,col_min_thre=0.01,col_max_thre=3,
                  blue_col=brewer.pal(9,'Set1')[2],
                  red_col=brewer.pal(9,'Set1')[1]){
  ## create vector for z-score, can change sig threshold
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

#' Generate a color vector for input character.
#'
#' \code{get.class.color} will generate a vector of colors for input character.
#'
#' This is a simple function to generate a vector of colors for input characters.
#' Users could define some of the colors for part of the inputs.
#'
#' @param x a vector of characters.
#' @param use_color a vector of characters, color bar used to generate the color vector for the input.Default is brewer.pal(12, 'Set3').
#' @param pre_define a vector of characters, pre-defined color code for some of the input characters. Default is NULL.
#' @param use_inner logical, indicating whether to use inner pre-defined color code for some characters. Default is TRUE.
#'
#' @return
#' Will return a vector of colors with names the input vector characters.

#' @examples
#' get.class.color(c('ClassA','ClassB','ClassC','ClassA','ClassC','ClassC'))
#' get.class.color(c('ClassA','ClassB','ClassC','SHH','WNT','Group3','Group4'),
#'                 use_inner=FALSE)
#' get.class.color(c('ClassA','ClassB','ClassC','SHH','WNT','Group3','Group4'),
#'                 use_inner=FALSE,use_color=brewer.pal(8, 'Set1'))
#'
#' pre_define <- c('blue', 'red', 'yellow', 'green','yellow', 'green')
#'                 ## pre-defined colors for MB
#' names(pre_define) <- c('WNT', 'SHH', 'Group3', 'Group4','GroupC', 'GroupD')
#'                 ##pre-defined color name for MB
#' get.class.color(c('ClassA','ClassB','ClassC','SHH','WNT','Group3','Group4'),
#'                 pre_define=pre_define,use_inner=FALSE)
#'
#' \dontrun{
#'}
#' @export
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

#' Draw 2D dimension plot for sample cluster visualization.
#'
#' \code{draw.2D} is a function to draw 2D dimension plot for sample cluster visualization.
#'
#' @param X a vector of numeric values, the first dimension values.
#' @param Y a vector of numeric values, the second dimension values.
#' @param class_label a vector of characters, with the names are the samples (optional) and the values are the sample cluster label.
#' The function will treat that the order are the same for X,Y and class_label.
#' @param xlab character, the label for x-axis.
#' @param ylab character, the label for y-axis.
#' @param legend_cex numeric, the cex for the legend.
#' @param main character, the title for the plot.
#' @param point_cex numeric, the cex for the points.
#'
#' @return logical value indicating whether the plot has been successfully generated
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.2D(X=pc[,1],Y=pc[,2],class_label=pred_label)
#' @export
draw.2D <- function(X,Y,class_label,xlab='PC1',ylab='PC2',legend_cex=0.8,main="",point_cex=1){
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
         horiz = FALSE,xpd = TRUE,border = NA,bty = 'n',cex=legend_cex)
  return(TRUE)
}

#' Draw 2D dimension plot for sample cluster visualization with user-defined text on each point.
#'
#' \code{draw.2D.text} is a function to draw 2D dimension plot for sample cluster visualization.
#'
#' @param X a vector of numeric values, the first dimension values.
#' @param Y a vector of numeric values, the second dimension values.
#' @param class_label a vector of characters, with the names are the samples (optional) and the values are the sample cluster label.
#' The function will treat that the order are the same for X,Y and class_label.
#' @param class_text a vector of characters, the user-defined text on each point.
#' If NULL, will use the names of class_label. Default is NULL.
#' @param xlab character, the label for x-axis.
#' @param ylab character, the label for y-axis.
#' @param legend_cex numeric, the cex for the legend.
#' @param main character, the title for the plot.
#' @param point_cex numeric, the cex for the points.
#' @param text_cex numeric, the cex for the points.
#'
#' @return logical value indicating whether the plot has been successfully generated
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.2D.text(X=pc[,1],Y=pc[,2],class_label=pred_label,
#'              point_cex=5,text_cex=0.5)
#' @export
draw.2D.text <- function(X,Y,class_label,class_text=NULL,xlab='PC1',ylab='PC2',legend_cex=0.8,main="",
                         point_cex=1,text_cex=NULL){
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
  if(is.null(class_text)==TRUE){
    class_text <- names(class_label)
  }
  cc <- 10/length(class_label)
  if(cc<0.05) cc<-0.05
  if(cc>1) cc<-1
  if(is.null(text_cex)==FALSE) cc <- text_cex
  class_label <- as.character(class_label)
  cls_cc <- get.class.color(class_label) ## get color for each label
  plot(Y ~ X,pch = 16,cex = point_cex,col = cls_cc,main=main,xlab=xlab,ylab=ylab)
  text(x=X,y=Y,labels=class_text,cex=cc,xpd=TRUE,adj=0.5)
  #print(cc);print(str(nn))
  legend(par()$usr[2],par()$usr[4],sort(unique(class_label)),fill = cls_cc[sort(unique(class_label))],
         horiz = FALSE,xpd = TRUE,border = NA,bty = 'n',cex=legend_cex)
  return(TRUE)
}

#' Draw 3D dimension plot for sample cluster visualization.
#'
#' \code{draw.3D} is a function to draw 3D dimension plot for sample cluster visualization.
#'
#' @param X a vector of numeric values, the first dimension values.
#' @param Y a vector of numeric values, the second dimension values.
#' @param Z a vector of numeric values, the third dimension values.
#' @param class_label a vector of characters, with the names are the samples (optional) and the values are the sample cluster label.
#' The function will treat that the order are the same for X,Y and class_label.
#' @param xlab character, the label for x-axis.
#' @param ylab character, the label for y-axis.
#' @param zlab character, the label for z-axis.
#' @param legend_cex numeric, the cex for the legend.
#' @param main character, the title for the plot.
#' @param point_cex numeric, the cex for the points.
#' @param legend_pos character, the position to put the legend. Default is 'topright'.
#' @param legend_ncol integer, number of columns used to display the legend. Default is 1.
#' @param ... other paramters used in \code{scatter3D}.
#'
#' @return logical value indicating whether the plot has been successfully generated
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.3D(X=pc[,1],Y=pc[,2],Z=pc[,3],class_label=pred_label)
#' @export
draw.3D <- function(X,Y,Z,class_label,xlab='PC1',ylab='PC2',zlab='PC3',
                    legend_cex=0.8,main="",point_cex=1,legend_pos='topright',legend_ncol=1,...){
  if(is.null(class_label)==TRUE){
    message('No class_label, please check and re-try !');return(FALSE);
  }
  if(length(X)!=length(Y)){
    message('Input two dimension vector with different length, please check and re-try !');return(FALSE);
  }
  if(length(X)!=length(class_label)){
    message('Input dimension vector has different length with the class_label, please check and re-try !');return(FALSE);
  }
  class_label <- as.character(class_label)
  c1 <- factor(class_label,levels=sort(unique(class_label)))
  cls_cc <- get.class.color(class_label) ## get color for each label
  #print(str(class_label))
  #print(str(sort(unique(class_label))))
  #print(str(cls_cc))
  #print(cls_cc[sort(unique(class_label))])
  scatter3D(X,Y,Z,pch = 16,xlab = xlab,ylab = ylab,zlab=zlab,bty='g',
            colvar=1:length(c1),
            col=cls_cc,colkey = FALSE,cex=point_cex,main=main,...)
  legend(legend_pos,sort(unique(class_label)),fill = cls_cc[sort(unique(class_label))],
         border = NA,bty = 'n',ncol = legend_ncol,cex = legend_cex)
  return(TRUE)
}

#' Draw 2D dimension plot with ellipse for sample cluster visualization.
#'
#' \code{draw.2D.ellipse} is a function to draw 2D dimension plot with ellipse to cover the samples in the sample cluster for visualization.
#'
#' @param X a vector of numeric values, the first dimension values.
#' @param Y a vector of numeric values, the second dimension values.
#' @param class_label a vector of characters, with the names are the samples (optional) and the values are the sample cluster label.
#' The function will treat that the order are the same for X,Y and class_label.
#' @param xlab character, the label for x-axis.
#' @param ylab character, the label for y-axis.
#' @param legend_cex numeric, the cex for the legend.
#' @param main character, the title for the plot.
#' @param point_cex numeric, the cex for the points.
#'
#' @return logical value indicating whether the plot has been successfully generated
#' @examples
#' mat1 <- matrix(rnorm(2000,mean=0,sd=1),nrow=100,ncol=20)
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' pc <- prcomp(t(mat1))$x
#' pred_label <- kmeans(pc,centers=4)$cluster ## this can use other cluster results
#' draw.2D.ellipse(X=pc[,1],Y=pc[,2],class_label=pred_label)
#' @export
draw.2D.ellipse <- function(X,Y,class_label,xlab='PC1',ylab='PC2',legend_cex=0.8,main="",point_cex=1){
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
      text(m1[1],m1[2],i,xpd=TRUE,adj=0,cex=legend_cex)
    }
    if(length(d1)==2){
      draw.ellipse(m1[1],m1[2],a=d1[1],b=d1[1],col=c1,border=NA,xpd=TRUE)
      text(m1[1],m1[2],i,xpd=TRUE,adj=0,cex=legend_cex)
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
        boxtext(m2[1]+d1,m2[2],labels=i,adj=0,cex=legend_cex,col.bg=c1)
        segments(x0=w1[s1,1],y0=w1[s1,2],x1=m2[1]+d1,y1=m2[2],col='dark grey')
      } else{
        boxtext(m2[1]-d1,m2[2],labels=i,adj=1,cex=legend_cex,col.bg=c1)
        segments(x0=w1[s1,1],y0=w1[s1,2],x1=m2[1]-d1,y1=m2[2],col='dark grey')
      }
    }
  }
  return(TRUE)
}

###
prepdata <- function (expressionset, intgroup, do.logtransform)
{
  conversions = c(RGList = "NChannelSet")
  for (i in seq_along(conversions)) {
    if (is(expressionset, names(conversions)[i])) {
      expressionset = try(as(expressionset, conversions[i]))
      if (is(expressionset, "try-error")) {
        stop(sprintf("The argument 'expressionset' is of class '%s', and its automatic conversion into '%s' failed. Please try to convert it manually, or contact the creator of that object.\n",
                     names(conversions)[i], conversions[i]))
      }
      else {
        break
      }
    }
  }
  x = platformspecific(expressionset, intgroup, do.logtransform)
  if (!all(intgroup %in% colnames(x$pData)))
    stop("all elements of 'intgroup' should match column names of 'pData(expressionset)'.")
  x = append(x, list(numArrays = ncol(x$M), intgroup = intgroup,
                     do.logtransform = do.logtransform))
  x = append(x, intgroupColors(x))
  return(x)
}

#' QC plot for ExpressionSet class object.
#'
#' \code{draw.eset.QC} is a function to draw the basic QC plots for an ExpressionSet class object.
#' The QC plots include heatmap, pca, density and meansd.
#'
#' @param eset ExpressionSet class, the input ExpressionSet class object to be plot.
#' @param outdir character, output directory to save the figures.
#' Suggest to set \code{network.par$out.dir.QC} or \code{analysis.par$out.dir.QC}
#' @param do.logtransform logical, whether to do log transformation before drawing the QC plots. Default is FALSE.
#' @param intgroup a vector of characters, the interested groups from the phenotype information of the eset to be used in plot.
#' If NULL, will automatcially extract all possible groups by \code{get_int_group}.
#' Default is NULL.
#' @param prefix character, the prefix for the QC figure name.Default is "".
#' @param choose_plot a vector of characters,
#' choose one or multiple from 'heatmap', 'pca', 'density', 'meansd.'
#' Default is 'heatmap','pca','density','meansd'.
#' @examples
#' \dontrun{
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' intgroups <- get_int_group(network.par$net.eset)
#' network.par$out.dir.QC <- getwd() ## set the output directory
#' draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=intgroups)
#' }
#' @export
draw.eset.QC <- function(eset,outdir = '.',do.logtransform = FALSE,intgroup=NULL,prefix = '',
                         choose_plot=c('heatmap','pca','density','meansd')) {
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
              legend_cex=0.8,main=sprintf('PCA/Kmeans plot for %s',intgroup[i]))
      draw.2D.ellipse(as.data.frame(pca$x)$PC1,as.data.frame(pca$x)$PC2,class_label=class_label,xlab='PC1',ylab='PC2',
                      legend_cex=0.8,main=sprintf('PCA/Kmeans plot for %s',intgroup[i]))
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
           type = 'l',xlab = "",ylab='Density',main = sprintf('Density plot for %s',intgroup[i]),
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

#' Draw the cluster plot by PCA (visulization algorithm) and Kmeans (cluster algorithm).
#'
#' \code{draw.pca.kmeans} is a visualization function to draw the cluster for the input data matrix.
#'
#' This function mainly aims for the visulization of sample clusters.
#' The input is the expression matrix for each row a gene/transcript/probe and each column a sample.
#' User need to input the real observation label for samples and this function will choose the best K by compared with predicted label and the observed label.
#' The output figure contains two sub-figures, left is labelled by the real observartion label and right is labelled by the predicted label
#' with comparison score (choose from ARI, NMI, Jaccard) shown above.
#'
#' @param mat a numeric data matrix, the column (e.g sample) will be clustered by the features (e.g genes) in rows.
#' @param all_k a vector of integers, the pre-defined K to evaluate.
#' If NULL, will use all possible K. Default is NULL.
#' @param obs_label a vector of characters, a vector for sample categories with names representing the sample name.
#' @param legend_pos character, position for the legend displayed on the plot. Default is 'topleft'.
#' @param legend_cex numeric, cex for the legend displayed on the plot. Default is 0.8.
#' @param plot_type character, the type for the plot, choose from '2D' or '2D.ellipse' or '3D'. Default is '2D.ellipse'.
#' @param point_cex numeric, cex for the point in the plot. Default is 1.
#' @param kmeans_strategy character, the strategy to run the kmeans algorith, choose from 'basic' and 'consensus',
#' here the consensus kmeans is performed by functions in \code{ConsensusClusterPlus}. Default is 'basic'.
#' @param choose_k_strategy character, the strategy to choose the best K,
#' choose from 'ARI (adjusted rand index)', 'NMI (normalized mutual information)', 'Jaccard'. Default is 'ARI'.
#' @param return_type character, the strategy to return the results, choose from 'optimal' and 'all'.
#' If choose 'all', cluster results from all k in \code{all_k} will be returned.
#' Default is 'optimal'.
#' @return a vector of predicted label (if \code{return_type} is 'optimal') and a list of all possible K (if \code{return_type} is 'all')
#' @examples
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' mat <- exprs(network.par$net.eset)
#' phe <- pData(network.par$net.eset)
#' intgroup <- get_int_group(network.par$net.eset)
#' for(i in 1:length(intgroup)){
#'  print(intgroup[i])
#'  pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]))
#'  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
#' }
#' pred_label <- draw.pca.kmeans(mat=mat,all_k = NULL,
#'                              obs_label=get_obs_label(phe,intgroup[i]),
#'                              kmeans_strategy='consensus')
#' @export
draw.pca.kmeans <- function(mat=NULL,all_k=NULL,obs_label=NULL,legend_pos = 'topleft',legend_cex = 0.8,
                               plot_type='2D.ellipse',point_cex=1,
                               kmeans_strategy='basic',choose_k_strategy='ARI',
                               return_type='optimal'){
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
  cluster_mat <- prcomp(t(mat))$x
  all_jac <- list()
  all_k_res <- list()
  if(kmeans_strategy=='basic'){
    for(k in all_k){
      tmp_k <- list()
      for(i in 1:10){
        tmp_k[[i]] <- kmeans(cluster_mat,centers=as.numeric(k))$cluster
      }
      pred_label <- tmp_k
      jac <- unlist(lapply(pred_label,function(x){get_clustComp(x, obs_label,strategy=choose_k_strategy)}))
      top_i <- which.max(jac)
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
  message('Optimal k is chosen byScore between predicted and observed label')
  print(all_jac)
  use_k <- all_k[which.max(all_jac)]
  pred_label <- all_k_res[[as.character(use_k)]]
  message(sprintf('Best Score occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
##
  d1 <- data.frame(id=colnames(mat),X=cluster_mat[,1],Y=cluster_mat[,2],Z=cluster_mat[,3],label=pred_label,stringsAsFactors=FALSE)
  layout(t(matrix(1:2)))
  if(plot_type=='2D.ellipse'){
    draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='PC1',ylab='PC2',legend_cex=legend_cex,point_cex=point_cex)
    draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,xlab='PC1',ylab='PC2',legend_cex=legend_cex,point_cex=point_cex,
                    main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  if(plot_type=='2D'){
    draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='PC1',ylab='PC2',legend_cex=legend_cex,point_cex=point_cex)
    draw.2D(d1$X,d1$Y,class_label=d1$label,xlab='PC1',ylab='PC2',legend_cex=legend_cex,point_cex=point_cex,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  if(plot_type=='3D'){
    draw.3D(d1$X,d1$Y,d1$Z,class_label=obs_label[d1$id],xlab='PC1',ylab='PC2',zlab='PC3',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos)
    draw.3D(d1$X,d1$Y,d1$Z,class_label=d1$label,xlab='PC1',ylab='PC2',zlab='PC3',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  layout(1);
  if(return_type=='optimal') return(pred_label)
  if(return_type=='all') return(all_k_res)
}

#' Draw the cluster plot by UMAP (visulization algorithm) and Kmeans (cluster algorithm).
#'
#' \code{draw.umap.kmeans} is a visualization function to draw the cluster for the input data matrix.
#'
#' This function mainly aims for the visulization of sample clusters.
#' The input is the expression matrix for each row a gene/transcript/probe and each column a sample.
#' User need to input the real observation label for samples and this function will choose the best K by compared with predicted label and the observed label.
#' The output figure contains two sub-figures, left is labelled by the real observartion label and right is labelled by the predicted label
#' with comparison score (choose from ARI, NMI, Jaccard) shown above.
#' Not suggested when sample size is small.
#'
#' @param mat a numeric data matrix, the column (e.g sample) will be clustered by the features (e.g genes) in rows.
#' @param all_k a vector of integers, the pre-defined K to evaluate.
#' If NULL, will use all possible K. Default is NULL.
#' @param obs_label a vector of characters, a vector for sample categories with names representing the sample name.
#' @param legend_pos character, position for the legend displayed on the plot. Default is 'topleft'.
#' @param legend_cex numeric, cex for the legend displayed on the plot. Default is 0.8.
#' @param plot_type character, the type for the plot, choose from '2D' or '2D.ellipse' or '3D'. Default is '2D.ellipse'.
#' @param point_cex numeric, cex for the point in the plot. Default is 1.
#' @param kmeans_strategy character, the strategy to run the kmeans algorith, choose from 'basic' and 'consensus',
#' here the consensus kmeans is performed by functions in \code{ConsensusClusterPlus}. Default is 'basic'.
#' @param choose_k_strategy character, the strategy to choose the best K,
#' choose from 'ARI (adjusted rand index)', 'NMI (normalized mutual information)', 'Jaccard'. Default is 'ARI'.
#' @param return_type character, the strategy to return the results, choose from 'optimal' and 'all'.
#' If choose 'all', cluster results from all k in \code{all_k} will be returned.
#' Default is 'optimal'.
#' @return a vector of predicted label (if \code{return_type} is 'optimal') and a list of all possible K (if \code{return_type} is 'all')
#' @examples
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' mat <- exprs(network.par$net.eset)
#' phe <- pData(network.par$net.eset)
#' intgroup <- get_int_group(network.par$net.eset)
#' for(i in 1:length(intgroup)){
#'  print(intgroup[i])
#'  pred_label <- draw.umap.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]))
#'  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
#' }
#' pred_label <- draw.umap.kmeans(mat=mat,all_k = NULL,
#'                              obs_label=get_obs_label(phe,intgroup[i]),
#'                              kmeans_strategy='consensus')
#' @export
draw.umap.kmeans <- function(mat=NULL,all_k=NULL,obs_label=NULL,
                             legend_pos = 'topleft',legend_cex = 0.8,
                             kmeans_strategy='basic',choose_k_strategy='ARI',
                             plot_type='2D.ellipse',point_cex=1,return_type='optimal'){
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
  ori_cc <- umap.defaults;
  ori_cc$n_epochs <- 2000;
  ori_cc$n_neighbors <- min(15,round(ncol(mat)/2));
  if(plot_type=='3D') ori_cc$n_components <- 3
  cluster_mat <- umap(as.matrix(t(mat)),config=ori_cc)
  #
  all_jac <- list()
  all_k_res <- list()
  if(kmeans_strategy=='basic'){
    for(k in all_k){
      tmp_k <- list()
      for(i in 1:10){
        tmp_k[[i]] <- kmeans(cluster_mat$layout,centers=as.numeric(k))$cluster
      }
      pred_label <- tmp_k
      jac <- unlist(lapply(pred_label,function(x){get_clustComp(x, obs_label,strategy = choose_k_strategy)}))
      top_i <- which.max(jac)
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
  print(all_jac)
  use_k <- all_k[which.max(all_jac)]
  message(sprintf('Best Score occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  pred_label <- all_k_res[[as.character(use_k)]]

  d1 <- data.frame(id=colnames(mat),X=cluster_mat$layout[,1],Y=cluster_mat$layout[,2],label=pred_label,stringsAsFactors=FALSE)
  if(plot_type=='3D')   d1 <- data.frame(id=colnames(mat),X=cluster_mat$layout[,1],Y=cluster_mat$layout[,2],Y=cluster_mat$layout[,3],label=pred_label,stringsAsFactors=FALSE)
  layout(t(matrix(1:2)))
  if(plot_type=='2D.ellipse'){
    draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],xlab="",ylab="",legend_cex=legend_cex,point_cex=point_cex)
    draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,xlab="",ylab="",legend_cex=legend_cex,point_cex=point_cex,
                    main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  if(plot_type=='2D'){
    draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],xlab="",ylab="",legend_cex=legend_cex,point_cex=point_cex)
    draw.2D(d1$X,d1$Y,class_label=d1$label,xlab="",ylab="",legend_cex=legend_cex,point_cex=point_cex,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  if(plot_type=='3D'){
    draw.3D(d1$X,d1$Y,d1$Z,class_label=obs_label[d1$id],xlab='',ylab='',zlab='',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos)
    draw.3D(d1$X,d1$Y,d1$Z,class_label=d1$label,xlab='',ylab='',zlab='',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  layout(1);
  if(return_type=='optimal') return(pred_label)
  if(return_type=='all') return(all_k_res)
}


## get consensus Kmeans results, using M3C(toooo slow),ConsensusClusterPlus
# return cluster results, RCSI score, P-value
# plot==TRUE, plot RSCI+BETA_P
# get_consensus_cluster
get_consensus_cluster <-function(mat,all_k=2:12,clusterAlg="km",plot=FALSE,...)
{
  maxK <- max(all_k)
  res1 <- ConsensusClusterPlus(mat,maxK=maxK,clusterAlg=clusterAlg,plot=plot,...)
  cls_res <- lapply(all_k,function(x){
    x1 <- res1[[x]]$consensusClass
    x1
  })
  names(cls_res) <- as.character(all_k)
  cluster_res <- cls_res
  return(cluster_res)
}


#' Draw the cluster plot by MICA (cluster algorithm).
#'
#' \code{draw.MICA} is a visualization function to draw the cluster results for the MICA output.
#'
#' This function mainly aims for the visulization of sample clusters.
#' The input is the MICA project information.
#' User need to input the real observation label for samples and this function will choose the best K by compared with predicted label and the observed label.
#' The output figure contains two sub-figures, left is labelled by the real observartion label and right is labelled by the predicted label
#' with comparison score (choose from ARI, NMI, Jaccard) shown above.
#' Not suggested when sample size is small.
#'
#' @param outdir character, the output directory for running MICA.
#' @param prjname charater, the project name set for running MICA.
#' @param all_k a vector of integers, the pre-defined K to evaluate.
#' If NULL, will use all possible K. Default is NULL.
#' @param obs_label a vector of characters, a vector for sample categories with names representing the sample name.
#' @param legend_pos character, position for the legend displayed on the plot. Default is 'topleft'.
#' @param legend_cex numeric, cex for the legend displayed on the plot. Default is 0.8.
#' @param plot_type character, the type for the plot, choose from '2D' or '2D.ellipse' or '3D'. Default is '2D.ellipse'.
#' @param point_cex numeric, cex for the point in the plot. Default is 1.
#' @param choose_k_strategy character, the strategy to choose the best K,
#' choose from 'ARI (adjusted rand index)', 'NMI (normalized mutual information)', 'Jaccard'. Default is 'ARI'.
#' @param visualization_type character,choose from 'tsne', 'umap' or 'mds'. Default is 'tsne'.
#' @param return_type character, the strategy to return the results, choose from 'optimal' and 'all'.
#' If choose 'all', cluster results from all k in \code{all_k} will be returned.
#' Default is 'optimal'.
#' @return a vector of predicted label (if \code{return_type} is 'optimal') and a list of all possible K (if \code{return_type} is 'all')
#' @export
draw.MICA <- function(outdir=NULL,prjname=NULL,all_k=NULL,obs_label=NULL,
                         legend_pos = 'topleft',legend_cex = 0.8,
                         point_cex=1,plot_type='2D.ellipse',
                         choose_k_strategy='ARI',
                         visualization_type='tsne',return_type='optimal') {
  if(plot_type=='3D' & visualization_type=='tsne'){
    message('Current tsne not support for 3D');return(FALSE)
  }
  # choose best k
  res1 <- get_clustComp_MICA(outdir=outdir, all_k=all_k, obs_label=obs_label, prjname = prjname,strategy = choose_k_strategy)
  all_k_res <- res1$all_k_res
  all_jac <- res1$all_jac
  use_k <- all_k[which.max(all_jac)]
  message(sprintf('Best Score occurs when k=%s, with value=%s',use_k,all_jac[as.character(use_k)]))
  #
  use_file <- sprintf('%s/scMINER_%s/scMINER_%s_MDS_%s/scMINER_MICA_out/%s.ggplot.txt',
                      outdir,prjname,prjname,use_k,prjname)
  d1 <- read.delim(use_file, stringsAsFactors = FALSE) ## get cluster results
  if(visualization_type=='umap' | visualization_type=='mds'){
    use_file <- sprintf('%s/scMINER_%s/scMINER_%s_MDS_%s/scMINER_MICA_out/%s_clust.h5',
                        outdir,prjname,prjname,use_k,prjname)
    fid <- H5Fopen(use_file)
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
        use_mat_umap <- umap(t(dist_mat),config=ori_cc)
        X <- use_mat_umap$layout[,1];Y <- use_mat_umap$layout[,2]; Z <- use_mat_umap$layout[,3];
        d1$X <- X; d1$Y <- Y;d1$Z <- Z
      }else{
        use_mat_umap <- umap(t(dist_mat),config=ori_cc)
        X <- use_mat_umap$layout[,1];Y <- use_mat_umap$layout[,2]
        d1$X <- X; d1$Y <- Y;
      }
    }
    H5Fclose(fid)
  }
  layout(t(matrix(1:2)))
  if(plot_type=='2D.ellipse'){
    draw.2D.ellipse(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex)
    draw.2D.ellipse(d1$X,d1$Y,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
                    main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  if(plot_type=='2D'){
    draw.2D(d1$X,d1$Y,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex)
    draw.2D(d1$X,d1$Y,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  if(plot_type=='2D.text'){
    draw.2D.text(d1$X,d1$Y,class_label=obs_label[d1$id],class_text=d1$id,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex)
    draw.2D.text(d1$X,d1$Y,class_label=d1$label,class_text=d1$id,xlab='MICA-1',ylab='MICA-2',legend_cex=legend_cex,point_cex=point_cex,
                 main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
  }
  if(plot_type=='3D'){
    draw.3D(d1$X,d1$Y,d1$Z,class_label=obs_label[d1$id],xlab='MICA-1',ylab='MICA-2',zlab='MICA-3',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos)
    draw.3D(d1$X,d1$Y,d1$Z,class_label=d1$label,xlab='MICA-1',ylab='MICA-2',zlab='MICA-3',legend_cex=legend_cex,point_cex=point_cex,legend_pos=legend_pos,
            main=sprintf('%s:%s',choose_k_strategy,format(all_jac[as.character(use_k)],digits=4)))
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
    use_file <-
      sprintf(
        '%s/scMINER_%s/scMINER_%s_MDS_%s/scMINER_MICA_out/%s.ggplot.txt',
        outdir,prjname,prjname,k,prjname
      )
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

#' Volcano plot for top DE (differentiated expressed) genes or DA (differentiated activity) drivers
#'
#' \code{draw.volcanoPlot} draw the volcano plot for top DE genes or DA drivers, could display the name of the top items in figures and will retrun the list of top items.
#'
#' This plot function input the master table and draw the volcano plot by setting significant threshold for logFC and P-value.
#'
#' @param dat data.frame, prefer the formatted master table generated by \code{generate.masterTable}, if not, must contain the columns for the following required parameters.
#' @param label_col character, the name of the column in \code{dat}, which contains the gene/driver label for display
#' @param logFC_col character, the name of the column in \code{dat}, which contains the logFC value
#' @param Pv_col character, the name of the column in \code{dat}, which contains the P-value
#' @param logFC_thre numeric, cutoff value for the logFC. Genes or drivers with absolute logFC value higher than the cutoff are remained.Default is 1.5.
#' @param Pv_thre numeric, cutoff value for the p-values. Genes or drivers with lower p-values are remained.Default is 0.01.
#' @param xlab character, title for the X axis
#' @param ylab character, title for the Y axis
#' @param show_label logical, whether or not to display the genes or drivers passed the cutoff on the plot. Default is FALSE
#' @param label_cex numeric, \code{cex} for the label displayed on the plot. Default is 0.5
#' @param legend_cex numeric, \code{cex} for the legend displayed on the plot. Default is 0.7
#' @param label_type character, the strategy for label display on the plot, by choosing 'origin' or 'distribute'. Default is 'distribute'.
#' @param main character, \code{main} for the title on the plot.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return data.frame containing the significant genes or drivers with the following components:
#'
#' \item{label_col}{}
#' \item{logFC_col}{}
#' \item{Pv_col}{}
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
draw.volcanoPlot <- function(dat=NULL,label_col=NULL,logFC_col=NULL,Pv_col=NULL,logFC_thre=1.5, Pv_thre=0.01,
                             xlab='log2 Fold Change',ylab='P-value',show_label=FALSE,label_cex=0.5,legend_cex=0.7,
                             label_type='distribute',main="",pdf_file=NULL){
  dat <- unique(dat[,c(label_col,logFC_col,Pv_col)])
  dat <- dat[order(dat[,3],decreasing=TRUE),]
  dat <- dat[which(is.na(dat[,2])==FALSE),]
  x <- as.numeric(dat[,logFC_col])
  y <- as.numeric(dat[,Pv_col]);
  y <- -log10(y)
  s1 <- which(abs(x)>=logFC_thre & y>= -log10(Pv_thre))
  geneWidth <- 0
  geneHeight <- 0
  if(length(s1)>0){
    s11 <- s1[which(x[s1]>=0)]
    s12 <- s1[which(x[s1]<0)]
    geneWidth  <- max(strwidth(dat[s1,label_col],'inches',cex=label_cex))
    geneHeight <- max(strwidth(toupper(letters),'inches',cex=label_cex))*max(length(s11),length(s12))*1.2
    if(is.null(pdf_file)==FALSE){
      if(show_label==TRUE & label_type=='distribute'){
        pdf(pdf_file,width=10+geneWidth*2,height=max(10,geneHeight))
      }else{
        pdf(pdf_file,width=10,height=10)
      }
    }
  }
  par(mai=c(1.5,2,1.5,1))
  mm <- max(abs(x))
  if(show_label==TRUE & label_type=='distribute'){
    plot(y~x,pch=16,col=get_transparent('grey',0.7),xlab=xlab,ylab="",
         xlim=c(-3*mm/7*geneWidth-1.5*mm,1.5*mm+3*mm/7*geneWidth),ylim=c(0,max(y)*1.5),yaxt='n',main=main,cex.lab=1.2,cex.main=1.6)
  }else{
    plot(y~x,pch=16,col=get_transparent('grey',0.7),xlab=xlab,ylab="",
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
###
#' Heatmap plot for displaying expression level or activity level for genes and drivers
#'
#' \code{draw.heatmap} draw the heatmap for expression level or activity level for genes and drivers across selected samples.
#'
#' This plot function input the expression/activity matrix, with each row as one gene or driver, each column as one sample.
#'
#' @param mat numeric matrix, each row as one gene or driver, each column as one sample
#' @param use_genes a vector of characters, the list of genes for display. Default is the rownames(mat).
#' @param use_gene_label a vector of characters, label for the list of genes for display. Default is the use_genes.
#' @param use_samples a vector of characters, the list of samples for display. Default is the colnames(mat).
#' @param use_sample_label a vector of characters, label for the list of samples for display. Default is the use_samples.
#' @param phenotype_info data.frame,contain the sample phenotype information, can be generated by \code{pData(eset)}.
#' The rownames should match the colnames of mat. Default is NULL.
#' @param use_phe a list of characters, selected phenotype for display,must be the subset of colnames of phenotype_info.Default is NULL.
#' @param main character, title for the draw. Default is "".
#' @param scale character, indicating if the values should be centered and scaled in either the row direction or the column direction, or none.
#' The default is "none".
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#' @param cluster_rows,cluster_columns parameters used in \code{Heatmap}, please check for details. Default is TRUE.
#' @param clustering_distance_rows,clustering_distance_columns parameters used in \code{Heatmap}, please check for details. Default is 'pearson'.
#' @param show_row_names,show_column_names parameters used in \code{Heatmap}, please check for details. Default is TRUE.
#' @param row_names_gp,column_names_gp parameters used in \code{Heatmap}, please check for details. Default is gpar(fontsize = 12).
#' @param ..., more parameters used in \code{Heatmap}
#'
#' @return logical value indicating whether the plot has been successfully generated
#'
#' @examples
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1/driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames matches originalID in ms_tab
#' ac_mat <- exprs(analysis.par$merge.ac.eset) ## ac,the rownames matches originalID_label in ms_tab
#' phe_info <- pData(analysis.par$cal.eset) ## phenotype information
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
#'             row_names_gp = gpar(fontsize = 12))
#' draw.heatmap(mat=ac_mat,use_genes=ms_tab[rownames(sig_driver),'originalID_label'],
#'              use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
#'              use_samples=colnames(exp_mat),
#'              use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
#'              phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup'),
#'              main='Activity for Top drivers',scale='row',
#'              cluster_rows=TRUE,cluster_columns=TRUE,
#'              clustering_distance_rows='pearson',
#'              clustering_distance_columns='pearson',
#'              row_names_gp = gpar(fontsize = 6))
#'
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1/driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames matches originalID in ms_tab
#' ac_mat <- exprs(analysis.par$merge.ac.eset) ## ac,the rownames matches originalID_label in ms_tab
#' phe_info <- pData(analysis.par$cal.eset) ## phenotype information
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
#'             analysis.par$out.dir.PLOT))
#' draw.heatmap(mat=ac_mat,use_genes=ms_tab[rownames(sig_driver),'originalID_label'],
#'              use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
#'              use_samples=colnames(exp_mat),
#'              use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
#'              phenotype_info=phe_info,
#'              use_phe=c('gender','pathology','subgroup'),
#'              main='Activity for Top drivers',scale='row',
#'              cluster_rows=TRUE,cluster_columns=TRUE,
#'              clustering_distance_rows='pearson',
#'              clustering_distance_columns='pearson',
#'              row_names_gp = gpar(fontsize = 6),
#'              pdf_file=sprintf('%s/heatmap_demo2.pdf',
#'              analysis.par$out.dir.PLOT))
#'}
#' @export
draw.heatmap <- function(mat=NULL,use_genes=rownames(mat),use_gene_label=use_genes,use_samples=colnames(mat),use_sample_label=use_samples,
                         phenotype_info=NULL,use_phe=NULL,main="",scale='none',pdf_file=NULL,
                         cluster_rows=TRUE,cluster_columns=TRUE,
                         clustering_distance_rows='pearson',clustering_distance_columns='pearson',
                         row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12),
                         show_row_names=TRUE,show_column_names=TRUE,...){
  names(use_gene_label) <- use_genes
  names(use_sample_label) <- use_samples
  if(is.null(rownames(phenotype_info))==FALSE){
    phenotype_info <- phenotype_info[colnames(mat),]
  }
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
#' find.gsByGene
#' @param gene character
#' @param use_gs a vector of characters
#' @export
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
#' Merge GeneSets by choosing several categories.
#'
#' \code{merge_gs} is a simple function to merge the gene set list.
#'
#' @param all_gs2gene list, which could be obtained by running \code{gs.preload()}.
#' @param use_gs a vector of characters, could check \code{all_gs2gene_info} for the cateogory description.
#' Default is 'H', 'CP:BIOCARTA', 'CP:REACTOME', 'CP:KEGG'.
#' @return a one-level list for geneset to genes.
#' @examples
#' gs.preload(use_spe='Homo sapiens',update=FALSE)
#' use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,
#'                        use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
#'
#' @export
merge_gs <- function(all_gs2gene=all_gs2gene,use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5')){
  if(is.null(use_gs)==TRUE){
    nn <- unlist(lapply(all_gs2gene,names))
    use_gs2gene <- unlist(all_gs2gene,recursive = FALSE)
    names(use_gs2gene)<-nn
  }else{
    nn <- unlist(lapply(all_gs2gene[use_gs],names))
    use_gs2gene <- unlist(all_gs2gene[use_gs],recursive = FALSE)
    names(use_gs2gene)<-nn
  }
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

#' Gene set enrichment analysis by Fisher's Exact Test.
#'
#' \code{funcEnrich.Fisher} will perform gene set enrichment analysis by Fisher's Exact Test.
#'
#' This is a function to find significant enriched gene sets for input gene list. Users could prepare gs2gene or use all_gs2gene preloaded by using \code{gs.preload}.
#' Background gene list is accepeted.
#'
#' @param input_list a vector of characters, the list of genes for analysis. Only accept gene symbols, and gene ID conversion could be done by preparing a transfer table
#' by using \code{get_IDtransfer} and using \code{get_name_transfertab} to transfer the gene IDs.
#' @param bg_list a vector of characters, the background list of genes for analysis. Only accept gene symbols.
#' Default is NULL, will use all genes in the gs2gene as the background list.
#' @param gs2gene a list for geneset to genes, the name for the list is the gene set name and the content in each list is the vector for genes belong to that gene set.
#' If NULL, will use all_gs2gene loaded by using \code{gs.preload}. Default is NULL.
#' @param use_gs a vector of characters, the name for gene set category used for anlaysis.
#' If gs2gene is set to NULL, use_gs must be the subset of \code{names(all_gs2gene)}.
#' Could check \code{all_gs2gene_info} for the cateogory description.
#' If set to 'all', all gene sets in gs2gene will be used.
#' Default is 'H', 'CP:BIOCARTA', 'CP:REACTOME', 'CP:KEGG'
#' if gs2gene is set to NULL (use all_gs2gene).
#' If user input own gs2gene list, use_gs will be set to 'all' as default.
#' @param min_gs_size numeric, minimum gene set size for analysis, default is 5.
#' @param max_gs_size numeric, maximum gene set size for analysis, default is 500.
#' @param Pv_adj character, p-value adjustment method, could check \code{p.adjust.methods} for the available options. Default is 'fdr'.
#' @param Pv_thre numeric, cutoff value for the adjusted p-values for significance. Default is 0.1.
#'
#' @return The function will return a list of gene sets with significant statistics, detailed as follows,
#'
#' \item{#Name}{Name for the enriched gene set}
#' \item{Total_item}{Number of background size}
#' \item{Num_item}{Number of genes in the gene set (filtered by the background list)}
#' \item{Num_list}{Number of input genes for testing (filtered by the background list)}
#' \item{Num_list_item}{Number input genes annotated by the gene set (filtered by the background list)}
#' \item{Ori_P}{Original P-value from Fisher's Exact Test}
#' \item{Adj_p}{Adjusted P-value}
#' \item{Odds_Ratio}{Odds ratio by the 2*2 matrix used for Fisher's Exact Test}
#' \item{Intersected_items}{List of the intersected genes, collapsed by ';', the number is equal to Num_list_item}
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
  if(is.null(gs2gene)==TRUE){ ## use inner gs2gene
    if(is.null(use_gs)==TRUE){
      use_gs <- c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG')
    }else{
      if(use_gs[1] == 'all'){
        use_gs <- c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)
      }
    }
    if(length(setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)))>0){
      message(sprintf('Input %s not in all_gs2gene, please check all_gs2gene_info (items in Category or Sub-Category) and re-try!',
                      paste(setdiff(use_gs,c(all_gs2gene_info$Category,all_gs2gene_info$`Sub-Category`)),collapse=';')));
      return(FALSE)
    }
    if(length(use_gs)>1){
      gs2gene <- merge_gs(all_gs2gene,use_gs = use_gs)
    }else{
      gs2gene <- all_gs2gene[[use_gs]]
    }
  }else{
    if(is.null(use_gs)==TRUE){
      use_gs <- 'all'
    }
    if(length(use_gs)>1){
        gs2gene <- merge_gs(gs2gene,use_gs = use_gs)
    }else{
      if(use_gs == 'all'){
        gs2gene <- merge_gs(gs2gene,use_gs = NULL)
      }else{
        gs2gene <- gs2gene[[use_gs]]
      }
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

#' Bar plot for the result of gene set enrichment analysis.
#'
#' \code{draw.funcEnrich.bar} will draw the barplot for the result of gene set enrichment analysis.
#'
#' This is a function to draw the barplot for the result of gene set enrichment analysis.
#' Two modes, one just display the P-value, and the other could show the top intersected genes for each gene set.
#'
#' @param funcEnrich_res data.frame, containing the result for the function enrichment analysis. Prefer the format generated by using \code{funcEnrich.Fisher}.
#' If not, users could prepare the required columns and indicate the column names in the following parameters.
#' @param top_number numeric, number for the top enriched gene sets to be displayed on the plot. Default is 30.
#' @param Pv_col character, the name of the column in \code{funcEnrich_res}, which contains the P-value. Default is 'Ori_P'.
#' @param name_col character, the name of the column in \code{funcEnrich_res}, which contains the gene set name. Default is '#Name'.
#' @param item_col character, the name of the column in \code{funcEnrich_res}, which contains the detailed intersected gene list, collapsed by ';'.
#' Default is 'Intersected_items'.
#' @param Pv_thre numeric, cutoff value for the p-values. Genes or drivers with lower p-values are remained. Default is 0.1.
#' @param display_genes logical, whether or not to display the intersected genes on the plot. Default is FALSE
#' @param gs_cex numeric, \code{cex} for the gene sets displayed on the plot. Default is 0.5
#' @param gene_cex numeric, \code{cex} for the genes displayed on the plot. Default is 0.5
#' @param main character, \code{main} for the title on the plot.
#' @param bar_col character, color code for the bar on the plot. Default is brewer.pal(8,'RdBu')[7].
#' @param eg_num numeric, example number of intersected genes shown on the plot. Default is 5.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return logical value indicating whether the plot has been successfully generated
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
#' draw.funcEnrich.bar(funcEnrich_res=res1,top_number=5,
#'                    main='Function Enrichment for Top drivers',
#'                    display_genes = TRUE,eg_num=3,
#'                    gs_cex=0.6,gene_cex=0.5)
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
#' Cluster plot for the result of gene set enrichment analysis.
#'
#' \code{draw.funcEnrich.cluster} will draw the cluster for the result of gene set enrichment analysis.
#'
#' This is a function to draw the cluster (genes and gene sets) for the result of gene set enrichment analysis.
#' The cluster is based on the binary matrix representing the gene's existence in the enriched gene sets.
#' Detailed matrix for the cluster, enriched P-value will be displayed on the plot.
#'
#' @param funcEnrich_res data.frame, containing the result for the function enrichment analysis. Prefer the format generated by using \code{funcEnrich.Fisher}.
#' If not, users could prepare the required columns and indicate the column names in the following parameters.
#' @param top_number numeric, number for the top enriched gene sets to be displayed on the plot. Default is 30.
#' @param Pv_col character, the name of the column in \code{funcEnrich_res}, which contains the P-value. Default is 'Ori_P'.
#' @param name_col character, the name of the column in \code{funcEnrich_res}, which contains the gene set name. Default is '#Name'.
#' @param item_col character, the name of the column in \code{funcEnrich_res}, which contains the detailed intersected gene list, collapsed by ';'.
#' Default is 'Intersected_items'.
#' @param Pv_thre numeric, cutoff value for the p-values. Genes or drivers with lower p-values are remained. Default is 0.1.
#' @param gs_cex numeric, \code{cex} for the gene sets displayed on the plot. Default is 0.7.
#' @param gene_cex numeric, \code{cex} for the genes displayed on the plot. Default is 0.8.
#' @param pv_cex numeric, \code{cex} for the P-value displayed on the plot. Default is 0.7.
#' @param main character, \code{main} for the title on the plot.
#' @param h numeric, cutoff for the cluster. This parameter will be used in the \code{cutree} function. Default is 0.95
#' @param cluster_gs logical, whether or not to cluster gene sets. Default is TRUE.
#' @param cluster_gene logical, whether or not to cluster genes. Default is TRUE.
#' @param use_genes a vector of characters, gene list used for display in plot,
#' if NULL will display all genes in the top enriched gene sets.Default is NULL.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#' @param return_mat logical, whether or not to return the matrix used for display. Default is FALSE.
#'
#' @return if return_mat==FALSE, will return logical value indicating whether the plot has been successfully generated,
#' otherwise will return the matrix used for cluster.
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
draw.funcEnrich.cluster <- function(funcEnrich_res=NULL,top_number=30,Pv_col='Ori_P',name_col='#Name',item_col='Intersected_items',Pv_thre=0.1,
                                    gs_cex=0.7,gene_cex=0.8,pv_cex=0.7,main="",h=0.95,cluster_gs=TRUE,cluster_gene=TRUE,
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
  plot(1,xaxt='n',yaxt='n',bty='n',xlim=c(pp[1],pp[2]),ylim=c(pp[3],pp[4]),col='white',xlab="",ylab="")
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
  plot(1,xaxt='n',yaxt='n',bty='n',xlim=c(pp[1],pp[2]),ylim=c(pp[3],pp[4]),col='white',xlab="",ylab="")
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
#' Bubble plot for the top drivers in NetBID2 analysis.
#'
#' \code{draw.bubblePlot} will draw the buble plot for the top drivers and the enriched gene sets for the targets of each driver.
#'
#' This is a function to draw the bubble plot for the top significant drivers. Each row is a gene set, and each column is a driver.
#' Each bubble represents the enrichment for each driver's target gene in the corresponding gene set.
#' The size for each bubble shows the intersected size for the target gene and the gene set.
#' The color for each bubble shows the significance of enrichment performed by Fisher's Exact Test.
#' Color bar and size bar are shown in the draw.
#' Besides, the target size and the driver gene/transcript bio-type is shown at the bottom of the draw.
#'
#' @param driver_list a vector of characters, the name for the top drivers.
#' @param show_label a vector of characters, the name for the top drivers to display on the plot.
#' If NULL, will display the name in driver_list. Default is NULL.
#' @param Z_val a vector of numeric values, the Z statistics for the driver_list.
#' Better to give name to the vector, otherwise will automatically use driver_list as the name.
#' @param driver_type a vector of characters, the bio-type or other characteristics for the driver.
#' If not NULL, will display the type on the plot.
#' Better to give name to the vector, otherwise will automatically use driver_list as the name.
#' @param target_list a list for the target gene information for the drivers. The names for the list must contain the driver in driver_list.
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param transfer2symbol2type data.frame, the transfer table for the original ID to the gene symbol and gene biotype (at gene level)
#' or transcript symbol and transcript biotype (at transcript level). strongly suggest to use \code{get_IDtransfer2symbol2type} to generate the transfer table.
#' @param gs2gene a list for geneset to genes, the name for the list is the gene set name and the content in each list is the vector for genes belong to that gene set.
#' If NULL, will use all_gs2gene loaded by using \code{gs.preload}. Default is NULL.
#' @param use_gs a vector of characters, the name for gene set category used for anlaysis.
#' If gs2gene is set to NULL, use_gs must be the subset of \code{names(all_gs2gene)}.
#' Could check \code{all_gs2gene_info} for the cateogory description.
#' Default is 'H', 'CP:BIOCARTA', 'CP:REACTOME', 'CP:KEGG'.
#' @param display_gs_list a vector of characters, the list of gene set names to display on the plot.
#' If NULL, the gene sets are shown according to the significant ranking.
#' Default is NULL.
#' @param bg_list a vector of characters, the background list of genes for analysis. Only accept gene symbols.
#' Default is NULL, will use all genes in the gs2gene as the background list.
#' @param min_gs_size numeric, minimum gene set size for analysis, default is 5.
#' @param max_gs_size numeric, maximum gene set size for analysis, default is 500.
#' @param Pv_adj character, p-value adjustment method, could check \code{p.adjust.methods} for the available options. Default is 'none'.
#' @param Pv_thre numeric, cutoff value for the adjusted p-values for significance.Default is 0.1.
#' @param top_geneset_number number for the top enriched gene sets to be displayed on the plot. Default is 30.
#' @param top_driver_number number for the top significant drivers to be displayed on the plot. Default is 30.
#' @param main character, \code{main} for the title on the plot.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#' @param mark_gene a vector of characters, if not NULL, the drivers in the mark_gene will be labelled red in the draw. Default is NULL.
#' @param driver_cex numeric, \code{cex} for the driver displayed on the plot. Default is 1.
#' @param gs_cex numeric, \code{cex} for the gene sets displayed on the plot. Default is 1.
#'
#' @return
#' Will return logical value indicating whether the plot has been successfully generated
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
#' use_genes <- unique(analysis.par$merge.network$network_dat$target.symbol)
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
#'                bg_list=ms_tab[,'geneSymbol'],min_gs_size=5,
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
#' use_genes <- unique(analysis.par$merge.network$network_dat$target.symbol)
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
#'                bg_list=ms_tab[,'geneSymbol'],
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
                            pdf_file=NULL,main="",mark_gene=NULL,driver_cex=1,gs_cex=1){
  ## check NULL

  if(is.null(names(show_label))==TRUE){names(show_label) <- driver_list}
  if(is.null(names(Z_val))==TRUE){names(Z_val) <- driver_list}
  if(is.null(driver_type)==FALSE){
    if(is.null(names(driver_type))==TRUE){names(driver_type) <- driver_list}
  }
  driver_list <- driver_list[order(abs(Z_val),decreasing = TRUE)]
  if(length(driver_list)>top_driver_number){
    driver_list <- driver_list[1:top_driver_number]
  }
  driver_list <- driver_list[order(Z_val[driver_list],decreasing=TRUE)]
  ## get target gene for driver_list
  transfer_tab <- transfer2symbol2type
  #rownames(transfer_tab) <- transfer_tab[,1]
  target_gene <- lapply(driver_list,function(x){
    x1 <- target_list[[x]]$target
    x1 <- x1[which(x1 %in% transfer_tab[,1])]
    x2 <- transfer_tab[which(transfer_tab[,1] %in% x1),]
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
  if(is.null(display_gs_list)==FALSE){
    all_path <- intersect(display_gs_list,all_path)
    if(length(all_path)<3){
      message('The number for passed gene sets is smaller than 3, please check the display_gs_list and re-try!')
      return(FALSE)
    }
  }
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
  plot(1,bty='n',col='white',xlim=c(0,nc+1),ylim=c(-2,nr+1),xaxt='n',yaxt='n',xlab="",ylab="",main=main)
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

#' GSEA (gene set enrichment analysis) plot for a gene set or a driver.
#'
#' \code{draw.GSEA} will generate a GSEA plot for a gene set (with annotated gene list) or a driver (with target gene list).
#'
#' This is a plot function to draw GSEA for a gene set or a driver by input differentiated expression profile.
#' User could input the annotation text for the significance or the function could display the significance calculated by Kolmogorov-Smirnov tests.
#' ATTENTION: when user input the \code{use_direction}, the rank profile will be duplicated (here Zero cross at the middle) and
#' genes with negative direction (-1 in the use_direction list) will be put in the mirror position of the original place.
#'
#' @param rank_profile a vector of numerics. The ranking profile for the differentiated values in a specific sample condition comparison.
#' The names of the vector must be the gene names. The value in the vector could be the logFC or t-statistics.
#' @param use_genes a vector of characters, the list of genes for analysis. The ID must be the subset of the names of \code{rank_profile}.
#' This could either be the annotated gene list for the gene set or the target gene list for the driver.
#' @param use_direction a vector of numerics, indicate the direction for the driver and the target gene list.
#' 1 indicates positive regulation and -1 indicates negative regulation.
#' If NULL, will not consider direction information. Default is NULL.
#' @param main character, title for the draw. Default is "".
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#' @param annotation character, annotation for the significance to be displayed on the plot.
#' If NULL, will perform Kolmogorov-Smirnov tests to get significance. Default is NULL.
#' @param annotation_cex numeric, \code{cex} for the annotation displayed on the plot. Default is 1.2
#' @param left_annotation character, annotation displayed on the left of the figure representing left condition of the rank_profile. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure representing right condition of the rank_profile. Default is "".
#'
#' @return logical value indicating whether the plot has been successfully generated

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
  #### start plot
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]

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
  plot(rank_profile,col='grey',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(-mm,mm))
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
  if(is.null(use_direction)==FALSE){
    axis(side=1,at=x1,labels=get_label_manual(x1/2))
  }else{
    axis(side=1,at=x1,labels=get_label_manual(x1))
  }
  # get zero cross
  w1 <- which.min(abs(rank_profile))
  abline(v=w1,lty=2,col='grey')
  if(is.null(use_direction)==FALSE){
    text(w1,-mm/4,sprintf('Zero cross at %s',round(w1/2)),adj=0.5)
  }else{
    text(w1,-mm/4,sprintf('Zero cross at %d',w1),adj=0.5)
  }
  if(is.null(use_direction)==FALSE){
    legend(w1,pp[3]-(pp[4]-pp[3])/4,lty=1,lwd=2,c('Enrichment profile','Hits_positive_direction','Hits_negative_direction','Ranking metric scores'),
           col=c('green',pos_col,neg_col,'grey'),xpd=TRUE,horiz=TRUE,xjust=0.5,cex=1.2)
  }else{
    legend(w1,pp[3]-(pp[4]-pp[3])/4,lty=1,lwd=2,c('Enrichment profile','Hits','Ranking metric scores'),
           col=c('green','black','grey'),xpd=TRUE,horiz=TRUE,xjust=0.5,cex=1.2)
  }

  pm <- par()$usr
  ## get image bar
  par(mar=c(0,6,0,2))
  use_col <- z2col(rank_profile,sig_thre = 0,n_len = 30,blue_col='blue',red_col='red')
  image(x=as.matrix(1:r_len),col=use_col,bty='n',xaxt='n',yaxt='n',xlim=c(pm[1],pm[2])/r_len)
  abline(v=use_pos/r_len,col='grey')
  ## mark gene position;
  par(mar=c(0,6,0,2))
  plot(1,col='white',xlab="",ylab="",bty='n',xlim=c(1,r_len),xaxt='n',yaxt='n')
  if(is.null(use_direction)==FALSE){
    use_pos_P <- which(names(rank_profile) %in% use_genes[grep('POS',use_genes)])
    use_pos_N <- which(names(rank_profile) %in% use_genes[grep('NEG',use_genes)])
    abline(v=use_pos_P,col=pos_col)
    abline(v=use_pos_N,col=neg_col)
  }else{
    abline(v=use_pos)
  }
  ## GSEA ES
  par(mar=c(0,6,5,2))
  # get ES score
  es_res <- get_ES(rank_profile,use_genes)
  y2 <- seq(min(es_res$RES),max(es_res$RES),length.out=7); y2 <- round(y2,1)
  if(is.null(use_direction)==FALSE){
    plot(es_res$RES,col='green',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',
         xlim=c(1,r_len),type='l',lwd=3,ylim=c(min(es_res$RES),max(y2)),main=main,xpd=TRUE)
  }else{
    plot(es_res$RES,col='green',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',
         xlim=c(1,r_len),type='l',lwd=3,ylim=c(min(es_res$RES),max(y2)),main=main,xpd=TRUE)
  }

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

#' GSEA (gene set enrichment analysis) plot for the list of drivers.
#'
#' \code{draw.GSEA.NetBID} will generate a GSEA plot for the list of drivers, including the target genes' position on the differentiated expression profile, with
#' statistics of differentiated expression (DE) and differentiated activity (DA) for each driver.
#'
#' This is a plot function to draw GSEA for the list of drivers by the input of differentiated expression profile.
#' User could choose to display the target genes in one row or two rows, by selecting black color or red to blue color bar.
#'
#' @param DE data.frame,the differentiated expression results.
#' This data.frame could be generated by using \code{getDE.limma.2G} or \code{getDE.BID.2G}.
#' If user want to generate this data.frame by other strategies, the rownames must be the gene names or need one column to be the gene name
#' (set in \code{name_col}) and must contain the columns indicating the differentiated expression profile.
#' @param name_col character, the name of the column in \code{DE}, which contains the gene name. If NULL, will use the rownames of \code{DE}.
#' Default is NULL.
#' @param profile_col character, the name of the column in \code{DE}, which will be used as the differentiated expression profile.
#' If DE is created by \code{getDE.limma.2G} or \code{getDE.BID.2G}, this parameter could be 'logFC' or 't'.
#' @param profile_trend character, the choice of how to display the profile, from high/positive to low/negative ('pos2neg')
#' or low/negative to high/positive ('neg2pos').Default is 'pos2neg'.
#' @param driver_list a vector of characters, the name for the top drivers.
#' @param show_label a vector of characters, the name for the top drivers to display on the plot.
#' If NULL, will display the name in driver_list. Default is NULL.
#' @param driver_DA_Z a vector of numeric values, the Z statistics of differentiated activity (DA) for the driver_list.
#' Better to give name to the vector, otherwise will automatically use driver_list as the name.
#' @param driver_DE_Z a vector of numeric values, the Z statistics of differentiated expression (DE) for the driver_list.
#' Better to give name to the vector, otherwise will automatically use driver_list as the name.
#' @param target_list a list for the target gene information for the drivers. The names for the list must contain the driver in driver_list.
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param top_driver_number numeric, number for the top significant drivers to be displayed on the plot. Default is 30.
#' @param target_nrow numeric, number of rows for each driver display on the plot. Two options, 1 or 2.
#' If set to 1, the target genes' position on the profile will be displayed in one row.
#' If set to 2, the target genes' position on the profile will be displayed in two rows,
#' with positive regulated genes displayed on the first row and negative regulated genes displayed on the second row.
#' Default is 2.
#' @param target_col character, choice of color pattern used to display the targets. Two options,'black' or 'RdBu'.
#' If set to 'black', the lines will be colored in black.
#' If set to 'RdBu', the lines will be colored into Red to Blue color bar.
#' If \code{target_col_type} is set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If \code{target_col_type} is set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' with significant high set for red and low for blue. The significant threshold is set by \code{profile_sig_thre}.
#' Default is 'RdBu'.
#' @param target_col_type character, choice of the pattern used to display the color for target genes, only work when \code{target_col} is set as 'RdBu'.
#' Two options,'PN' or 'DE'.
#' If set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' Default is 'PN'.
#' @param left_annotation character, annotation displayed on the left of the figure representing left condition of the rank_profile. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure representing right condition of the rank_profile. Default is "".
#' @param main character, title for the plot. Default is "".
#' @param profile_sig_thre numeric, threshold for the absolute values in profile to be treated as significance.
#' Target genes without signifcant values in the profile will be colored in grey. Only work when \code{target_col_type} is set as "DE" and \code{target_col} is set as "RdBu".
#' Default is 0.
#' @param Z_sig_thre numeric, threshold for the Z statistics in \code{driver_DA_Z} and \code{driver_DE_Z} to be treated as signifcance.
#' Only signifcant values will have background color. Default is 1.64.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return logical value indicating whether the plot has been successfully generated

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
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
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
  if(is.null(name_col)==TRUE){
    DE <- cbind(DE[,setdiff(colnames(DE),'name')],name=rownames(DE),stringsAsFactors=FALSE)
    name_col <- 'name'
  }
  w1 <- which(is.na(DE[,profile_col])==FALSE)
  DE <- DE[w1,]
  if(profile_trend=='pos2neg') DE <- DE[order(DE[,profile_col],decreasing = TRUE),] else DE <- DE[order(DE[,profile_col],decreasing = FALSE),]
  DE_profile <- DE[,profile_col]
  DE_profile_name <- DE[,name_col]
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
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,n_gene),xaxt='n',yaxt='n')
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
    #names(cc) <- DE_profile_name
    cc[which(cc=='white')] <- 'light grey'
  }
  if(target_nrow==1){
    for(i in 1:length(driver_list)){
      t1 <- use_target_list[[driver_list[[i]]]]
      w0 <- which(DE_profile_name %in% t1$target)
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
                             red_col = pos_col,blue_col=neg_col))
        }
      }
    }
  }
  if(target_nrow==2){
    for(i in 1:length(driver_list)){
      t1 <- use_target_list[[driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(DE_profile_name %in% t11)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy1[2*i+1],col=pos_col,lwd=1.5)
          }
        }
      }
      w0 <- which(DE_profile_name %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy1[2*i],col=neg_col,lwd=1.5)
          }
        }
      }
    }
  }
  ## plot 3
  par(mar=c(2,0.5,2,2))
  plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
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
           ybottom=yy22[i],ytop=yy22[i]+dyy*0.35,col=pos_col,border=NA)
      rect(xleft=(pp[1]+pp[2])*0.55,xright=(pp[1]+pp[2])*0.55+target_size[i,2]/mm*tt,
           ytop=yy22[i],ybottom=yy22[i]-dyy*0.35,col=neg_col,border=NA)
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
#' GSEA (gene set enrichment analysis) plot for the list of gene sets.
#'
#' \code{draw.GSEA.NetBID.GS} will generate a GSEA plot for the list of gene sets, including the annotated genes' position on the differentiated expression profile, with
#' statistics of differentiated activity (DA) for each gene set.
#'
#' This is a plot function to draw GSEA for the list of gene sets by the input of differentiated expression profile.
#' User could choose to display the annotated genes by selecting black color or red to blue color bar.
#' @param DE data.frame,the differentiated expression results.
#' This data.frame could be generated by using \code{getDE.limma.2G} or \code{getDE.BID.2G}.
#' If user want to generate this data.frame by other strategies, the rownames must be the gene names and
#' must contain the columns indicating the differentiated expression profile.
#' @param name_col character, the name of the column in \code{DE}, which contains the gene name. If NULL, will use the rownames of \code{DE}.
#' Default is NULL.
#' @param profile_col character, the name of the column in \code{DE}, which will be used as the differentiated expression profile.
#' If DE is created by \code{getDE.limma.2G} or \code{getDE.BID.2G}, this parameter could be 'logFC' or 't'.
#' @param profile_trend character, the choice of how to display the profile, from high/positive to low/negative ('pos2neg')
#' or low/negative to high/positive ('neg2pos').Default is 'pos2neg'.
#' @param use_gs2gene a list for geneset to genes, the name for the list is the gene set name and the content in each list is the vector for genes belong to that gene set.
#' This parameter could be obtained by choosing one from \code{all_gs2gene[['CP:KEGG']]}, or merge several categories by \code{merge_gs}.
#' @param sig_gs_list a vector of characters, the name for the top gene sets.
#' @param gs_DA_Z a vector of numeric values, the Z statistics of differentiated activity (DA) for the sig_gs_list.
#' Better to give name to the vector, otherwise will automatically use sig_gs_list as the name.
#' @param top_gs_number numeric, number for the top significant gene sets to be displayed on the plot. Default is 30.
#' @param target_col character, choice of color pattern used to display the targets. Two options,'black' or 'RdBu'.
#' If set to 'black', the lines will be colored in black.
#' If set to 'RdBu', the lines will be colored into Red to Blue color bar. The color for the annotated genes is set according
#' to its value in the differentiated expression profile, with significant high set for red and low for blue.
#' The significant threshold is set by \code{profile_sig_thre}.
#' Default is 'RdBu'.
#' @param left_annotation character, annotation displayed on the left of the figure representing left condition of the rank_profile. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure representing right condition of the rank_profile. Default is "".
#' @param main character, title for the plot. Default is "".
#' @param profile_sig_thre numeric, threshold for the absolute values in profile to be treated as significance.
#' annotated genes without signifcant values in the profile will be colored in grey. Only work when \code{target_col_type} is set as "DE" and \code{target_col} is set as "RdBu".
#' Default is 0.
#' @param Z_sig_thre numeric, threshold for the Z statistics in \code{driver_DA_Z} and \code{driver_DE_Z} to be treated as signifcance.
#' Only signifcant values will have background color. Default is 1.64.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return logical value indicating whether the plot has been successfully generated

#' @examples
#' \dontrun{
#' db.preload(use_level='transcript',use_spe='human',update=FALSE)
#'
#' ## get all_gs2gene
#'
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo','gene set/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#'
#' ms_tab <- analysis.par$final_ms_tab
#' comp <- 'G4.Vs.others'
#' DE <- analysis.par$DE[[comp]]
#' analysis.par$out.dir.PLOT <- getwd()
#'
#' ## directory for saving the pdf files
#' exp_mat_gene <- exprs(analysis.par$cal.eset)
#'
#' ## calculate activity for all genesets
#' use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,
#'                        use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
#' ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,cal_mat = exp_mat_gene)
#'
#' ## get DA for the gene set
#' phe_info <- pData(analysis.par$cal.eset)
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
#'                     sig_gs_list = sig_gs$names,
#'                     gs_DA_Z=DA_gs[sig_gs$names,'Z-statistics'],
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
    DE <- cbind(DE[,setdiff(colnames(DE),'name')],name=rownames(DE),stringsAsFactors=FALSE)
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
  w1 <- which(is.na(DE[,profile_col])==FALSE)
  DE <- DE[w1,]
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
  rr1 <- ceiling(gswidth/2)
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
  print(matrix(c(rep(0,length.out=rr1),rep(1,length.out=ratio2),rep(0,length.out=rr2*1),
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
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,n_gene),xaxt='n',yaxt='n')
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,1),xaxt='n',yaxt='n',bty='n')
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
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

#' GSEA (gene set enrichment analysis) plot for the Synergy Inference by data-driven Network-based Bayesian Analysis (SINBA) analysis results.
#'
#' \code{draw.GSEA.NetBID.SINBA} will generate a GSEA plot for Synergy Inference by data-driven Network-based Bayesian Analysis (SINBA)
#' analysis results. SINBA calculates the synergistic effect between a seed driver and a partner driver.
#' The plot includes the GSEA plot for the seed driver and the partner driver independently and
#' the GSEA plot for the combination for the seed driver to each partner driver.
#' The statistics on the plot include the differentiated expression (DE), differentiated activity (DA) for each driver,
#' and the different Z (deltaZ) showing the difference between the combination of the seed and the partner driver to the sum of the original Z statistics.
#'
#' This is a plot function to draw GSEA for synergistic effect prediction between the seed driver and a list of partner drivers.
#' User need to input the differentiated expression information, and choose to display the target genes in one row or two rows,
#' by selecting black color or red to blue color bar.
#'
#' @param DE data.frame,the differentiated expression results.
#' This data.frame could be generated by using \code{getDE.limma.2G} or \code{getDE.BID.2G}.
#' If user want to generate this data.frame by other strategies, the rownames must be the gene names or need one column to be the gene name
#' (set in \code{name_col}) and must contain the columns indicating the differentiated expression profile.
#' @param name_col character, the name of the column in \code{DE}, which contains the gene name. If NULL, will use the rownames of \code{DE}.
#' Default is NULL.
#' @param profile_col character, the name of the column in \code{DE}, which will be used as the differentiated expression profile.
#' If DE is created by \code{getDE.limma.2G} or \code{getDE.BID.2G}, this parameter could be 'logFC' or 't'.
#' @param profile_trend character, the choice of how to display the profile, from high/positive to low/negative ('pos2neg')
#' or low/negative to high/positive ('neg2pos').Default is 'pos2neg'.
#' @param seed_driver character, name for the seed driver.
#' @param partner_driver_list a vector of characters, name for the partner driver list.
#' @param seed_driver_label character, label for the seed driver displayed on the plot. Default is seed_driver.
#' @param partner_driver_label a vector of characters, label for the partner driver list displayed on the plot. Default is partner_driver_list
#' @param driver_DA_Z a vector of numeric values, the Z statistics of differentiated activity (DA) for the driver list.
#' Better to give name to the vector, otherwise will automatically use driver list (seed + partner) as the name.
#' @param driver_DE_Z a vector of numeric values, the Z statistics of differentiated expression (DE) for the driver list.
#' Better to give name to the vector, otherwise will automatically use driver list (seed + partner) as the name.
#' @param target_list a list for the target gene information for the drivers. The names for the list must contain the driver in driver list (seed + partner)
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param DA_Z_merge a vector of numeric values, the Z statistics of differentiated activity (DA) for the combination of the seed driver to partner drivers saperately.
#' Better to give name to the vector, otherwise will automatically use partner driver list as the name.
#' @param target_list_merge a list for the target gene information for thecombination of the seed driver to partner drivers saperately.
#' The names for the list must contain the driver in partner_driver_list
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list_merge} could be automatically generated by \code{merge_target_list}.
#' @param top_driver_number numeric, number for the top significant partner drivers to be displayed on the plot. Default is 10.
#' @param top_order character, choice of order pattern used to display the partner drivers. Two options,'merge' or 'diff'.
#' 'merge' means the partner drivers will be sorted by the combined Z statistics.
#' 'diff' means the partner drivers will be sorted by the delta Z statistics.
#' Default is 'merge'.
#' @param target_nrow numeric, number of rows for each driver display on the plot. Two options, 1 or 2.
#' If set to 1, the target genes' position on the profile will be displayed in one row.
#' If set to 2, the target genes' position on the profile will be displayed in two rows,
#' with positive regulated genes displayed on the first row and negative regulated genes displayed on the second row.
#' Default is 2.
#' @param target_col character, choice of color pattern used to display the targets. Two options,'black' or 'RdBu'.
#' If set to 'black', the lines will be colored in black.
#' If set to 'RdBu', the lines will be colored into Red to Blue color bar.
#' If \code{target_col_type} is set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If \code{target_col_type} is set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' with significant high set for red and low for blue. The significant threshold is set by \code{profile_sig_thre}.
#' Default is 'RdBu'.
#' @param target_col_type character, choice of the pattern used to display the color for target genes, only work when \code{target_col} is set as 'RdBu'.
#' Two options,'PN' or 'DE'.
#' If set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' Default is 'PN'.
#' @param left_annotation character, annotation displayed on the left of the figure representing left condition of the rank_profile. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure representing right condition of the rank_profile. Default is "".
#' @param main character, title for the plot. Default is "".
#' @param profile_sig_thre numeric, threshold for the absolute values in profile to be treated as significance.
#' Target genes without signifcant values in the profile will be colored in grey. Only work when \code{target_col_type} is set as "DE" and \code{target_col} is set as "RdBu".
#' Default is 0.
#' @param Z_sig_thre numeric, threshold for the Z statistics in \code{driver_DA_Z} and \code{driver_DE_Z} to be treated as signifcance.
#' Only signifcant values will have background color. Default is 1.64.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return logical value indicating whether the plot has been successfully generated

#' @examples
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
#' ## choose seed driver and partner driver list
#' seed_driver <- driver_list[1]
#' part_driver <- ms_tab$originalID_label
#' ## get merge target
#' merge_target <- lapply(part_driver,function(x){
#'   m1 <- merge_target_list(driver1=seed_driver,driver2=x,
#'                           target_list=analysis.par$merge.network$target_list)
#' })
#' names(merge_target) <- part_driver
#' ## get activity matrix for the merge target network
#' ac_combine_mat <- cal.Activity(all_target=merge_target,
#'                                cal_mat=exprs(analysis.par$cal.eset),
#'                                es.method='weightedmean')
#' ## get DA for the combined drivers
#' comp_name <- 'G4.Vs.others'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')]
#' # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')]
#' # get sample list for G1
#' DA_driver_combine <- getDE.limma.2G(eset=generate.eset(ac_combine_mat),
#'                                     G1=G1,G0=G0,
#'                                     G1_name='G4',G0_name='others')
#' ## or use: DA_driver_combine <- getDE.BID.2G(eset=generate.eset(ac_combine_mat),
#'                                     G1=G1,G0=G0,
#'                                     G1_name='G4',G0_name='others')
#' ## prepare for SINBA input
#' ori_part_Z <- analysis.par$DA[[comp_name]][part_driver,'Z-statistics']
#' ori_seed_Z <- analysis.par$DA[[comp_name]][seed_driver,'Z-statistics']
#' DE <- analysis.par$DE[[comp_name]]
#' driver_DA_Z <- analysis.par$DA[[comp_name]][,'Z-statistics']
#' names(driver_DA_Z) <- rownames(analysis.par$DA[[comp_name]])
#' driver_DE_Z <- analysis.par$DE[[comp_name]][,'Z-statistics']
#' names(driver_DE_Z) <- rownames(analysis.par$DE[[comp_name]])
#' DA_Z_merge <- DA_driver_combine[,'Z-statistics']
#' names(DA_Z_merge) <- rownames(DA_driver_combine)
#' target_list_merge <- merge_target
#' seed_driver_label <- ms_tab[seed_driver,'gene_label']
#' partner_driver_list <- part_driver
#' profile_col <- 't'
#' partner_driver_label <- ms_tab[partner_driver_list,'gene_label']
#' target_list <- analysis.par$merge.network$target_list
##
#' draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,
#'                        seed_driver=seed_driver,
#'                        partner_driver_list=partner_driver_list,
#'                        seed_driver_label=seed_driver_label,
#'                        partner_driver_label=partner_driver_label,
#'                        driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,
#'                        target_list=target_list,
#'                        DA_Z_merge=DA_Z_merge,
#'                        target_list_merge=target_list_merge,
#'                        top_driver_number=20,profile_trend='pos2neg',
#'                        top_order='merge',Z_sig_thre = 1.64,
#'                        target_nrow=1,target_col='RdBu',target_col_type='PN',
#'                        pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo1.pdf',
#'                        analysis.par$out.dir.PLOT))
#'}
#' @export
draw.GSEA.NetBID.SINBA <- function(DE=NULL,name_col=NULL,profile_col=NULL,profile_trend='pos2neg',
                                   seed_driver=NULL,partner_driver_list=NULL,
                                   seed_driver_label=seed_driver,partner_driver_label=partner_driver_list,
                                   driver_DA_Z=NULL,driver_DE_Z=NULL,target_list=NULL,
                                   DA_Z_merge=NULL,target_list_merge=NULL,
                                   top_driver_number=10,top_order='merge',target_nrow=2,target_col='RdBu',target_col_type='PN',
                                   left_annotation="",right_annotation="",main="",
                                   profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
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
      #partner_driver_list <- partner_driver_list[order(abs(DA_Z_merge[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)][1:top_driver_number] ## only consider positive part
    }
    if(profile_trend=='pos2neg')
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = FALSE)]
    else
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)]
  }else{
    if(length(partner_driver_list)>top_driver_number){
      #partner_driver_list <- partner_driver_list[order(abs(diff_Z[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = TRUE)][1:top_driver_number] ## only consider positive increase
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
  if(is.null(name_col)==TRUE){
    DE <- cbind(DE[,setdiff(colnames(DE),'name')],name=rownames(DE),stringsAsFactors=FALSE)
    name_col <- 'name'
  }
  w1 <- which(is.na(DE[,profile_col])==FALSE)
  DE <- DE[w1,]
  if(profile_trend=='pos2neg') DE <- DE[order(DE[,profile_col],decreasing = TRUE),] else DE <- DE[order(DE[,profile_col],decreasing = FALSE),]
  DE_profile <- DE[,profile_col]
  DE_profile_name <- DE[,name_col]
  ##############################################
  ## calculate layout

  target_list <- lapply(target_list,function(x)x[which(x$target %in% DE_profile_name),])
  target_list_merge <- lapply(target_list_merge,function(x)x[which(x$target %in% DE_profile_name),])

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
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,n_gene),xaxt='n',yaxt='n')

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow

  rect(xleft = pp[1],xright=pp[2],ybottom = yy1[length(yy1)-target_nrow],ytop=yy1[length(yy1)],border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.3)) ## for seed rows
  rect(xleft = pp[1],xright=pp[2],ybottom = yy2,ytop=yy4,border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.2)) ## for combine rows

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
    #names(cc) <- names(DE_profile)
    cc[which(cc=='white')] <- 'light grey'
  }
  if(target_nrow==1){
    # for seed driver
    t1 <- use_target_list[[seed_driver]]
    w0 <- which(DE_profile_name %in% t1$target)
    w1 <- w0*rr+pp[1]
    if(target_col=='black'){
      segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
    }else{
      if(target_col_type=='DE'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=cc[w0])
      }else{
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      w0 <- which(DE_profile_name %in% t1$target)
      w1 <- w0*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
        }
      }
    }
    # for combine
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      w0 <- which(DE_profile_name %in% t1$target)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(DE_profile_name %in% t_over)
      w1_over <- w0_over*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
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
    w0 <- which(DE_profile_name %in% t11)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=pos_col,lwd=1.5)
        }
      }
    }
    w0 <- which(DE_profile_name %in% t12)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=neg_col,lwd=1.5)
        }
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(DE_profile_name %in% t11)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=pos_col,lwd=1.5)
          }
        }
      }
      w0 <- which(DE_profile_name %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=neg_col,lwd=1.5)
          }
        }
      }
    }
    # for each partner driver + seed combination
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(DE_profile_name %in% t11)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(DE_profile_name %in% t_over) ## setdiff(t_over,names(DE_profile)) !!!
      w1_over <- w0_over*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=pos_col,lwd=1.5)
          }
        }
      }
      points(intersect(w1_over,w1),rep((yy11[4*i-1]+yy1[4*i-1])/2,length.out=length(intersect(w1_over,w1))),pch='*',col='black')
      w0 <- which(DE_profile_name %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=neg_col,lwd=1.5)
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,3),xaxt='n',yaxt='n',bty='n')
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
  plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')

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
         ybottom=yy1[i],ytop=yy1[i]+dyy22/2*0.35,col=pos_col,border=NA)
    rect(xleft=xleft,xright=xleft+target_size[1,2]/mm*tt,
         ytop=yy1[i],ybottom=yy1[i]-dyy22/2*0.35,col=neg_col,border=NA)
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
           ybottom=yy33[2*i],ytop=yy33[2*i]+dyy33*0.35,col=pos_col,border=NA)
      rect(xleft=xleft,xright=xleft+target_size[i+1,2]/mm*tt,
           ytop=yy33[2*i],ybottom=yy33[2*i]-dyy33*0.35,col=neg_col,border=NA)
      # merge
      rect(xleft=xleft,xright=xleft+merge_target_size[i,1]/mm*tt,
           ybottom=yy33[2*i-1],ytop=yy33[2*i-1]+dyy33*0.35,col=pos_col,border=NA)
      rect(xleft=xleft,xright=xleft+merge_target_size[i,2]/mm*tt,
           ytop=yy33[2*i-1],ybottom=yy33[2*i-1]-dyy33*0.35,col=neg_col,border=NA)
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
                                   target1=intersect(use_target_list[[seed_driver]]$target,DE_profile_name),
                                   target2=intersect(use_target_list[[partner_driver_list[i]]]$target,DE_profile_name),
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

#' Merge target list for two drivers.
#'
#' \code{merge_target_list} is a function to merge the target list for two drivers.
#' Higher MI statistics for the shared target genes by the two drivers will be kept in the final target list.
#'
#' @param driver1 character, the name for the first driver to merge.
#' @param driver2 character, the name for the second driver to merge.
#' @param target_list a list for the target gene information for the drivers. The names for the list must contain the driver1 and driver2.
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @return a list for the target gene information.
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
#' Box plot and stripchart for the gene's expression levels and driver's activity values in samples with different categories.
#'
#' \code{draw.categoryValue} will draw the box plot with stripchart for the gene's expression levels and/or driver's activity values
#'  in samples with different phenotype categories.
#'
#' This is a function to draw the gene's expression level and driver's activity values at the same time in one plot
#' across samples with different phenotype categories. Also, only draw of expression level or activity level is accepted.
#'
#' @param ac_val a vector of numeric values, the activity level for the interested driver across all samples.
#' @param exp_val a vector of numeric values, the expression level for the interested gene across all samples.
#' @param use_obs_class a vector of characters, the cateogory class for all samples.
#' The order of samples in \code{use_obs_class} must be the same with \code{ac_val} or \code{exp_val}.
#' This vector could be generated by the function \code{get_obs_label} to extract this vector from the dataframe of \code{pData(eset)}
#' by selecting the column name.
#' @param class_order a vector of characters, the order of category class displayed on the figure.
#' If NULL, will use the alphabetical order of the category class name. Default is NULL.
#' @param category_color a vector of characters, each item is the color for the class in \code{class_order}.
#' If NULL, will automatically use function \code{get.class.color} generate the color bar. Default is NULL.
#' @param stripchart_color character, the color for the stripchart. Default is 'black' with transparent alpha set at 0.7.
#' @param strip_cex numeric, \code{cex} for points on the plot. Default is 1.
#' @param class_srt numeric, the displayed category class label rotation in degrees. Default is 90.
#' @param class_cex numeric, \code{cex} for the category class label displayed on the plot. Default is 1.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#' @param main_ac character, main for the sub-plot of activity level. Default is "".
#' @param main_exp character,main for the sub-plot of expression level. Default is "".
#' @param main_cex numeric, \code{cex} for the main title displayed on the plot. Default is 1.
#'
#' @return
#' Will return logical value indicating whether the plot has been successfully generated

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
#' exp_mat <- exprs(analysis.par$cal.eset)
#' ## expression,the rownames could match originalID
#' ac_mat  <- exprs(analysis.par$merge.ac.eset)
#' ## activity,the rownames could match originalID_label
#' phe_info <- pData(analysis.par$cal.eset)
#' use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
#' draw.categoryValue(ac_val=ac_mat[use_driver,],
#'                    exp_val=exp_mat[ms_tab[use_driver,'originalID'],],
#'                    use_obs_class=use_obs_class,
#'                    class_order=c('WNT','SHH','G4'),
#'                    class_srt=30,
#'                    main_ac = ms_tab[use_driver,'gene_label'],
#'                    main_exp=ms_tab[use_driver,'geneSymbol'])
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
#' exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames could match originalID
#' ac_mat  <- exprs(analysis.par$merge.ac.eset) ## activity,the rownames could match originalID_label
#' phe_info <- pData(analysis.par$cal.eset)
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
                               main_ac="",main_exp="",main_cex=1){
  if(is.null(class_order)){
    class_order <- sort(unique(use_obs_class))
  }
  if(is.null(category_color)==TRUE){
    class_col <- get.class.color(class_order)
    class_col1 <- get_transparent(class_col,0.5)
  }else{
    class_col1 <- category_color
  }
  c1 <- 0
  if(is.null(ac_val)==FALSE){c1 <- c1+1}
  if(is.null(exp_val)==FALSE){c1 <- c1+1}
  labelWidth <- max(strwidth(class_order,'inches',cex=class_cex)*sin(class_srt*pi/180))
  if(is.null(pdf_file)==FALSE){
    hh <- 5+1.5+labelWidth
    if(c1==1) pdf(pdf_file,width=1.5+3,height=hh)
    if(c1==2) pdf(pdf_file,width=1.5+3*2,height=hh)
  }
  if(c1>1) layout(t(matrix(1:c1)))
  par(mai=c(labelWidth+0.5,1,1,0.5))
  if(is.null(ac_val)==FALSE){
    ddf <- data.frame(data=ac_val,class=factor(use_obs_class,levels=class_order))
    a <- boxplot(data~class,data=ddf,ylab='Activity Value',col=class_col1,outline=FALSE,border='dark grey',cex.lab=1.2,names=NA,bty='n',
                 ylim=c(min(ddf$data),max(ddf$data)),main=main_ac,cex.main=main_cex)
    text(1:length(class_order),par()$usr[3]-(par()$usr[4]-par()$usr[3])/20,adj=0.5+class_srt/180,class_order,srt=class_srt,xpd=TRUE,cex=class_cex)
    stripchart(data~class,data=ddf,add=TRUE,pch=16,method='jitter',vertical=TRUE,col=stripchart_color,cex=strip_cex)
  }
  if(is.null(exp_val)==FALSE){
    ddf <- data.frame(data=exp_val,class=factor(use_obs_class,levels=class_order))
    a <- boxplot(data~class,data=ddf,col=class_col1,ylab='Expression Value',outline=FALSE,border='dark grey',cex.lab=1.2,names=NA,bty='n',
                 ylim=c(min(ddf$data),max(ddf$data)),main=main_exp,cex.main=main_cex)
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
    x3 <- unlist(strsplit(as.character(x2),""))
    x4 <- length(x3)%/%3 ## add number
    if(x4>0){
      pp <- length(x3)-seq(1,x4)*3; x3[pp] <- paste0(x3[pp],','); paste(x3,collapse="")
    }else{
      x2
    }
  })
  unlist(x1)
}

#' Target network structure plot for the driver.
#'
#' \code{draw.targetNet} will draw the network structure for the selected driver and its target genes.
#'
#' This is a function to draw target network structure for the selected driver.
#' The color bar represents the positive (red) or negative (blue) regulation with line width showing the strength.
#'
#' @param source_label the label for the driver displayed on the plot.
#' @param source_z numeric, the Z statistic for the driver, used to color the driver point.
#' If NULL, the driver will be colored in grey. Default is NULL.
#' @param edge_score a vector of numeric values, indicating the correlation between the driver and the target genes.
#' The value ranges from -1 to 1, with positive value indicating postivie regulation and negative value indicating negative correlation.
#' The names for the vector is the gene labels displayed on the plot.
#' @param label_cex numeric, \code{cex} for the target genes displayed on the plot. Default is 0.7.
#' @param source_cex numeric, \code{cex} for the source genes displayed on the plot. Default is 1.
#' @param arrow_direction character, choose from 'in' or 'out'. Default is 'out'.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return Will return logical value indicating whether the plot has been successfully generated
#'
#' @examples
#' source_label <- 'test1'
#' source_z <- 1.96
#' edge_score <- (sample(1:200,size=100,replace=TRUE)-100)/100
#' names(edge_score) <- paste0('G',1:100)
#' draw.targetNet(source_label=source_label,source_z=source_z,
#'                edge_score=edge_score)
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
                           pdf_file=NULL,arrow_direction='out'){
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  edge_score<- sort(edge_score)
  tmp1 <- sapply(unique(names(edge_score)),function(x){
    x1 <- edge_score[which(names(edge_score)==x)]
    x1[which.max(abs(x1))]
  })
  names(tmp1) <- unique(names(edge_score))
  edge_score <- tmp1
  edge_score<- sort(edge_score)
  g1 <- names(edge_score)
  ec <- z2col(edge_score*100,sig_thre=0,n_len=length(edge_score),red_col=pos_col,blue_col=neg_col);names(ec) <- names(edge_score)
  ec <- get_transparent(ec,alpha=0.8)
  ew <- 2*label_cex*(abs(edge_score)-min(abs(edge_score)))/(max(abs(edge_score))-min(abs(edge_score)))+label_cex/2; names(ew) <- names(edge_score)
  t2xy <- function(tt,radius=1) {
    t2p <- pi*2 * tt
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  geneWidth <- max(strwidth(g1,'inches',cex=label_cex))
  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=6+2*geneWidth,height=6+2*geneWidth)
  par(mai=c(1,1,1,1))
  plot(1,xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white',xlab="",ylab="",xaxt='n',yaxt='n')
  pp <- par()$usr
  tt <- seq(-0.5,0.5,length.out=length(g1)+1)[-1];
  p1<-t2xy(tt,radius=0.8);
  for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g1[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
  geneWidth <- strwidth(source_label,'user',cex=source_cex)
  if(arrow_direction=='out'){
    p2<-t2xy(tt,radius=0.8-label_cex/36);
    p3<-t2xy(tt,radius=0.8-label_cex/48);
    arrows(x0=0,y0=0,x1=p2$x,y1=p2$y,col=ec,lwd=1,angle=10,length=0.1*label_cex);
  }else{
    p2<-t2xy(tt,radius=0.8-label_cex/36);
    p3<-t2xy(tt,radius=0.8-label_cex/36);
    p4<-t2xy(tt,radius=geneWidth/2);
    arrows(x0=p2$x,y0=p2$y,x1=p4$x,y1=p4$y,col=ec,lwd=1,angle=5,length=0.1*label_cex);
  }
  points(p3$x,p3$y,pch=16,col='dark grey',cex=label_cex)
  if(is.null(source_z)==TRUE){
    #points(0,0,col='light grey',cex=geneWidth*36,pch=16)
    draw.ellipse(0,0,a=geneWidth/2,b=geneWidth/2,col='light grey',border=NA)
  }else{
    #points(0,0,col=z2col(source_z),cex=geneWidth*36,pch=16)
    draw.ellipse(0,0,a=geneWidth/2,b=geneWidth/2,col=z2col(source_z),border=NA)
  }
  #points(0,0,col='light grey',cex=14,pch=16)
  text(0,0,source_label,adj=0.5,xpd=TRUE,cex=source_cex)
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}

#' Target network structure plot for two drivers.
#'
#' \code{draw.targetNet.TWO} will draw the network structure for the selected two drivers and their target genes.
#'
#' This is a function to draw target network structure for the selected two drivers.
#' The color bar represents the positive (red) or negative (blue) regulation with line width showing the strength.
#'
#' @param source1_label the label for the first(left) driver displayed on the plot.
#' @param source2_label the label for the second(right) driver displayed on the plot.
#' @param source1_z numeric, the Z statistic for the first driver, used to color the driver point.
#' If NULL, the driver will be colored in grey. Default is NULL.
#' @param source2_z numeric, the Z statistic for the second driver, used to color the driver point.
#' If NULL, the driver will be colored in grey. Default is NULL.
#' @param edge_score1 a vector of numeric values, indicating the correlation between the first driver and the target genes.
#' The value ranges from -1 to 1, with positive value indicating postivie regulation and negative value indicating negative correlation.
#' The names for the vector is the gene labels displayed on the plot.
#' @param edge_score2 a vector of numeric values, indicating the correlation between the second driver and the target genes.
#' Similar with \code{edge_score1}
#' @param arrow1_direction character, the arrow direction for source1, choose from 'in' or 'out'. Default is 'out'.
#' @param arrow2_direction character, the arrow direction for source2, choose from 'in' or 'out'. Default is 'out'.
#' @param label_cex numeric, \code{cex} for the target genes displayed on the plot. Default is 0.7.
#' @param source_cex numeric, \code{cex} for the source genes displayed on the plot. Default is 1.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#' @param total_possible_target numeric or a vector of characters. If input numeric, will be the total number of possible targets.
#' If input a vector of characters, will be the background list of all possible target genes.
#' This parameter will be passed to function \code{test.targetNet.overlap} to test whether the target genes of the two drivers are significantly intersected.
#' If NULL, will do not perform this test. Default is NULL.
#' @param show_test logical, indicating whether the testing results will be printed and returned. Default is FALSE.
#'
#' @return if \code{show_test}==FALSE, will return logical value indicating whether the plot has been successfully generated, otherwise will return the statistics of testing.
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
#'
#' \dontrun{
#' source1_label <- 'test1'
#' source1_z <- 1.96
#' edge_score1 <- (sample(1:160,size=100,replace=TRUE)-80)/80
#' names(edge_score1) <- sample(paste0('G',1:1000),size=80)
#' source2_label <- 'test2'
#' source2_z <- -2.36
#' edge_score2 <- (sample(1:240,size=100,replace=TRUE)-120)/120
#' names(edge_score2) <- sample(paste0('G',1:1000),size=120)
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
                               arrow1_direction='out',arrow2_direction='out',
                               label_cex=0.7,source_cex=1,pdf_file=NULL,
                               total_possible_target=NULL,show_test=FALSE){
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
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
  ec1 <- z2col(edge_score1*100,sig_thre=0,n_len=length(edge_score1),red_col=pos_col,blue_col=neg_col);names(ec1) <- names(edge_score1)
  ec2 <- z2col(edge_score2*100,sig_thre=0,n_len=length(edge_score2),red_col=pos_col,blue_col=neg_col);names(ec2) <- names(edge_score2)
  ew1 <- 2*label_cex*(abs(edge_score1)-min(abs(edge_score1)))/(max(abs(edge_score1))-min(abs(edge_score1)))+label_cex/2; names(ew1) <- names(edge_score1)
  ew2 <- 2*label_cex*(abs(edge_score2)-min(abs(edge_score2)))/(max(abs(edge_score2))-min(abs(edge_score2)))+label_cex/2; names(ew2) <- names(edge_score2)
  t2xy <- function(tt,radius=1) {
    t2p <- pi*2 * tt + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  geneWidth <- max(strwidth(g1,'inches',cex=label_cex))

  if(is.null(pdf_file)==FALSE) pdf(pdf_file,width=10+4*geneWidth,height=8+2*geneWidth)
  plot(1,xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white',xlab="",ylab="",xaxt='n',yaxt='n')
  par(mai=c(1,1,1,1))

  geneWidth1 <- strwidth(source1_label,'user',cex=source_cex)
  geneWidth2 <- strwidth(source2_label,'user',cex=source_cex)

  if(length(g1)>0){
    tt <- seq(-0.225,0.225,length.out=length(g1));init.angle <- -180;p1<-t2xy(tt,radius=0.8);
    for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g1[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    p1<-t2xy(tt,radius=0.8-label_cex/36);
    if(arrow1_direction=='out'){
      p2<-t2xy(tt,radius=0.8-label_cex/36);
      p3<-t2xy(tt,radius=0.8-label_cex/48);
      arrows(x0=-0.2,y0=0,x1=p1$x,y1=p1$y,col=ec1[g1],lwd=ew1[g1],angle=10,length=0.1*label_cex);
    }else{
      p2<-t2xy(tt,radius=0.8-label_cex/36);
      p3<-t2xy(tt,radius=0.8-label_cex/36);
      p4<-t2xy(tt,radius=geneWidth1/2);
      arrows(x0=p2$x,y0=p2$y,x1=p4$x-0.2,y1=p4$y,col=ec1[g1],lwd=ew1[g1],angle=5,length=0.1*label_cex);
    }
    points(p3$x,p3$y,pch=16,col='dark grey',cex=label_cex)
  }
  if(length(g2)>0){
    tt <- seq(-0.225,0.225,length.out=length(g2));init.angle <- 0;
    p1<-t2xy(tt,radius=0.8);
    for(i in 1:length(p1$x)) text(p1$x[i],p1$y[i],g2[i],cex=label_cex,srt=180*atan(p1$y[i]/p1$x[i])/pi,adj=ifelse(p1$x[i]>0,0,1),xpd=TRUE)
    p1<-t2xy(tt,radius=0.8-label_cex/36);
    if(arrow2_direction=='out'){
      p2<-t2xy(tt,radius=0.8-label_cex/36);
      p3<-t2xy(tt,radius=0.8-label_cex/48);
      arrows(x0=0.2,y0=0,x1=p1$x,y1=p1$y,col=ec2[g2],lwd=ew2[g2],angle=10,length=0.1*label_cex);
    }else{
      p2<-t2xy(tt,radius=0.8-label_cex/36);
      p3<-t2xy(tt,radius=0.8-label_cex/36);
      p4<-t2xy(tt,radius=geneWidth2/2);
      arrows(x0=p2$x,y0=p2$y,x1=p4$x+0.2,y1=p4$y,col=ec2[g2],lwd=ew2[g2],angle=5,length=0.1*label_cex);
    }
    points(p3$x,p3$y,pch=16,col='dark grey',cex=label_cex)
  }
  if(length(g12)>0){
    tt <- seq(min(0.1*length(g12),0.7),-min(0.1*length(g12),0.7),length.out=length(g12));
    if(arrow1_direction=='out'){
      arrows(x0=-0.2,y0=0,x1=0,y1=tt,col=ec1[g12],lwd=ew1[g12],angle=10,length=0.1*label_cex)
    }else{
      p4<-t2xy(tt,radius=geneWidth1/2);
      arrows(x0=0,y0=tt,x1=-0.2,y1=0,col=ec1[g12],lwd=ew1[g12],angle=5,length=0.1*label_cex);
      #arrows(x0=-0.2,y0=0,x1=0,y1=tt,col=ec1[g12],lwd=ew1[g12],angle=10,length=0.1*label_cex)
    }
    if(arrow2_direction=='out'){
      arrows(x0=0.2,y0=0,x1=0,y1=tt,col=ec2[g12],lwd=ew2[g12],angle=10,length=0.1*label_cex)
    }else{
      p4<-t2xy(tt,radius=geneWidth2/2);
      arrows(x0=0,y0=tt,x1=0.2,y1=0,col=ec1[g12],lwd=ew1[g12],angle=5,length=0.1*label_cex);
      #arrows(x0=0.2,y0=0,x1=0,y1=tt,col=ec2[g12],lwd=ew2[g12],angle=10,length=0.1*label_cex)
    }
    #segments(x0=-0.2,y0=0,x1=0,y1=tt,col=ec1[g12],lwd=ew1[g12])
    #segments(x0=0.2,y0=0,x1=0,y1=tt,col=ec2[g12],lwd=ew2[g12])
    boxtext(0,tt,labels=g12,col.bg=get_transparent('light grey',0.3),cex=label_cex)
  }

  if(is.null(source2_z)==TRUE)
    draw.ellipse(0.2,0,a=geneWidth1/2,b=geneWidth1/2,col='light grey',border=NA)
  else
    draw.ellipse(0.2,0,a=geneWidth1/2,b=geneWidth1/2,col=z2col(source2_z),border=NA)

  text(0.2,0,source2_label,adj=0.5,cex=source_cex)

  if(is.null(source1_z)==TRUE)
    draw.ellipse(-0.2,0,a=geneWidth1/2,b=geneWidth1/2,col='light grey',border=NA)
  else
    draw.ellipse(-0.2,0,a=geneWidth1/2,b=geneWidth1/2,col=z2col(source1_z),border=NA)

  text(-0.2,0,source1_label,adj=0.5,cex=source_cex)
  # fisher test for target
  if(is.null(total_possible_target)==FALSE & show_test==TRUE){
    res <- test.targetNet.overlap(source1_label,source2_label,names(edge_score1),names(edge_score2),total_possible_target)
    if(is.null(pdf_file)==FALSE) dev.off()
    return(res)
  }
  if(is.null(pdf_file)==FALSE) dev.off()
  return(TRUE)
}

#' Test for the target genes' intersection between two drivers.
#'
#' \code{test.targetNet.overlap} will test whether the target genes of two drivers are significantly intersected.
#'
#' This is a function to perform Fisher's Exact Test for the intersection of the target genes from two drivers.
#'
#' @param source1_label character, the label for the first driver.
#' @param source2_label character, the label for the second driver.
#' @param target1 a vector of characters, the list of target genes for the first driver.
#' @param target2 a vector of characters, the list of target genes for the second driver.
#' @param total_possible_target numeric or a vector of characters. If input numeric, will be the total number of possible targets.
#' If input a vector of characters, will be the background list of all possible target genes.
#'
#' @return Return the statistics of testing, including the \code{P.Value}, \code{Odds_Ratio} and \code{Intersected_Nuumber}.
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
    draw.spByGene(net,driver_list=ori_driver_list,target_list=ori_target_list,transfer_tab=transfer_tab,...)
  }
  return(list(all_sp=all_r1,use_sp=net))
}
### plot igraph
draw.spByGene <- function(net,driver_list=NULL,target_list=NULL,transfer_tab=NULL,vertex.size=25,vertex.label.cex=1,
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

#' Get target list by input network information from data.frame.
#'
#' \code{get_net2target_list} is included in the \code{get.SJAracne.network},
#' the reason to make it invokable just for user to read in network files prepared by themselves.
#'
#' @param net_dat data.frame, must contain columns named with "source" and "target",
#' "MI" and "spearman" are strongly suggested but not required.
#' If these two columns missed, some options may be error in some of the functions in the following steps
#' (such as es.method='weightedmean' in \code{cal.Activity}).
#' @return a list for the target gene information for the drivers. Each object in the list is a data.frame to save the target genes.
#' @examples
#' tf.network.file <- sprintf('%s/demo1/network/SJAR/project_2019-02-14/%s/%s',
#'                    system.file(package = "NetBID2"),
#'                    'output_tf_sjaracne_project_2019-02-14_out_.final',
#'                    'consensus_network_ncol_.txt')
#' net_dat      <- read.delim(file=tf.network.file,stringsAsFactors = FALSE)
#' target_list  <- get_net2target_list(net_dat)
#' @export
get_net2target_list <- function(net_dat=NULL) {
  all_source <- unique(net_dat$source)
  all_target <- lapply(all_source, function(x) {
    n1 <- net_dat[which(net_dat$source == x), intersect(c('target', 'MI', 'spearman'),colnames(net_dat))]
    rownames(n1) <- n1$target
    return(n1)
  })
  names(all_target) <- all_source
  return(all_target)
}

#' Generate data structure to save network information by input network file generated by SJAracne
#'
#' \code{get.SJAracne.network} is a function to read in network file generated by SJAracne
#' (consensus_network_ncol_.txt file in the result directory)
#'
#' This function aims to read in network files generated by SJAracne and save the network information into three lists,
#' \code{network_dat} is a data.frame to save all the information in the network file;
#' \code{target_list} is a list for containing the target genes' information for the drivers. The names for the list is the driver name
#' and each object in the list is a data.frame to save the target genes.
#' \code{igraph_obj} is an igraph object to save the network, it is a directed, weighted network.
#' The function will set two edge attributes to the igraph_obj, \code{weight} is the MI values and \code{sign} is the sign for the spearman value
#' to indicate positive regulation (1) or negative regulation (-1).
#'
#' @param network_file character, file path for the network file. Must use the file (consensus_network_ncol_.txt) in the result directory.
#' The directory may like:
#' <SJAR_main_dir>/<project_name>/output_[tf|sig]_sjaracne_<project_name>_out_.final/
#' <SJAR_main_dir>/SJARACNE_<project_name>/SJARACNE_out.final
#' consensus_network_ncol_.txt
#'
#' @return This function will return a list containing three items, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#'
#' @examples
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             prject_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
#'
#' \dontrun{
#' }
#' @export
get.SJAracne.network <- function(network_file=NULL){
  if(is.null(network_file)){
    message('No input network file, please check and re-try !');return(FALSE)
  }
  net_dat      <- read.delim(file=network_file,stringsAsFactors = FALSE)
  target_list  <- get_net2target_list(net_dat)
  igraph_obj   <- graph_from_data_frame(net_dat[,c('source','target')],directed=TRUE) ## add edge weight ???
  igraph_obj   <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  igraph_obj   <- set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}

#' Update the network information in the structured network list dataset
#'
#' \code{update_SJAracne.network} is a function to update the network information by input threshold for statistics or input gene list for use.
#'
#' This function aims to update the network list dataset generated by \code{get.SJAracne.network}
#' and return the list dataset passed the filtration with the same data structure.
#'
#' @param network_list list,the network list dataset generated by \code{get.SJAracne.network},
#' contains \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#' @param all_possible_drivers a vector of characters,all possible driver list used to filter the network file.
#' If NULL, will set to the possible drivers from \code{network_list}. Default is NULL.
#' @param all_possible_targets a vector of characters,all possible target list used to filter the network file.
#' If NULL, will set to the possible targets \code{network_list}. Default is NULL.
#' @param force_all_drivers logical, whether or not to include all genes in the \code{all_possible_drivers} in the final network.
#' For \code{network_dat} and \code{target_list}, if \code{force_all_drivers} is set to TRUE, genes in \code{all_possible_drivers}
#' will not be filtered by the following statistical thresholds.
#' For \code{igraph_obj}, if \code{force_all_drivers} is set to TRUE, all genes in \code{all_possible_drivers},
#' even not exist in the original network file will be include in the vertice of the igraph object.
#' Default is TRUE.
#' @param force_all_targets logical, whether or not to include all genes in the \code{all_possible_targets} in the final network.
#' For \code{network_dat} and \code{target_list}, if \code{all_possible_targets} is set to TRUE, genes in \code{all_possible_targets}
#' will not be filtered by the following statistical thresholds.
#' For \code{igraph_obj}, if \code{force_all_targets} is set to TRUE, all genes in \code{all_possible_targets},
#' even not exist in the original network file will be include in the vertice of the igraph object.
#' Default is TRUE.
#' @param min_MI numeric, minimum threshold for MI. Default is 0.
#' @param max_p.value numeric, maximum threshold for p.value. Default is 1.
#' @param min_spearman_value numeric, minimum threshold for spearman absolute value. Default is 0.
#' @param min_pearson_value numeric, minimum threshold for pearson absolute value. Default is 0.
#' @param spearman_sign_use a vector of numeric value, 1 indicates positve values in spearman will be used. -1 indicates negative values will be used.
#' If only want to include positive values, set \code{spearman_sign_use} to 1.
#' Default is c(1,-1).
#' @param pearson_sign_use a vector of numeric value, 1 indicates positve values in pearson will be used. -1 indicates negative values will be used.
#' If only want to include positive values, set \code{pearson_sign_use} to 1.
#' Default is c(1,-1).
#' @param directed logical, whether the network in igraph is directed or not. Default is TRUE.
#' @param weighted logical, whether to add the edge weight in the igraph object. Default is TRUE.
#'
#' @return This function will return a list containing three items, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#'
#' @examples
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             prject_name=project_name,
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
#' print(intersect(c('addition_driver_1','addition_driver_2'),
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
  n1 <- names(network_list)
  n2 <- c('network_dat','target_list','igraph_obj')
  if(length(setdiff(n2,n1))>0){
    message(sprintf('%s not included in the network_list, please chech and re-try !',paste(setdiff(n2,n1),collapse=';')));
    return(FALSE)
  }
  ori_net_dat <- network_list$network_dat
  net_dat <- ori_net_dat
  if(is.null(all_possible_drivers)==TRUE) all_possible_drivers <- unique(net_dat$source)
  if(is.null(all_possible_targets)==TRUE) all_possible_targets <- unique(net_dat$target)
  # basic statistics filter
  w1 <- which(net_dat$MI>=min_MI & net_dat$p.value<=max_p.value
              & abs(net_dat$pearson)>=min_pearson_value &
                abs(net_dat$spearman)>=min_spearman_value)
  all_w1 <- w1 ## use rows
  # sign choose
  if(!1 %in% spearman_sign_use){ ## do not use the positive ones
    w1 <- which(net_dat$spearman<=0)
    all_w1 <- setdiff(all_w1,w1)
  }
  if(!-1 %in% spearman_sign_use){ ## do not use the negative ones
    w1 <- which(net_dat$spearman>=0)
    all_w1 <- setdiff(all_w1,w1)
  }
  if(!1 %in% pearson_sign_use){ ## do not use the positive ones
    w1 <- which(net_dat$pearson<=0)
    all_w1 <- setdiff(all_w1,w1)
  }
  if(!-1 %in% pearson_sign_use){ ## do not use the negative ones
    w1 <- which(net_dat$pearson>=0)
    all_w1 <- setdiff(all_w1,w1)
  }
  # nodes filter
  all_possible_nodes <- unique(c(unique(net_dat$source),unique(net_dat$target))) ## original all possible
  if(force_all_drivers==TRUE){
    w1 <- which(net_dat$source %in% all_possible_drivers)
    all_w1 <- unique(c(all_w1,w1))
    all_possible_nodes <- unique(c(all_possible_nodes,all_possible_drivers))
  }
  if(force_all_targets==TRUE){
    w1 <- which(net_dat$target %in% all_possible_targets)
    all_w1 <- unique(c(all_w1,w1))
    all_possible_nodes <- unique(c(all_possible_nodes,all_possible_targets))
  }
  message(sprintf('%d from %d edges are kept in the network !',length(all_w1),nrow(ori_net_dat)))
  message(sprintf('%d nodes will be used to generate the igraph!',length(all_possible_nodes)))
  # keep all genes in all* to be in the igraph
  net_dat <- net_dat[which(net_dat$source %in% all_possible_drivers & net_dat$target %in% all_possible_targets),] ## filter by all
  igraph_obj   <- graph_from_data_frame(net_dat[,c('source','target')],directed=directed,vertices = all_possible_nodes) ##
  if(weighted==TRUE) igraph_obj <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  target_list  <- get_net2target_list(net_dat)
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}


#' Generate QC plot for the network object
#'
#' \code{draw.network.QC} is a function to draw the QC plot for the network object,
#' mainly the density plot for the degree and the check for scale-free feature.
#'
#' @param igraph_obj igraph object, this could be generated by \code{get.SJAracne.network}
#' @param outdir character, the output directory to save the QC figures.
#' @param prefix character, the prefix for the QC figure name.Default is "".
#'
#' @return logical value indicating whether the plot has been successfully generated
#'
#' @examples
#' \dontrun{
#' if(exists('analysis.par')==TRUE) rm(analysis.par)
#' network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo
#' network.project.name <- 'project_2019-02-14' # demo project name
#' project_main_dir <- 'test/'
#' project_name <- 'test_driver'
#' analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir,
#'                                             prject_name=project_name,
#'                                             network_dir=network.dir,
#'                                             network_project_name=network.project.name)
#' analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
#' analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
#' draw.network.QC(analysis.par$tf.network$igraph_obj,
#'                 outdir=analysis.par$out.dir.QC,prefix='TF_net_')
#' }
#' @export
draw.network.QC <- function(igraph_obj,outdir=NULL,prefix=""){
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
  return(TRUE)
}

## functions to check the scale free feature of the network
check_scalefree <- function(igraph_obj) {
  gr1 <- igraph_obj
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

#' Inner use: prepare for MICA
#' @param mat matrix
#' @param outdir character
#' @param prjname character
#' @param all_k a vector of integers
#' @param retransformation character
#' @param perplexity numeric
#' @export
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

#' Inner use: prepare for MICA with file
#' @param input_exp character
#' @param outdir character
#' @param prjname character
#' @param all_k a vector of integers
#' @param retransformation character
#' @param perplexity numeric
#' @export
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


#' Prepare for running SJaracne
#'
#' \code{SJAracne.prepare} is a function to prepare files for running SJAracne.
#'
#' Detailed description to run SJAracne could be found in Github.
#' Check \url{https://github.com/jyyulab/SJARACNe/} for detail.
#'
#' @param eset an ExpressionSet class object, which contains the expression matrix.
#' @param use.samples a vector of characters, the sample list used for running SJAracne.
#' @param TF_list a vector of characters, the TF list used for analysis.
#' @param SIG_list a vector of characters, the SIG list used for analysis.
#' @param SJAR.main_dir character, the file path to save the results for SJAracne
#' @param SJAR.project_name character, the name of the project used to label the output directory.
#' @param IQR.thre numeric, threshold for IQR filter for all non-driver genes.
#' @param IQR.loose_thre numeric, threshold for IQR filter for all driver(TF/SIG) genes.
#' @param add_options additional options used to run sjaracne.
#' @examples
#' \dontrun{
#' network.par <- list()
#' network.par$out.dir.DATA <- system.file('demo1','network/DATA/',package = "NetBID2")
#' NetBID.loadRData(network.par=network.par,step='exp-QC')
#' db.preload(use_level='gene',use_spe='human',update=FALSE)
#' use_gene_type <- 'external_gene_name' ## this should user-defined !!!
#' use_genes <- rownames(fData(network.par$net.eset))
#' use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
#' #select sample for analysis
#' phe <- pData(network.par$net.eset)
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
  function(eset,use.samples = rownames(pData(eset)),
           TF_list=NULL,SIG_list=NULL,
           SJAR.main_dir='',
           SJAR.project_name = "",
           IQR.thre=0.5,IQR.loose_thre=0.1,add_options='') {
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
    choose1 <- IQR.filter(d, rownames(d),thre = IQR.thre,loose_thre=IQR.loose_thre,loose_gene=unique(c(TF_list,SIG_list)))
    d <- d[choose1, ]
    use.genes <- rownames(d)
    use.genes <- use.genes[which(is.na(use.genes)==FALSE)]
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
    network_project_name <- SJAR.project_name
    tf_out_path <- SJAR.main_dir
    tf_pj <- sprintf('%s_TF',network_project_name)
    sig_out_path <- SJAR.main_dir
    sig_pj <- sprintf('%s_SIG',network_project_name)
    cmd_tf <- sprintf('sjaracne %s %s %s %s %s',tf_pj,SJAR.expression_matrix,SJAR.hub_genes.tf,tf_out_path,add_options)
    cmd_sig <- sprintf('sjaracne %s %s %s %s %s',sig_pj,SJAR.expression_matrix,SJAR.hub_genes.sig,sig_out_path,add_options)
    cat(cmd_tf,file = SJAR.bash_file.tf,sep = '\n')
    cat(cmd_sig,file = SJAR.bash_file.sig,sep = '\n')
    message(sprintf('The running command is:  %s  for TF, and the bash file is generated in %s',cmd_tf,SJAR.bash_file.tf))
    message(sprintf('The running command is:  %s  for SIG, and the bash file is generated in %s',cmd_sig,SJAR.bash_file.sig))
    return(TRUE)
}


#' Inner use: prepare for SJaracne run on SJ server
#' prepare SJAracne dataset for net-dataset
#' project_name expression_matrix hub_genes outdir
#' @param eset eSet object
#' @param use.samples a vector of characters
#' @param TF_list a vector of characters
#' @param SIG_list a vector of characters
#' @param SJAR.project_name character
#' @param SJAR.main_dir character
#' @param SJAR.bash.main_dir character
#' @param mem integer
#' @param IQR.thre numeric
#' @param IQR.loose_thre numeric
#' @export
SJ.SJAracne.prepare <-
  function(eset,use.samples = rownames(pData(eset)),TF_list=NULL,SIG_list=NULL,
           SJAR.project_name = "",
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
    choose1 <- IQR.filter(d, rownames(d),thre = IQR.thre,loose_thre=IQR.loose_thre,loose_gene=unique(c(TF_list,SIG_list)))
    d <- d[choose1, ]
    use.genes <- rownames(d)
    use.genes <- use.genes[which(is.na(use.genes)==FALSE)]
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

#' Inner use: auto generate bash files for Step1, Step2, Step3, Step4 for all bash files under one directory
#' @param out.dir.SJAR character
#' @export
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


#combRowEvid.2grps
#d: dataframe as below: (first column is the annotation)
#			geneSymbol  DMSO_1 DMSO_2 DMSO_3    DEX_1    DEX_2    DEX_3
#P0112072    PRKAR2B 0.50099 1.2108 1.0524 -0.34881 -0.13441 -0.87112
#P0112073    PRKAR2B 1.84579 2.0356 2.6025  1.62954  1.88281  1.29604
#comp as below:
# [1] 0 0 0 1 1 1
# Levels: 0 1
#d<-d[1,]

#' combRowEvid.2grps is a inner function for getDE.BID.2G.
#' @param d input expression matrix
#' @param comp an indicator vector
#' @param family family Objects for Models, default is 'gaussian'
#' @param method choose from 'MLE' or 'Bayesian'
#' @param pooling choose from 'full','no','partial'
#' @param n.iter number of interaction, default is 1000
#' @param prior.V.scale default is 0.02
#' @param prior.R.nu default is 1
#' @param prior.G.nu default is 2
#' @param nitt default is 13000
#' @param burnin default is 3000
#' @param thin default is 10
#' @param restand default is TRUE
#' @param logTransformed default is TRUE
#' @param log.base default is 2
#' @param average.method choose from c('geometric','arithmetic')
#' @param pseudoCount default is 0
#' @export
combRowEvid.2grps<-function(d,comp,family=gaussian,method=c('MLE','Bayesian'),pooling=c('full','no','partial'),n.iter=1000,
                            prior.V.scale=0.02,prior.R.nu=1,prior.G.nu=2,nitt = 13000,burnin =3000,thin=10,
                            restand=TRUE,logTransformed=TRUE,log.base=2,average.method=c('geometric','arithmetic'),pseudoCount=0
){

  if(!all(comp %in% c(1,0))){
    stop('comp only takes 1 or 0 !!! \n')
  }

  if(missing(pooling))
    pooling<-c('full','no','partial')
  else
    pooling<-tolower(pooling)

  pooling<-match.arg(pooling,several.ok=T)


  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (!family$family %in% c('gaussian','binomial', 'poisson')) {
    print(family)
    stop("Only Gaussian Poisson, and Binomial are supported!!! \n")
  }

  if(missing(method))
    method<-'Bayesian'
  if(grepl('Bayes',method,ignore.case = T)) method<-'Bayesian'
  if(grepl('MLE|MaxLikelihood',method,ignore.case = T)) method<-'MLE'
  method<-match.arg(method)


  #remove NAs
  nna<-apply(!is.na(d),2,all)
  d<-d[,nna]
  #comp<-comp[nna[-1]]

  #d<-d[1,]
  #cat(dim(d),'\n')
  #cat(d[1,1],'\n')
  d<-data.frame(t(d[,-1]))
  #cat(dim(d),'\n')




  dat<-data.frame(
    response=
      c(unlist(d[comp==levels(comp)[1],]),unlist(d[comp==levels(comp)[2],]))
    ,
    treatment=
      factor(c(rep(levels(comp)[1],sum(comp==levels(comp)[1])*ncol(d)),rep(levels(comp)[2],sum(comp==levels(comp)[2])*ncol(d))))
    ,
    probe=
      factor(c(rep(colnames(d),each=sum(comp==levels(comp)[1])),rep(colnames(d),each=sum(comp==levels(comp)[2]))))
  )


  #calculate FC
  FC.val<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='geometric',pseudoCount=pseudoCount)

  AveExp<-mean(dat$response)
  n.levels<-nlevels(dat$probe)

  rs<-c(FC=FC.val,AveExp=AveExp,n.levels=as.integer(n.levels))

  #cat(rs,'\n')

  if('arithmetic' %in% average.method){
    FC.ari<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='arithmetic',pseudoCount=0)
    rs<-c(FC.ari_raw=FC.ari,rs)
  }


  #re-standarize the input data
  if(restand & sd(dat$response)>0)
    dat$response<-0.5*(dat$response-mean(dat$response))/sd(dat$response)

  #MLE approach to estimate parameters
  if(method=='MLE'){
    if(family$family=='gaussian'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat)
        else
          M.no<-glm(response ~ treatment, data=dat)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }


      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          #consider sd==0 case
          if(sd(dat$response)==0)
          {
            dat$response<-rnorm(nrow(dat),mean(dat$response),sd(dat$response)+0.001)
          }
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(	sum.tmp$devcomp$dims['n']-	sum.tmp$devcomp$dims['p']-	sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          #########################################
          #Another way to calculate pvalue
          #########################################
          #M.null<-glmer(response ~ (treatment + 1 | probe), data=dat)
          #anova(M.partial, M.null)
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='binomial'){
      if('full'%in%pooling){
        M.full<-glm(treatment ~ response, data=dat, family=family)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(treatment ~ response + probe, data=dat, family=family)
        else
          M.no<-glm(treatment ~ response, data=dat, family=family)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(treatment ~ response + (response + 1 | probe), data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            Dev.partial=as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(treatment ~ response, data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }


    }else if(family$family=='poisson'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat,family='poisson')
        else
          M.no<-glm(response ~ treatment, data=dat,family='poisson')

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            #REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            #REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }

    }
    else{
      stop('Only liner model with gaussian or poisson distrn and binomial family model are supported !!! \n')
    }
  }

  #Bayesian approach
  else if(method=='Bayesian'){
    if(family$family=='gaussian'){
      if('full'%in%pooling){
        #comlete pooling
        #				M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter)
        if(sd(dat$response)>0){
          M.full<-bayesglm(response ~ treatment, data=dat)
          sum.tmp<-summary(M.full)

          rs.full<-c(
            coef.full=sum.tmp$coef[2,1],
            se.full=sum.tmp$coef[2,2],
            t.full=sum.tmp$coef[2,3],
            pval.full=sum.tmp$coef[2,4],
            z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
            df.full=sum.tmp$df.residual-sum.tmp$df[1],
            AIC.full=sum.tmp$aic,
            BIC.full=AIC(M.full,k=log(nrow(dat))),
            Dev.full=sum.tmp$deviance
          )
        }else{
          rs.full<-c(
            coef.full=0,
            se.full=0,
            t.full=0,
            pval.full=1,
            z.full=0,
            df.full=NA,
            AIC.full=NA,
            BIC.full=NA,
            Dev.full=NA
          )
        }
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        #no pooling
        if(n.levels>1)
          #					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment + probe, data=dat)
        else
          #M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.no)

        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)

      }

      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))

        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
        }else{
          prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))

          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)

          df.partial<-nrow(dat)-2

          #cat(n.levels,df.partial,'\n')
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        if(is.na(pval.partial))
          pval.partial<-sum.tmp$sol[2,5]
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }

    }else if(family$family=='binomial'){

      if(family$link=='logit')
        prior.scale<-2.5
      else if(family$link=='probit')
        prior.scale<-2.5*1.6

      if('full'%in%pooling){
        #M.full<-bayesglm(treatment ~ response, data=dat, family=family, n.iter = n.iter, prior.scale = prior.scale)
        M.full<-bayesglm(treatment ~ response, data=dat, family=family, prior.scale = prior.scale)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #a bug in bayeglm for summary, fix: total obs - rank of the model
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        if(n.levels>1)
          #M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family,prior.scale = prior.scale)
        else
          #M.no<-bayesglm(treatment ~ response, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response, data=dat, family=family, prior.scale = prior.scale)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      if('partial'%in%pooling){
        #categorial=logit, ordinal=probit
        if(family$link=='logit'){
          glmm.family<-'categorial'
          if(n.levels>1){
            prior<-list(R = list(V = prior.V.scale, nu=n.levels), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
          }else{
            prior<-list(R = list(V = prior.V.scale, nu=n.levels))
          }
        }
        else if(family$link=='probit'){
          glmm.family<-'ordinal'
          if(n.levels>1){
            prior<-list(R = list(V = prior.V.scale, nu=n.levels+1), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
          }else{
            prior<-list(R = list(V = prior.V.scale, nu=n.levels+1))
          }
        }
        else
          stop('For multilevel model with Binomial family and Bayeisan method, only logit and probit model are supported!!! \n')

        #			coef.scale<-ifelse(family$link=='logit', 2.5^2, (2.5*1.6)^2)
        #			intercept.scale<-ifelse(family$link=='logit', 10^2, (10*1.6)^2)
        #			scale<-var(dat$treatment)/var(dat$response)
        ################################# prior to be modified #################################
        #alpha.mu=rep(0,2),alpha.V=diag(2)*25^2)))
        #########################################################################################
        sd.threshold<-prior.scale*4
        #coef.sign<- -sign(logFC)
        #				while(sd.threshold>prior.scale*2){
        if(n.levels>1){
          M.partial<-MCMCglmm(treatment~response, random=~idh(1+response):probe,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
        }else{
          M.partial<-MCMCglmm(treatment~response,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #						z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        #coef.sign<-sign(rs.partial['coef.partial'])
        sd.threshold<-as.double(rs.partial['sd.partial'])
        #				}
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='poisson'){
      if('full'%in%pooling){
        #comlete pooling
        #M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter,family='poisson')
        M.full<-bayesglm(response ~ treatment, data=dat,family='poisson')

        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        #no pooling
        if(n.levels>1)
          #M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment + probe, data=dat,family='poisson')
        else
          #M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))

        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-(n.levels+1)*2
        }
        else{
          prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #					z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
    }else{
      stop('Only liner model and binomial family model are supported !!! \n')
    }
  }

  #taking care of extreme cases where se or sd is zero
  w1 <- grep('^t|^z',names(rs))
  w2 <- grep('^se|^sd',names(rs))
  w3 <- grep('^pval',names(rs))
  if(rs[w2]==0 & is.na(rs[w1[1]])==TRUE){
    rs[w1]<-0
    rs[w3]<-1
  }
  rs
}



#comparsion of mult-igroups
#only support linear Gausian model so far
combRowEvid.multgrps<-function(d,comp,family=gaussian,method=c('MLE','Bayesian'),pooling=c('full','no','partial'), n.iter=1000,
                               prior.V.scale=0.02, prior.R.nu=1, prior.G.nu=2, nitt = 13000, burnin = 3000,thin=10,restand=TRUE,logTransformed=TRUE,log.base=2,pseudoCount=0){

  if(missing(pooling))
    pooling<-c('full','no','partial')
  else
    pooling<-tolower(pooling)

  pooling<-match.arg(pooling,several.ok=T)

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  #	if (!family$family %in% c('gaussian','binomial', 'poisson')) {
  if (!family$family %in% c('gaussian')) {
    print(family)
    #stop("Only Gaussian Poisson, and Binomial are supported!!! \n")
    stop("Only Gaussian model is supported!!! \n")
  }

  if(missing(method))
    method<-'Bayesian'
  if(grepl('Bayes',method,ignore.case = T)) method<-'Bayesian'
  if(grepl('MLE|MaxLikelihood',method,ignore.case = T)) method<-'MLE'
  method<-match.arg(method)

  d<-melt(t(d[,-1]))

  dat<-data.frame(response=d$value,treatment=rep(comp,nlevels(d$X2)),probe=d$X2)

  AveExp<-mean(dat$response)
  n.levels<-nlevels(dat$probe)
  n.treatments<-nlevels(dat$treatment)
  rs<-c(AveExp=AveExp,n.levels=as.integer(n.levels))

  #re-standarize the input data
  if(restand & sd(dat$response)>0)
    dat$response<-0.5*(dat$response-mean(dat$response))/sd(dat$response)

  #MLE approach to estimate parameters
  if(method=='MLE'){
    if(family$family=='gaussian'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ as.ordered(treatment), data=dat)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[-1,1],
          se.full=sum.tmp$coef[-1,2],
          t.full=sum.tmp$coef[-1,3],
          pval.full=sum.tmp$coef[-1,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat)
        else
          M.no<-glm(response ~ treatment, data=dat)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2:n.treatments,1],
          se.no=sum.tmp$coef[2:n.treatments,2],
          t.no=sum.tmp$coef[2:n.treatments,3],
          pval.no=sum.tmp$coef[2:n.treatments,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }


      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          #consider sd==0 case
          if(sd(dat$response)==0)
          {
            dat$response<-rnorm(nrow(dat),mean(dat$response),sd(dat$response)+0.001)
          }
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat)
          sum.tmp<-summary(M.partial)

          t.partial<-sum.tmp$coef[-1,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          #########################################
          #Another way to calculate pvalue
          #########################################
          #M.null<-glmer(response ~ (treatment + 1 | probe), data=dat)
          #anova(M.partial, M.null)
          rs.partial<-c(
            coef.partial=sum.tmp$coef[-1,1],
            se.partial=sum.tmp$coef[-1,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[-1,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[-1,1],
            se.partial=sum.tmp$coef[-1,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }

    }else if(family$family=='binomial'){
      if('full'%in%pooling){
        M.full<-glm(treatment ~ response, data=dat, family=family)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(treatment ~ response + probe, data=dat, family=family)
        else
          M.no<-glm(treatment ~ response, data=dat, family=family)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(treatment ~ response + (response + 1 | probe), data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            Dev.partial=as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(treatment ~ response, data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }


    }else if(family$family=='poisson'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat,family='poisson')
        else
          M.no<-glm(response ~ treatment, data=dat,family='poisson')

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            #REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            #REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }

    }
    else{
      stop('Only liner model with gaussian or poisson distrn and binomial family model are supported !!! \n')
    }
  }

  #Bayesian approach
  else if(method=='Bayesian'){
    require(arm)
    require(MCMCglmm)
    if(family$family=='gaussian'){
      if('full'%in%pooling){
        #comlete pooling
        #M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter)
        M.full<-bayesglm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.full)

        rs.full<-c(
          coef.full=sum.tmp$coef[-1,1],
          se.full=sum.tmp$coef[-1,2],
          t.full=sum.tmp$coef[-1,3],
          pval.full=sum.tmp$coef[-1,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        #no pooling
        if(n.levels>1)
          #					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment + probe, data=dat)
        else
          #					M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2:n.treatments,1],
          se.no=sum.tmp$coef[2:n.treatments,2],
          t.no=sum.tmp$coef[2:n.treatments,3],
          pval.no=sum.tmp$coef[2:n.treatments,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(nlevels(dat$treatment))*prior.V.scale, nu = prior.G.nu )))
        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*n.treatments
        }else{
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-n.treatments
        }
        sum.tmp<-summary(M.partial)
        sd.partial<-apply(M.partial$Sol[,-1],2,sd)
        t.partial<-sum.tmp$sol[-1,1]/sd.partial
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))

        rs.partial<-c(
          coef.partial=sum.tmp$sol[-1,1],
          'l-95% CI.partial'=sum.tmp$sol[-1,2],
          'u-95% CI.partial'=sum.tmp$sol[-1,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[-1,5],
          #effective sample size
          n.effet.partial=sum.tmp$sol[-1,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd.partial,
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='binomial'){

      if(family$link=='logit')
        prior.scale<-2.5
      else if(family$link=='probit')
        prior.scale<-2.5*1.6

      if('full'%in%pooling){
        #M.full<-bayesglm(treatment ~ response, data=dat, family=family, n.iter = n.iter, prior.scale = prior.scale)
        M.full<-bayesglm(treatment ~ response, data=dat, family=family, prior.scale = prior.scale)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #a bug in bayeglm for summary, fix: total obs - rank of the model
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        if(n.levels>1)
          #M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family,prior.scale = prior.scale)
        else
          #					M.no<-bayesglm(treatment ~ response, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response, data=dat, family=family,prior.scale = prior.scale)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      if('partial'%in%pooling){
        #categorial=logit, ordinal=probit
        if(family$link=='logit'){
          glmm.family<-'categorial'
          prior<-list(R = list(V = prior.V.scale, nu=n.levels), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
        }
        else if(family$link=='probit'){
          glmm.family<-'ordinal'
          prior<-list(R = list(V = prior.V.scale, nu=n.levels+1), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
        }
        else
          stop('For multilevel model with Binomial family and Bayeisan method, only logit and probit model are supported!!! \n')

        #			coef.scale<-ifelse(family$link=='logit', 2.5^2, (2.5*1.6)^2)
        #			intercept.scale<-ifelse(family$link=='logit', 10^2, (10*1.6)^2)
        #			scale<-var(dat$treatment)/var(dat$response)
        ################################# prior to be modified #################################
        #alpha.mu=rep(0,2),alpha.V=diag(2)*25^2)))
        #########################################################################################
        sd.threshold<-prior.scale*4
        #coef.sign<- -sign(logFC)
        #				while(sd.threshold>prior.scale*2){
        if(n.levels>1){
          M.partial<-MCMCglmm(treatment~response, random=~idh(1+response):probe,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
        }else{
          M.partial<-MCMCglmm(treatment~response,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #						z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        #coef.sign<-sign(rs.partial['coef.partial'])
        sd.threshold<-as.double(rs.partial['sd.partial'])
        #				}
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='poisson'){
      if('full'%in%pooling){
        #comlete pooling
        #M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter,family='poisson')
        M.full<-bayesglm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        #no pooling
        if(n.levels>1)
          #					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment + probe, data=dat, family='poisson')
        else
          #					M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))
        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-(n.levels+1)*2
        }
        else{
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #					z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
    }else{
      stop('Only liner model and binomial family model are supported !!! \n')
    }
  }
  rs
}
##
#fold change function
#positive: class1/class0
#negative: class0/class1
FC <- function(x,cl,logTransformed = TRUE,
           log.base = 2,average.method = c('geometric', 'arithmetic'),
           pseudoCount = 0) {
    x.class0 <- x[(cl == 0)] + pseudoCount
    x.class1 <- x[(cl == 1)] + pseudoCount
    if (missing(average.method))
      average.method <- 'geometric'
    if (logTransformed) {
      if (is.na(log.base) | log.base < 0)
        stop('You must specify log.bsae !\n')
      logFC <- mean(x.class1) - mean(x.class0)
      FC.val <- sign(logFC) * log.base ^ abs(logFC)
    } else{
      logFC <-
        ifelse(average.method == 'arithmetic',
               log(mean(x.class1)) - log(mean(x.class0)),
               mean(log(x.class1) - mean(log(x.class0))))
      FC.val <- sign(logFC) * exp(abs(logFC))
    }
    FC.val[FC.val == 0 | is.na(FC.val)] <- 1
    FC.val
  }
