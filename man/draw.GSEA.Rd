% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_functions.R
\name{draw.GSEA}
\alias{draw.GSEA}
\title{GSEA (Gene Set Enrichment Analysis) Plot for one Gene Set or one Driver}
\usage{
draw.GSEA(
  rank_profile = NULL,
  use_genes = NULL,
  use_direction = NULL,
  main = "",
  pdf_file = NULL,
  annotation = NULL,
  annotation_cex = 1.2,
  left_annotation = NULL,
  right_annotation = NULL
)
}
\arguments{
\item{rank_profile}{a named vector of numerics, the differential values (DE or DA) calculated from a sample comparison (e.g. "G4 vs. Others").
Names of the vector must be gene names.
For the DA, user could use `processDriverProfile()` to convert the DA profile into gene-name based profile.
The differential values can be "logFC" or "t-statistics".}

\item{use_genes}{a vector of characters, a vector of genes to display. The genes can either be annotated genes in gene set or the targe genes from a specific driver.
The gene names must be a subset of \code{names(rank_profile)}.}

\item{use_direction}{a vector of numeric 1s and -1s, 1 is positive regulation from driver, -1 is negative regulation from driver.
Users can get this vector by converting the signs of "spearman". If NULL, no regulation direction will be displayed. Default is NULL.}

\item{main}{character, an overall title for the plot. Default is "".}

\item{pdf_file}{character, the file path to save as PDF file. If NULL, no PDF file will be saved. Default is NULL.}

\item{annotation}{character, the annotation set by users for easier reference.
Normally the annotation is the P-value or other statistics to show the significance of the interested gene set or driver.
If NULL, will perform a Kolmogorov-Smirnov test to get the significance value.
If want to get the statistics from GSEA test, could use `funcEnrich.GSEA()` to get the statistics first.
Default is NULL.}

\item{annotation_cex}{numeric, giving the amount by which the text of annotation should be magnified relative to the default. Default is 1.2.}

\item{left_annotation}{character, annotation displayed on the left of the figure, representing left condition of the \code{rank_profile}. Default is "".}

\item{right_annotation}{character, annotation displayed on the right of the figure, representing right condition of the \code{rank_profile}. Default is "".}
}
\value{
Return a logical value. If TRUE, the plot has been created successfully.
}
\description{
\code{draw.GSEA} draws a GSEA plot to analyze one gene set (with gene list annotated) or one driver (with list of target genes).
}
\examples{
gs.preload(use_spe='Homo sapiens',update=FALSE)
analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
ms_tab <- analysis.par$final_ms_tab

## draw for the most significant gene set
# by driver's DA profile
DA_profile <- processDriverProfile(Driver_name=ms_tab$gene_label,
                                    Driver_profile=ms_tab$logFC.G4.Vs.others_DA,
                                    choose_strategy='absmax',
                                    return_type ='gene_statistics')
res1 <- funcEnrich.GSEA(rank_profile=DA_profile,
                         use_gs=c('H'),
                         Pv_thre=0.1,Pv_adj = 'none')
top_gs <- res1[1,'#Name'] ## draw for the top 1
annot <- sprintf('NES: \%s \nAdjusted P-value: \%s',
          signif(res1[1,'NES'],2),
          signif(res1[1,'Adj_P'],2))
draw.GSEA(rank_profile=DA_profile,
          use_genes=all_gs2gene$H[[top_gs]],
          main=sprintf('GSEA plot for gene set \%s',
          top_gs),
          annotation=annot,annotation_cex=1.2,
          left_annotation='high in G4',
          right_annotation='high in others')

## draw for the most significant driver
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
                               logFC_col='logFC.G4.Vs.others_DA',
                               Pv_col='P.Value.G4.Vs.others_DA',
                               logFC_thre=0.4,
                               Pv_thre=1e-7,
                               main='Volcano Plot for G4.Vs.others_DA',
                               show_label=FALSE,
                               label_type = 'origin',
                               label_cex = 0.5)
driver_list <- rownames(sig_driver)
DE_profile <- analysis.par$DE[[1]]$`Z-statistics`;
names(DE_profile) <- rownames(analysis.par$DE[[1]])
use_driver <- driver_list[1]
use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
use_target_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman) ## 1/-1
annot <- sprintf('P-value: \%s',signif(ms_tab[use_driver,'P.Value.G4.Vs.others_DA'],2))

## draw for the driver
draw.GSEA(rank_profile=DE_profile,
          use_genes=use_target_genes,
          use_direction=use_target_direction,
          main=sprintf('GSEA plot for driver \%s',
          ms_tab[use_driver,'gene_label']),
          annotation=annot,annotation_cex=1.2,
          left_annotation='high in G4',
          right_annotation='high in others')
draw.GSEA(rank_profile=DE_profile,
          use_genes=use_target_genes,
          use_direction=NULL,
          main=sprintf('GSEA plot for driver \%s',
          ms_tab[use_driver,'gene_label']),
          annotation=NULL,annotation_cex=1.2,
          left_annotation='high in G4',
          right_annotation='high in others')

## draw for the gene set
gs.preload(use_spe='Homo sapiens',update=FALSE)
use_target_genes <- all_gs2gene[[1]][[1]]
draw.GSEA(rank_profile=DE_profile,
         use_genes=use_target_genes,
         main=sprintf('GSEA plot for \%s',names(all_gs2gene[[1]][1])),
         left_annotation='high in G4',
         right_annotation='high in others')

\dontrun{
#' analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
ms_tab <- analysis.par$final_ms_tab
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
                               logFC_col='logFC.G4.Vs.others_DA',
                               Pv_col='P.Value.G4.Vs.others_DA',
                               logFC_thre=0.4,
                               Pv_thre=1e-7,
                               main='Volcano Plot for G4.Vs.others_DA',
                               show_label=FALSE,
                               label_type = 'origin',
                               label_cex = 0.5)
driver_list <- rownames(sig_driver)
DE_profile <- analysis.par$DE[[1]]$`Z-statistics`;
names(DE_profile) <- rownames(analysis.par$DE[[1]])
use_driver <- driver_list[1]
use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
use_target_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman) ## 1/-1
annot <- sprintf('P-value: \%s',signif(ms_tab[use_driver,'P.Value.G4.Vs.others_DA'],2))
analysis.par$out.dir.PLOT <- getwd() ## directory for saving the pdf files
draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,
          use_direction=use_target_direction,
          main=sprintf('GSEA plot for driver \%s',ms_tab[use_driver,'gene_label']),
          pdf_file = sprintf('\%s/GSEA_driver.pdf',
          analysis.par$out.dir.PLOT),
          annotation=annot,annotation_cex=1.2,
          left_annotation='high in G4',
          right_annotation='high in others')

## draw for the gene set
gs.preload(use_spe='Homo sapiens',update=FALSE)
use_target_genes <- all_gs2gene[[1]][[1]]
draw.GSEA(rank_profile=DE_profile,
         use_genes=use_target_genes,
         main=sprintf('GSEA plot for \%s',names(all_gs2gene[[1]][1])),
         pdf_file = sprintf('\%s/GSEA_GS_each.pdf',analysis.par$out.dir.PLOT),
         left_annotation='high in G4',
         right_annotation='high in others')
}
}
