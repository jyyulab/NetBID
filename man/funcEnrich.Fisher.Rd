% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_functions.R
\name{funcEnrich.Fisher}
\alias{funcEnrich.Fisher}
\title{Gene Set Enrichment Analysis by Fisher's Exact Test}
\usage{
funcEnrich.Fisher(
  input_list = NULL,
  bg_list = NULL,
  use_gs = NULL,
  gs2gene = NULL,
  min_gs_size = 5,
  max_gs_size = 500,
  Pv_adj = "fdr",
  Pv_thre = 0.1
)
}
\arguments{
\item{input_list}{a vector of characters, a vector of gene symbols. If gene symbols are not available, users can call \code{get_IDtransfer}
and \code{get_name_transfertab} for ID conversion.}

\item{bg_list}{a vector of characters, a vector of background gene symbols. If NULL, genes in \code{gs2gene} will be used as background.
Default is NULL.}

\item{use_gs}{a vector of characters, the names of gene sets.
If \code{gs2gene} is NULL, \code{all_gs2gene} will be used. The \code{use_gs} must be the subset of \code{names(all_gs2gene)}.
If "all", all the gene sets in \code{gs2gene} will be used.
If user input his own \code{gs2gene} list, \code{use_gs} will be set to "all" as default.
Default is c("H", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG").}

\item{gs2gene}{list, a list contains elements of gene sets.
The name of the element is gene set, each element contains a vector of genes in that gene set.
If NULL, will use \code{all_gs2gene}, which is created by function \code{gs.preload}. Default is NULL.}

\item{min_gs_size}{numeric, the minimum size of gene set to analysis. Default is 5.}

\item{max_gs_size}{numeric, the maximum size of gene set to analysis. Default is 500.}

\item{Pv_adj}{character, method to adjust P-value. Default is "fdr".
For details, please check \code{p.adjust.methods}.}

\item{Pv_thre}{numeric, threshold for the adjusted P-values. Default is 0.1.}
}
\value{
Return a data.frame, contains gene sets with significant enrichment statistics. Column details are as follows,

\item{#Name}{Name of the enriched gene set}
\item{Total_item}{Background size}
\item{Num_item}{Number of genes in the gene set (filtered by the background list)}
\item{Num_list}{Number of input genes for testing (filtered by the background list)}
\item{Num_list_item}{Number of input genes annotated by the gene set (filtered by the background list)}
\item{Ori_P}{Original P-value from Fisher's Exact Test}
\item{Adj_P}{Adjusted P-value}
\item{Odds_Ratio}{Odds ratio from the 2*2 matrix used for Fisher's Exact Test}
\item{Intersected_items}{A vector of the intersected genes, collapsed by ';'. Number is equal to Num_list_item}
}
\description{
\code{funcEnrich.Fisher} performs gene set enrichment analysis to the input gene list, by using the Fisher's Exact Test.
Background gene list is accepeted.
}
\examples{
analysis.par <- list()
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
gs.preload(use_spe='Homo sapiens',update=FALSE)
res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
                               bg_list=ms_tab[,'geneSymbol'],
                               use_gs=c('H','C5'),
                               Pv_thre=0.1,Pv_adj = 'none')
\dontrun{
}
}
