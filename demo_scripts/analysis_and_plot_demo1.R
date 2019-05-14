## This script aims to do following analysis ##
############# preparations ####################

library(NetBID2)

# set the directory to get the pre-saved RData
analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
print(analysis.par$out.dir.PLOT)
#analysis.par$out.dir.PLOT <- 'test/driver/PLOT' ## directory to save the plot figures

## Q1: How to get the top list of drivers with significantly different activity in G4 Vs. other subtypes ?
ms_tab <- analysis.par$final_ms_tab ## get the master table data frame
ms_tab <- ms_tab[which(ms_tab$Size>=30 & ms_tab$Size <=1000),] ## target size filteration
comp_name <- 'G4.Vs.others' ## get the comparison name
##
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                               Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.4,Pv_thre=1e-8,
                               main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                               pdf_file=sprintf('%s/vocalno_label_DA.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.4,Pv_thre=1e-8,
                                main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=FALSE,
                                pdf_file=sprintf('%s/vocalno_nolabel_DA.pdf',analysis.par$out.dir.PLOT))
sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                             Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-4,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                             pdf_file=sprintf('%s/vocalno_label_DE.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                             Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-4,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=FALSE,
                             pdf_file=sprintf('%s/vocalno_nolabel_DE.pdf',analysis.par$out.dir.PLOT))

### QI.2: How to understand the significance of those top DA drivers ?
DE <- analysis.par$DE[[comp_name]]
driver_list <- rownames(sig_driver) ## the rownames is the originalID_label
##
draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',name_col='ID',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in others',right_annotation = 'high in G4',
                 main= comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo1.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main= comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo2.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='black',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main= comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo3.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='RdBu',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main=comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo4.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='black',
                 left_annotation = 'high in G4',right_annotation = 'high in others',
                 main= comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo5.pdf',analysis.par$out.dir.PLOT))
                 
                 
### QI.3: What is the expression/activity pattern of those top DA drivers ?
exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames must be the originalID
ac_mat <- exprs(analysis.par$merge.ac.eset) ## ac,the rownames must be the originalID_label
phe_info <- pData(analysis.par$cal.eset)
# could use addtional paramters in Heatmap()
draw.heatmap(mat=exp_mat,use_genes=ms_tab[driver_list,'originalID'],use_gene_label=ms_tab[driver_list,'geneSymbol'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo1.pdf',analysis.par$out.dir.PLOT))

draw.heatmap(mat=ac_mat,use_genes=driver_list,use_gene_label=ms_tab[driver_list,'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo2.pdf',analysis.par$out.dir.PLOT))

# try for the top DE genes             
gene_list <- rownames(sig_gene) 
draw.heatmap(mat=exp_mat,use_genes=ms_tab[gene_list,'originalID'],use_gene_label=ms_tab[gene_list,'geneSymbol'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo3.pdf',analysis.par$out.dir.PLOT))

draw.heatmap(mat=ac_mat,use_genes= gene_list,use_gene_label=ms_tab[gene_list,'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup','age'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo4.pdf',analysis.par$out.dir.PLOT))        
             
             
### QI.4: What is the biological function of those top DA drivers ?
# load RData for all_gs2gene
gs.preload(use_spe='Homo sapiens',update=FALSE)

# get function enrichment results
driver_list_up <- rownames(sig_driver)[which(sig_driver[,2]>0)] # up
driver_list_down <- rownames(sig_driver)[which(sig_driver[,2]<0)] # down

res_up <- funcEnrich.Fisher(input_list=ms_tab[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP','CGP','CP:KEGG'),
                           Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
                           
res_down <- funcEnrich.Fisher(input_list=ms_tab[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP','CGP','CP:KEGG'),
                           Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

out2excel(list(up=res_up,down=res_down),out.xlsx=sprintf('%s/fisher_res.xlsx',analysis.par$out.dir.PLOT))
                           
# draw barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers',pdf_file=sprintf('%s/funcEnrich_bar_nogene.pdf',analysis.par$out.dir.PLOT))
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_withgene.pdf',analysis.par$out.dir.PLOT))


draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95)
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGS.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = FALSE)
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGENE.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = TRUE)
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 1.5,gene_cex=1.4,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterNO.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = FALSE)

### QI.5: What is the biological function of the target genes of those top DA drivers ?
transfer_tab <- analysis.par$transfer_tab
## draw
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[driver_list,'gene_label'],
                Z_val=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=5,max_gs_size=500,use_gs=c('CP:KEGG','H'),
                top_geneset_number=10,top_driver_number=10,
                pdf_file = sprintf('%s/bubblePlot.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')

draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[driver_list,'gene_label'],
                Z_val=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=10,max_gs_size=300,use_gs=c('CP:KEGG'),
                top_geneset_number=30,top_driver_number=30,
                pdf_file = sprintf('%s/bubblePlot_KEGG.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')
# add marker gene
mark_gene <- c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D') ## marker for Group4
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[rownames(sig_driver),'gene_label'],
                Z_val=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=10,max_gs_size=300,use_gs=c('CP:KEGG','CP:BIOCARTA','H'),
                top_geneset_number=30,top_driver_number=80,
                pdf_file = sprintf('%s/bubblePlot_combine.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets',
                mark_gene=ms_tab[which(ms_tab$geneSymbol %in% mark_gene),'originalID_label'],gs_cex = 1,driver_cex=1.2)

## Part II: for a selected interested driver

### QII.1: How to understand the significance of the selected driver ?
DE <- analysis.par$DE[[comp_name]]
DE_profile <- DE$`Z-statistics`; 
names(DE_profile) <- rownames(DE)

use_driver <- driver_list[1]
use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
use_target_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman) ## 1/-1
annot <- sprintf('P-value: %s',signif(ms_tab[use_driver,sprintf('P.Value.%s_DA',comp_name)],2))

## with direction
draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,use_direction=use_target_direction,
          main=sprintf('GSEA plot for driver %s',ms_tab[use_driver,'gene_label']),
          pdf_file = sprintf('%s/GSEA_each_direction.pdf',analysis.par$out.dir.PLOT),
          annotation=annot,annotation_cex=1.2,
          left_annotation='high in G4',right_annotation='high in others')

## GSEA plot without direction, without annotation(use ks)
draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,use_direction=NULL,
          main=sprintf('GSEA plot for driver %s',ms_tab[use_driver,'gene_label']),
          pdf_file = sprintf('%s/GSEA_each.pdf',analysis.par$out.dir.PLOT),
          annotation='',annotation_cex=1.2,
          left_annotation='high in G4',right_annotation='high in others')


### QII.2: How to visualize the network structure of the selected driver ?
use_driver <- driver_list[1]
edge_score <- analysis.par$merge.network$target_list[[use_driver]]$MI*sign(analysis.par$merge.network$target_list[[use_driver]]$spearman)
names(edge_score) <- analysis.par$merge.network$target_list[[use_driver]]$target

draw.targetNet(source_label=ms_tab[use_driver,'gene_label'],source_z=ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],
               edge_score = edge_score,pdf_file=sprintf('%s/targetNet_out.pdf',analysis.par$out.dir.PLOT),label_cex = 0.4,n_layer=4, alphabetical_order=TRUE)
               
draw.targetNet(source_label=ms_tab[use_driver,'gene_label'],source_z=ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],
               edge_score = edge_score,pdf_file=sprintf('%s/targetNet_in.pdf',analysis.par$out.dir.PLOT),label_cex = 0.35,arrow_direction = 'in',n_layer=6)

# for two target
use_driver2 <- driver_list[2]
edge_score2 <- analysis.par$merge.network$target_list[[use_driver2]]$MI*sign(analysis.par$merge.network$target_list[[use_driver2]]$spearman)
names(edge_score2) <- analysis.par$merge.network$target_list[[use_driver2]]$target
use_genes <- unique(analysis.par$merge.network$network_dat$target.symbol)
draw.targetNet.TWO(source1_label=ms_tab[use_driver,'gene_label'],edge_score1 = edge_score,
                   source2_label=ms_tab[use_driver2,'gene_label'],edge_score2 = edge_score2,
                   source1_z=ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],source2_z=ms_tab[use_driver2,sprintf('Z.%s_DA',comp_name)],
                   pdf_file=sprintf('%s/targetNetTWO.pdf',analysis.par$out.dir.PLOT),total_possible_target=use_genes,show_test=TRUE,
                   label_cex = 0.4,n_layer=1,source_cex=0.6,alphabetical_order=FALSE)

# or directly get the test result
test.targetNet.overlap(source1_label=ms_tab[use_driver,'gene_label'],source2_label=ms_tab[use_driver2,'gene_label'],
                       target1 = names(edge_score),target2 = names(edge_score2),total_possible_target=use_genes)

### QII.3: What is the expression/activity for this selected driver in samples with different subtypes ?
use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
#
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=exp_mat[ms_tab[use_driver,'originalID'],],use_obs_class=use_obs_class,
                   class_order=c('WNT','SHH','G4'),class_srt=30,main_ac = ms_tab[use_driver,'gene_label'],main_exp=ms_tab[use_driver,'geneSymbol'],
                   pdf_file=sprintf('%s/categoryValue_demo1.pdf',analysis.par$out.dir.PLOT))
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=NULL,use_obs_class=use_obs_class,class_order=c('WNT','SHH','G4'),
                   pdf_file=sprintf('%s/categoryValue_demo2.pdf',analysis.par$out.dir.PLOT))


#
use_obs_class <- get_obs_label(phe_info = phe_info,c('subgroup','gender'))
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=exp_mat[ms_tab[use_driver,'originalID'],],use_obs_class=use_obs_class,
                   class_srt=30,main_ac = ms_tab[use_driver,'gene_label'],main_exp=ms_tab[use_driver,'geneSymbol'],
                   pdf_file=sprintf('%s/categoryValue_demo3.pdf',analysis.par$out.dir.PLOT))


### QII.4: What is the function of the target genes for this selected driver ?

use_driver <- driver_list[1]
use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
use_target_genes <- get_name_transfertab(use_genes= use_target_genes,transfer_tab=transfer_tab,ignore_version=TRUE)	
bg_list <- get_name_transfertab(use_genes= unique(analysis.par$merge.network$network_dat$target),transfer_tab=transfer_tab,ignore_version=TRUE)	
res <- funcEnrich.Fisher(input_list= use_target_genes,bg_list= bg_list,use_gs=c('H','CP:REACTOME','BP','CGP','CP:KEGG'),
                           Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
draw.funcEnrich.cluster(funcEnrich_res= res,top_number=20,gs_cex = 1.2,gene_cex=1,pv_cex=1,Pv_thre=0.1,
                        pdf_file = sprintf('%s/funcEnrich_clusterBOTH_%s.pdf',analysis.par$out.dir.PLOT,use_driver),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95)


## Part III: further analysis

### QIII.1: What are the activities of the curated gene sets in all samples and what are the top significantly differed gene sets ?
db.preload(use_level='gene',use_spe='human',update=FALSE)
exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames must be the originalID
## get expression matrix for the transfered gene name, if original is gene-based expression matrix, just use the exp_mat
exp_mat_gene <- exp_mat

## calculate activity for all genesets
use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,cal_mat = exp_mat_gene)

## get DA
phe_info <- pData(analysis.par$cal.eset)
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # get sample list for G0
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # get sample list for G1
DA_gs_bid <- getDE.BID.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,G1_name='G4',G0_name='others')

## draw vocalno plot for top sig-GS
sig_gs <- draw.volcanoPlot(dat= DA_gs_bid,label_col='ID',logFC_col='logFC',
                               Pv_col='P.Value',logFC_thre=0.2,Pv_thre=1e-3,
                               main='Volcano Plot for gene sets',show_label=TRUE,label_type = 'distribute',label_cex = 0.5,
                               pdf_file=sprintf('%s/vocalno_GS_DA.pdf',analysis.par$out.dir.PLOT))
## draw heatmap for top sig-GS
draw.heatmap(mat=ac_gs[sig_gs$ID,],pdf_file=sprintf('%s/heatmap_GS.pdf',analysis.par$out.dir.PLOT),scale='row',
             phenotype_info=phe_info,use_phe=c('gender','subgroup'))

## draw GSEA plot for top sig-GS
DE <- analysis.par$DE[[comp_name]]
#
draw.GSEA.NetBID.GS(DE=DE,name_col='ID',profile_col='t',profile_trend='pos2neg',
                 sig_gs_list = sig_gs$ID,
                 gs_DA_Z= DA_gs_bid[sig_gs$ID,'Z-statistics'],
                 use_gs2gene = use_gs2gene,
                 top_gs_number=20,target_col='RdBu',
                 left_annotation = 'high in others',right_annotation = 'high in G4',
                 main= comp_name,Z_sig_thre=1.64,profile_sig_thre = 0,
                 pdf_file=sprintf('%s/NetBID_GSEA_GS_demo1.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.GS(DE=DE,name_col='ID',profile_col='t',profile_trend='pos2neg',
                    sig_gs_list = sig_gs$ID,
                    gs_DA_Z= DA_gs_bid[sig_gs$ID,'Z-statistics'],
                    use_gs2gene = use_gs2gene,
                    top_gs_number=20,target_col='black',
                    left_annotation = 'high in others',right_annotation = 'high in G4',
                    main= comp_name,Z_sig_thre=1.64,profile_sig_thre = 0,
                    pdf_file=sprintf('%s/NetBID_GSEA_GS_demo2.pdf',analysis.par$out.dir.PLOT))

## draw GSEA plot for each sig-GS
DE_profile <- DE$`Z-statistics`; names(DE_profile) <- rownames(DE)
use_gs <- sig_gs$ID[1]
use_target_genes <- rownames(DE)[which(DE$ID %in% use_gs2gene[[use_gs]])]
draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,
          main=sprintf('GSEA plot for %s',use_gs),
          pdf_file = sprintf('%s/GSEA_GS_each.pdf',analysis.par$out.dir.PLOT),
          left_annotation='high in G4',right_annotation='high in others',
          annotation=sprintf('P-value: %s',signif(sig_gs[use_gs,'P.Value'],2)))          
## draw category plot for each sig-GS
use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
draw.categoryValue(ac_val=ac_gs[use_gs,],use_obs_class=use_obs_class,
                   class_order=c('WNT','SHH','G4'),class_srt=30,pdf_file=sprintf('%s/categoryValue_GS_demo1.pdf',analysis.par$out.dir.PLOT),
                   main_ac= use_gs,main_cex=0.8)



### QIII.2: How to find drivers with significantly overlapped target genes ?
gs2gene_target <- analysis.par$merge.network$target_list[driver_list]
gs2gene_target <- lapply(gs2gene_target,function(x)x$target)
transfer_tab_fake <- data.frame(from=transfer_tab[,1],to=transfer_tab[,1],type=transfer_tab[,3],stringsAsFactors=F)
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[driver_list,'gene_label'],
                Z_val=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab_fake,
                bg_list=ms_tab[,'geneSymbol'],gs2gene=gs2gene_target,
                min_gs_size=10,max_gs_size=1000,
                use_gs='all',
                top_geneset_number=10,top_driver_number=10,
                pdf_file = sprintf('%s/bubblePlot_overlap.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')


