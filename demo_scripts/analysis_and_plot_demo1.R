## This script aims to do following analysis ##
############# preparations ####################

library(NetBID2)

# set the directory to get the pre-saved RData
analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")

############ Step0: load
# RData for ms table
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
analysis.par$out.dir.PLOT <- 'test/PLOT' ## directory to save the plot figures

################# PartI: overview plots
############ I.1: volcano plot for DE/DA, p-value Vs. fold-change
ms_tab <- analysis.par$final_ms_tab
sig_driver1 <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col='logFC.G4.Vs.others_DA',
                               Pv_col='P.Value.G4.Vs.others_DA',logFC_thre=0.4,Pv_thre=1e-8,
                               main='Volcano Plot for G4.Vs.others_DA',show_label=TRUE,
                               pdf_file=sprintf('%s/vocalno_showlabel_distribute.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
sig_driver2 <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col='logFC.G4.Vs.others_DA',
                               Pv_col='P.Value.G4.Vs.others_DA',logFC_thre=0.4,Pv_thre=1e-9,
                               main='Volcano Plot for G4.Vs.others_DA',show_label=TRUE,label_type = 'origin',label_cex = 0.5,
                               pdf_file=sprintf('%s/vocalno_showlabel_origin.pdf',analysis.par$out.dir.PLOT))
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col='logFC.G4.Vs.others_DA',
                                Pv_col='P.Value.G4.Vs.others_DA',logFC_thre=0.4,Pv_thre=1e-7,
                                main='Volcano Plot for G4.Vs.others_DA',show_label=FALSE,label_type = 'origin',label_cex = 0.5,
                                pdf_file=sprintf('%s/vocalno_nolabel_DA.pdf',analysis.par$out.dir.PLOT))
sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col='logFC.G4.Vs.others_DE',
                             Pv_col='P.Value.G4.Vs.others_DE',logFC_thre=2,Pv_thre=1e-3,
                             main='Volcano Plot for G4.Vs.others_DE',show_label=FALSE,
                             pdf_file=sprintf('%s/vocalno_nolabel_DE.pdf',analysis.par$out.dir.PLOT))

############ I.2: Heatmap for top drivers
ms_tab  <- analysis.par$final_ms_tab
exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames must be the originalID
ac_mat <- exprs(analysis.par$merge.ac.eset) ## ac,the rownames must be the originalID_label
phe_info <- pData(analysis.par$cal.eset)
# could use addtional paramters in Heatmap()
draw.heatmap(mat=exp_mat,use_genes=ms_tab[rownames(sig_driver),'originalID'],use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','pathology','subgroup'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_demo1.pdf',analysis.par$out.dir.PLOT))

draw.heatmap(mat=ac_mat,use_genes=ms_tab[rownames(sig_driver),'originalID_label'],use_gene_label=ms_tab[rownames(sig_driver),'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=phe_info[colnames(exp_mat),'geo_accession'],
             phenotype_info=phe_info,use_phe=c('gender','subgroup'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 6),pdf_file=sprintf('%s/heatmap_demo2.pdf',analysis.par$out.dir.PLOT))

############ I.3: Function enrichment plot for top drivers
# load RData for all_gs2gene
gs.preload(use_spe='Homo sapiens',update=FALSE)
# get function enrichment results
res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],bg_list=ms_tab[,'geneSymbol'],use_gs=c('H','C5'),Pv_thre=0.1,Pv_adj = 'none')
# draw barplot
draw.funcEnrich.bar(funcEnrich_res=res1,top_number=30,main='Function Enrichment for Top drivers',pdf_file=sprintf('%s/funcEnrich_bar_nogene.pdf',analysis.par$out.dir.PLOT))
draw.funcEnrich.bar(funcEnrich_res=res1,top_number=30,main='Function Enrichment for Top drivers',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_withgene.pdf',analysis.par$out.dir.PLOT))

# draw cluster results for function enrichment (here provide pdf_file as input for directly output into pdf file, need for other functions?)
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_cluster.pdf',analysis.par$out.dir.PLOT))
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE)
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGS.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = FALSE)
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGENE.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = TRUE)
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 1.5,gene_cex=1.4,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterNO.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = FALSE)

############ I.4: Bubble plot for top drivers
# need to prepare ID transfer tab first
# load db
db.preload(use_level='gene',use_spe='human',update=FALSE)
# get transfer table for all target genes
ms_tab  <- analysis.par$final_ms_tab
use_genes <- unique(analysis.par$merge.network$network_dat$target.symbol)
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes) ## get transfer table !!!
## draw
draw.bubblePlot(driver_list=rownames(sig_driver),show_label=ms_tab[rownames(sig_driver),'gene_label'],
                Z_val=ms_tab[rownames(sig_driver),'Z.G4.Vs.others_DA'],
                driver_type=ms_tab[rownames(sig_driver),'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=5,max_gs_size=500,use_gs=c('H'),
                top_geneset_number=30,top_driver_number=50,
                pdf_file = sprintf('%s/bubblePlot.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')

draw.bubblePlot(driver_list=rownames(sig_driver),show_label=ms_tab[rownames(sig_driver),'gene_label'],
                Z_val=ms_tab[rownames(sig_driver),'Z.G4.Vs.others_DA'],
                driver_type=ms_tab[rownames(sig_driver),'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=10,max_gs_size=300,use_gs=c('CP:KEGG'),
                top_geneset_number=30,top_driver_number=50,
                pdf_file = sprintf('%s/bubblePlot_KEGG.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')
# add marker gene
mark_gene <- c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D') ## marker for Group4
draw.bubblePlot(driver_list=rownames(sig_driver),show_label=ms_tab[rownames(sig_driver),'gene_label'],
                Z_val=ms_tab[rownames(sig_driver),'Z.G4.Vs.others_DA'],
                driver_type=ms_tab[rownames(sig_driver),'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=10,max_gs_size=300,use_gs=c('CP:KEGG','CP:BIOCARTA','H'),
                top_geneset_number=50,top_driver_number=80,
                pdf_file = sprintf('%s/bubblePlot_combine.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets',mark_gene=ms_tab[which(ms_tab$geneSymbol %in% mark_gene),'originalID_label'],gs_cex = 1,driver_cex=1.2)

############ I.5: GSEA plot for top driver
## NetBID GSEA plot
db.preload(use_level='gene',use_spe='human',update=FALSE)
## prepare
# get DE profile for all genes
ms_tab <- analysis.par$final_ms_tab
ms_tab <- ms_tab[which(ms_tab$Size>50),]
comp <- 'G4.Vs.others'
DE <- analysis.par$DE[[comp]]
driver_list <- rownames(sig_driver)
driver_list <- driver_list[which(driver_list %in% ms_tab$originalID_label)]
###
draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',name_col='ID',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'test_left',right_annotation = 'test_right',
                 main='test',target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo1.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'test_left',right_annotation = 'test_right',
                 main='test',target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo2.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='black',
                 left_annotation = 'test_left',right_annotation = 'test_right',
                 main='test',target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo3.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='RdBu',
                 left_annotation = 'test_left',right_annotation = 'test_right',
                 main='test',target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo4.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,'Z.G4.Vs.others_DA'],
                 driver_DE_Z=ms_tab[driver_list,'Z.G4.Vs.others_DE'],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=1,target_col='black',
                 left_annotation = 'test_left',right_annotation = 'test_right',
                 main='test',target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_demo5.pdf',analysis.par$out.dir.PLOT))

################# PartII: for each driver
############ II.1: GSEA plot for the driver
ms_tab <- analysis.par$final_ms_tab
DE_profile <- analysis.par$DE[[1]]$`Z-statistics`; names(DE_profile) <- rownames(analysis.par$DE[[1]])
driver_list <- rownames(sig_driver)
use_driver <- driver_list[order(sig_driver[,3])[1]]
use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
use_target_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman) ## 1/-1
annot <- sprintf('P-value: %s',signif(ms_tab[use_driver,'P.Value.G4.Vs.others_DA'],2))
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
          annotation=NULL,annotation_cex=1.2,
          left_annotation='high in G4',right_annotation='high in others')


############ II.2: target network structure for each driver
use_driver <- driver_list[1]
#transfer_tab <- get_IDtransfer2symbol2type(from_type = 'ensembl_transcript_id',use_genes=use_genes) ## get transfer table !!!
edge_score <- analysis.par$merge.network$target_list[[use_driver]]$MI*sign(analysis.par$merge.network$target_list[[use_driver]]$spearman)
names(edge_score) <- analysis.par$merge.network$target_list[[use_driver]]$target
#
use_driver2 <- driver_list[3]
edge_score2 <- analysis.par$merge.network$target_list[[use_driver2]]$MI*sign(analysis.par$merge.network$target_list[[use_driver2]]$spearman)
names(edge_score2) <- analysis.par$merge.network$target_list[[use_driver2]]$target
#
draw.targetNet(source_label=ms_tab[use_driver,'gene_label'],source_z=ms_tab[use_driver,'Z.G4.Vs.others_DA'],
               edge_score = edge_score,pdf_file=sprintf('%s/targetNet_out.pdf',analysis.par$out.dir.PLOT),label_cex = 0.2)

draw.targetNet(source_label=ms_tab[use_driver,'gene_label'],source_z=ms_tab[use_driver,'Z.G4.Vs.others_DA'],
               edge_score = edge_score,pdf_file=sprintf('%s/targetNet_in.pdf',analysis.par$out.dir.PLOT),label_cex = 0.2,arrow_direction = 'in')

# for two target
use_genes <- unique(analysis.par$merge.network$network_dat$target.symbol)
draw.targetNet.TWO(source1_label=ms_tab[use_driver,'gene_label'],edge_score1 = edge_score,
                   source2_label=ms_tab[use_driver2,'gene_label'],edge_score2 = edge_score2,
                   source1_z=ms_tab[use_driver,'Z.G4.Vs.others_DA'],source2_z=ms_tab[use_driver2,'Z.G4.Vs.others_DA'],
                   pdf_file=sprintf('%s/targetNetTWO.pdf',analysis.par$out.dir.PLOT),total_possible_target=use_genes,show_test=TRUE,label_cex = 0.2)

# or directly get the test result
test.targetNet.overlap(source1_label=ms_tab[use_driver,'gene_label'],source2_label=ms_tab[use_driver2,'gene_label'],
                       target1 = names(edge_score),target2 = names(edge_score2),total_possible_target=use_genes)

############ II.3: category plot for the expression/activity for each driver
use_driver <- driver_list[1]
exp_mat <- exprs(analysis.par$cal.eset) ## expression,the rownames must be the originalID
ac_mat <- exprs(analysis.par$merge.ac.eset) ## ac,the rownames must be the originalID_label
phe_info <- pData(analysis.par$cal.eset)
use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
#
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=exp_mat[ms_tab[use_driver,'originalID'],],use_obs_class=use_obs_class,
                   class_order=c('WNT','SHH','G4'),class_srt=30,main_ac = ms_tab[use_driver,'gene_label'],main_exp=ms_tab[use_driver,'geneSymbol'],
                   pdf_file=sprintf('%s/categoryValue_demo1.pdf',analysis.par$out.dir.PLOT))
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=NULL,use_obs_class=use_obs_class,class_order=c('WNT','SHH','G4'),
                   pdf_file=sprintf('%s/categoryValue_demo2.pdf',analysis.par$out.dir.PLOT))

## function enrichment plot for the target genes for the driver, similar as above

################# PartIII: advanced part

############## III.1: gene set-based activity analysis, including vocalno, heatmap, category and GSEA plot
db.preload(use_level='gene',use_spe='human',update=FALSE)
gs.preload()
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
DA_gs <- getDE.limma.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,G1_name='G4',G0_name='others')
#DA_gs_bid <- getDE.BID.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,G1_name='G4',G0_name='others')

## draw vocalno plot for top sig-GS
sig_gs <- draw.volcanoPlot(dat=DA_gs,label_col='ID',logFC_col='logFC',
                               Pv_col='P.Value',logFC_thre=0.25,Pv_thre=1e-4,
                               main='Volcano Plot for gene sets',show_label=TRUE,label_type = 'distribute',label_cex = 0.5,
                               pdf_file=sprintf('%s/vocalno_GS_DA.pdf',analysis.par$out.dir.PLOT))

## draw heatmap for top sig-GS
draw.heatmap(mat=ac_gs[sig_gs$ID,],pdf_file=sprintf('%s/heatmap_GS.pdf',analysis.par$out.dir.PLOT),scale='row',
             phenotype_info=phe_info,use_phe=c('gender','subgroup'))

## draw GSEA plot for top sig-GS
comp <- 'G4.Vs.others'
DE <- analysis.par$DE[[comp]]
#
draw.GSEA.NetBID.GS(DE=DE,name_col='ID',profile_col='t',profile_trend='pos2neg',
                 sig_gs_list = sig_gs$ID,
                 gs_DA_Z=DA_gs[sig_gs$ID,'Z-statistics'],
                 use_gs2gene = use_gs2gene,
                 top_gs_number=20,target_col='RdBu',
                 left_annotation = 'test_left',right_annotation = 'test_right',
                 main='test',Z_sig_thre=1.64,profile_sig_thre = 0,
                 pdf_file=sprintf('%s/NetBID_GSEA_GS_demo1.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.GS(DE=DE,name_col='ID',profile_col='t',profile_trend='pos2neg',
                    sig_gs_list = sig_gs$ID,
                    gs_DA_Z=DA_gs[sig_gs$ID,'Z-statistics'],
                    use_gs2gene = use_gs2gene,
                    top_gs_number=20,target_col='black',
                    left_annotation = 'test_left',right_annotation = 'test_right',
                    main='test',Z_sig_thre=1.64,profile_sig_thre = 0,
                    pdf_file=sprintf('%s/NetBID_GSEA_GS_demo2.pdf',analysis.par$out.dir.PLOT))

## draw category plot for each sig-GS
use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
draw.categoryValue(ac_val=ac_gs[sig_gs$ID[1],],use_obs_class=use_obs_class,
                   class_order=c('WNT','SHH','G4'),class_srt=30,pdf_file=sprintf('%s/categoryValue_GS_demo1.pdf',analysis.par$out.dir.PLOT),
                   main_ac=sig_gs$ID[1],main_cex=0.8)

## draw GSEA plot for each sig-GS
DE_profile <- DE$`Z-statistics`; names(DE_profile) <- rownames(DE)
use_target_genes <- rownames(DE)[which(DE$ID %in% use_gs2gene[[sig_gs$ID[1]]])]
draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,
          main=sprintf('GSEA plot for %s',sig_gs$ID[1]),
          pdf_file = sprintf('%s/GSEA_GS_each.pdf',analysis.par$out.dir.PLOT),
          left_annotation='high in G4',right_annotation='high in others')


############## III.2: SINBA plot for synergistic effect

seed_driver <- sig_driver[which(sig_driver$logFC.G4.Vs.others_DA>0)[1],1]
part_driver <- ms_tab$originalID_label

######
comp_name <- 'G4.Vs.others' ## each comparison must give a name !!!
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # get sample list for G0
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # get sample list for G1
##
merge_target <- lapply(part_driver,function(x){
  m1 <- merge_target_list(driver1=seed_driver,driver2=x,target_list=analysis.par$merge.network$target_list)
})
names(merge_target) <- part_driver
ac_combine_mat <- cal.Activity(target_list=merge_target,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
DA_driver_combine <- getDE.BID.2G(eset=generate.eset(ac_combine_mat),G1=G1,G0=G0,G1_name='G4',G0_name='others')
ori_part_Z <- analysis.par$DA[[comp_name]][part_driver,'Z-statistics']
ori_seed_Z <- analysis.par$DA[[comp_name]][seed_driver,'Z-statistics']
#diff_Z <- 2*DA_driver_combine[part_driver,'Z-statistics']-(ori_part_Z+ori_seed_Z)
#names(diff_Z) <- part_driver
comp <- 'G4.Vs.others'
DE <- analysis.par$DE[[comp]]
driver_DA_Z <- analysis.par$DA[[comp_name]][,'Z-statistics']
names(driver_DA_Z) <- rownames(analysis.par$DA[[comp_name]])
driver_DE_Z <- analysis.par$DE[[comp_name]][,'Z-statistics']
names(driver_DE_Z) <- rownames(analysis.par$DE[[comp_name]])
DA_Z_merge <- DA_driver_combine[,'Z-statistics']
names(DA_Z_merge) <- rownames(DA_driver_combine)
target_list_merge <- merge_target
seed_driver_label <- ms_tab[seed_driver,'gene_label']
partner_driver_list <- part_driver
profile_col <- 't'
partner_driver_label <- ms_tab[partner_driver_list,'gene_label']
target_list <- analysis.par$merge.network$target_list
##
draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=20,profile_trend='pos2neg',top_order='merge',Z_sig_thre = 1.64,
                       target_nrow=1,target_col='RdBu',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo1.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=15,profile_trend='pos2neg',top_order='merge',Z_sig_thre = 1.64,
                       target_nrow=2,target_col='RdBu',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo2.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=15,profile_trend='pos2neg',top_order='diff',Z_sig_thre = 1.64,
                       target_nrow=2,target_col='RdBu',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo3.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.SINBA(DE=DE,profile_col = profile_col,seed_driver=seed_driver,partner_driver_list=partner_driver_list,
                       seed_driver_label=seed_driver_label,partner_driver_label=partner_driver_label,
                       driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,target_list=target_list,
                       DA_Z_merge=DA_Z_merge,target_list_merge=target_list_merge,
                       top_driver_number=5,profile_trend='pos2neg',top_order='merge',Z_sig_thre = 1.64,
                       target_nrow=1,target_col='black',target_col_type='PN',
                       pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo4.pdf',analysis.par$out.dir.PLOT))




