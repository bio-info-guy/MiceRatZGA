require(GenomicFeatures)
mouse_gtf <- makeTxDbFromGFF('./dataset/mouse_data/mouse_spike.gtf', 'gtf')
rat_gtf <- makeTxDbFromGFF('./dataset/rat_data/rat_spike.gtf', 'gtf')
tx2gene_mouse <- AnnotationDbi::select(mouse_gtf, keys(mouse_gtf, keytype = "TXNAME"), "GENEID", "TXNAME")
tx2gene_rat <- AnnotationDbi::select(rat_gtf, keys(rat_gtf, keytype = "TXNAME"), "GENEID", "TXNAME")
write.table(tx2gene_mouse, './dataset/mouse_data/tx2gene.tsv', sep = '\t', quote = F)
write.table(tx2gene_rat, './dataset/rat_data/tx2gene.tsv', sep = '\t', quote = F)

mouse_star_lenient <- read.csv('./dataset/mouse_data/mouse_lenient.star.ct.txt', sep ='\t', header = T, row.names = 1)
mouse_star_lenient_3 <- mouse_star_lenient[,grepl('_3_', colnames(mouse_star_lenient))]
colnames(mouse_star_lenient_3 ) <- sapply(strsplit(colnames(mouse_star_lenient_3 ), '_'), FUN= function(x){x[1]})
mouse_star_lenient_r3 <- mouse_star_lenient[,grepl('_r3_', colnames(mouse_star_lenient))]
colnames(mouse_star_lenient_r3 ) <- sapply(strsplit(colnames(mouse_star_lenient_r3 ), '_'), FUN= function(x){x[1]})

mouse_0 <- prepareMeta(list(), './mouse_data/', QC = F)
mouse_0 <- prepareCount(mouse_0, dirc = './mouse_data/')
mouse_0 <- prepareGeneFeatures(mouse_0, './dataset/mouse_data/mouse_spike.gtf')
mouse_0$meta$mito_rate <- colSums(mouse_0$tpm$bio[row.names(subset(mouse_0$features_tpm, seqnames == 'NC_005089.1')),row.names(mouse_0$meta)])/colSums(mouse_0$tpm$bio[,row.names(mouse_0$meta)])
mouse_0$meta$rRNA_rate <- colSums(mouse_0$tpm$bio[row.names(subset(mouse_0$features_tpm, gene_biotype == 'rRNA')),row.names(mouse_0$meta)])/colSums(mouse_0$tpm$bio)[row.names(mouse_0$meta)]
rat_0 <- prepareMeta(list(), './rat_data/', QC = F)
rat_0 <- prepareCount(rat_0, dirc = './rat_data/')
rat_0 <- prepareGeneFeatures(rat_0, './dataset/rat_data/rat_spike.gtf')
rat_0$meta$mito_rate <- colSums(rat_0$tpm$bio[row.names(subset(rat_0$features_tpm, seqnames == 'NC_001665.2')),row.names(rat_0$meta)])/colSums(rat_0$tpm$bio)[row.names(rat_0$meta)]
rat_0$meta$rRNA_rate <- colSums(rat_0$tpm$bio[row.names(subset(rat_0$features_tpm, gene_biotype == 'rRNA')),row.names(rat_0$meta)])/colSums(rat_0$tpm$bio)[row.names(rat_0$meta)]

tableS1 <- rbind(mouse_0$meta, rat_0$meta)
mouse_map_info <- read.csv('./dataset/mouse_data/mapping.info', sep = '\t', row.names = 1)
rat_map_info <- read.csv('./dataset/rat_data/mapping.info', sep = '\t', row.names = 1)
map_info <- rbind(mouse_map_info, rat_map_info)
tableS1$HISAT_Unique_Mapping_Rate <- map_info[row.names(tableS1), "hisat_uniq_rate"]
tableS1$SALMON_Unique_Mapping_Rate <- map_info[row.names(tableS1), "salmon_rate"]
tableS1$TOTAL_Reads <- map_info[row.names(tableS1), "total_reads"]
write.csv(tableS1, './table_S1.csv')
options(ucscChromosomeNames=FALSE)

# Seperate Conditions
mouse_0 <- prepareMeta(list(), './mouse_data/', QC = F, filter_samples = c('MU18', 'MU28', 'MU32', 'MU27', 'MU23', 'MU34', 'MU26', 'MU3'))

mouse_0$meta[c('MU26', 'MU3'),]$cellType <- 'mouseZygote'
mouse_0 <- prepareCount(mouse_0, dirc = './mouse_data/')
mouse_0 <- prepareGeneFeatures(mouse_0, './dataset/mouse_data/mouse_spike.gtf')
mouse_0$meta$rRNA_rate <- colSums(mouse_0$tpm$bio[row.names(subset(mouse_0$features_tpm, gene_biotype == 'rRNA')),row.names(mouse_0$meta)])/colSums(mouse_0$tpm$bio)[row.names(mouse_0$meta)]
mouse_0$meta$rRNA_rate2 <- colSums(mouse_0$tpm$bio[c('Gm26917', 'Lars2'),row.names(mouse_0$meta)])*2/colSums(mouse_0$tpm$bio[,row.names(mouse_0$meta)])
mouse_0$meta$mito_rate <- colSums(mouse_0$tpm$bio[row.names(subset(mouse_0$features_tpm, seqnames == 'NC_005089.1')),row.names(mouse_0$meta)])/colSums(mouse_0$tpm$bio[,row.names(mouse_0$meta)])

mouse_0$meta$rRNA_rate_ct <- colSums(mouse_0$ct$bio[row.names(subset(mouse_0$features_ct, gene_biotype == 'rRNA')),row.names(mouse_0$meta)])/colSums(mouse_0$ct$bio)[row.names(mouse_0$meta)]
mouse_0$meta$rRNA_rate_ct_star <- colSums(mouse_0$ct$bio_star[row.names(subset(mouse_0$features_ct, gene_biotype == 'rRNA')),row.names(mouse_0$meta)])/colSums(mouse_0$ct$bio_star)[row.names(mouse_0$meta)]

mouse_0$meta$rRNA_rate2_ct <- colSums(mouse_0$ct$bio[c('Gm26917', 'Lars2'),row.names(mouse_0$meta)])*2/colSums(mouse_0$ct$bio[,row.names(mouse_0$meta)])
mouse_0$meta$rRNA_rate2_ct_star <- colSums(mouse_0$ct$bio_star[c('Gm26917', 'Lars2'),row.names(mouse_0$meta)])*2/colSums(mouse_0$ct$bio_star[,row.names(mouse_0$meta)])

mouse_0$meta$mito_rate_ct <- colSums(mouse_0$ct$bio[row.names(subset(mouse_0$features_ct, seqnames == 'NC_005089.1')),row.names(mouse_0$meta)])/colSums(mouse_0$ct$bio[,row.names(mouse_0$meta)])
mouse_0$meta$mito_rate_ct_star <- colSums(mouse_0$ct$bio_star[row.names(subset(mouse_0$features_ct, seqnames == 'NC_005089.1')),row.names(mouse_0$meta)])/colSums(mouse_0$ct$bio_star[,row.names(mouse_0$meta)])




mouse_0$tpm$bio <- mouse_0$tpm$bio[row.names(subset(mouse_0$features_tpm, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1' & !gene_id %in% c('Gm26917', 'Lars2'))),]
mouse_0$tpm$bio <- t(1e6*t(mouse_0$tpm$bio)/colSums(mouse_0$tpm$bio))
mouse_0$ct$bio <- mouse_0$ct$bio[row.names(subset(mouse_0$features_ct, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1' & !gene_id %in% c('Gm26917', 'Lars2'))),]
mouse_0$ct$bio_star <- mouse_0$ct$bio_star[row.names(subset(mouse_0$features_ct, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1' & !gene_id %in% c('Gm26917', 'Lars2'))),]


sample_PCA(DESeq2::varianceStabilizingTransformation(as.matrix(mouse_0$ct$bio)), mouse_0$meta, umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = T, point_size = 2, color_by = 'cellType', umap.config = umap_config)
sample_PCA(log(mouse_0$tpm$bio+1), mouse_0$meta, umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = T, point_size = 2, color_by = 'cellType', umap.config = umap_config)




mouse <- prepareMeta(list(), './mouse_data/', QC = F, filter_samples = c('MU18', 'MU28', 'MU32', 'MU27', 'MU23', 'MU34'))
mouse$meta[c('MU26', 'MU3'),]$cellType <- 'mouseZygote'
mouse$meta$mito_rate <- colSums(mouse$tpm$bio[row.names(subset(mouse$features_tpm, seqnames == 'NC_005089.1')),row.names(mouse$meta)])/colSums(mouse$tpm$bio[,row.names(mouse$meta)])
mouse$meta$rRNA_rate <- colSums(mouse$tpm$bio[c('Gm26917', 'Lars2'),row.names(mouse$meta)])/colSums(mouse$tpm$bio[,row.names(mouse$meta)])
mouse_pc_lncRNA_tpseudo_nonmito_nonrRNA <- row.names(subset(mouse$features_tpm, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1'))
mouse$tpm$bio <- mouse$tpm$bio[mouse_pc_lncRNA_tpseudo_nonmito_nonrRNA,]
mouse$tpm$bio <- t(1e6*t(mouse$tpm$bio)/colSums(mouse$tpm$bio))
mouse$ct$bio <- mouse$ct$bio[row.names(subset(mouse$features_ct, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1')),]
mouse$genes_subset <- intersect(row.names(mouse$tpm$bio[rowMeans(as.matrix(mouse$tpm$bio)) > 2,]), row.names(mouse$ct$bio[rowMeans(as.matrix(mouse$ct$bio)) > 10,]))

#mouse_filtered <- mouse
#mouse_filtered$tpm$bio <- mouse_filtered$tpm$bio[row.names(mouse$tpm$bio[rowMeans(as.matrix(mouse$tpm$bio)) > 2,]),]
#mouse_filtered$ct$bio <- mouse_filtered$ct$bio[mouse_filtered$genes_subset,]


rat <- list()
rat <- prepareMeta(rat, './rat_data/', QC = F)
rat <- prepareMeta(rat, './rat_data/', QC = F,filter_samples = c('RU14','RF11','RF13', 'RF7', 'RF16', 'RF8'))
#rat$meta <- MeanInGroupCorrelation(rat$meta, mat =  rat$bio.reads[,row.names(rat$meta)], log = T)
rat <- prepareCount(rat, dirc = './rat_data/')
rat <- prepareGeneFeatures(rat, gtf_file = './dataset/rat_data/rat_spike.gtf')
rat$meta$mito_rate <- colSums(rat$tpm$bio[row.names(subset(rat$features_tpm, seqnames == 'NC_001665.2')),row.names(rat$meta)])/colSums(rat$tpm$bio)[row.names(rat$meta)]
rat$meta$rRNA_rate <- colSums(rat$tpm$bio[row.names(subset(rat$features_tpm, gene_biotype == 'rRNA')),row.names(rat$meta)])/colSums(rat$tpm$bio)[row.names(rat$meta)]
rat$tpm$bio <- rat$tpm$bio[row.names(subset(rat$features_tpm, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1')),]
rat$tpm$bio <- t(1e6*t(rat$tpm$bio)/colSums(rat$tpm$bio))
rat$ct$bio <- rat$ct$bio[row.names(subset(rat$features_ct, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1')),]
rat$genes_subset <- intersect(row.names(rat$tpm$bio[rowMeans(as.matrix(rat$tpm$bio)) > 2,]), row.names(rat$ct$bio[rowMeans(as.matrix(rat$ct$bio)) > 10,]))


rat$meta$rRNA_rate_ct <- colSums(rat$ct$bio[row.names(subset(rat$features_ct, gene_biotype == 'rRNA')),row.names(rat$meta)])/colSums(rat$ct$bio)[row.names(rat$meta)]
rat$meta$rRNA_rate_ct_star <- colSums(rat$ct$bio_star[row.names(subset(rat$features_ct, gene_biotype == 'rRNA')),row.names(rat$meta)])/colSums(rat$ct$bio_star)[row.names(rat$meta)]

rat$meta$rRNA_rate2_ct <- colSums(rat$ct$bio[c('Gm26917', 'Lars2'),row.names(rat$meta)])*2/colSums(rat$ct$bio[,row.names(rat$meta)])
rat$meta$rRNA_rate2_ct_star <- colSums(rat$ct$bio_star[c('Gm26917', 'Lars2'),row.names(rat$meta)])*2/colSums(rat$ct$bio_star[,row.names(rat$meta)])

rat$meta$mito_rate_ct <- colSums(rat$ct$bio[row.names(subset(rat$features_ct, seqnames == 'NC_001665.2')),row.names(rat$meta)])/colSums(rat$ct$bio[,row.names(rat$meta)])
rat$meta$mito_rate_ct_star <- colSums(rat$ct$bio_star[row.names(subset(rat$features_ct, seqnames == 'NC_005089.1')),row.names(rat$meta)])/colSums(rat$ct$bio_star[,row.names(rat$meta)])



expr_genes_summary <- data.frame(num_of_Genes = c(colSums(mouse$ct$bio > 0), colSums(rat$ct$bio > 0)),
                          Species = c(rep('Mouse', 28), rep('Rat', 25)),
                          cellType = c(mouse$meta$cellType, rat$meta$cellType)) %>% group_by(cellType) %>% summarise(Num_of_Genes = mean(num_of_Genes), sd_num = sd(num_of_Genes), Species = unique(Species))

ggplot(expr_genes_summary, aes(Species, Num_of_Genes, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Num_of_Genes - sd_num, ymax = Num_of_Genes + sd_num), color = '#ffa500',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+theme(plot.title = element_text(hjust = 0.5, size = 12),axis.title.y = element_text(size = 14) ,axis.title.x = element_text(size = 14))+ylab('Number of Genes Detected')
ggsave('number_gene_detected.png', width = 5,height = 4)




#rat_kb <- rat
#rat_filtered$tpm$bio <- rat_filtered$tpm$bio[row.names(rat$tpm$bio[rowMeans(as.matrix(rat$tpm$bio)) > 2,]),]
#rat_filtered$ct$bio <- rat_filtered$ct$bio[rat_filtered$genes_subset,]



# Plots
setwd('./results//Exploratory/')

## Make Embedding Plots
embed_vis <- list()
embed_vis[['mouse']] <- makePCA_UMAP(mouse, cellTypes =unique(mouse$meta$cellType), output_dir = './exploratory/mouse', counts = T)
embed_vis[['rat']] <- makePCA_UMAP(rat, cellTypes=unique(rat$meta$cellType), output_dir = './exploratory/rat', counts = T)



embed_vis[['mouse']][['hmap']] <- makeHmap(mouse, output_dir = './exploratory/mouse_', legend = F)
embed_vis[['rat']][['hmap']] <- makeHmap(rat, output_dir = './exploratory/rat_', legend = F)

dev.off()
vis <- list()
vis[['hmap']][[1]] <- NULL

for(i in 1:length(embed_vis)){
  vis[['pca']][[i]] <- embed_vis[[i]][['pca']]+theme(
    plot.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(), 
    legend.position = "none", 
    axis.ticks.length=unit(-0.1, "cm"), 
    axis.title.x = element_text(size=17, vjust=-0, color="black"), 
    axis.title.y = element_text(size=17, vjust=2, hjust = -0.1, color="black"))
  vis[['umap']][[i]] <- embed_vis[[i]][['umap']]+theme(
    plot.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(), 
    legend.position = "none", 
    axis.ticks.length=unit(-0.1, "cm"), 
    axis.title.x = element_text(size=16, vjust=-0, color="black"), 
    axis.title.y = element_text(size=16, vjust=2, color="black"))
  if(i < length(embed_vis)){
  vis[['hmap']][[i+1]] <- remove_heatmap_element(embed_vis[[i]][['hmap']])$gtable
  }else{
    vis[['hmap']][[i+1]] <- embed_vis[[i]][['hmap']]$gtable
  }
}

cowplot::plot_grid(cowplot::plot_grid(plotlist = vis[['pca']], ncol = 3), cowplot::plot_grid(plotlist = vis[['umap']], ncol = 3), cowplot::plot_grid(plotlist = vis[['hmap']], ncol = 4, rel_widths = c(0.12,1,1,1)), nrow = 3, rel_heights = c( 1, 1, 1.3), ncol =1, align = 'hv', hjust = -0.3 , vjust = c(1.3, 1, 2), labels = LETTERS[1:3])
ggsave(filename = './results/Exploratory/mouse_rat_vis.png', width = 8, height = 8, units = 'in', dpi = 300)


"
mast_res_10 <- mast_diff(mouse, mast_res, plot = T, tpm = T, nbins = 0, min_per_bin = 50, freq = 0.99, max_thres = 10)
mast_res_8 <- mast_diff(mouse, mast_res, plot = T, tpm = T, nbins = 0, min_per_bin = 50, freq = 0.8, max_thres = 10)
mast_res_5 <- mast_diff(mouse, mast_res, plot = T, tpm = T, nbins = 0, min_per_bin = 50, freq = 0.5, max_thres = 6)
mast_res_2 <- mast_diff(mouse, mast_res, plot = T, tpm = T, nbins = 0, min_per_bin = 50, freq = 0.2, max_thres = 6)
mast_res_1 <- mast_diff(mouse, mast_res, plot = T, tpm = T, nbins = 30, min_per_bin = 50, freq = 0.1, max_thres = 6, min_cell_grp = 0, min_cell = 0)
mast_gse_10 <- gse_CP(row.names(subset(mast_res_10$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
logFC = sapply(row.names(mast_res_10$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res_10$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)
mast_gse_8 <- gse_CP(row.names(subset(mast_res_8$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                      logFC = sapply(row.names(mast_res_8$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res_8$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)
mast_gse_5 <- gse_CP(row.names(subset(mast_res_5$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                      logFC = sapply(row.names(mast_res_5$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res_5$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)
mast_gse_2 <- gse_CP(row.names(subset(mast_res_2$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                      logFC = sapply(row.names(mast_res_2$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res_2$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)
mast_gse_1 <- gse_CP(row.names(subset(mast_res_1$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                      logFC = sapply(row.names(mast_res_1$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res_1$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)
scde_res <- deg_scde(ct = mouse_spliced, meta = mouse$meta, control = 'mouseEgg', scde.model = T, min_cell_grp = 5, min_cell = 8, min.pct = 0.1)
scde_gse <- gse_CP(row.names(scde_res),
                   logFC = sapply(row.names(scde_res), FUN =function(x){scde_res[x,'Log2FC']}), organisms = 'mouse', GSE = T)
mast_res_kb <- mast_diff(ct= mouse_spliced, meta = mouse$meta, plot = T, tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, max_thres = log2(1000), min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F)
mast_gse_kb <- gse_CP(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                   logFC = sapply(row.names(mast_res_kb$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res_kb$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)

"
mast_res_MU26 <- mast_diff(mouse_test, plot = T, tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, max_thres = log2(1000), min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F, correct_wild_coef = T)
mast_ora_MU26_up <- enrich_CP(row.names(subset(mast_res_MU26$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 'mouse', universe = row.names(mast_res_MU26$mouseEgg_v_mouseZygote$DESig))
mast_ora_MU26_down <- enrich_CP(row.names(subset(mast_res_MU26$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)), 'mouse', universe = row.names(mast_res_MU26$mouseEgg_v_mouseZygote$DESig))


mast_res_star <- mast_diff(ct = mouse_0$ct$bio_star, meta = mouse_0$meta, plot = T, tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, max_thres = log2(1000), min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F, correct_wild_coef = T)
mast_ora_star_up <- enrich_CP(row.names(subset(mast_res_star$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 'mouse', universe = row.names(mast_res_star$mouseEgg_v_mouseZygote$DESig))
mast_ora_star_down <- enrich_CP(row.names(subset(mast_res_star$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)), 'mouse', universe = row.names(mast_res_star$mouseEgg_v_mouseZygote$DESig))

mast_gse <- gse_CP(row.names(subset(mast_res_star$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                   logFC = sapply(row.names(mast_res_star$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res_star$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)


mast_res <- mast_diff(mouse, plot = T, tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, max_thres = log2(1000), min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F, correct_wild_coef = T)

mast_gse <- gse_CP(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                     logFC = sapply(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)

mast_ora_mouse_up <- enrich_CP(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 'mouse', universe = row.names(mast_res$mouseEgg_v_mouseZygote$DESig))
mast_ora_mouse_down <- enrich_CP(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)), 'mouse', universe = row.names(mast_res$mouseEgg_v_mouseZygote$DESig))



# representative pathways from mast_gse:
# 'R-MMU-176408', 'GO:0140013', 'GO:0051028', 'R-MMU-453276', "R-MMU-194441","GO:0051123","	R-MMU-2980766","R-MMU-5693532","GO:0006086","R-MMU-69306"
#'R-MMU-3247509', 'GO:0002181', 'R-MMU-927802','GO:0006090','R-MMU-6791226','GO:0050684','GO:0043433', 'GO:0060148','R-MMU-111933','GO:0042448','GO:0050821','R-MMU-1428517','WP403', 'GO:0043547', 'GO:0043406','GO:0060070','GO:0061614','GO:0007565','mmu04064','WP265','WP113','GO:0051017'


"
mast_rat_10 <- mast_diff(rat, mast_res, control = 'ratEgg',  plot = T ,tpm = T, nbins = 0, min_per_bin = 50, freq = 0.99, max_thres = 10)
mast_rat_8 <- mast_diff(rat, mast_res, control = 'ratEgg',  plot = T ,tpm = T, nbins = 0, min_per_bin = 50, freq = 0.8, max_thres = 10)
mast_rat_5 <- mast_diff(rat, mast_res, control = 'ratEgg',  plot = T ,tpm = T, nbins = 0, min_per_bin = 50, freq = 0.5, max_thres = 10)
mast_rat_2 <- mast_diff(rat, mast_res, control = 'ratEgg',  plot = T ,tpm = T, nbins = 0, min_per_bin = 50, freq = 0.2, max_thres = 10)
mast_rat_1 <- mast_diff(rat, mast_res, control = 'ratEgg',  plot = T ,tpm = T, nbins = 0, min_per_bin = 50, freq = 0.1, max_thres = 10)
mast_gse_rat_10 <- gse_CP(row.names(subset(mast_rat_10$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                         logFC = sapply(row.names(mast_rat_10$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat_10$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)
mast_gse_rat_8 <- gse_CP(row.names(subset(mast_rat_8$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                     logFC = sapply(row.names(mast_rat_8$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat_8$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)
mast_gse_rat_5 <- gse_CP(row.names(subset(mast_rat_5$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                         logFC = sapply(row.names(mast_rat_5$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat_5$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)
mast_gse_rat_2 <- gse_CP(row.names(subset(mast_rat_2$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                         logFC = sapply(row.names(mast_rat_2$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat_2$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)
mast_gse_rat_1 <- gse_CP(row.names(subset(mast_rat_1$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                         logFC = sapply(row.names(mast_rat_1$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat_1$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)
scde_rat_kb <- deg_scde(ct = round(rat_spliced), meta = rat$meta, control = 'ratEgg', scde.model = T, min_cell_grp = 5, min_cell = 8, min.pct = 0.1)
scde_rat_hisat <- deg_scde(rat, control = 'ratEgg', scde.model = T, min_cell_grp = 5, min_cell = 8, min.pct = 0.1)
scde_gse_rat <- gse_CP(row.names(scde_rat),
                       logFC = sapply(row.names(scde_rat), FUN =function(x){scde_rat[x,'Log2FC']}), organisms = 'rat', GSE = T)
"
#mast_rat_kb <- mast_diff(ct = rat_spliced, meta=rat$meta, control = 'ratEgg',  plot = T ,tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, max_thres = log2(1000), min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F, correct_wild_coef = T)

mast_rat <- mast_diff(obj = rat, control = 'ratEgg',  plot = T ,tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, max_thres = log2(1000), min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F, correct_wild_coef = T)

mast_gse_rat <- gse_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                         logFC = sapply(row.names(mast_rat$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)

#mast_gse_rat_kb <- gse_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                       #logFC = sapply(row.names(mast_rat_kb$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat_kb$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)


mast_ora_rat_up <- enrich_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 0.5)), 'rat', universe = row.names(mast_rat$ratEgg_v_ratZygote$DESig))
mast_ora_rat_down <- enrich_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -0.5)), 'rat', universe = row.names(mast_rat$ratEgg_v_ratZygote$DESig))



custom_ridgeplot(mast_gse$combined, terms = c('R-MMU-3247509', 'GO:0002181', 'R-MMU-927802','GO:0006090',
                                            'GO:0050684', 'GO:0016441','R-MMU-72312',
                                              'R-MMU-111933','GO:0042448','R-MMU-1428517',
                                              'WP403', 'GO:0043547', 'GO:0043405','GO:0060070',
                                              'GO:0061614','GO:0007565','GO:0038061','WP265','R-MMU-75953',
                                              'WP113','GO:0051017','R-MMU-69620', 'GO:0140013', 
                                              'GO:0051028', 'R-MMU-453276', "R-MMU-194441",
                                              "R-MMU-2980766","R-MMU-5693532","GO:0006086","R-MMU-69306"), top_n = 0)+theme(axis.text.y = element_text(face="bold", color="black", size=8), plot.title = element_text(hjust = 1, size = 12))+ggtitle('GSEA of Mouse DEGs')
ggsave('./deg/mouse_gsea_ridge.png', height = 8, width =4, units = 'in')
                                                                                                                                                                                                                                              
custom_ridgeplot(mast_gse_rat$combined,  terms = c('GO:0016581', 'GO:0090545', 'R-RNO-73728', 'R-RNO-212300', 'GO:0030527', 
     'rno03010', 'GO:0090092', 'rno00020','rno00280','R-RNO-927802',
    'GO:0004930','rno01200','R-RNO-204005'), top_n = 0)+theme(plot.title = element_text(hjust = 1, size = 12), axis.text.y = element_text(face="bold", color="black", size=8))+ggtitle('GSEA of Rat DEGs')
                                                                                                                                                                                                                                          
ggsave('./deg/rat_gsea_ridge.png', height = 8, width =4, units = 'in')

DEG_full <- list()

DEG_full <- performAllDEG(deg_res = DEG_full, obj = rat, control = 'ratEgg', METHODS = c( 'monocle', 'monocle_ct', 'limma', 'edgeR'), min_pct = 0.1)
DEG_full <- performAllDEG(deg_res = DEG_full, obj = mouse, METHODS = c( 'monocle', 'monocle_ct', 'limma', 'edgeR', 'wilcox' ), min_pct = 0.1)
DEG_GSE_full <- ClusterProfilerORAResult(DEG_full, mode = 'all', method = 'GSE')
#DEG_expr_filt <- list()

#DEG_expr_filt <- performAllDEG(deg_res = DEG_expr_filt, obj = rat_filtered, control = 'ratEgg', METHODS = c('scde', 'scdd',  'monocle', 'monocle_ct', 'limma', 'edgeR', 'basics', 'bpsc','wilcox', 'roc', 'LR', 'mast', 't', 'DESeq2', 'bimod'), FCThresh = 1)
#DEG_expr_filt <- performAllDEG(deg_res = DEG_expr_filt, obj = mouse_filtered, METHODS = c('scde', 'scdd',  'monocle', 'monocle_ct', 'limma', 'edgeR', 'basics', 'bpsc','wilcox', 'roc', 'LR', 'mast', 't', 'DESeq2', 'bimod'), FCThresh = 1)

#DEG_expr_filt <- performAllDEG(deg_res = DEG_expr_filt, obj = mouse_filtered, res_list = DEG_expr_filt$mouseEgg_v_mouseZygote, METHODS = c( 'mast'), FCThresh = 1)
#DEG_expr_filt <- performAllDEG(deg_res = DEG_expr_filt, obj = rat_filtered, res_list = DEG_expr_filt$ratEgg_v_ratZygote, control = 'ratEgg', METHODS = c( 'mast'), FCThresh = 1)



volcanoplot(mast_res$mouseEgg_v_mouseZygote$DESig, fc = log2(2))



volcanoplot(mast_res$mouseEgg_v_mouseZygote$DESig, pval_col = 'Pr..Chisq.',fc = log2(2), top_genes = c('Nup37', 'Nup54', 'Obox1', 'Obox2', 'Obox5', 'Obox7', 'Nup35', 'Rpl9', 
                                                                                                       'Rpl3', 'Rpl26', 'Rpl12', 'Rpl19', 'Rpl18', 'Rpl3','Mastl','Cdc20', 
                                                                                                       'Cnot7','Mad2l1', 'Ube2e1','Taf6', 'Taf9','Gtf2b', 'Gtf2a2', 'Ndufc1',
                                                                                                       'Ndufa2', 'Ndufb7','Atp5h','Atp1a1','Ndufa8','Atp5j','Ndufa3', 'Dlat', ))
ggsave('./deg/mouse_volcano_mast.png', width = 4, height = 4)
volcanoplot(mast_rat$ratEgg_v_ratZygote$DESig, pval_col = 'Pr..Chisq.',fc = log2(2), top_genes = c('Rps20', 'Rps2', 'Rps15', 'Rplp1','Hist2h3c2','H2ac1', 'H2ax', 'Gata3',
                                                                                                   'Sdhc','Aco2','Idh1','Sdhd','Pdha1','Dlst','Ndufab1','Ndufc1'))
ggsave('./deg/rat_volcano_mast.png', width = 4, height = 4)

#makeDEGmaps(mouse, res = DEG$mouseEgg_v_mouseZygote, main = 'Mouse Egg vs Zygote', output_dir = './deg/mouse_deg')
#makeDEGmaps(rat, res = DEG$ratEgg_v_ratZygote, main = 'Rat Egg vs Zygote', output_dir = './deg/rat_deg')

#makeDEGmaps(mouse, res = DEG$mouseEgg_v_mouseZygote, main = 'Mouse Egg vs Zygote', output_dir = './deg/mouse_deg', individual = F)
#makeDEGmaps(rat, res = DEG$ratEgg_v_ratZygote, main = 'Rat Egg vs Zygote', output_dir = './deg/rat_deg', individual = F)


#### As of Jan 16th, we use DEGs from MAST, input is CPM produced by HISAT2
#### Only protein-coding, transcribed_pseudogene and lncRNA are used
#### genes must be expressed in 10% of either zygote or oocyte samples
#### No other filtering is done through MAST
#### Genes are recorded in the vectors below


mouse_deg <- mast_res$mouseEgg_v_mouseZygote$DESig
rat_deg <- mast_rat$ratEgg_v_ratZygote$DESig

mouse_deg_up <- row.names(subset(mouse_deg, Log2FC > log2(2) & fdr < 0.05))
mouse_deg_down <- row.names(subset(mouse_deg, Log2FC < -log2(2) & fdr < 0.05))

rat_deg_up <- row.names(subset(rat_deg, Log2FC > log2(2) & fdr < 0.05))
rat_deg_down <- row.names(subset(rat_deg, Log2FC < -log2(2) & fdr < 0.05))


{" Emrichment of All Methods
#DEG_ORA <- list()
#DEG_ORA <- ClusterProfilerORAResult(DEG, mode = 'all', method = 'ORA')
#DEG_GSE <- ClusterProfilerORAResult(DEG_new, mode = 'all', method = 'GSE')

#DEG_ORA_new <- list()
#DEG_ORA_new <- ClusterProfilerORAResult(DEG_new, mode = 'all', method = 'ORA')


#DEG_ORA_expr_filt <- list()
#DEG_ORA_expr_filt <- ClusterProfilerORAResult(DEG_expr_filt, mode = 'all', method = 'ORA')
#DEG_ORA_expr_filt_pos <- list()
#DEG_ORA_expr_filt_pos <- ClusterProfilerORAResult(DEG_expr_filt, mode = 'positive', method = 'ORA')
#DEG_ORA_expr_filt_neg <- list()
#DEG_ORA_expr_filt_neg <- ClusterProfilerORAResult(DEG_expr_filt, mode = 'negative', method = 'ORA')


#DEG_GSE_full <- ClusterProfilerORAResult(DEG_full, mode = 'all', method = 'GSE')
#DEG_GSE_full_pval_sign <- ClusterProfilerORAResult(DEG_full, mode = 'all', method = 'GSE', gsea_rank = 'pval_sign')
#DEG_GSE_full_pval_logfc <- ClusterProfilerORAResult(DEG_full, mode = 'all', method = 'GSE', gsea_rank = 'pval_logfc')

#DEG_GSE_filt <- ClusterProfilerORAResult(DEG_expr_filt, mode = 'all', method = 'GSE')
#DEG_GSE_filt_pval_sign <- ClusterProfilerORAResult(DEG_expr_filt, mode = 'all', method = 'GSE', gsea_rank = 'pval_sign')
#DEG_GSE_filt_pval_logfc <- ClusterProfilerORAResult(DEG_expr_filt, mode = 'all', method = 'GSE', gsea_rank = 'pval_logfc')



### 
#mouse_deg_up_gs <- DEG_ORA_expr_filt_pos$mast_tpm$mouseEgg_v_mouseZygote
#mouse_deg_down_gs <- DEG_ORA_expr_filt_neg$mast_tpm$mouseEgg_v_mouseZygote


#rat_deg_up_gs <- DEG_ORA_expr_filt_pos$mast_tpm$ratEgg_v_ratZygote
#rat_deg_down_gs <- DEG_ORA_expr_filt_neg$mast_tpm$ratEgg_v_ratZygote




#CP_Plots <- clusterProfilerPlots(DEG_ORA_expr_filt, dir = './enrichment/')

#enrich_plot_grid <- function(plts, method = 'mast'){
  #cowplot::plot_grid(cowplot::plot_grid(plts$pcnet$isotonic_sorbitol_v_AA_starvation$WKP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()),
                #                        plts$pcnet$isotonic_sorbitol_v_Glucose_starvation$WKP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), ncol = 2),
                #  cowplot::plot_grid(plts$pcnet$isotonic_sorbitol_v_AA_starvation$KEGG_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), 
                #                     plts$pcnet$isotonic_sorbitol_v_Glucose_starvation$KEGG_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), ncol = 2),
                #  cowplot::plot_grid(plts$pcnet$isotonic_sorbitol_v_AA_starvation$GO_BP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), 
                #                   plts$pcnet$isotonic_sorbitol_v_Glucose_starvation$GO_BP_ora[[method]]+theme(legend.position = "none", plot.title = element_blank()), ncol = 2), nrow = 3, labels = c('A','B','C'), hjust = c(-0.2,-0.2, -0.2), label_fontface = 'bold', label_size = 18)
# ggsave(filename = paste('./results/Enrichment/yeast_',method,'_res.png', sep=''), width = 10, height = 15, units = 'in', dpi = 300)
#
"}

### kalisto-bustool estimated mature/nascent transcript counts


mouse_spliced <-  t(read.csv('./velocity/mouse_spliced.csv', header = T, row.names = 1))[,colnames(mouse$ct$bio)]

mouse_spliced <- mouse_spliced[intersect(row.names(mouse_spliced), row.names(mouse$ct$bio)),]


rat_spliced <-  t(read.csv('./velocity/rat_spliced.csv', header = T, row.names = 1))[,colnames(rat$ct$bio)]
rat_spliced <- rat_spliced[intersect(row.names(rat_spliced), row.names(rat$ct$bio)),]

#rat_spliced_tpm <- rat_spliced/rat$salmon$gene$length[row.names(rat_spliced), colnames(rat_spliced)]
#rat_spliced_tpm <- t(t(rat_spliced_tpm)/colSums(rat_spliced_tpm))*1e6


mouse_intron <- read_filter_splice('./velocity/mouse_spliced.csv', './velocity/mouse_unspliced.csv')
mouse_intron$genes <- intersect(row.names(mouse$ct$bio), mouse_intron$genes)
rat_intron <- read_filter_splice('./velocity/rat_spliced.csv', './velocity/rat_unspliced.csv')
rat_intron$genes <- intersect(row.names(rat$ct$bio), rat_intron$genes)

perc_mouse_genes_unspliced <- colSums(mouse_intron$unspliced[mouse_intron$genes,])/(colSums(mouse_intron$spliced)+colSums(mouse_intron$unspliced[mouse_intron$genes,]))

perc_rat_genes_unspliced <- colSums(rat_intron$unspliced[rat_intron$genes,])/(colSums(rat_intron$spliced )+colSums(rat_intron$unspliced[rat_intron$genes,]))


percentage_summary <- data.frame(percentages = c(perc_mouse_genes_unspliced, perc_rat_genes_unspliced),
           Species = c(rep('Mouse', 28), rep('Rat', 25)),
           cellType = c(mouse$meta$cellType, rat$meta$cellType)) %>% group_by(cellType) %>% summarise(Percentage = mean(percentages*100), sd_Perc = sd(percentages*100), Species = unique(Species))

ggplot(percentage_summary, aes(Species, Percentage, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Percentage - sd_Perc, ymax = Percentage + sd_Perc), color = '#ffa500',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+
  annotate("text", x = 1, y = 2.3, label = "ns", size = 8) +annotate("text", x = 2, y = 2.5, label = "**", size = 8)+theme(axis.title = element_text(size = 16), axis.text =element_text(size = 12), legend.position="none")
ggsave('percent_gene_unspliced.png', width = 3.5,height = 4)

num_mouse_genes_unspliced <- colSums(mouse_intron$unspliced[mouse_intron$genes,] > 0)#/(colSums(mouse_intron$spliced >0 ))#+colSums(mouse_intron$unspliced[mouse_intron$genes,]))
num_rat_genes_unspliced <- colSums(rat_intron$unspliced[rat_intron$genes,] > 0)#/(colSums(mouse_intron$spliced >0 ))#+colSums(mouse_intron$unspliced[mouse_intron$genes,]))


num_summary <- data.frame(num_of_Genes = c(num_mouse_genes_unspliced, num_rat_genes_unspliced),
                                      Species = c(rep('Mouse', 28), rep('Rat', 25)),
                                      cellType = c(mouse$meta$cellType, rat$meta$cellType)) %>% group_by(cellType) %>% summarise(Num_of_Genes = mean(num_of_Genes), sd_num = sd(num_of_Genes), Species = unique(Species))

ggplot(num_summary, aes(Species, Num_of_Genes, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Num_of_Genes - sd_num, ymax = Num_of_Genes + sd_num), color = '#ffa500',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+
  annotate("text", x = 1, y = 2900, label = "*", size = 5) +annotate("text", x = 2, y = 2900, label = "**", size = 5)+ggtitle('# Genes w Nascent Reads')+theme(plot.title = element_text(hjust = 0.5, size = 12))
ggsave('number_gene_unspliced.png', width = 3,height = 2.5)


#mouse_unsplic_full_mast <- mast_diff(ct = mouse_intron$unspliced[mouse_intron$genes,], meta = mouse$meta, control = 'mouseEgg', tpm = F, nbins = 10, min_per_bin = 50, freq = 0.1, max_thres = 6, plot = T, min_cell_grp = 5, min_cell = 8)[[1]]$DESig

#rat_unsplic_full_mast <- mast_diff(ct = rat_intron$unspliced[rat_intron$genes,], meta = rat$meta, control = 'ratEgg', tpm = F, nbins = 10, min_per_bin = 50, freq = 0.1, max_thres = 6, plot = T)[[1]]$DESig





#length(intersect(row.names(subset(rat_unsplic_mast, fdr < 0.05 & Log2FC >= 1)), row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC >= log2(2)))))

#length(intersect(row.names(subset(rat_unsplic_mast, fdr < 0.05 & Log2FC <= -1)), row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC <= -log2(2)))))


mouse_unsplic_mast <- mast_diff(ct = mouse_intron$unspliced[mouse_intron$genes,], meta = mouse$meta, normFactor = colMeans(mouse_intron$spliced/edgeR::cpm(mouse_intron$spliced), na.rm = T),control = 'mouseEgg', tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, min_cell_grp = 2, min_cell = 5, max_thres = 6, plot = T)[[1]]$DESig
mouse_unsplic_mast1 <- mast_diff(ct = mouse_intron$unspliced[mouse_intron$genes,], meta = mouse$meta,control = 'mouseEgg', tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, min_cell_grp = 2, min_cell = 5, max_thres = 6, plot = T)[[1]]$DESig

mouse_unsplic_mast_ora_up <- enrich_CP(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & Log2FC > 1)), 'mouse', universe = mouse_unsplic_mast$features)
mouse_unsplic_mast_ora_down <- enrich_CP(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & Log2FC < -1)), 'mouse', universe = mouse_unsplic_mast$features)
mouse_unsplic_mast_gsea <- gse_CP(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & abs(Log2FC) > 1)), logFC = sapply(row.names(mouse_unsplic_mast), FUN =function(x){mouse_unsplic_mast[x,'Log2FC']}), organisms = 'mouse', GSE = T)


mouse_unsplic_mast_ora_up <- enrich_CP(setdiff(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & Log2FC > 1)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05))), 'mouse', universe = mouse_unsplic_mast$features)

rat_unsplic_mast <- mast_diff(ct = rat_intron$unspliced[rat_intron$genes,], meta = rat$meta, normFactor = colMeans(rat_intron$spliced/edgeR::cpm(rat_intron$spliced), na.rm = T), control = 'ratEgg', tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, min_cell_grp = 2, min_cell = 5, max_thres = 6, plot = T)[[1]]$DESig

mast_gse_spliced <- gse_CP(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                   logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), row.names(subset(mouse_unsplic_mast, fdr < 0.05))), FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)


rat_unsplic_mast_ora_up <- enrich_CP(row.names(subset(rat_unsplic_mast, fdr < 0.05 & Log2FC > 1)), 'rat', universe = rat_unsplic_mast$features)
rat_unsplic_mast_ora_down <- enrich_CP(row.names(subset(rat_unsplic_mast, fdr < 0.05 & Log2FC < -1)), 'rat', universe = rat_unsplic_mast$features)
rat_unsplic_mast_gsea <- gse_CP(row.names(subset(rat_unsplic_mast, fdr < 0.05 & abs(Log2FC) > 1)), logFC = sapply(row.names(rat_unsplic_mast), FUN =function(x){rat_unsplic_mast[x,'Log2FC']}), organisms = 'rat', GSE = T)


length(intersect(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & Log2FC >= 1)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC >= log2(2)))))

length(intersect(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & Log2FC <= -1)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC <= -log2(2)))))


unsplic_splic_log2FC <- data.frame(x=c(mouse_unsplic_mast$Log2FC, rat_unsplic_mast$Log2FC), 
                                   y=c(mast_res$mouseEgg_v_mouseZygote$DESig[row.names(mouse_unsplic_mast),'Log2FC'], mast_rat$ratEgg_v_ratZygote$DESig[row.names(rat_unsplic_mast),'Log2FC']), 
                                   Species = c(rep('Mouse', length(mouse_unsplic_mast$Log2FC)), rep('Rat', length(rat_unsplic_mast$Log2FC))))


ggplot(data = unsplic_splic_log2FC, aes(x = x, y = y)) + 
  geom_point()+
  facet_wrap(~Species)+
  theme_classic()+
  ylab(expression('Mature Log'[2]*'FC'))+
  xlab(expression('Nascent Log'[2]*'FC'))+
  geom_text(data = data.frame(x = c(-2, -2), y = c(4,4), lab = paste('PCC:',c(round(cor(mouse_unsplic_mast$Log2FC, mast_res$mouseEgg_v_mouseZygote$DESig[row.names(mouse_unsplic_mast),'Log2FC'],use = 'na.or.complete', method = 'pearson'),3),
                                                              round(cor(rat_unsplic_mast$Log2FC, mast_rat$ratEgg_v_ratZygote$DESig[row.names(rat_unsplic_mast),'Log2FC'], use = 'na.or.complete',method = 'pearson'),3)), sep = ' '),
                                Species = c('Mouse', 'Rat')),mapping = aes(x = x, y = y, label = lab))
ggsave('nascent_mature_log2FC.png', height = 3, width =3.5)


  

mouse_intron_prop <- fisher_proportion_test(mouse_intron$unspliced[mouse_intron$genes,], mouse_intron$spliced[mouse_intron$genes,], mouse$meta$cellType)
rat_intron_prop <- fisher_proportion_test(rat_intron$unspliced[rat_intron$genes,], rat_intron$spliced[rat_intron$genes,], rat$meta$cellType, 'ratEgg')

length(intersect(row.names(subset(test, prop_diff > 0 & qvalue < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2)))))



mouse_unsplic_prop_ora_up <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'mouse', universe = row.names(mouse_intron_prop))
mouse_unsplic_prop_ora_down <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff < 0)), 'mouse', universe = row.names(mouse_intron_prop))
mouse_unsplic_prop_gsea <- gse_CP(row.names(subset(mouse_intron_prop, fdr < 0.05 & abs(Log2FC) > 1)), logFC = sapply(row.names(mouse_intron_prop), FUN =function(x){mouse_intron_prop[x,'prop_diff']}), organisms = 'mouse', GSE = T)


rat_unsplic_prop_ora_up <- enrich_CP(row.names(subset(rat_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'rat', universe = row.names(rat_intron_prop))
rat_unsplic_prop_ora_down <- enrich_CP(row.names(subset(rat_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'rat', universe = row.names(rat_intron_prop))
rat_unsplic_prop_gsea <- gse_CP(row.names(subset(rat_intron_prop, qvalue < 0.05 & prop_diff > 0)), logFC = sapply(row.names(rat_intron_prop), FUN =function(x){rat_intron_prop[x,'prop_diff']}), organisms = 'rat', GSE = T)

{"
unsplic_splic_log2FC_v_PropFC <- data.frame(x=c(mouse_intron_prop$prop_diff, rat_intron_prop$prop_diff), 
                                   y=c(mast_res$mouseEgg_v_mouseZygote$DESig[row.names(mouse_intron_prop),'Log2FC'], mast_rat$ratEgg_v_ratZygote$DESig[row.names(rat_intron_prop),'Log2FC']), 
                                   Species = c(rep('Mouse', length(mouse_intron_prop$prop_diff)), rep('Rat', length(rat_intron_prop$prop_diff))))


ggplot(data = unsplic_splic_log2FC_v_PropFC, aes(x = x, y = y)) + 
  geom_point()+
  facet_wrap(~Species)+
  theme_classic()+
  ylab(expression('Mature Log'[2]*'FC'))+
  xlab(expression('Nascent Log'[2]*'PropFC'))+
  geom_text(data = data.frame(x = c(-2, -2), y = c(4,4), lab = paste('PCC:',c(round(cor(mouse_intron_prop$prop_diff, mast_res$mouseEgg_v_mouseZygote$DESig[row.names(mouse_intron_prop),'Log2FC'],use = 'na.or.complete', method = 'pearson'),3),
                                                                              round(cor(rat_intron_prop$prop_diff, mast_rat$ratEgg_v_ratZygote$DESig[row.names(rat_intron_prop),'Log2FC'], use = 'na.or.complete',method = 'pearson'),3)), sep = ' '),
                              Species = c('Mouse', 'Rat')),mapping = aes(x = x, y = y, label = lab))

mast_gse_wo_unspliced <- gse_CP(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                   logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), mouse_intron$genes), FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)
"}

dnp_deg_mouse <- list('DEG Up'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'DEG Down'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'DNP Up'= row.names(subset(mouse_intron_prop, prop_diff > 0 & qvalue < 0.05)),
                      'DNP Down' =row.names(subset(mouse_intron_prop, prop_diff < 0 & qvalue < 0.05)))
make_comb_mat(dnp_deg_mouse)
png('dnp_deg_mouse.png', width = 3, height = 3, res = 300, units = 'in')
UpSet(make_comb_mat(dnp_deg_mouse)[1:4], comb_col = c('black'))
dev.off()


dne_deg_mouse <- list('Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'Mouse Down DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'Mouse Nascent Genes'= mouse_intron$genes)



dnp_deg_rat <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                    'Rat Down DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                    'Rat Increase DNP'= row.names(subset(rat_intron_prop, prop_diff > 0 & qvalue < 0.05)),
                    'Rat Decrease DNP' =row.names(subset(rat_intron_prop, prop_diff < 0 & qvalue < 0.05)))
png('dnp_deg_rat.png', width = 3, height = 3, res = 300, units = 'in')
UpSet(make_comb_mat(dnp_deg_rat)[1:4], comb_col = c('black'))
dev.off()

dne_deg_rat <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                    'Rat Down DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                    'Rat Nascent Genes'= rat_intron$genes)


sample_PCA(log(t(t(mouse_intron$unspliced[mouse_intron$genes,])/estimateSizeFactorsForMatrix(mouse_intron$spliced))+1), mouse$meta, umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = F, point_size = 2, color_by = 'cellType', umap.config = list(n_neighbors = 15, min_dist = 0.5, metric='pearson'))
ggsave(paste(output_dir,"_PCA_label.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')
sample_PCA(log(t(t(rat_intron$unspliced[rat_intron$genes,])/estimateSizeFactorsForMatrix(rat_intron$spliced))+1), rat$meta[colnames(rat_intron$spliced),], umap =  T, dimension=2 ,main = "UMAP of Log Expression", labeling = F, point_size = 2, color_by = 'cellType', umap.config = list(n_neighbors = 15, min_dist = 0.5, metric='pearson'))
ggsave(paste(output_dir,"_PCA_label.tiff", sep=''), units="in", width=5, height=4, dpi=300, compression = 'lzw')



### SUPPA2
{
## SUPPA2 results

## diffsplice parameters are abs(dPSI) > 0.2 and SUPPA2 pval < 0.05
### AS of Nov 28th, total of 5430 genes have diff splice events in rat data and 6275 genes have splice events in mouse
### total 358 sig diff splice in mouse, 82 genes sig diffsplice in rat
### events by type:
# mouse rat
#SE 220 36 
#MX 6 1
#AL 4 4
#AF 81  20
#RI 11  2
#A5 59  11
#A3 64  12

suppa_mouse_eve_res <- read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/suppa/diffsplice/per_local_event/MF-MU_local_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_mouse_eve_res$gene_short_name <- sapply(row.names(suppa_mouse_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_mouse_eve_res$event_type <- sapply(sapply(row.names(suppa_mouse_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][2]}), FUN = function(x){strsplit(x,"\\:")[[1]][1]})
suppa_mouse_eve_res <- suppa_mouse_eve_res[suppa_mouse_eve_res[,1] != 'NaN',]
suppa_mouse_eve_res$diff <- abs(suppa_mouse_eve_res[,1]) > 0.2 & suppa_mouse_eve_res[,2] < 0.05

suppa_mouse_iso_res <- read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/suppa/diffsplice/per_isoform/MF-MU_transcript_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_mouse_iso_res$gene_short_name <- sapply(row.names(suppa_mouse_iso_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_mouse_iso_res <- suppa_mouse_iso_res[suppa_mouse_iso_res[,1] != 'NaN',]
suppa_mouse_iso_res$diff <- abs(suppa_mouse_iso_res[,1]) > 0.2 & suppa_mouse_iso_res[,2] < 0.05

suppa_mouse_iso_ora <- enrich_CP(subset(suppa_mouse_iso_res, diff)$gene_short_name, universe = suppa_mouse_iso_res$gene_short_name, organisms = 'mouse')

suppa_mouse_eve_ora <- lapply(unique(suppa_mouse_eve_res$event_type), FUN = function(x){
  sub_df <- subset(suppa_mouse_eve_res, event_type == x)
  diff_g = subset(sub_df, diff)$gene_short_name
  if(length(diff_g) <= 2){
    return(NULL)
  }else{
    enrich_CP(diff_g, universe = sub_df$gene_short_name, organisms = 'mouse')
  }
})
names(suppa_mouse_eve_ora) <- unique(suppa_mouse_eve_res$event_type)
suppa_mouse_eve_ora$all <- enrich_CP(do.call(c, lapply(unique(subset(suppa_mouse_eve_res, diff)$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})), 
                                     universe = do.call(c, lapply(unique(suppa_mouse_eve_res$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})), 
                                     organisms = 'mouse')


suppa_rat_eve_res <- read.csv('./splicing/rat/suppa/rat_suppa_drimseq/suppa/diffsplice/per_local_event/RF-RU_local_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_rat_eve_res$gene_short_name <- sapply(row.names(suppa_rat_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_rat_eve_res$event_type <- sapply(sapply(row.names(suppa_rat_eve_res), FUN = function(x){strsplit(x,"\\;")[[1]][2]}), FUN = function(x){strsplit(x,"\\:")[[1]][1]})
suppa_rat_eve_res <- suppa_rat_eve_res[suppa_rat_eve_res[,1] != 'NaN',]
suppa_rat_eve_res$diff <- abs(suppa_rat_eve_res[,1]) > 0.2 & suppa_rat_eve_res[,2] < 0.05

suppa_rat_iso_res <- read.csv('./splicing/rat/suppa/rat_suppa_drimseq/suppa/diffsplice/per_isoform/RF-RU_transcript_diffsplice.dpsi', header = T, row.names = 1, sep = '\t')
suppa_rat_iso_res$gene_short_name <- sapply(row.names(suppa_rat_iso_res), FUN = function(x){strsplit(x,"\\;")[[1]][1]})
suppa_rat_iso_res <- suppa_rat_eve_res[suppa_rat_iso_res[,1] != 'NaN',]
suppa_rat_iso_res$diff <- abs(suppa_rat_iso_res[,1]) > 0.2 & suppa_rat_iso_res[,2] < 0.05

suppa_rat_iso_ora <- enrich_CP(subset(suppa_rat_iso_res, diff)$gene_short_name, universe = suppa_rat_iso_res$gene_short_name, organisms = 'rat')
suppa_rat_eve_ora <- lapply(unique(suppa_rat_eve_res$event_type), FUN = function(x){
  sub_df <- subset(suppa_rat_eve_res, event_type == x)
  diff_g = subset(sub_df, diff)$gene_short_name
  if(length(diff_g) <= 2){
    return(NULL)
  }else{
    enrich_CP(diff_g, universe = sub_df$gene_short_name, organisms = 'rat')
  }
})

names(suppa_rat_eve_ora) <- unique(suppa_rat_eve_res$event_type)
suppa_rat_eve_ora$all <- enrich_CP(do.call(c, lapply(unique(subset(suppa_rat_eve_res, diff)$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]}))
                                   , universe = do.call(c, lapply(unique(suppa_mouse_eve_res$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})), 
                                   organisms = 'rat')

  #Count genes
length(unique(suppa_mouse_eve_res$gene_short_name))
length(do.call(c, lapply(unique(subset(suppa_mouse_eve_res, local_MF.local_MU_p.val < 0.05 )$gene_short_name), FUN =function(x){strsplit(x, '_and_')[[1]]})))
tapply(row.names(suppa_mouse_eve_res), suppa_mouse_eve_res$event_type, FUN =function(x){length(unique(subset(suppa_mouse_eve_res[x,], diff)$gene_short_name))})

length(unique(suppa_rat_eve_res$gene_short_name))
length(unique(subset(suppa_rat_eve_res, diff)$gene_short_name))
tapply(row.names(suppa_rat_eve_res), suppa_rat_eve_res$event_type, FUN =function(x){length(unique(subset(suppa_rat_eve_res[x,], diff)$gene_short_name))})

}



### DRIMSEQ/STAGER
stager_rat <- read.csv('./splicing/rat/suppa/rat_suppa_drimseq/dexseq_dtu/results/stager/getAdjustedPValues.RF-RU.tsv', header = T,  sep = '\t')

#stager_rat_ora1 <- enrich_CP(subset(stager_rat, gene < 0.05)$geneID, universe = stager_rat$geneID, organisms = 'rat')

stager_mouse <- read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/results/stager/getAdjustedPValues.MF-MU.tsv', header = T,  sep = '\t')

#stager_mouse_ora1 <- enrich_CP(subset(stager_mouse, gdene < 0.05)$geneID, universe = stager_mouse$geneID, organisms = 'mouse')


mouse_dtu <- read_rnasplice_dex_dtu('./splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqDataSet.MF-MU.rds', './splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/filter/drimseq/dmDSdata.rds')
rat_dtu <- read_rnasplice_dex_dtu('./splicing/rat/suppa/rat_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqDataSet.RF-RU.rds', './splicing/rat/suppa/rat_suppa_drimseq/dexseq_dtu/filter/drimseq/dmDSdata.rds')
mouse_dtu$prop <- get_dexseq_prop(mouse_dtu)
rat_dtu$prop <- get_dexseq_prop(rat_dtu)

system.time({
  mouse_dtu$dexseq = estimateSizeFactors(mouse_dtu$dexseq)
  mouse_dtu$dexseq = estimateDispersions(mouse_dtu$dexseq)
  mouse_dtu$dexseq = testForDEU(mouse_dtu$dexseq, reducedModel = ~sample + exon)
})


mouse_dtu <- read_rnasplice_dex_dtu('./splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqDataSet.MF-MU.rds', './splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/filter/drimseq/dmDSdata.rds')
rat_dtu <- read_rnasplice_dex_dtu('./splicing/rat/suppa/rat_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqDataSet.RF-RU.rds', './splicing/rat/suppa/rat_suppa_drimseq/dexseq_dtu/filter/drimseq/dmDSdata.rds')
mouse_dtu$prop <- get_dexseq_prop(mouse_dtu)
rat_dtu$prop <- get_dexseq_prop(rat_dtu)

mouse_dtu$drimseq@samples$condition <- c('RF'='ratZygote', 'MF'='mouseZygote', 'RU'='ratEgg', 'MU'='mouseEgg')[as.character(mouse_dtu$drimseq@samples$condition)]
rat_dtu$drimseq@samples$condition <- c('RF'='ratZygote', 'MF'='mouseZygote', 'RU'='ratEgg', 'MU'='mouseEgg')[as.character(rat_dtu$drimseq@samples$condition)]

{"
mouse_drim <- drim_stageR(mouse_dtu, correct_det = F, filter = F)
mouse_drim_ora <- enrich_CP(subset(data.frame(mouse_drim$stageR_res), gene < 0.05)$gene_id, universe = mouse_drim$res.g$gene_id, organisms = 'mouse')
mouse_drim_prop <- get_drim_prop(mouse_drim$drim)
mouse_drim_diff <- data.frame(gene = mouse_drim_prop$gene_id, transcript=mouse_drim_prop$feature_id, diff = prop_diff(subset(mouse_drim_prop, select = -c(gene_id, feature_id)), group = mouse_drim$samps$group))
mouse_drim_prop_gene <- tapply(mouse_drim_diff, mouse_drim_diff$gene, FUN = function(x){sum(x$diff) > 0.4})
mouse_drim_prop_gene <- names(mouse_drim_prop_gene[mouse_drim_prop_gene])
print(length(intersect(mouse_drim_prop_gene, subset(data.frame(mouse_drim$stageR_res), gene < 0.05)$gene_id)))
mouse_drim_ora1 <- enrich_CP(intersect(mouse_drim_prop_gene, subset(data.frame(mouse_drim$stageR_res), gene < 0.05)$gene_id), universe = mouse_drim$res.g$gene_id, organisms = 'mouse')




rat_drim <- drim_stageR(rat_dtu, correct_det = F, filter = F, cond_names = c('RF' = 'ratZygote', 'RU' = 'ratEgg'), otherCond = 'ratZygote', seed = 12345)
rat_drim_ora <- enrich_CP(subset(data.frame(rat_drim$stageR_res), gene < 0.05)$gene_id, universe = rat_drim$res.g$gene_id, organisms = 'rat')
rat_drim_prop <- get_drim_prop(rat_drim$drim)
rat_drim_diff <- data.frame(gene = rat_drim_prop$gene_id, transcript=rat_drim_prop$feature_id, diff = prop_diff(subset(rat_drim_prop, select = -c(gene_id, feature_id)), group = rat_drim$samps$group))
rat_drim_prop_gene <- tapply(rat_drim_diff, rat_drim_diff$gene, FUN = function(x){sum(x$diff) > 0.4})
rat_drim_prop_gene <- names(rat_drim_prop_gene[rat_drim_prop_gene])
print(length(intersect(rat_drim_prop_gene, subset(data.frame(rat_drim$stageR_res), gene < 0.05)$gene_id)))
rat_drim_ora1 <- enrich_CP(intersect(rat_drim_prop_gene, subset(data.frame(rat_drim$stageR_res), gene < 0.05)$gene_id), universe = rat_drim$res.g$gene_id, organisms = 'rat')
"}
dexseq_mouse2 <- DEXSeqResults(mouse_dtu$dexseq, independentFiltering = T)
dexseq_mouse_stageR2 <- stageR_dexseqRes(dexseq_mouse2)

dexseq_mouse <- DEXSeqResults(mouse_dtu$dexseq, independentFiltering = T) #read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqResults.MF-MU.tsv', header = T,  sep = '\t')
dexseq_rat <- DEXSeqResults(rat_dtu$dexseq, independentFiltering = T) #read.csv('./splicing/rat/suppa/rat_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqResults.RF-RU.tsv', header = T,  sep = '\t')
dexseq_mouse_stageR <- stageR_dexseqRes(dexseq_mouse)
dexseq_rat_stageR <- stageR_dexseqRes(dexseq_rat)


#Number of transcripts DE causing DTU 
mouse_det_bygene <- tapply(subset(dexseq_mouse_stageR, gene < 0.05), subset(dexseq_mouse_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)})

stager_mouse_ora <- enrich_CP(subset(dexseq_mouse_stageR, gene < 0.05)$geneID, universe = dexseq_mouse_stageR$geneID, organisms = 'mouse')
stager_rat_ora <- enrich_CP(subset(dexseq_rat_stageR, gene < 0.05)$geneID, universe = dexseq_rat_stageR$geneID, organisms = 'rat')


stager_mouse_ora2 <- enrich_CP(subset(dexseq_mouse_stageR2, gene < 0.05)$geneID, universe = dexseq_mouse_stageR2$geneID, organisms = 'mouse')


{""
mouse_prop_diff <- data.frame(gene = mouse_dtu$prop$groupID, transcript=mouse_dtu$prop$featureID, diff = prop_diff(mouse_dtu$prop[,-c(1,2)], mouse_dtu$drimseq@samples$condition))
mouse_prop_gene <- tapply(mouse_prop_diff, mouse_prop_diff$gene, FUN = function(x){sum(x$diff) > 0.3})
mouse_prop_gene <- names(mouse_prop_gene[mouse_prop_gene])
print(length(intersect(mouse_prop_gene, subset(dexseq_mouse_stageR, gene < 0.05)$geneID)))
stager_mouse_ora1 <- enrich_CP(intersect(mouse_prop_gene, subset(dexseq_mouse_stageR, gene < 0.05)$geneID), universe = dexseq_mouse_stageR$geneID, organisms = 'mouse')

rat_prop_diff <- data.frame(gene = rat_dtu$prop$groupID, transcript=rat_dtu$prop$featureID, diff = prop_diff(rat_dtu$prop[,-c(1,2)], rat_drim$samps$group))
rat_prop_gene <- tapply(rat_prop_diff, rat_prop_diff$gene, FUN = function(x){sum(x$diff) > 0.3})
rat_prop_gene <- names(rat_prop_gene[rat_prop_gene])
print(length(intersect(rat_prop_gene, subset(dexseq_rat_stageR, gene < 0.05)$geneID)))
stager_rat_ora1 <- enrich_CP(unique(intersect(rat_prop_gene, subset(dexseq_rat_stageR, gene < 0.05)$geneID)), universe = dexseq_rat_stageR$geneID, organisms = 'rat')
""}

dotplot(stager_rat_ora$GO_BP_ora, showCategory=10, color='pvalue') +theme(axis.text.y = element_text(angle = 0, vjust = 0, hjust=0.5, size = 8), 
                                                                          axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                                                          legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
                                                                          legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
                                                                          legend.title = element_text(size=8))+ggtitle('Rat DTU GO: BP')+ scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.002))+
                                                                              scale_size(range=c(3,8), limits = c(5,240))+xlim(c(0, 0.12))
ggsave('Rat_dex_dtu_GOBP.png', width = 4, height = 4.5)

dotplot(clusterProfiler::simplify(stager_mouse_ora$GO_BP_ora), showCategory=10, color='pvalue') +theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8, colour = c(rep('black', 2), rep('red', 8))), 
                                                                    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                                                    legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
                                                                    legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
                                                                    legend.title = element_text(size=8))+ggtitle('Mouse DTU GO: BP')+
                                                                    scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.002))+
                                                                    scale_size(range=c(3,8), limits = c(5,240))+xlim(c(0, 0.12))
ggsave('Mouse_dex_dtu_GOBP.png', width = 4, height = 4.5)


View(subset(dexseq_mouse_stageR, gene < 0.05))

table(tapply(subset(dexseq_mouse_stageR, gene < 0.05), subset(dexseq_mouse_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)}))

table(tapply(subset(dexseq_rat_stageR, gene < 0.05), subset(dexseq_rat_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)}))


test_enrich <- enricher(
     bitr(unique( subset(dexseq_rat_stageR, gene < 0.05)$geneID), fromType = 'ALIAS', 'ENTREZID',org.Rn.eg.db)$ENTREZID,
     pvalueCutoff = 0.05,
     pAdjustMethod = "BH",
     universe = bitr(dexseq_rat_stageR$geneID, fromType = 'ALIAS', 'ENTREZID',org.Rn.eg.db)$ENTREZID,
     minGSSize = 10,
     maxGSSize = 500,
     qvalueCutoff = 0.2,
     gson = NULL,
     TERM2GENE = test$T2G,
     TERM2NAME = data.frame(TERM=names(test$T2N), NAME=test$T2N)
  )

dtu_deg_mouse <- list('Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                'Mouse Down DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                'Mouse DTU'= unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID))

UpSet(make_comb_mat(dtu_deg_mouse)[1:2], comb_col = c('black'))
UpSet(make_comb_mat(dtu_deg_rat)[1:2], comb_col = c('black'))

"#CD534CFF"
ggvenn(
  dtu_deg_mouse, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3, text_size = 5, show_percentage = F
)

dtu_deg_mouse <- list('Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'Mouse Down DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'Mouse DTU'= unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID))
dtu_deg_rat <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'Rat Down DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'Rat DTU'= unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID))

"#CD534CFF"
ggvenn(
  dtu_deg_rat, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3, text_size = 5, show_percentage = F
)

dtu_deg_mouse_ora_up <- enrich_CP(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID)),
                                'mouse', universe = intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), unique(dexseq_mouse_stageR$geneID)))
dtu_deg_mouse_ora_down <- enrich_CP(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)), unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID)),
                                  'mouse', universe = intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), unique(dexseq_mouse_stageR$geneID)))

dtu_deg_rat_ora_up <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'mouse', universe = row.names(mouse_intron_prop))
dtu_deg_rat_ora_down <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff < 0)), 'mouse', universe = row.names(mouse_intron_prop))




plotProportions(mouse_drim$drim, gene_id = 'Anapc5', group_variable = "condition", plot_type = "ribbonplot")+theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8))
plotProportions(rat_drim$drim, gene_id = 'Msl3', group_variable = "condition", plot_type = "boxplot1")+theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8))


# plot dexseq proportions

plotDEXSeqDTU(rat_dtu$prop, 'Foxm1', rat_dtu$drimseq@samples, isProportion = T)
plotDEXSeqDTU(mouse_dtu$prop, 'Abi3bp', mouse_drim$samps, isProportion = T)







print(length(unique(subset(stager_mouse, gene < 0.05)$geneID)))
print(length(unique(stager_mouse$geneID)))
print(length(unique(subset(stager_rat, gene < 0.05)$geneID)))
print(length(unique(stager_rat$geneID)))

##Rmats results
rmats_read_enrich <- function(res_dir, plot_dir, org = 'mouse'){
  rmats <- list()
  rmats[['a3ss']] <- read.csv(paste(res_dir, 'A3SS.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['a5ss']] <- read.csv(paste(res_dir, 'A5SS.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['se']] <- read.csv(paste(res_dir, 'SE.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['ri']] <- read.csv(paste(res_dir, 'RI.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats[['mxe']] <- read.csv(paste(res_dir, 'MXE.MATS.JC.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats_gs <- list()
  for(n in names(rmats)){
    rmats_gs[[n]] <- enrich_CP(subset(rmats[[n]], FDR < 0.05 & abs(IncLevelDifference) > 0.2)$GeneID, universe = rmats[[n]]$GeneID , organisms = org)
  }
  rmats_plot <- list()
  rmats_plot[['a3ss']] <- read.csv(paste(plot_dir, 'A3SS/event.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats_plot[['a5ss']] <- read.csv(paste(plot_dir, 'A5SS/event.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats_plot[['se']] <- read.csv(paste(plot_dir, 'SE/event.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats_plot[['ri']] <- read.csv(paste(plot_dir, 'RI/event.txt', sep = '/'), sep = '\t', header = T, row.names = 1)
  rmats_plot[['mxe']] <- read.csv(paste(plot_dir, 'MXE/event.txt', sep = '/'), sep = '\t', header = T, row.names = 1)

  return(list(rmats = rmats, rmats_plot = rmats_plot, rmats_gs = rmats_gs))
}




#### DAPARS APA
## The 

mouse_dapars_sf_st25 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)

mouse_dapars_ov_sf_st25 <- deg_utr2('./apa/mouse.overlap.sf.st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)

mouse_dapars_nr_sf_st25 <- deg_utr2('./apa/mouse.nonredundant.sf.st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta, combine_p = 'simes', diff_thresh = 0.2)
utr_mouse_nr_sf_st25_gs <- enrich_CP(unique(subset(mouse_dapars_nr_sf_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(mouse_dapars_nr_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

utr_mouse_ov_sf_st25_gs <- enrich_CP(unique(subset(mouse_dapars_ov_sf_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(mouse_dapars_nr_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')


mouse_dapars_sf_st25$gene_res[mouse_dapars_sf_st25$gene_res$APA_dist <= 75,]$diff = F
utr_mouse_sf_st25_up <- enrich_CP(unique(subset(mouse_dapars_sf_st25$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_st25_down <- enrich_CP(unique(subset(mouse_dapars_sf_st25$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_st25_gs <- enrich_CP(unique(subset(mouse_dapars_sf_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

utr_mouse_sf_st25_gs_strict <- enrich_CP(unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_names), universe = unique(subset(mouse_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'mouse')




mouse_dapars_sf_long_300_st25 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_300_long_st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_300_st25$gene_res[mouse_dapars_sf_long_300_st25$gene_res$APA_dist <= 75,]$diff = F
utr_mouse_sf_long_300_st25_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_long_300_st25$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_300_st25_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_300_st25$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_300_st25_gs <- enrich_CP(unique(subset(mouse_dapars_sf_long_300_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')


{"
mouse_dapars_sf_nodiff <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_short_long_nodiff.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_min_long <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_min_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_utr_100 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_utr_100.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_utr_100_med <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_utr_100_median.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_utr_100_min_long_median <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_utr_100_min_long_median.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_min_long_median <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_min_long_median.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)

mouse_dapars_sf_long_100 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_100_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_200 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_200_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_300 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_300_long.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)

mouse_dapars_sf_long_300$gene_res[mouse_dapars_sf_long_300$gene_res$APA_dist == 100,]$diff = F
mouse_dapars_sf_long_200$gene_res[mouse_dapars_sf_long_200$gene_res$APA_dist == 100,]$diff = F
mouse_dapars_sf_long_100$gene_res[mouse_dapars_sf_long_100$gene_res$APA_dist == 100,]$diff = F

mouse_dapars_sf$gene_res[mouse_dapars_sf$gene_res$APA_dist == 100,]$diff = F
utr_mouse_sf_gs_long_100_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_100$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_100$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_200_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_200$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_200$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_long_300_down <- enrich_CP(unique(subset(mouse_dapars_sf_long_300$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_long_300$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')



mouse_dapars_sf <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf$gene_res[mouse_dapars_sf$gene_res$APA_dist == 100,]$diff = F
utr_mouse_sf_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')


mouse_dapars_sf_long_100_st25 <- deg_utr('./apa/mouse_uniq_sfnorm_cov_0_100_long_st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), mouse$meta)
mouse_dapars_sf_long_100_st25$gene_res[mouse_dapars_sf_long_100_st25$gene_res$APA_dist <= 75,]$diff = F



utr_mouse_sf_utr_100_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_utr_100$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_utr_100_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_utr_100$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

utr_mouse_sf_utr_100_min_long_median_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100_min_long_median$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_utr_100_min_long_median$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_utr_100_min_long_median_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf_utr_100_min_long_median$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_utr_100_min_long_median$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_min_long_gs_up <- enrich_CP(unique(subset(mouse_dapars_sf_min_long$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_min_long$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')
utr_mouse_sf_min_long_gs_down <- enrich_CP(unique(subset(mouse_dapars_sf_min_long$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_min_long$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

  
"}


rat_dapars_sf_st25 <- deg_utr('./apa/rat_uniq_sfnorm_cov_0_st25.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)
rat_dapars_nr_sf_st25 <- deg_utr('./apa/rat.nonredundant.sf.st25.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)
rat_dapars_sf_st25$gene_res[rat_dapars_sf_st25$gene_res$APA_dist <= 75,]$diff = F
utr_rat_sf_st25_gs_up <- enrich_CP(unique(subset(rat_dapars_sf_st25$gene_res, diff & mean.diff > 0)$gene_short_names), universe = unique(rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_st25_gs_down <- enrich_CP(unique(subset(rat_dapars_sf_st25$gene_res, diff & mean.diff < 0)$gene_short_names), universe = unique(rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_st25_gs <- enrich_CP(unique(subset(rat_dapars_sf_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')

rat_dapars_nr_sf_st25 <- deg_utr('./apa/rat.nonredundant.sf.st25.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)
utr_rat_nr_sf_st25_gs <- enrich_CP(unique(subset(rat_dapars_nr_sf_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(rat_dapars_nr_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')


rat_dapars_sf_long_300_st25 <- deg_utr('./apa/rat_uniq_sfnorm_cov_0_300_long_st25.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)
rat_dapars_sf_long_300_st25$gene_res[rat_dapars_sf_long_300_st25$gene_res$APA_dist <= 75,]$diff = F
utr_rat_sf_long_300_st25_gs_up <- enrich_CP(unique(subset(rat_dapars_sf_long_300_st25$gene_res, diff & mean.diff > 0)$gene_short_names), universe = unique(rat_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_long_300_st25_gs_down <- enrich_CP(unique(subset(rat_dapars_sf_long_300_st25$gene_res, diff & mean.diff < 0)$gene_short_names), universe = unique(rat_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_long_300_st25_gs <- enrich_CP(unique(subset(rat_dapars_sf_long_300_st25$gene_res, diff & abs(mean.diff) >= 0.2)$gene_short_names), universe = unique(rat_dapars_sf_long_300_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')


utr_rat_sf_st25_gs_strict <- enrich_CP(unique(subset(rat_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_names), universe = unique(subset(rat_PAS_dist$PAS_motif, num_motif > 0)$gene_short_names), logFC = NULL, organisms = 'rat')


{"
rat_dapars_sf <- deg_utr('./apa/rat_uniq_sfnorm_cov_0.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)
rat_dapars_sf_long_100 <- deg_utr('./apa/rat_uniq_sfnorm_cov_0_100_long.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)
rat_dapars_sf_long_200 <- deg_utr('./apa/rat_uniq_sfnorm_cov_0_200_long.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)
rat_dapars_sf_long_300 <- deg_utr('./apa/rat_uniq_sfnorm_cov_0_300_long.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)

rat_dapars_sf$gene_res[rat_dapars_sf$gene_res$APA_dist == 100,]$diff = F
rat_dapars_sf_long_300$gene_res[rat_dapars_sf_long_300$gene_res$APA_dist == 100,]$diff = F
rat_dapars_sf_long_200$gene_res[rat_dapars_sf_long_200$gene_res$APA_dist == 100,]$diff = F
rat_dapars_sf_long_100$gene_res[rat_dapars_sf_long_100$gene_res$APA_dist == 100,]$diff = F

rat_dapars_sf_long_100_st25 <- deg_utr('./apa/rat_uniq_sfnorm_cov_0_100_long_st25.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), rat$meta)

print(length(subset(rat_dapars_sf$gene_res, diff & mean.diff > 0)$gene_short_names))
print(length(subset(rat_dapars_sf$gene_res, diff & mean.diff < 0)$gene_short_names))
utr_rat_sf_gs_up <- enrich_CP(unique(subset(rat_dapars_sf$gene_res, diff & mean.diff > 0)$gene_short_names), universe = unique(rat_dapars_sf$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_gs_down <- enrich_CP(unique(subset(rat_dapars_sf$gene_res, diff & mean.diff < 0)$gene_short_names), universe = unique(rat_dapars_sf$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')

utr_rat_sf_long_300_gs_up <- enrich_CP(unique(subset(rat_dapars_sf_long_300$gene_res, diff & mean.diff > 0.2)$gene_short_names), universe = unique(rat_dapars_sf_long_300$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_long_300_gs_down <- enrich_CP(unique(subset(rat_dapars_sf_long_300$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(rat_dapars_sf_long_300$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_long_200_gs_down <- enrich_CP(unique(subset(rat_dapars_sf_long_200$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(rat_dapars_sf_long_100$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
utr_rat_sf_long_100_gs_down <- enrich_CP(unique(subset(rat_dapars_sf_long_100$gene_res, diff & mean.diff < -0.2)$gene_short_names), universe = unique(rat_dapars_sf_long_200$gene_res$gene_short_names), logFC = NULL, organisms = 'rat')
"}


intersect(subset(mouse_dapars_sf_st25$gene_res, diff)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, diff)$gene_short_names)
intersect(subset(mouse_dapars_sf_long_300_st25$gene_res, diff)$gene_short_names, subset(rat_dapars_sf_long_300_st25$gene_res, diff)$gene_short_names)

both_utr_gs_down <- enrich_CP(intersect(subset(mouse_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff < -0.1)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff < -0.1)$gene_short_names),
          universe = intersect(mouse_dapars_sf_st25$gene_res$gene_short_names,rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

both_utr_gs_up <- enrich_CP(intersect(subset(mouse_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff > 0.2)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, fdr < 0.05 & mean.diff > 0.2)$gene_short_names),
                              universe = intersect(mouse_dapars_sf_st25$gene_res$gene_short_names,rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')

both_utr_gs <- enrich_CP(intersect(subset(mouse_dapars_sf_st25$gene_res, fdr < 0.05 & abs(mean.diff) > 0.1)$gene_short_names, subset(rat_dapars_sf_st25$gene_res, fdr < 0.05 & abs(mean.diff) > 0.1)$gene_short_names),
                              universe = intersect(mouse_dapars_sf_st25$gene_res$gene_short_names,rat_dapars_sf_st25$gene_res$gene_short_names), logFC = NULL, organisms = 'mouse')


utr_v_deg_plot <- ggplot(utr_v_deg2, aes(color=deg_sig, y=deg, x=utr, alpha = utr_sig, fill = deg_sig)) + 
  geom_point()+theme_classic()+
  theme(legend.position = c(0.01, .01),
        legend.justification = c("left", "bottom"),
        legend.box.just = "right",
        legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        legend.margin = margin(6, 10, 6, 6), 
        axis.text.x = element_text(size=14, color = 'black'),
        axis.text.y = element_text(size=14, color = 'black'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  guides(color=guide_legend(title=expression('Sig DEG')))+scale_alpha_manual(values = c('DEG up'=1, 'DEG down'=1, 'Non-DEG' = 0.3), guide = 'none')+
  scale_color_manual(values = c('blue','black', 'red'))+xlab('UTR PDUI difference')+ylab(expression('Log'[2]*'FC'))+geom_hline(yintercept=0, linetype="dashed", color = "blue", size = 2)

ggExtra::ggMarginal(utr_v_deg_plot, type = "density", 
                    yparams=list(mapping = aes(fill = utr_sig, x = deg), data = utr_v_deg2, alpha = 0.3), 
                    xparams=list(mapping = aes(fill = deg_sig, x = utr, alpha = 0.3), data = utr_v_deg2))


utr_deg_mouse <- list('DEG Up'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))),
                      'DEG Down'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -log2(2))),
                      'Long UTR'= unique(subset(mouse_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)$gene_short_name),
                      'Short UTR' =unique(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)$gene_short_name))




utr_deg_mouse <- list('DEG Up'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))),
                      'DEG Down'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -log2(2))),
                        'Long UTR'= unique(subset(mouse_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                      'Short UTR' =unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name))

utr_deg_rat <- list('DEG Up'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'DEG Down'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'Long UTR'= unique(subset(rat_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                      'Short UTR' =unique(subset(rat_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name))

png('utr_deg_mouse.png', width = 3, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_mouse)[1:4], comb_col = c('black', 'red', 'black','black'))
dev.off()

png('utr_deg_rat.png', width = 2.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_rat)[2:4], comb_col = c('black', 'black','black'))
dev.off()

setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names)


setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names)



mast_go_noployA_up <- enrich_CP(setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 
                                subset(mouse_dapars_sf$gene_res, mean.diff < -0.2 & diff)$gene_short_names),
                                universe = row.names(mast_res$mouseEgg_v_mouseZygote$DESig),
                                logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names), 
                                FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                                organisms = 'mouse',
                                GSE = F)

mast_go_short_utr_up <- enrich_CP(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 0.05)), 
                                        subset(mouse_dapars_sf$gene_res, mean.diff < -0.05 & diff)$gene_short_names),
                                universe = intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), 
                                                    mouse_dapars_sf$gene_res$gene_short_names),
                                logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), subset(mouse_dapars_sf$gene_res, mean.diff < -0.1)$gene_short_names), 
                                               FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                                organisms = 'mouse',
                                GSE = F)


mast_gse_noployA <- gse_CP(setdiff(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 
                                         subset(mouse_dapars_sf$gene_res, mean.diff < -0.05)$gene_short_names),
                                 logFC = sapply(setdiff(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), subset(mouse_dapars_sf$gene_res, mean.diff < -0.05)$gene_short_names), 
                                                FUN =function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                                 organisms = 'mouse',
                                 GSE = T)

mast_gse_rat <- gse_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & abs(Log2FC) > 1)),
                       logFC = sapply(row.names(mast_rat$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat', GSE = T)


rat_dapars_sf$pdui <- rat_dapars_sf$df[,grepl('PDUI', colnames(rat_dapars_sf$df))]
colnames(rat_dapars_sf$pdui) <- colnames(rat_dapars_sf$long)
rat_dapars_sf$pdui <- rat_dapars_sf$pdui[row.names(rat_dapars_sf$deg),]
rat_dapars_sf$pdui <- t(apply(rat_dapars_sf$pdui, 1, FUN = function(x){x[is.na(x)] = mean(x, na.rm =T); x}))


mouse_dapars_sf$pdui <- mouse_dapars_sf$df[,grepl('PDUI', colnames(mouse_dapars_sf$df))]
colnames(mouse_dapars_sf$pdui) <- colnames(mouse_dapars_sf$long)
mouse_dapars_sf$pdui <- mouse_dapars_sf$pdui[row.names(mouse_dapars_sf$deg),]
mouse_dapars_sf_st25$pdui <- t(apply(mouse_dapars_sf_st25$pdui, 1, FUN = function(x){x[is.na(x)] = mean(x, na.rm =T); x}))



sample_PCA(log(rat_dapars_sf_st25$pdui_imp+1), rat$meta, umap =  T, dimension=2 ,
           main = "Rat", labeling = F, point_size = 2, color_by = 'cellType',
           umap.config = list(n_neighbors = 25, min_dist = 0.5, metric='cosine', seed = 12345))+
            theme(axis.text.y = element_text( size = 8), axis.title.x = element_text(size = 8),
          axis.text.x = element_text( size = 8) , axis.title.y = element_text(size = 8))+guides(color=FALSE)
ggsave('rat_umap_DAP.png', units = 'in', dpi = 300, height = 2.5, width = 2.5)
sample_PCA(log(mouse_dapars_sf_st25$pdui_imp+1), mouse$meta, umap =  T, dimension=2 ,
           main = "Mouse", labeling = F, point_size = 2, color_by = 'cellType', 
           umap.config = list(n_neighbors = 25, min_dist = 0.5, metric='cosine', seed = 12345))+
          theme(axis.text.y = element_text( size = 8), axis.title.x = element_text(size = 8),
        axis.text.x = element_text( size = 8) , axis.title.y = element_text(size = 8))+guides(color=FALSE)
ggsave('mouse_umap_DAP.png', units = 'in', dpi = 300, height = 2.5, width = 2.5)



mouse_dapars_sf_st25$gene_res = mouse_dapars_sf_st25$gene_res[unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),]

test_cols =  rep('black', length(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS)$mean.diff))
names(test_cols) <- row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS))
test_alpha = rep(0.3, length(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS)$mean.diff))
names(test_alpha) <- row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS))

#test_cols[unique(do.call(c, lapply(utr_mouse_sf_st25_down$GO_BP_ora@result$geneID[1:40], FUN = function(x){strsplit(x,'/')[[1]]})))] <- 'red'
#test_alpha[unique(do.call(c, lapply(utr_mouse_sf_st25_down$GO_BP_ora@result$geneID[1:40], FUN = function(x){strsplit(x,'/')[[1]]})))] <- 1

test_cols[intersect(row.names(subset(mouse_dapars_sf_st25$gene_res,mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))))] <- 'red'
test_alpha[intersect(row.names(subset(mouse_dapars_sf_st25$gene_res,mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))))] <- 1

sig_PAS <- subset(mouse_PAS_dist$PAS_motif, num_motif > 0 & fdr < 0.05)$gene_short_names


png('short_UTR_vs_Expr_mouse.png', width = 3, height = 3.6, units = 'in', res = 300)
plot(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS)$mean.diff, 
     mast_res$mouseEgg_v_mouseZygote$DESig[row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05 & gene_short_names %in% sig_PAS)),]$Log2FC, 
     pch = 20, col = alpha(test_cols,  test_alpha), ylab='', xlab='')
title(line =0.8, main = expression('Shortened UTR vs Log'[2]*'FC'))
title( line=1.8, cex.lab=1.2, family="Arial", xlab = 'UTR PDUI Diff', ylab = expression('Expr Log'[2]*'FC'))
dev.off()
plot(mouse_dapars_sf_st25$gene_res[intersect(mouse_dapars_sf_st25$gene_res$gene_short_names, row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > 1 & fdr < 0.05))),]$mean.diff, 
     mast_res$mouseEgg_v_mouseZygote$DESig[intersect(mouse_dapars_sf_st25$gene_res$gene_short_names, row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > 1 & fdr < 0.05))),]$Log2FC, 
     pch = 20, col = alpha(test_cols,  test_alpha))



dotplot(clusterProfiler::simplify(utr_mouse_sf_st25_gs$GO_BP_ora, 0.7), showCategory=10, color='pvalue') +theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8, colour =rep('red', 10)), 
                                                                            axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                                                            legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
                                                                            legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
                                                                            legend.title = element_text(size=8))+ xlab('Gene ratio')+
  scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.002))+
  scale_size(range=c(3,8), limits = c(5,220))+ theme_classic()+
  theme(legend.position = c(0.99, .01),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        legend.margin = margin(6, 10, 6, 6), 
        axis.text.x = element_text(size=14, color = 'black'),
        axis.text.y = element_text(size=14, color = 'red'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave('Mouse_dap_GOBP.png', width = 5, height = 6)



dotplot(clusterProfiler::simplify(utr_rat_sf_st25_gs_strict$GO_BP_ora, 0.5), showCategory=10, color='pvalue') +theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8, colour =rep('black', 10)), 
                                                                                                                         axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                                                                                                         legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
                                                                                                                         legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
                                                                                                                         legend.title = element_text(size=8))+ggtitle('Rat DAP GO: BP')+
  scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.002))+
  scale_size(range=c(3,8), limits = c(3,50))


cnetplot(clusterProfiler::simplify(utr_mouse_sf_st25_gs$GO_BP_ora, 0.6), 10, foldChange = sapply(row.names(mouse_dapars_sf_st25$gene_res), FUN = function(x){
  mouse_dapars_sf_st25$gene_res[x,]$mean.diff
}), cex_label_gene = 0.3, cex_label_category = 1)
ggsave('Rat_dap_GOBP.png', width = 4, height = 4.5)
#### Orthology
utr_ortho <- list('Rat Long UTR'= row.names(subset(rat_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)),
                      'Rat Short UTR' =row.names(subset(rat_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)),
                      'Mouse Long UTR'= row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)),
                      'Mouse Short UTR' =row.names(subset(mouse_dapars_sf_st25$gene_res,mean.diff < -0.2 & fdr < 0.05)))

utr_ortho <- list('Rat Long UTR'= unique(subset(rat_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                  'Rat Short UTR' =unique(subset(rat_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                  'Mouse Long UTR'= unique(subset(mouse_PAS_dist$PAS_motif, mean.diff > 0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name),
                  'Mouse Short UTR' =unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name))

unique(subset(mouse_PAS_dist$PAS_motif, mean.diff < -0.2 & fdr < 0.05 & num_motif > 0)$gene_short_name)

png('utr_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_ortho)[1:4], comb_col = c('black'))
dev.off()
length(intersect(row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)), row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05))))

deg_ortho <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)),
                  'Rat Down DEG' =row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25)  & fdr < 0.05)),
                  'Mouse Up DEG'= row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)),
                  'Mouse Down DEG' =row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25)  & fdr < 0.05)))

mu_ru <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)))
md_ru <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)))
mu_rd <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)))
md_rd <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)))
m_r <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)), row.names(subset(mast_res$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)))

m_r_state <- c(rep('Mouse Up / Rat Up', length(mu_ru)),rep('Mouse down / Rat Up', length(md_ru)),rep('Mouse up / Rat down', length(mu_rd)),rep('Mouse down / Rat down', length(md_rd)))

names(m_r_state) <- c(mu_ru, md_ru,mu_rd,md_rd)

mouse_rat_all_deg <- intersect(row.names(mast_rat$ratEgg_v_ratZygote$DESig), row.names(mast_res$mouseEgg_v_mouseZygote$DESig))

deg_ortho_ora <- list('mu_ru' <- enrich_CP(mu_ru,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'mu_rd'<- enrich_CP(mu_rd,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'md_ru' <- enrich_CP(md_ru,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'md_rd'<- enrich_CP(md_rd,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse'),
                      'm_r'<- enrich_CP(md_rd,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse')
)
m_r_ora <- enrich_CP(m_r,universe = mouse_rat_all_deg, logFC = NULL, organisms = 'mouse')
m_r_ora$combined <- m_r_ora$GO_BP_ora
m_r_ora$combined@result <- rbind(m_r_ora$GO_BP_ora@result, m_r_ora$WKP_ora@result, m_r_ora$KEGG_ora@result, m_r_ora$REACT_ora@result)
m_r_ora$combined@geneSets <- c(m_r_ora$GO_BP_ora@geneSets, m_r_ora$REACT_ora@geneSets,m_r_ora$KEGG_ora@geneSets, m_r_ora$WKP_ora@geneSets)
m_r_ora_cnetplots <- custom_cnet_plot(m_r_ora$combined, category = c('R-MMU-1428517', 'R-MMU-927802', 'R-MMU-72172', 'R-MMU-72202', 'R-MMU-69620', 'R-MMU-72312', 'GO:0008380', 'GO:0002181', 'GO:0016574','R-MMU-4551638'),
                         gene_color = sapply(mast_rat$ratEgg_v_ratZygote$DESig$features, FUN = function(x){mast_rat$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), 
                         gene_color2 = sapply(mast_res$mouseEgg_v_mouseZygote$DESig$features, FUN = function(x){mast_res$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                         layout = 'fr', color_cat_pval = F)

ggsave(plot = m_r_ora_cnetplots$plot1, 'm_r_ora_mouse.png', width = 5, height = 3)
ggsave(plot = m_r_ora_cnetplots$plot2, 'm_r_ora_rat.png', width = 5, height = 3)



png('deg_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(deg_ortho)[1:4], comb_col = c('red', 'black',  'black', 'blue'))
dev.off()


dnp_ortho <- list('Mouse Increase DNP'= row.names(subset(test, prop_diff > 0 & qvalue < 0.05)),
                    'Mouse Decrease DNP' =row.names(subset(test, prop_diff < 0 & qvalue < 0.05)),
                    'Rat Increase DNP'= row.names(subset(test1, prop_diff > 0 & qvalue < 0.05)),
                    'Rat Decrease DNP' =row.names(subset(test1, prop_diff < 0 & qvalue < 0.05)))
png('dnp_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(dnp_ortho)[1:4], comb_col = c('black'))
dev.off()

dtu_ortho <- list( 'Rat DTU'= unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID),
                      'Mouse DTU'= unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID))

dtu_ortho_ora <- enrich_CP(intersect(unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID), unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID)), 
                                            universe = intersect(dexseq_rat_stageR$geneID, dexseq_mouse_stageR$geneID), organisms = 'mouse')



png('dtu_ortho.png', width = 2.6, height = 2.6, res = 300, units = 'in')

ggvenn(
  dtu_ortho, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 0, text_size = 5, show_percentage = F, auto_scale = T
)

dev.off()


png('utr_deg_rat.png', width = 2.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_rat)[1:3], comb_col = c('black', 'black','black'))
dev.off()

ggvenn(
  utr_ortho, 
  fill_color = c("#0073C2FF", "#EFC000FF","#CD534CFF",  "#868686FF"),
  stroke_size = 0.5, set_name_size = 5, text_size = 5, show_percentage = F
)


### Plotting Genes Regions From R

test_vp <- plot_range_coverage('NC_000082.7:87,251,316-87,281,194', strand = '+', mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), log = F, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('mouseZygote', 'mouseEgg'), file_name_suffix = 'mouse', gene_name = 'Usp16')
test_vp <- plot_range_coverage('NC_051349.1:35,194,767-35,197,558', strand = '-', rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), log = T, cols = DOT_COLOR[c('ratZygote', 'ratEgg')], cell_names = c('ratZygote', 'ratEgg'), file_name_suffix = 'rat', 'Oog1')
test_vp <- plot_range_coverage('NC_000073.7:15,130,624-15,132,520', strand = '+', mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('mouseZygote', 'mouseEgg'), file_name_suffix = 'mouse')





test_vp <- plot_utr_coverage('Cdc25a', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat Zygote', 'Rat Egg'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Cnot6l', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/igv/RF.sf.average.sorted.bw', './dataset/igv/RU.sf.average.sorted.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat Zygote', 'Rat Egg'), file_name_suffix = 'rat', loci = 'NC_051349.1:13495326-13502511')
test_vp <- plot_utr_coverage('Nek7', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse Zygote', 'Mouse Egg'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Btg4', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/igv/MF.sf.average.sorted.bw', './dataset/igv/MU.sf.average.sorted.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse Zygote', 'Mouse Egg'), file_name_suffix = 'mouse')




# Comparisons with previous studies on ZGA in mice
## Park, 2013, 2015 DBMTEE
park_data <- do.call(rbind, lapply(1:25,FUN = function(x){
  df <- read.csv(paste('~/Downloads/tables/hclust_result/CLS_norm_',x,'.dat', sep = ''), sep = '\t', skip = 1);
  df$clust <- x
  if(x %in% c(1,2,4,12,13,17,18,19,20)){
    df$pattern <- 'maternal'
  }else if(x %in% c(5,6,7,8,9,10,11,14,22,23)){
    df$pattern <- 'mZGA'
  }else{
    df$pattern <- 'Other'
  }
  df}))

common_genes_park <- intersect(row.names(mast_res$mouseEgg_v_mouseZygote$DESig), row.names(park_data))
park_data_sub <- park_data[common_genes_park,]
mast_data_sub <- mast_res$mouseEgg_v_mouseZygote$DESig[common_genes_park,]

dnp_ortho <- list('Our data mZGA'= row.names(subset(mast_data_sub, Log2FC > 1 )),
                  'Our data maternal' =row.names(subset(mast_data_sub, Log2FC < -1)),
                  'DMTEE Maternal'= row.names(subset(park_data_sub, pattern == 'maternal')),
                  'DMTEE mZGA ' =row.names(subset(park_data_sub, pattern == 'mZGA')))
UpSet(make_comb_mat(dnp_ortho )[1:4], comb_col = c('black'))



######Intergenic and close-to-gene reads from QoRTs
mouse_intergenic <- lapply(list.files('./dataset/mouse_data/QCData/hisat/'), function(x){
  df <- read.csv(paste('./dataset/mouse_data/QCData/hisat/', x, 'QC.geneCounts.detailed.txt.gz', sep ='/'), header = T, sep = '\t', row.names = 1)
  df$far_AMBIG.. <- as.numeric(str_remove_all(df$far_AMBIG..,  "[()]"))
  df
})
names(mouse_intergenic) <- list.files('./dataset/mouse_data/QCData/hisat/')


mouse_qorts_gene <- do.call(cbind, lapply(mouse_intergenic, function(x){x[,'COUNT']}))
mouse_close_ig <- do.call(cbind, lapply(mouse_intergenic, function(x){x[,'nearby']}))
row.names(mouse_close_ig) <- row.names(mouse_intergenic$MF1)
mouse_far_ig <- do.call(cbind, lapply(mouse_intergenic,function(x){x[,'far']}))
row.names(mouse_far_ig) <- row.names(mouse_intergenic$MF1)
mouse_far_ig_amb <- do.call(cbind, lapply(mouse_intergenic,function(x){x[,'far_AMBIG..']}))
mouse_qorts_intron <- do.call(cbind, lapply(mouse_intergenic, function(x){x[,'intronic']}))
row.names(mouse_qorts_intron) <- row.names(mouse_intergenic$MF1)

intersect(row.names(mouse_qorts_intron), row.names(mouse_intron$unspliced))



rat_intergenic <- lapply(list.files('./dataset/rat_data/QCData/hisat/'), function(x){
  df <- read.csv(paste('./dataset/rat_data/QCData/hisat/', x, 'QC.geneCounts.detailed.txt.gz', sep ='/'), header = T, sep = '\t', row.names = 1)
  df$far_AMBIG.. <- as.numeric(str_remove_all(df$far_AMBIG..,  "[()]"))
  df
})
names(rat_intergenic) <- list.files('./dataset/rat_data/QCData/hisat/')

rat_close_ig <- do.call(cbind, lapply(rat_intergenic, function(x){x[,'nearby']}))
rat_far_ig <- do.call(cbind, lapply(rat_intergenic,function(x){x[,'far']}))
rat_far_ig_amb <- do.call(cbind, lapply(rat_intergenic,function(x){x[,'far_AMBIG..']}))
rat_qorts_intron <- do.call(cbind, lapply(rat_intergenic, function(x){x[,'intronic']}))

## Intergenic from mostly 1000bp segments in intergenic regions


rat_hisat_intergenic <- read.csv('./dataset/rat_intergenic.hisat.ct.txt', row.names = 1, header =T, sep = '\t')
colnames(rat_hisat_intergenic) <- sapply(strsplit(colnames(rat_hisat_intergenic), '_'), function(x){x[1]})
mouse_hisat_intergenic <- read.csv('./dataset/mouse_intergenic.hisat.ct.txt', row.names = 1, header =T, sep = '\t')
colnames(mouse_hisat_intergenic) <- sapply(strsplit(colnames(mouse_hisat_intergenic), '_'), function(x){x[1]})


rat_hisat_intergenic1 <- read.csv('./dataset/mouse_data/rat_intergenic_1000.hisat.ct.txt', row.names = 1, header =T, sep = '\t')
colnames(rat_hisat_intergenic1) <- sapply(strsplit(colnames(rat_hisat_intergenic1), '_'), function(x){x[1]})
mouse_intergenic1 <- read_htseq_intergenic(samples = colnames(mouse$ct$bio), cov_thresh = 0)
rat_intergenic1 <- read_htseq_intergenic('./dataset/mouse_data/rat_intergenic_1000.hisat.ct.txt', samples = colnames(rat$ct$bio), 0)


rat_intergenic1$df %>% melt('cellType') %>% ggbarplot( x = "variable", y = "value",color = "cellType",fill='cellType', add.params = list(color = 'black'), palette = 'rickandmorty', add = "mean_se",position = position_dodge(0.8))+ 
  stat_compare_means(aes(group = cellType), method = 't.test', label = "p.signif", label.y = c(5000,4000,8000,8000,3000,20000))+xlab('Location')+ylab('Number of Regions')


rat_intergenic1$df %>% melt('cellType') %>% ggbarplot( x = "variable", y = "value",color = "cellType",fill='cellType', add.params = list(color = 'black'), palette = 'rickandmorty', add = "mean_se",position = position_dodge(0.8))+ 
  stat_compare_means(aes(group = cellType), method = 't.test', label = "p.signif", label.y = c(5000,4000,8000,8000,3000,20000))+xlab('Location')+ylab('Number of Regions')








mouse_intergenic1 <- read_htseq_intergenic('./dataset/mouse_data/mouse_intergenic_1000.hisat.ct.txt', ct = mouse$ct$bio, meta = mouse$meta, cov_thresh = 1)
rat_intergenic1 <- read_htseq_intergenic('./dataset/mouse_data/rat_intergenic_1000.hisat.ct.txt', ct = rat$ct$bio, meta = rat$meta, 1)
mouse_rat_int_df <- rbind(rat_intergenic1$small_df, mouse_intergenic1$small_df)
mouse_rat_int_df$organism <- c(rep('Rat', 25), rep('Mouse', 28))

mouse_rat_int_df %>% melt(c('cellType', 'organism')) %>% ggbarplot( x = "variable", y = "value", color = "cellType",fill='cellType', facet.by = 'organism', add.params = list(color = '#ffa500'), 
                                                                    palette = DOT_COLOR, add = "mean_se",position = position_dodge(0.8))+
  stat_compare_means(aes(group = cellType), method = 't.test', label = "p.signif", label.y.npc = c(0.76, 0.76, 0.76 ,0.76), size = 8)+
  xlab('Location')+ylab('Number of Regions')+theme_classic2(base_size = 15)+ theme(legend.position="none", axis.text.x = element_text(size = 14))+scale_y_continuous(name="Number of Regions", labels = label_number(scale_cut = cut_short_scale()))
ggsave('intergenic_regions.png', width = 4,height = 3.5)


only_up_1000_regions = grepl('_u_0', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3

only_down_1000_regions = grepl('_d_0', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3

shared_1000_regions = grepl('_0', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) > 3

close_up_regions = grepl('_u', row.names(mouse_hisat_intergenic1)) & !grepl('_0', row.names(mouse_hisat_intergenic1)) & !grepl('far', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3
close_down_regions = grepl('_d', row.names(mouse_hisat_intergenic1)) & !grepl('_0', row.names(mouse_hisat_intergenic1)) & !grepl('_far', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) == 3

close_shared_regions = !grepl('_0', row.names(mouse_hisat_intergenic1)) & !grepl('far', row.names(mouse_hisat_intergenic1)) & sapply(strsplit(row.names(mouse_hisat_intergenic1), '_'), length) > 3

far_regions = grepl('_far', row.names(mouse_hisat_intergenic1)) 

sum(only_up_1000_regions)

grepl('int_s', colnames(mouse_hisat_intergenic1)), grepl('int_n', colnames(mouse_hisat_intergenic1)), grepl('.u.', colnames(mouse_hisat_intergenic1))


num_mouse_genes_intergenic <- colSums(mouse_hisat_intergenic1[grepl('_0', row.names(mouse_hisat_intergenic1)),] > 0)[row.names(mouse$meta)]
num_rat_genes_intergenic <- colSums(rat_hisat_intergenic1[grepl('_0', row.names(rat_hisat_intergenic1)),] > 0)[row.names(rat$meta)]






intergenic_num_summary <- data.frame(percentages = c(num_mouse_genes_intergenic, num_rat_genes_intergenic),
                                 Species = c(rep('Mouse', 28), rep('Rat', 25)),
                                 cellType = c(mouse$meta$cellType, rat$meta$cellType)) %>% group_by(cellType) %>% summarise(Num_of_genes = mean(percentages), sd_Perc = sd(percentages), Species = unique(Species))



ggplot(intergenic_num_summary, aes(Species, Num_of_genes, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Num_of_genes - sd_Perc, ymax = Num_of_genes + sd_Perc), color = '#ffa500',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+
  annotate("text", x = 1, y = 0.023, label = "", size = 5) +annotate("text", x = 2, y = 0.025, label = "**", size = 5)+ggtitle('Upstream Intergenic Read Number of Genes')+theme(plot.title = element_text(hjust = 0.5, size = 12))
ggsave('up_intergenic.png', width = 3,height = 2.5)









