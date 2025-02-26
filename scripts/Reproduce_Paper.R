source('./scripts/Reproduce_functions.R')
source('./scripts/Misc.R')
# read in functions
mouse_gtf <- makeTxDbFromGFF('./dataset/mouse_spike.gtf', 'gtf')
rat_gtf <- makeTxDbFromGFF('./dataset/rat_spike.gtf', 'gtf')
tx2gene_mouse <- AnnotationDbi::select(mouse_gtf, keys(mouse_gtf, keytype = "TXNAME"), "GENEID", "TXNAME")
tx2gene_rat <- AnnotationDbi::select(rat_gtf, keys(rat_gtf, keytype = "TXNAME"), "GENEID", "TXNAME")
write.table(tx2gene_mouse, './dataset/output_mouse/tx2gene.tsv', sep = '\t', quote = F)
write.table(tx2gene_rat, './dataset/output_rat/tx2gene.tsv', sep = '\t', quote = F)

# Table S1 
mouse_0 <- list()
mouse_0$meta <- read.table("./dataset/output_mouse/meta.tsv", sep = '\t', header=TRUE, row.names = 1)
mouse_0 <- prepareCount(mouse_0, dirc = './output_mouse/')
mouse_0 <- prepareGeneFeatures(mouse_0, './dataset/mouse_spike.gtf')
mouse_0$meta$mito_rate <- colSums(mouse_0$ct$bio[row.names(subset(mouse_0$features_ct, seqnames == 'NC_005089.1')),row.names(mouse_0$meta)])/colSums(mouse_0$ct$bio[,row.names(mouse_0$meta)])
mouse_0$meta$rRNA_rate <- colSums(mouse_0$ct$bio[row.names(subset(mouse_0$features_ct, gene_biotype == 'rRNA')),row.names(mouse_0$meta)])/colSums(mouse_0$ct$bio)[row.names(mouse_0$meta)]
mouse_0$meta$sensitivity <- colSums(mouse_0$ct$bio > 0)

rat_0 <- list()
rat_0$meta <- read.table("./dataset/output_rat/meta.tsv", sep = '\t', header=TRUE, row.names = 1)
rat_0 <- prepareCount(rat_0, dirc = './output_rat/')
rat_0 <- prepareGeneFeatures(rat_0, './dataset/rat_spike.gtf')
rat_0$meta$mito_rate <- colSums(rat_0$ct$bio[row.names(subset(rat_0$features_ct, seqnames == 'NC_001665.2')),row.names(rat_0$meta)])/colSums(rat_0$ct$bio)[row.names(rat_0$meta)]
rat_0$meta$rRNA_rate <- colSums(rat_0$ct$bio[row.names(subset(rat_0$features_ct, gene_biotype == 'rRNA')),row.names(rat_0$meta)])/colSums(rat_0$ct$bio)[row.names(rat_0$meta)]
rat_0$meta$sensitivity <- colSums(rat_0$ct$bio > 0)

tableS1 <- rbind(mouse_0$meta, rat_0$meta)
mouse_map_info <- read.csv('./dataset/output_mouse/multiqc_data/multiqc_star.txt', sep = '\t', row.names = 1)
rat_map_info <- read.csv('./dataset/output_rat/multiqc_data/multiqc_star.txt', sep = '\t', row.names = 1)
map_info <- rbind(mouse_map_info, rat_map_info)
tableS1$star_unique_Mapping_Rate <- map_info[row.names(tableS1), 'uniquely_mapped_percent']
tableS1$input_reads_to_star <- map_info[row.names(tableS1), "total_reads"]
write.csv(tableS1, './table_S1.csv')
options(ucscChromosomeNames=FALSE)



mouse <- list()
mouse$meta <- read.table("./dataset/output_mouse/meta.tsv", sep = '\t', header=TRUE, row.names = 1)
mouse$meta <- subset(mouse$meta, !unique.ID %in% c('MU28', 'MU32', 'MU18'))
mouse <- prepareCount(mouse, dirc = './output_mouse/')
mouse <- prepareGeneFeatures(mouse, './dataset/mouse_spike.gtf')
mouse$meta$mito_rate <- colSums(mouse$ct$bio[row.names(subset(mouse$features_ct, seqnames == 'NC_005089.1')),row.names(mouse$meta)])/colSums(mouse$ct$bio[,row.names(mouse$meta)])
mouse$meta$rRNA_rate <- colSums(mouse$ct$bio[row.names(subset(mouse$features_ct, gene_biotype == 'rRNA')),row.names(mouse$meta)])/colSums(mouse$ct$bio)[row.names(mouse$meta)]
mouse$meta$sensitivity <- colSums(mouse$ct$bio > 0)
mouse_pc_lncRNA_tpseudo_nonmito_nonrRNA <- row.names(subset(mouse$features_tpm, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1'))
mouse$tpm$bio <- mouse$tpm$bio[mouse_pc_lncRNA_tpseudo_nonmito_nonrRNA,]
#mouse$tpm$bio <- t(1e6*t(mouse$tpm$bio)/colSums(mouse$tpm$bio))
mouse$ct$bio <- mouse$ct$bio[row.names(subset(mouse$features_ct, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_005089.1' )),]



rat <- list()
rat$meta <- read.table("./dataset/output_rat/meta.tsv", sep = '\t', header=TRUE, row.names = 1)
rat$meta <- subset(rat$meta, !unique.ID %in% c('RF7', 'RF8','RF16', 'RF13', 'RU14', 'RF11'))
rat <- prepareCount(rat, dirc = './output_rat/')
rat <- prepareGeneFeatures(rat, './dataset/rat_spike.gtf')
rat$meta$mito_rate <- colSums(rat$ct$bio[row.names(subset(rat$features_ct, seqnames == 'NC_001665.2')),row.names(rat$meta)])/colSums(rat$ct$bio)[row.names(rat$meta)]
rat$meta$rRNA_rate <- colSums(rat$ct$bio[row.names(subset(rat$features_ct, gene_biotype == 'rRNA')),row.names(rat$meta)])/colSums(rat$ct$bio)[row.names(rat$meta)]
rat$meta$sensitivity <- colSums(rat$ct$bio > 0)
rat_pc_lncRNA_tpseudo_nonmito_nonrRNA <- row.names(subset(rat$features_tpm, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_001665.2'))
rat$tpm$bio <- rat$tpm$bio[row.names(subset(rat$features_tpm, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_001665.2')),]
#rat$tpm$bio <- t(1e6*t(rat$tpm$bio)/colSums(rat$tpm$bio))
rat$ct$bio <- rat$ct$bio[row.names(subset(rat$features_ct, gene_biotype %in% c('protein_coding', 'lncRNA', 'transcribed_pseudogene') & seqnames != 'NC_001665.2')),]


#make gene coverage plots with QoRTs
mouse.qorts <- subset(read.csv('./dataset/output_mouse/meta.tsv', row.names = 1, header = T, sep = '\t'), organism == 'mouse')[row.names(mouse$meta),]
rat.qorts <- subset(read.csv('./dataset/output_rat/meta.tsv', row.names = 1, header = T, sep = '\t'), organism == 'rat')[row.names(rat$meta),]
mouse.qorts$input.read.pair.count <- get_qorts_summary('./dataset/output_mouse/QCData/')['PREALIGNMENT_READ_CT', row.names(mouse.qorts)]
rat.qorts$input.read.pair.count <- get_qorts_summary('./dataset/output_rat/QCData/')['PREALIGNMENT_READ_CT', row.names(rat.qorts)]

write.table(mouse.qorts,'./dataset/mouse_data/meta.tsv', sep = '\t', quote = F, row.names = T, col.names = T )
write.table(rat.qorts,'./dataset/rat_data/meta.tsv', sep = '\t', quote = F, row.names = T, col.names = T )



mouse.res.pair <- read.qc.results.data("./dataset/output_mouse/", decoder=mouse.qorts, autodetectMissingSamples = TRUE,  debugMode = F, calc.DESeq2 = T, calc.edgeR = T)
rat.res.pair <- read.qc.results.data("./dataset/output_rat/", decoder=rat.qorts, autodetectMissingSamples = TRUE, debugMode = F,  calc.DESeq2 = T, calc.edgeR = T)

plot_qorts(res = mouse.res.pair, sep = 'group.ID', plot_path = './mouse_qort_plot.pdf', colorby_colors = DOT_COLOR)
plot_qorts(res = rat.res.pair, sep = 'group.ID', plot_path = './rat_qort_plot.pdf', colorby_colors = DOT_COLOR)

#genes detected plots with a wilcoxon rank-sum test for difference
genes_detected_S1B <- data.frame(detected = c(mouse$meta$sensitivity, rat$meta$sensitivity),
                                   organism = c(rep('Mouse', nrow(mouse$meta)), rep('Rat', nrow(rat$meta))),
                                   cellType = c(mouse$meta$cellType, rat$meta$cellType), 
                                   tissue = as.factor(sapply(c(mouse$meta[,'cellType'], rat$meta[,'cellType']), function(x) ifelse(grepl('Egg', x), 'Oocyte', 'Zygote')))
)

genes_detected_S1B_plot <- genes_detected_S1B  %>% reshape2::melt(c('cellType', 'organism', 'tissue'))   %>% ggbarplot( x = 'tissue', y = "value", color = "cellType",fill='cellType', facet.by = 'organism', add.params = list(color = '#ffa500'), 
                                                                                                                   palette = DOT_COLOR,position = position_dodge(0.8))+
  stat_compare_means( method = 'wilcox.test', label = "p.format", label.y.npc = c(0.8, 0.8), label.x.npc = c(0.05, 0.4), size = 3.5)+
  xlab('Cell Type')+theme_classic2(base_size = 15)+ 
  theme(legend.position="none", 
        axis.text.x = element_text(size = 14, color = 'black'),
        strip.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 14, color = 'black'))+scale_y_continuous(name='Number of Genes Detected', labels = label_number(scale_cut = cut_short_scale()))
ggsave('genes_detected_S1B_plot.png', plot = genes_detected_S1B_plot, width = 4,height = 3.5)





# plot UMAPs 
umap_config= list(n_neighbors = 10, min_dist = 0.5, metric='pearson', seed = 123456)
umap_axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)
mouse_umap <- sample_PCA(DESeq2::varianceStabilizingTransformation(as.matrix(mouse$ct$bio)), mouse$meta, umap =  T, dimension=2 ,
           main = "Mouse", labeling = F, point_size = 2, color_by = 'cellType', umap.config = umap_config, legend.position = c(0.1, 0.9))+theme(aspect.ratio=1)
ggsave('./mouse_umap.png', plot = mouse_umap,height = 2.5, width = 2.5)
rat_umap <- sample_PCA(DESeq2::varianceStabilizingTransformation(as.matrix(rat$ct$bio)), rat$meta, umap =  T, dimension=2 ,
           main = "Rat", labeling = F, point_size = 2, color_by = 'cellType', umap.config = umap_config, legend.position = c(0.8, 0.9))+theme(aspect.ratio=1)
ggsave('./rat_umap.png',plot = rat_umap, height = 2.5, width = 2.5)





# DGE with MAST
mast_mouse <- mast_diff(mouse, plot = T, tpm = F, nbins = 0, freq = 0.1,  min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F, correct_wild_coef = T)
mast_gse_mouse <- gse_CP('mouse', n_type = 'ALIAS', logFC = sapply(row.names(mast_mouse$mouseEgg_v_mouseZygote$DESig), FUN =function(x){mast_mouse$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), full_combine = T)
mast_ora_mouse_up <- enrich_CP(row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)), 'mouse', universe = row.names(mast_mouse$mouseEgg_v_mouseZygote$DESig), n_type = 'ALIAS', full_combine = T)
mast_ora_mouse_down <- enrich_CP(row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)), 'mouse', universe = row.names(mast_mouse$mouseEgg_v_mouseZygote$DESig), n_type = 'ALIAS', full_combine = T)


mast_rat <- mast_diff(obj = rat, control = 'ratEgg',  plot = T ,tpm = F, nbins = 0, freq = 0.1, min_cell_grp = 2, min_cell = 5, include_filt_as_NA = F, correct_wild_coef = T)
mast_gse_rat <- gse_CP(logFC = sapply(row.names(mast_rat$ratEgg_v_ratZygote$DESig), FUN =function(x){mast_rat$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), organisms = 'rat',n_type = 'ALIAS',  full_combine = T)
mast_ora_rat_up <- enrich_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)), 'rat', universe = row.names(mast_rat$ratEgg_v_ratZygote$DESig), n_type = 'ALIAS', full_combine = T)
mast_ora_rat_down <- enrich_CP(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)), 'rat', universe = row.names(mast_rat$ratEgg_v_ratZygote$DESig), n_type = 'ALIAS',  full_combine = T)

# plot ridge plots

tree_ridge_mouse <- cp_tree_ridge_plot(mast_gse_mouse$combined_full, alpha = 0.05, geneSet = c('R-MMU-3247509', 'GO:0002181', 'R-MMU-927802','GO:0006090',
                                                                                   'GO:0050684', 'GO:0016441','R-MMU-72312',
                                                                                   'R-MMU-111933','GO:0042448','R-MMU-1428517',
                                                                                   'WP403', 'GO:0043547', 'GO:0043405','GO:0060070',
                                                                                   'GO:0061614','GO:0007565','GO:0038061','WP265','R-MMU-75953',
                                                                                   'WP113','GO:0051017','R-MMU-69620', 'GO:0140013', 
                                                                                   'GO:0051028', 'R-MMU-453276', "R-MMU-194441",
                                                                                   "R-MMU-2980766","R-MMU-5693532","GO:0006086","R-MMU-69306"))
print(tree_ridge_mouse$both) & theme(legend.position='bottom', legend.text = element_text(angle = 20), legend.title.position = 'top')
ggsave( './mouse_gsea_ridge.png', height = 9, width =6, units = 'in')

tree_ridge_rat <- cp_tree_ridge_plot(mast_gse_rat$combined_full, alpha = 0.05, nclust = 5)
print(tree_ridge_rat$both) & theme(legend.position='bottom', legend.text = element_text(angle = 20), legend.title.position = 'top')
ggsave('./rat_gsea_ridge.png', height = 9, width =4, units = 'in')

mouse_vol_plot <- volcano_plot(mast_mouse$mouseEgg_v_mouseZygote$DESig, pval_col = 'fdr',fc = log2(2), top_genes = c('Nup37', 'Nup54', 'Obox1', 'Obox2', 'Obox5', 'Obox7', 'Nup35', 'Rpl9', 
                                                                                                       'Rpl3', 'Rpl26', 'Rpl12', 'Rpl19', 'Rpl18', 'Rpl3','Mastl','Cdc20', 
                                                                                                       'Cnot7','Mad2l1', 'Ube2e1','Taf6', 'Taf9','Gtf2b', 'Gtf2a2', 'Ndufc1',
                                                                                                       'Ndufa2', 'Ndufb7','Atp5h','Atp1a1','Ndufa8','Atp5j','Ndufa3', 'Dlat' ))+theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15))
ggsave('./mouse_volcano_mast.png', plot = mouse_vol_plot, width = 4, height = 4)
rat_vol_plot <- volcano_plot(mast_rat$ratEgg_v_ratZygote$DESig, pval_col = 'fdr',fc = log2(2), top_genes = c('Rps20', 'Rps2', 'Rps15', 'Rplp1','Hist2h3c2','H2ac1', 'H2ax', 'Gata3',
                                                                                                   'Sdhc','Aco2','Idh1','Sdhd','Pdha1','Dlst','Ndufab1','Ndufc1'))+theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15))
ggsave('./rat_volcano_mast.png', plot = rat_vol_plot, width = 4, height = 4)

mouse_deg <- mast_mouse$mouseEgg_v_mouseZygote$DESig
rat_deg <- mast_rat$ratEgg_v_ratZygote$DESig

#1360 up and 3002 down, total 4362 genes in mice
mouse_deg_up <- row.names(subset(mouse_deg, Log2FC > log2(2) & fdr < 0.05))
mouse_deg_down <- row.names(subset(mouse_deg, Log2FC < -log2(2) & fdr < 0.05))

#711 up and 212 down, total 923 genes in rat
rat_deg_up <- row.names(subset(rat_deg, Log2FC > log2(2) & fdr < 0.05))
rat_deg_down <- row.names(subset(rat_deg, Log2FC < -log2(2) & fdr < 0.05))

# Stacked Number of DEG barplot
deg_summary <- data.frame(Species =  c(rep("Mouse" , 2) , rep("Rat" , 2)),
                          deg = rep(c('> 1' , '< -1'), 2),
                          value = c(length(mouse_deg_up), length(mouse_deg_down), length(rat_deg_up), length(rat_deg_down)))
deg_summary_plot <- ggplot(deg_summary, aes(fill=deg, y=value, x=Species, label = value)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(legend.position = c(.98, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.margin = margin(6, 6, 6, 6), 
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        axis.title = element_text(size = 14)) +
  guides(fill=guide_legend(title=expression('| log'[2]*'FC |'))) +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) + 
  scale_y_continuous(name="Number of DEGs", labels = scales::label_number(scale_cut = scales::cut_short_scale()))+scale_fill_manual(values = c('#6EE2FF', '#FF410D'))

ggsave('./number_of_DEGs.png', plot = deg_summary_plot, width = 2.5, height = 3.5)




##Nascent Transcript Analysis with kb estimated spliced/unspliced reads
#mouse_intron <- read_filter_splice('./velocity/mouse_spliced.csv', './velocity/mouse_unspliced.csv')
#rat_intron <- read_filter_splice('./velocity/rat_spliced.csv', './velocity/rat_unspliced.csv')

# Nascent Transcript Analysis with velocyto estimated spliced/unspliced reads
mouse_intron <- read_filter_splice_loom('./dataset/output_mouse/kb/mouse_trimmed_kb/', row.names(mouse$meta), f_format = 'kb_loom', filt_samp_by_cor = 0.8)
rat_intron <- read_filter_splice_loom('./dataset/output_rat/kb/rat_trimmed_kb/', row.names(rat$meta), f_format = 'kb_loom', filt_samp_by_cor = 0.8)



mouse_intron$genes <- intersect(row.names(mouse$ct$bio), mouse_intron$genes)
rat_intron$genes <- intersect(row.names(rat$ct$bio), rat_intron$genes)

for(n in colnames(rat_intron$spliced_raw)){
  new_df <- data.frame(cbind(rat_intron$spliced_raw[,n], rat_intron$unspliced_raw[,n]))
  colnames(new_df) <- c('kb_estimated_spliced_counts', 'kb_estimated_nascent_counts')
  write.table(x = new_df, file = paste('kb_nascent_geo/',n, '.kb.txt', sep = ''), sep = '\t', row.names = T, quote = F)
}





# Stacked number of genes with nascent expression barplot
nascent_summary <- data.frame(Species =  c(rep("Mouse" , 2) , rep("Rat" , 2)),
                          deg = rep(c('w Nascent reads' , 'wo Nascent reads'), 2),
                          value = c(length(mouse_intron$genes), sum(rowSums(mouse$ct$bio > 0) > 0)-length(mouse_intron$genes), length(rat_intron$genes), sum(rowSums(rat$ct$bio > 0) > 0)-length(rat_intron$genes)))
nascent_summary_plot <- ggplot(nascent_summary, aes(fill=deg, y=value, x=Species, label = value)) + 
  geom_bar(position="stack", stat="identity")+theme_classic()+
  theme(legend.position="bottom",legend.direction="vertical",
        legend.text=element_text(size=10),
        legend.title = element_text(size = 1),
        legend.margin = margin(6, 10, 6, 6), 
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        axis.title = element_text(size = 14))+
  guides(fill=guide_legend(title=expression('')))+
  geom_text(size = 4, position = position_stack(vjust = 0.5))+ 
  scale_y_continuous(name="Number of genes", labels = label_number(scale_cut = cut_short_scale()))+scale_fill_manual(values = c('#6EE2FF', '#FF410D'))

ggsave('./number_of_NascentGenes.png', plot = nascent_summary_plot, width = 2.5, height = 4)




perc_mouse_genes_unspliced <- colSums(mouse_intron$unspliced[mouse_intron$genes,])/(colSums(mouse_intron$spliced_raw[, colnames(mouse_intron$unspliced)])+colSums(mouse_intron$unspliced_raw[mouse_intron$genes,colnames(mouse_intron$unspliced)]))
perc_rat_genes_unspliced <- colSums(rat_intron$unspliced[rat_intron$genes,])/(colSums(rat_intron$spliced_raw[,colnames(rat_intron$unspliced)] )+colSums(rat_intron$unspliced_raw[rat_intron$genes,colnames(rat_intron$unspliced)]))



mouse_rat_nascent_df <- data.frame(percentages = c(perc_mouse_genes_unspliced*100, perc_rat_genes_unspliced*100),
           organism = c(rep('Mouse', length(perc_mouse_genes_unspliced)), rep('Rat', length(perc_rat_genes_unspliced))),
           cellType = c(mouse$meta[names(perc_mouse_genes_unspliced),'cellType'], rat$meta[names(perc_rat_genes_unspliced),'cellType']), 
           tissue = as.factor(sapply(c(mouse$meta[names(perc_mouse_genes_unspliced),'cellType'], rat$meta[names(perc_rat_genes_unspliced),'cellType']), function(x) ifelse(grepl('Egg', x), 'Oocyte', 'Zygote')))
)





nascent_perc_plot <- mouse_rat_nascent_df %>% reshape2::melt(c('cellType', 'organism', 'tissue'))   %>% ggbarplot( x = 'tissue', y = "value", color = "cellType",fill='cellType', facet.by = 'organism', add.params = list(color = '#ffa500'), 
                                                                    palette = DOT_COLOR, add = "mean_se",position = position_dodge(0.8))+
  stat_compare_means( method = 'wilcox.test', label = "p.signif", label.y.npc = c(0.8, 0.8), label.x.npc = c(0.05, 0.4), size = 9)+
  xlab('Cell Type')+ylab('Perctange of Nascent Reads (%)')+theme_classic2(base_size = 15)+ 
  theme(legend.position="none", 
        axis.text.x = element_text(size = 14, color = 'black'),
        strip.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 14, color = 'black'))+scale_y_continuous(name="Perctange of Nascent Reads (%)", labels = label_number(scale_cut = cut_short_scale()))
ggsave('percent_gene_unspliced.png', plot = nascent_perc_plot, width = 4,height = 3.5)





# Percentage of intergenic reads in both species
mouse_intergenic1 <- read_htseq_intergenic('./dataset/output_mouse/mouse_intergenic_1000.star.ct', ct = mouse$ct$bio, meta = mouse$meta, cov_thresh = 1)
rat_intergenic1 <- read_htseq_intergenic('./dataset/output_rat/rat_intergenic_1000.star.ct', ct = rat$ct$bio, meta = rat$meta, 1)
mouse_rat_int_df <- rbind(rat_intergenic1$small_df, mouse_intergenic1$small_df)[c(colnames(rat_intron$spliced), colnames(mouse_intron$spliced)),]
mouse_rat_int_df$organism <- c(rep('Rat', 25), rep('Mouse', 28))
colnames(mouse_rat_int_df)[c(1,2)] <- c('> 10Kb from Gene', 'All Regions')
intergenic_perc_plot <- mouse_rat_int_df %>% melt(c('cellType', 'organism')) %>% ggbarplot( x = "variable", y = "value", color = "cellType",fill='cellType', facet.by = 'organism', add.params = list(color = '#ffa500'), 
                                                                    palette = DOT_COLOR, add = "mean_se",position = position_dodge(0.8))+
  stat_compare_means(aes(group = cellType), method = 'wilcox.test', label = "p.signif", label.y.npc = c(0.6, 0.6, 0.6 ,0.6), size = 9)+
  xlab('Intergenic Location')+ylab('Number of Intergenic Regions w Coverage')+theme_classic2(base_size = 15)+ 
  theme(legend.position="none", axis.text.x = element_text(size = 12, color = 'black', angle = 15, vjust = 0.5),
        strip.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 14, color = 'black'))+
  scale_y_continuous(name="No. Regions w Coverage", labels = label_number(scale_cut = cut_short_scale()))
ggsave('intergenic_regions.png', plot=intergenic_perc_plot, width = 4,height = 3.5)

mouse_unsplic_mast <- mast_diff(ct = mouse_intron$unspliced[mouse_intron$genes,], meta = mouse$meta[colnames(mouse_intron$unspliced),], normFactor = colMeans(mouse_intron$spliced/edgeR::cpm(mouse_intron$spliced), na.rm = T),control = 'mouseEgg', tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, min_cell_grp = 2, min_cell = 5, max_thres = 6, plot = T)[[1]]$DESig
mouse_unsplic_mast_ora_up <- enrich_CP(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & Log2FC > 1)), 'mouse', universe = mouse_unsplic_mast$features, n_type = 'ALIAS', full_combine = T)
mouse_unsplic_mast_ora_down <- enrich_CP(row.names(subset(mouse_unsplic_mast, fdr < 0.05 & Log2FC < -1)), 'mouse', universe = mouse_unsplic_mast$features, n_type = 'ALIAS',  full_combine = T)
mouse_unsplic_mast_gsea <- gse_CP( logFC = sapply(row.names(mouse_unsplic_mast), FUN =function(x){mouse_unsplic_mast[x,'Log2FC']}), organisms = 'mouse', n_type = 'ALIAS',  full_combine = T)

rat_unsplic_mast <- mast_diff(ct = rat_intron$unspliced[rat_intron$genes,], meta = rat$meta[colnames(rat_intron$unspliced),], normFactor = colMeans(rat_intron$spliced/edgeR::cpm(rat_intron$spliced), na.rm = T), control = 'ratEgg', tpm = F, nbins = 0, min_per_bin = 50, freq = 0.1, min_cell_grp = 2, min_cell = 5, max_thres = 6, plot = T)[[1]]$DESig
rat_unsplic_mast_ora_up <- enrich_CP(row.names(subset(rat_unsplic_mast, fdr < 0.05 & Log2FC > 1)), 'rat', universe = rat_unsplic_mast$features,n_type = 'ALIAS',full_combine = T)
rat_unsplic_mast_ora_down <- enrich_CP(row.names(subset(rat_unsplic_mast, fdr < 0.05 & Log2FC < -1)), 'rat', universe = rat_unsplic_mast$features,n_type = 'ALIAS', full_combine = T)
rat_unsplic_mast_gsea <- gse_CP(logFC = sapply(row.names(rat_unsplic_mast), FUN =function(x){rat_unsplic_mast[x,'Log2FC']}), organisms = 'rat', n_type = 'ALIAS', full_combine = T)


mast_gse_spliced <- gse_CP(logFC = sapply(setdiff(row.names(mast_mouse$mouseEgg_v_mouseZygote$DESig), row.names(subset(mouse_unsplic_mast, fdr < 0.05))), FUN =function(x){mast_mouse$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), organisms = 'mouse', GSE = T)



unsplic_splic_log2FC <- data.frame(x=c(mouse_unsplic_mast$Log2FC, rat_unsplic_mast$Log2FC), 
                                   y=c(mast_mouse$mouseEgg_v_mouseZygote$DESig[row.names(mouse_unsplic_mast),'Log2FC'], mast_rat$ratEgg_v_ratZygote$DESig[row.names(rat_unsplic_mast),'Log2FC']), 
                                   Species = c(rep('Mouse', length(mouse_unsplic_mast$Log2FC)), rep('Rat', length(rat_unsplic_mast$Log2FC))))


# LogFC comparison between nascent expression and gene (mature) expression
nascent_mature_fc_plot <- ggplot(data = unsplic_splic_log2FC, aes(x = x, y = y)) + 
  geom_point()+
  facet_wrap(~Species)+
  theme_classic()+
  ylab(expression('Mature log'[2]*'FC'))+
  xlab(expression('Nascent Log'[2]*'FC'))+
  geom_text(data = data.frame(x = c(-2, -2),  y = c(4,4), lab = paste('PCC:',c(round(cor(mouse_unsplic_mast$Log2FC, mast_mouse$mouseEgg_v_mouseZygote$DESig[row.names(mouse_unsplic_mast),'Log2FC'],use = 'na.or.complete', method = 'pearson'),3),
                                                                              round(cor(rat_unsplic_mast$Log2FC, mast_rat$ratEgg_v_ratZygote$DESig[row.names(rat_unsplic_mast),'Log2FC'], use = 'na.or.complete',method = 'pearson'),3)), sep = ' '),
                              Species = c('Mouse', 'Rat')),mapping = aes(x = x, y = y, label = lab), size = 4)+
  theme(axis.text.x = element_text(size = 14, color = 'black'), 
        axis.text.y = element_text(size = 14, color = 'black'), 
        axis.title.y = element_text(size = 16, color = 'black'),
        axis.title.x = element_text(size = 16, color = 'black'),
        strip.text.x = element_text(size = 16))
ggsave('nascent_mature_log2FC.png', plot = nascent_mature_fc_plot, height = 3.5, width =4)




mouse_intron_prop <- fisher_proportion_test(mouse_intron$unspliced[mouse_intron$genes,], mouse_intron$spliced[mouse_intron$genes,], mouse$meta[colnames(mouse_intron$unspliced),]$cellType)
rat_intron_prop <- fisher_proportion_test(rat_intron$unspliced[rat_intron$genes,], rat_intron$spliced[rat_intron$genes,], rat$meta$cellType, 'ratEgg')

mouse_unsplic_prop_ora_up <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'mouse', universe = row.names(mouse_intron_prop),  full_combine = T, n_type = 'ALIAS')
mouse_unsplic_prop_ora_down <- enrich_CP(row.names(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff < 0)), 'mouse', universe = row.names(mouse_intron_prop), full_combine = T, n_type = 'ALIAS')


rat_unsplic_prop_ora_up <- enrich_CP(row.names(subset(rat_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'rat', universe = row.names(rat_intron_prop), full_combine = T, n_type = 'ALIAS')
rat_unsplic_prop_ora_down <- enrich_CP(row.names(subset(rat_intron_prop, qvalue < 0.05 & prop_diff > 0)), 'rat', universe = row.names(rat_intron_prop),  full_combine = T, n_type = 'ALIAS')



# Stacked number of genes with nascent expression barplot
dnp_summary <- data.frame(Species =  c(rep("Mouse" , 2) , rep("Rat" , 2)),
                              deg = rep(c('Increase' , 'Decrease'), 2),
                              value = c(nrow(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff > 0)),
                                        nrow(subset(mouse_intron_prop, qvalue < 0.05 & prop_diff < 0)), 
                                        nrow(subset(rat_intron_prop, qvalue < 0.05 & prop_diff > 0)), 
                                        nrow(subset(rat_intron_prop, qvalue < 0.05 & prop_diff < 0))))
dnp_summary_plot <- ggplot(dnp_summary, aes(fill=deg, y=value, x=Species, label = value)) + 
  geom_bar(position="stack", stat="identity")+theme_classic()+
  theme(legend.position = c(1.06, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.margin = margin(6, 6, 6, 6), 
        axis.text.x = element_text(size=14, color = 'black'),
        axis.text.y = element_text(size=14, color = 'black'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  guides(fill=guide_legend(title=expression('Nascent proportion')))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))+ 
  scale_y_continuous(name="Number of DNP genes", labels = label_number(scale_cut = cut_short_scale()))+scale_fill_manual(values = c('#6EE2FF', '#FF410D'))
ggsave('./number_of_DNP.png', width = 3.5, height = 4)







dnp_deg_mouse <- list('DEG Up'= row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'DEG Down'= row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'DNP Increase'= row.names(subset(mouse_intron_prop, prop_diff > 0 & qvalue < 0.05)),
                      'DNP Decrease' =row.names(subset(mouse_intron_prop, prop_diff < 0 & qvalue < 0.05)))
make_comb_mat(dnp_deg_mouse)
png('dnp_deg_mouse.png', width = 3, height = 3, res = 300, units = 'in')
UpSet(make_comb_mat(dnp_deg_mouse)[1:4], comb_col = c('black'))
dev.off()

dnp_deg_rat <- list('DEG Up'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                    'DEG Down'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                    'DNP Increase'= row.names(subset(rat_intron_prop, prop_diff > 0 & qvalue < 0.05)),
                    'DNP Decrease' =row.names(subset(rat_intron_prop, prop_diff < 0 & qvalue < 0.05)))
png('dnp_deg_rat.png', width = 3, height = 3, res = 300, units = 'in')
UpSet(make_comb_mat(dnp_deg_rat)[1:2], comb_col = c('black'))
dev.off()

oog2_mouse <- plot_range_coverage(range = 'NC_000070.7:143,920,015-143,923,663',txdb=mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse', y_lim = c(0.08, 0.73))
obox2_mouse <- plot_range_coverage(range = 'NC_000073.7:15,130,780-15,132,520', txdb=mouse_gtf, bw_file_list= c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse', y_lim = c(0.08, 0.73))
zar1_rat <- plot_range_coverage(range = 'NC_051349.1:35,194,801-35,196,994', txdb=rat_gtf, bw_file_list= c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), log = T, cols = DOT_COLOR[c('ratZygote', 'ratEgg')], cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat', y_lim = c(0.08, 0.73))
obox5_mouse <- plot_range_coverage(range = 'NC_000073.7:15,491,480-15,491,980', txdb=mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('Mouse zygote', 'Mouse oocyte' ), file_name_suffix = 'mouse', y_lim = c(0.08, 0.73))
nlrp4g_mouse <- plot_range_coverage(range = 'NC_000075.7:124,117,990-124,121,100', txdb=mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse', y_lim = c(0.08, 0.73))

oog1_mouse <- plot_range_coverage(txdb=mouse_gtf,  bw_file_list= c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse', gene_name = 'Oog1', y_lim = c(0.08, 0.73))
oog1_rat <- plot_range_coverage( txdb=rat_gtf, bw_file_list= c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), log = T, cols = DOT_COLOR[c('ratZygote', 'ratEgg')], cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat', gene_name = 'Oog1', y_lim = c(0.08, 0.73))
test_vp <- plot_range_coverage(range = 'NC_000073.7:10448500-10464076', txdb=mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), log = T, cols = DOT_COLOR[c('mouseZygote', 'mouseEgg')], cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse', y_lim = c(0.08, 0.73), gene_name = 'Nlrp4b')
test_vp <- plot_range_coverage(range= 'NC_051336.1:70426100-70446800', txdb=rat_gtf, bw_file_list= c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), log = T, cols = DOT_COLOR[c('ratZygote', 'ratEgg')], cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat', y_lim = c(0.08, 0.73))

celltype_bar_legend <- get_legend(ggplot(data.frame(y=c(0.1, 0.1, 0.1, 0.1), 
                  x=c(1,2,3,4), 
                  cellType = c('mouseEgg', 'mouseZygote', 'ratEgg', 'ratZygote')), aes(x=x, y=y, fill=cellType))+geom_bar(stat='identity')+
  scale_fill_manual(values = DOT_COLOR, labels = c('Mouse zygote', 'Mouse oocyte', 'Rat zygote', 'Rat oocyte' ))+guides(fill=guide_legend(title=expression('')))+
    theme(legend.direction = 'horizontal')
)

ggdraw(celltype_bar_legend)

ggsave('./celltype_bar_legend.png', plot = ggdraw(celltype_bar_legend),width = 5, height = 2)

## DEXSeq DTU
mouse_dtu <- read_rnasplice_dex_dtu('./dataset/output_mouse/splice/salmon/dexseq_dtu/results/dexseq/DEXSeqDataSet.MF-MU.rds', './dataset/output_mouse/splice/salmon/dexseq_dtu/filter/drimseq/dmDSdata.rds')
rat_dtu <- read_rnasplice_dex_dtu('./dataset/output_rat/splice/salmon/dexseq_dtu/results/dexseq/DEXSeqDataSet.RF-RU.rds', './dataset/output_rat/splice/salmon/dexseq_dtu/filter/drimseq/dmDSdata.rds')
mouse_dtu$prop <- get_dexseq_prop(mouse_dtu)
rat_dtu$prop <- get_dexseq_prop(rat_dtu)

system.time({
  mouse_dtu$dexseq = estimateSizeFactors(mouse_dtu$dexseq)
  mouse_dtu$dexseq = estimateDispersions(mouse_dtu$dexseq)
  mouse_dtu$dexseq = testForDEU(mouse_dtu$dexseq, reducedModel = ~sample + exon)
})


system.time({
  rat_dtu$dexseq = estimateSizeFactors(rat_dtu$dexseq)
  rat_dtu$dexseq = estimateDispersions(rat_dtu$dexseq)
  rat_dtu$dexseq = testForDEU(rat_dtu$dexseq, reducedModel = ~sample + exon)
})


dexseq_mouse <- DEXSeqResults(mouse_dtu$dexseq, independentFiltering = T) #read.csv('./splicing/mouse/suppa/mouse_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqResults.MF-MU.tsv', header = T,  sep = '\t')
dexseq_rat <- DEXSeqResults(rat_dtu$dexseq, independentFiltering = T) #read.csv('./splicing/rat/suppa/rat_suppa_drimseq/dexseq_dtu/results/dexseq/DEXSeqResults.RF-RU.tsv', header = T,  sep = '\t')
dexseq_mouse_stageR <- stageR_dexseqRes(dexseq_mouse)
dexseq_rat_stageR <- stageR_dexseqRes(dexseq_rat)

##make venn diagrams
dtu_deg_mouse <- list('DEG up'= row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                      'DEG down'= row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                      'DTU'= unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID))
dtu_deg_rat <- list('DEG up'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                    'DEG down'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                    'DTU'= unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID))

mouse_venn_deg_dtu <- plotVenn(dtu_deg_mouse, outFile='./dtu_deg_venn_mouse.svg', nCycles = 10000, fontScale = 3, labelRegions = F)
rat_venn_deg_ddtu <- plotVenn(dtu_deg_rat, outFile='./dtu_deg_venn_rat.svg', nCycles = 10000,  fontScale = 3, labelRegions = F)


#DTU genes summary BARPLOT
DTU_summary <- data.frame(Species =  c(rep("Mouse" , 2) , rep("Rat" , 2)),
                                     deg = rep(c('DTU', 'Non-DTU'),2),
                                     value = c(length(unique(subset(dexseq_mouse_stageR, gene < 0.05)$geneID)), 
                                               length(unique(subset(dexseq_mouse_stageR, gene >= 0.05)$geneID)),
                                               length(unique(subset(dexseq_rat_stageR, gene < 0.05)$geneID)),
                                               length(unique(subset(dexseq_rat_stageR, gene >= 0.05)$geneID))))

DTU_summary_plot <- ggplot(DTU_summary, aes(fill=deg, y=value, x=Species, label = value)) + 
  geom_bar(position="stack", stat="identity")+theme_classic()+
  theme(legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.margin = margin(6, 10, 6, 6), 
        axis.text.x = element_text(size=14, color = 'black'),
        axis.text.y = element_text(size=14, color = 'black'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  guides(fill=guide_legend(title=expression('')))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))+ 
  scale_y_continuous(name="Number of genes", labels = label_number(scale_cut = cut_short_scale()))+scale_fill_manual(values = paletteer::paletteer_d("ggsci::springfield_simpsons"))

ggsave('./number_of_DTU.png', plot= DTU_summary_plot, width = 3.7, height = 4.2)




## Num Transcript DTU summary BARPLOT
mouse_det_bygene <- table(tapply(subset(dexseq_mouse_stageR, gene < 0.05), subset(dexseq_mouse_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)}))
rat_det_bygene <- table(tapply(subset(dexseq_rat_stageR, gene < 0.05), subset(dexseq_rat_stageR, gene < 0.05)$geneID, FUN = function(x){sum(x$transcript < 0.05)}))
transcript_DTU_summary <- data.frame(Species =  c(rep("Mouse" , length(mouse_det_bygene)) , rep("Rat" , length(rat_det_bygene))),
                              deg = c(names(mouse_det_bygene), names(rat_det_bygene)),
                              value = c(mouse_det_bygene, rat_det_bygene), 
                              label = c(mouse_det_bygene, rat_det_bygene))
transcript_DTU_summary$label[transcript_DTU_summary$label < 50] <- NA
num_sig_tx_DTU <- ggplot(transcript_DTU_summary, aes(fill=deg, y=value, x=Species, label = label)) + 
  geom_bar(position="stack", stat="identity")+theme_classic()+
  theme(legend.position = c(1.07, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.margin = margin(6, 10, 6, 6), 
        axis.text.x = element_text(size=14, color = 'black'),
        axis.text.y = element_text(size=14, color = 'black'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  guides(fill=guide_legend(title=expression('No. Sig Transcripts')))+
  geom_text(size = 4, position = position_stack(vjust = 0.5))+ 
  scale_y_continuous(name="Number of DTU genes", labels = label_number(scale_cut = cut_short_scale()))+scale_fill_manual(values = paletteer::paletteer_d("ggsci::springfield_simpsons"))

ggsave('./number_of_sig_tx_DTU.png', plot = num_sig_tx_DTU, width = 3.7, height = 4.2)

# Enrich GO BP 
stager_mouse_ora <- enrich_CP(subset(dexseq_mouse_stageR, gene < 0.05)$geneID, universe = dexseq_mouse_stageR$geneID, organisms = 'mouse', full_combine = T)
stager_rat_ora <- enrich_CP(subset(dexseq_rat_stageR, gene < 0.05)$geneID, universe = dexseq_rat_stageR$geneID, organisms = 'rat', full_combine = T)

#make dotplots
rat_dtu_dot <- test_dotplot(stager_rat_ora$combined_full, showCategory=10, color='pvalue') 
rat_dtu_dot <- rat_dtu_dot + theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 12, colour = ifelse(rat_dtu_dot$data[order(rat_dtu_dot$data$GeneRatio),'qvalue'] < 0.05, 'red', 'black')), 
                                   axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 10), 
                                   legend.text = element_text(size=8), 
                                   legend.key.size = unit(1, 'cm'),
                                   legend.key.height = unit(0.3, 'cm'),
                                   legend.key.width = unit(0.3, 'cm'),
                                   legend.title = element_text(size=8),
                                   legend.position = c(0.6, 0.2), legend.box = 'horizontal')+xlab('Gene ratio')+
  scale_fill_continuous(low="red", high="blue", name = 'pvalue', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.005))+
  scale_size(range=c(3,8), limits = c(5,300))+annotate("text", label = 'Rat', x=0.03, y=10, size = 8)
ggsave('Rat_dex_dtu_GOBP.png',plot = rat_dtu_dot, width = 5, height = 5.5)

mouse_dtu_dot <- test_dotplot(stager_mouse_ora$combined_full, showCategory=10, color ='qvalue') 
mouse_dtu_dot <- mouse_dtu_dot + theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 12, colour = ifelse(mouse_dtu_dot$data[order(mouse_dtu_dot$data$GeneRatio),'qvalue'] < 0.05, 'red', 'black')), 
                                       axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
                                       legend.text = element_text(size=8), 
                                       legend.key.size = unit(1, 'cm'),
                                       legend.key.height = unit(0.3, 'cm'),
                                       legend.key.width = unit(0.3, 'cm'),
                                       legend.title = element_text(size=8),
                                       legend.position = c(0.6, 0.2), legend.box = 'horizontal')+xlab('Gene ratio')+
  scale_fill_continuous(low="red", high="blue", name = 'FDR', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.1))+
  scale_size(range=c(3,8), limits = c(5,300))+annotate("text", label = 'Mouse', x=0.03, y=10, size = 8)
ggsave('Mouse_dex_dtu_GOBP.png', plot = mouse_dtu_dot,width = 5, height = 5.5)


mouse_dtu$drimseq@samples$condition <- c('MU' = 'mouseEgg', 'MF' = 'mouseZygote')[mouse_dtu$drimseq@samples$condition]
rat_dtu$drimseq@samples$condition <- c('RU' = 'ratEgg', 'RF' = 'ratZygote')[rat_dtu$drimseq@samples$condition]

# Plot specific Genes

plotDEXSeqDTU(rat_dtu$prop, 'Foxm1', rat_dtu$drimseq@samples, isProportion = T)
plotDEXSeqDTU(mouse_dtu$prop, 'Abi3bp', mouse_drim$samps, isProportion = T)


plotDEXSeqDTU(mouse_dtu$prop, 'Bap1', mouse_dtu$drimseq@samples, isProportion = T)
ggsave('mouse_Bap1.png', width =3, height = 3)
plotDEXSeqDTU(mouse_dtu$prop, 'Usp3', mouse_dtu$drimseq@samples, isProportion = T)
ggsave('mouse_Usp3.png', width =3, height = 3)
plotDEXSeqDTU(rat_dtu$prop, 'Usp16', rat_dtu$drimseq@samples, isProportion = T)
ggsave('rat_Usp16.png', width =3, height = 3)
plotDEXSeqDTU(rat_dtu$prop, 'Rnf2', rat_dtu$drimseq@samples, isProportion = T)
ggsave('rat_Rnf2.png', width =3, height = 3)




## Dapars Differential APA (DAP)

#DAPARS fishers exact test custom with filtering
mouse_dapars_sf_st25 <- deg_utr2('./dataset/mouse_apa_res/apa2_e0_q50_st25.txt', mouse$ct$bio, c('mouseEgg', 'mouseZygote'), 
                                 mouse$meta, filter_by_PAS_motif = T,combine_p = 'simes',
                                 fasta_file = '../mouse_rat_proposal/dataset/igv/mouse/mouse_spike.fa',
                                 utr_cov_filt = 1, min_samp_filt = -1, na_filt = 3, expr_filt = 10, old_apa_dist = F, filter_by_pas_dist = 0
                                  )


utr_mouse_sf_st25_up <- enrich_CP(unique(subset(mouse_dapars_sf_st25$gene_res,  diff & mean.diff > 0.2)$gene_short_names), universe = unique(mouse_dapars_sf_st25$gene_res$gene_short_names),organisms = 'mouse', full_combine = T)
utr_mouse_sf_st25_down <- enrich_CP(unique(subset(mouse_dapars_sf_st25$gene_res,  diff & mean.diff < -0.2)$gene_short_names), universe = unique(mouse_dapars_sf_st25$gene_res$gene_short_names),  organisms = 'mouse', full_combine = T)




rat_dapars_sf_st25 <- deg_utr2('./dataset/rat_apa_res/apa2_e0_q50_st25.txt', rat$ct$bio, c('ratEgg', 'ratZygote'), 
                                 rat$meta, filter_by_PAS_motif = T, combine_p = 'simes',
                                 fasta_file = '../mouse_rat_proposal/dataset/igv/rat/rat_spike.fa',
                                 utr_cov_filt = 1, min_samp_filt = -1, na_filt = 3, expr_filt = 10, old_apa_dist = F, filter_by_pas_dist = 0
)

utr_rat_sf_st25_up <- enrich_CP(unique(subset(rat_dapars_sf_st25$gene_res, diff & mean.diff >= 0.2)$gene_short_names), universe = unique(rat_dapars_sf_st25$gene_res$gene_short_names),  organisms = 'rat', full_combine = T)
utr_rat_sf_st25_down <- enrich_CP(unique(subset(rat_dapars_sf_st25$gene_res, diff & mean.diff <= -0.2)$gene_short_names), universe = unique(rat_dapars_sf_st25$gene_res$gene_short_names), organisms = 'rat', full_combine = T)


#UMAP plots
rat_utr_umap <- sample_PCA(log(rat_dapars_sf_st25$pdui_imp+1), rat$meta, umap =  T, dimension=2 ,
           main = "Rat", labeling = F, point_size = 2, color_by = 'cellType',
           umap.config = list(n_neighbors = 25, min_dist = 0.5, metric='cosine', seed = 12345), legend.position = c(0.2, 0.8))
ggsave('rat_umap_DAP.png', plot = rat_utr_umap, units = 'in', dpi = 300, height = 3, width = 2.2)
mouse_utr_umap <- sample_PCA(log(mouse_dapars_sf_st25$pdui_imp+1), mouse$meta, umap =  T, dimension=2 ,
           main = "Mouse", labeling = F, point_size = 2, color_by = 'cellType', 
           umap.config = list(n_neighbors = 25, min_dist = 0.5, metric='cosine', seed = 12345),legend.position = c(0.8, 0.8))#+theme(plot.title = element_text(vjust = -5))
ggsave('mouse_umap_DAP.png', plot = mouse_utr_umap, units = 'in', dpi = 300, height = 3, width = 2.2)

for(n in colnames(mouse$ct$bio)){
  new_df <- cbind(mouse_dapars_sf_st25$df[,grepl(paste(n, '\\.', sep = ''), colnames(mouse_dapars_sf_st25$df))], mouse_dapars_sf_st25$df[,1:3])
  colnames(new_df)[1:3] <- c('long_utr_coverage', 'short_utr_coverage', 'PDUI')
  write.table(new_df, paste('apa_geo_submission/', n , '.dapars2.tsv', sep = ''), sep = '\t', quote = F, row.names = T)
}


for(n in colnames(rat$ct$bio)){
  new_df <- cbind(rat_dapars_sf_st25$df[,grepl(paste(n, '\\.', sep = ''), colnames(rat_dapars_sf_st25$df))], rat_dapars_sf_st25$df[,1:3])
  colnames(new_df)[1:3] <- c('long_utr_coverage', 'short_utr_coverage', 'PDUI')
  write.table(new_df, paste('apa_geo_submission/', n , '.dapars2.tsv', sep = ''), sep = '\t', quote = F, row.names = T)
}


# Maybe include PAS Signal filter
#mouse_PAS_dist <- post_dapars_pas_filter(mouse_dapars_sf_st25$deg, './dataset/igv/mouse/mouse_spike.fa', up_range = 80, down_range = 120, offset = 0)
#rat_PAS_dist <- post_dapars_pas_filter(rat_dapars_sf_st25$deg, './dataset/igv/rat/rat_spike.fa', up_range = 80, down_range = 120, offset = 0)

#DAP summary barplot
DAP_summary <- data.frame(Species =  c(rep("Mouse" , 2) , rep("Rat" , 2)),
                          deg = rep(c('Lengthened UTR', 'Shortened UTR'),2),
                          value = c(length(unique(subset(mouse_dapars_sf_st25$gene_res, mean.diff > 0.2 & diff)$gene_short_names)), 
                                    length(unique(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & diff)$gene_short_names)),
                                    length(unique(subset(rat_dapars_sf_st25$gene_res, mean.diff > 0.2 & diff)$gene_short_names)),
                                    length(unique(subset(rat_dapars_sf_st25$gene_res, mean.diff < -0.2 & diff)$gene_short_names))))

DAP_summary_plot <- ggplot(DAP_summary, aes(fill=deg, y=value, x=Species, label = value)) + 
  geom_bar(position="stack", stat="identity")+theme_classic()+
  theme(legend.position = c(1.06, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        legend.margin = margin(6,10, 6, 6), 
        axis.text.x = element_text(size=14, color = 'black'),
        axis.text.y = element_text(size=14, color = 'black'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  guides(fill=guide_legend(title=expression('')))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))+ 
  scale_y_continuous(name="Number of genes", labels = label_number(scale_cut = cut_short_scale()))+scale_fill_manual(values = c('#1A9993FF', '#F05C3BFF'))

ggsave('./number_of_DTU.png', plot = DAP_summary_plot, width = 3.7, height = 4.2)



utr_deg_mouse <- list('DEG Up'= row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2))),
                      'DEG Down'= row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -log2(2))),
                      'Long UTR'= unique(subset(mouse_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)$gene_short_name),
                      'Short UTR' =unique(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)$gene_short_name))




utr_deg_rat <- list('DEG Up'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC > 1)),
                    'DEG Down'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, fdr < 0.05 & Log2FC < -1)),
                    'Long UTR'= unique(subset(rat_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05 )$gene_short_name),
                    'Short UTR' =unique(subset(rat_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)$gene_short_name))

png('utr_deg_mouse.png', width = 3, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_mouse)[1:4], comb_col = c('black', 'red', 'black','black'))
dev.off()

png('utr_deg_rat.png', width = 2.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_deg_rat)[2:4], comb_col = c('black', 'black','black'))
dev.off()

# utr_v_deg <- data.frame(row.names = row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)), utr=subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)$mean.diff, 
#                         deg = mast_mouse$mouseEgg_v_mouseZygote$DESig[row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)),]$Log2FC
#                         )
# utr_v_deg$deg_sig <- 'FALSE'
# utr_v_deg[intersect(row.names(subset(mouse_dapars_sf_st25$gene_res,mean.diff < -0.2 & fdr < 0.05)), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2)))), 'deg_sig'] <- 'TRUE'
# utr_v_deg$alpha <- c('TRUE'=1, 'FALSE'=0.3)[utr_v_deg$deg_sig]



## Compare all genes that have DAP analysis D+PDUI change vs Log2FC, no correlation but association
utr_v_deg2 <- data.frame(row.names = row.names(mouse_dapars_sf_st25$gene_res), utr=mouse_dapars_sf_st25$gene_res$mean.diff, 
                        deg = mast_mouse$mouseEgg_v_mouseZygote$DESig[row.names(mouse_dapars_sf_st25$gene_res),]$Log2FC
)

utr_v_deg2$deg_sig <- 'Non'
utr_v_deg2[intersect(row.names(mouse_dapars_sf_st25$gene_res), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC > log2(2)))), 'deg_sig'] <- 'DEG up'
utr_v_deg2[intersect(row.names(mouse_dapars_sf_st25$gene_res), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, fdr < 0.05 & Log2FC < -log2(2)))), 'deg_sig'] <- 'DEG down'

utr_v_deg2$utr_sig <- 'Non'
utr_v_deg2[row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)), 'utr_sig'] <- 'Shortened UTR'
utr_v_deg2[row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)), 'utr_sig'] <- 'Lengthened UTR'


utr_v_deg2$utr_sig2 <- 'Non'
utr_v_deg2[row.names(subset(mouse_dapars_sf_st25$gene_res, abs(mean.diff) > 0.2 & fdr < 0.05)), 'utr_sig2'] <- 'Sig DAP'


# code for producing combined density plots with different legends and color groupings
pmain <- ggplot(utr_v_deg2, aes(x = utr, y = deg, color = deg_sig, alpha = utr_sig, size = utr_sig2)) +
  geom_point(aes(fill = utr_sig)) + theme_classic()+
  scale_color_manual(values=c('DEG down'="blue",'Non'= "grey",'DEG up'="red"))+ 
  theme(legend.position = 'none',
        axis.text.x = element_text(size=19, color = 'black'),
        axis.text.y = element_text(size=19, color = 'black'),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))+
  xlab('UTR PDUI difference')+
  ylab(expression('Log'[2]*'FC'))+
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1, alpha = 0.8)+
  geom_vline(xintercept=c(-0.2, 0.2), linetype="dashed", color = "black", size = 1, alpha= 0.8)+
  scale_alpha_manual(values = c('Shortened UTR'=1, 'lengthened UTR'=1, 'Non' = 0.04))+
  scale_size_manual(values = c('Sig DAP'=1.5, 'Non'=0.5))

# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x") +
  geom_density(data = subset(utr_v_deg2, deg_sig != 'Non'), aes(x = utr, fill = deg_sig),
               alpha = 0.7, size = 0.2) +
  scale_fill_manual(values=c('DEG down'="blue",'DEG up'="red"))+ theme(legend.position = 'none')
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_density(data = subset(utr_v_deg2, utr_sig != 'Non'), aes(x = deg, fill = utr_sig),
               alpha = 0.7, size = 0.2) + theme(legend.position = 'none')+
  coord_flip() +
  scale_fill_manual(values=c('Lengthened UTR'="#008331", 'Shortened UTR'='violet'))
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

legend_y <- cowplot::get_legend( ggplot(subset(utr_v_deg2, !is.na(utr_sig)), aes(x = utr, fill=utr_sig))+
                                   geom_density() + theme(legend.margin = margin(1, 130, -30, 1), 
                                                          legend.title = element_text(size = 14), 
                                                          legend.text=element_text(size=14))+
                                   scale_fill_manual(values=c('Lengthened UTR'="#008331", 'Shortened UTR'='violet'))+
  guides(fill = guide_legend(ncol = 1, title = 'Significant DAP')))

legend_x <- cowplot::get_legend( ggplot(subset(utr_v_deg2, !is.na(deg_sig)), aes(x = utr, fill=deg_sig))+
                                   geom_density() +theme(legend.margin = margin(1, 150, 1, 1), 
                                                         legend.title = element_text(size = 14),
                                                         legend.text=element_text(size=14))+
                                   scale_fill_manual(values=c('DEG down'="blue",'DEG up'="red"))+
                      guides(fill = guide_legend(ncol = 1, title = 'Significant DEG')))


utr_deg_mix_plot <- plot_grid(p2, plot_grid(legend_x, legend_y, ncol = 1), nrow = 1, rel_widths = c(1,0.3))
ggsave('utr_v_deg_scatter.png', plot = utr_deg_mix_plot , width = 8, height = 5)


# GO Bp DAP mouse genes dotplot
DAP_ora_dotplot <- test_dotplot(utr_mouse_sf_st25_down$combined_full, showCategory=10, color='qvalue') +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 8, colour =rep('red', 10)), 
  axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 8), 
  legend.text = element_text(size=8), legend.key.size = unit(1, 'cm'),  
  legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
  legend.title = element_text(size=8))+ 
  xlab('Gene ratio')+
  scale_fill_continuous(low="red", high="blue", name = 'FDR', guide=guide_colorbar(reverse=TRUE), limits=c(0,0.05))+
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

ggsave('Mouse_dap_GOBP.png', plot = DAP_ora_dotplot, width = 5, height = 6)


# UTR coverage plots
test_vp <- plot_utr_coverage('Btg4', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Cnot7', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Nek7', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Dazl', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Mastl', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Cnot6l', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat', loci = 'NC_051349.1:13495326-13502511')




test_vp <- plot_utr_coverage('Eif4e', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Usp28', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Bmi1', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Arid1a', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')
test_vp <- plot_utr_coverage('Anapc1', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')



test_vp <- plot_utr_coverage('Eif4e', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Usp28', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Bmi1', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat')

test_vp <- plot_utr_coverage('Arid1a', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Cdc40', utr_res = rat_dapars_sf_st25$gene_res, rat_gtf, c('./dataset/rat_apa_res/RF.sf.bw', './dataset/rat_apa_res/RU.sf.bw'), cols= c(DOT_COLOR['ratZygote'], DOT_COLOR['ratEgg']),cell_names = c('Rat zygote', 'Rat oocyte'), file_name_suffix = 'rat')
test_vp <- plot_utr_coverage('Cdc40', utr_res = mouse_dapars_sf_st25$gene_res, mouse_gtf, c('./dataset/mouse_apa_res/MF.sf.bw', './dataset/mouse_apa_res/MU.sf.bw'), cols= c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']),cell_names = c('Mouse zygote', 'Mouse oocyte'), file_name_suffix = 'mouse')




## Orthology



deg_ortho <- list('Rat Up DEG'= row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)),
                  'Rat Down DEG' =row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25)  & fdr < 0.05)),
                  'Mouse Up DEG'= row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)),
                  'Mouse Down DEG' =row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25)  & fdr < 0.05)))

mu_ru <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)))
md_ru <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)))
mu_rd <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, Log2FC > log2(1.25) & fdr < 0.05)))
md_rd <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, Log2FC < -log2(1.25) & fdr < 0.05)))
m_r <- intersect(row.names(subset(mast_rat$ratEgg_v_ratZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)), row.names(subset(mast_mouse$mouseEgg_v_mouseZygote$DESig, abs(Log2FC) > log2(1.25) & fdr < 0.05)))

m_r_state <- c(rep('Mouse Up / Rat Up', length(mu_ru)),rep('Mouse down / Rat Up', length(md_ru)),rep('Mouse up / Rat down', length(mu_rd)),rep('Mouse down / Rat down', length(md_rd)))

names(m_r_state) <- c(mu_ru, md_ru,mu_rd,md_rd)

mouse_rat_all_deg <- intersect(row.names(mast_rat$ratEgg_v_ratZygote$DESig), row.names(mast_mouse$mouseEgg_v_mouseZygote$DESig))

deg_ortho_ora <- list('mu_ru' = enrich_CP(mu_ru,universe = mouse_rat_all_deg, organisms = 'mouse'),
                      'mu_rd' = enrich_CP(mu_rd,universe = mouse_rat_all_deg,  organisms = 'mouse'),
                      'md_ru' = enrich_CP(md_ru,universe = mouse_rat_all_deg, organisms = 'mouse'),
                      'md_rd'= enrich_CP(md_rd,universe = mouse_rat_all_deg,  organisms = 'mouse'),
                      'm_r'= enrich_CP(m_r,universe = mouse_rat_all_deg,  organisms = 'mouse')
                      )

m_r_ora <- enrich_CP(m_r,universe = mouse_rat_all_deg,organisms = 'mouse', n_type = 'ALIAS', full_combine = T)
m_r_ora$combined <- m_r_ora$GO_BP_ora
m_r_ora$combined@result <- rbind(m_r_ora$GO_BP_ora@result, m_r_ora$WKP_ora@result,  m_r_ora$REACT_ora@result)
m_r_ora$combined@geneSets <- c(m_r_ora$GO_BP_ora@geneSets, m_r_ora$REACT_ora@geneSets,m_r_ora$KEGG_ora@geneSets, m_r_ora$WKP_ora@geneSets)
m_r_ora_cnetplots <- custom_cnet_plot(m_r_ora$combined_full, top_n_cat = 13, seed = 777,
                                      category = c('R-MMU-1428517', 'GO:0006338', 'R-MMU-72172', 
                                                   'R-MMU-72202', 'GO:0044772', 'GO:0006364',
                                                   'GO:0008380', 'GO:0016574', 'R-MMU-157118',
                                                   'R-MMU-4551638', 'GO:0006413', 'R-MMU-2990846', 
                                                   'R-MMU-927802'),
                                      gene_color = sapply(mast_rat$ratEgg_v_ratZygote$DESig$features, FUN = function(x){mast_rat$ratEgg_v_ratZygote$DESig[x,'Log2FC']}), 
                                      gene_color2 = sapply(mast_mouse$mouseEgg_v_mouseZygote$DESig$features, FUN = function(x){mast_mouse$mouseEgg_v_mouseZygote$DESig[x,'Log2FC']}), 
                                      layout = 'fr', color_cat_pval = T)

ggsave(plot = m_r_ora_cnetplots$plot2+ theme(legend.position = 'none',legend.box = "horizontal"), 'm_r_ora_mouse.png', width = 6, height = 3)
ggsave(plot = m_r_ora_cnetplots$plot1+ theme(legend.position = c(0.85,0.7),legend.box = "horizontal"), 'm_r_ora_rat.png', width = 6, height = 3)



png('deg_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(deg_ortho)[1:4], comb_col = c('black'))
dev.off()





dnp_ortho <- list('Mouse Increase DNP'= row.names(subset(mouse_intron_prop, prop_diff > 0 & qvalue < 0.05)),
                  'Mouse Decrease DNP' =row.names(subset(mouse_intron_prop, prop_diff < 0 & qvalue < 0.05)),
                  'Rat Increase DNP'= row.names(subset(rat_intron_prop, prop_diff > 0 & qvalue < 0.05)),
                  'Rat Decrease DNP' =row.names(subset(rat_intron_prop, prop_diff < 0 & qvalue < 0.05)))
png('dnp_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(dnp_ortho)[1:3], comb_col = c('black'))
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




utr_ortho <- list('Rat Long UTR'= row.names(subset(rat_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)),
                  'Rat Short UTR' =row.names(subset(rat_dapars_sf_st25$gene_res, mean.diff < -0.2 & fdr < 0.05)),
                  'Mouse Long UTR'= row.names(subset(mouse_dapars_sf_st25$gene_res, mean.diff > 0.2 & fdr < 0.05)),
                  'Mouse Short UTR' =row.names(subset(mouse_dapars_sf_st25$gene_res,mean.diff < -0.2 & fdr < 0.05)))


png('utr_ortho.png', width = 3.6, height = 2.6, res = 300, units = 'in')
UpSet(make_comb_mat(utr_ortho)[1:4], comb_col = c('black'))
dev.off()

