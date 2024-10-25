mouse$ct$bio

mouse_motif2tf <- read.csv('./dataset/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl', sep = '\t')

mouse_scenic_exprs <- mouse$ct$bio[rowMeans(mouse$ct$bio) > 5 & rowSums(mouse$ct$bio > 0) > 15,]
mouse_scenic_tfs <- intersect(mouse_motif2tf$gene_name, row.names(mouse_scenic_exprs))

rat_scenic_exprs <- rat$ct$bio[rowMeans(rat$ct$bio) > 5 & rowSums(rat$ct$bio > 0) > 15,]
rat_scenic_tfs <- intersect(mouse_motif2tf$gene_name, row.names(rat_scenic_exprs))

write.table(mouse_scenic_exprs, './pySCENIC_notebooks/mouse_exprs.txt', quote = F,  row.names = T, sep = '\t')
write.table(rat_scenic_exprs, './pySCENIC_notebooks/rat_exprs.txt', quote = F,  row.names = T, sep = '\t')
write(mouse_scenic_tfs, './pySCENIC_notebooks/mouse.tfs.txt', ncolumns = 1)
write(rat_scenic_tfs, './pySCENIC_notebooks/rat.tfs.txt', ncolumns = 1)
