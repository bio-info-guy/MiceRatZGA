# Load packages
suppressPackageStartupMessages({
  library(aplot)
  library(betareg)
  library(biomaRt)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(cowplot)
  library(data.table)
  library(DEXSeq)
  library(doParallel)
  library(dplyr)
  library(DRIMSeq)
  library(edgeR)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GGally)
  library(ggbeeswarm)
  library(ggfortify)
  library(ggplot2)
  library(ggpubr)
  library(ggraph)
  library(ggrepel)
  library(ggtree)
  library(grid)
  library(gridExtra)
  library(heatmap3)
  library(knitr)
  library(lmtest)
  library(MAST)
  library(meshes)
  library(MKmisc)
  library(msigdbr)
  library(NMF)
  library(nVennR)
  library(pathview)
  library(pbapply)
  library(pcaMethods)
  library(qvalue)
  library(QoRTs)
  library(RColorBrewer)
  library(ReactomePA)
  library(reshape2)
  library(reticulate)
  library(rsvd)
  library(scales)
  library(stageR)
  library(stringr)
  library(tidyr)
  library(trackViewer)
  library(tximport)
  library(umap)
  library(xlsx)
})
reticulate::use_condaenv("scvelo")
scv <- reticulate::import("scvelo")

FILES <- paste(getwd(), "/dataset/", sep = '')
# Define Global colors and conditions
DOT_COLOR <- c(
  "#293acc",
  "#B20000",
  "#008331",
  "#000000"
)
names(DOT_COLOR) <- c(
  "mouseEgg",
  "mouseZygote",
  "ratEgg",
  "ratZygote"
)
CONDITIONS_NAMES <- c(
  "Mouse Egg",
  "Mouse Zygote",
  "Rat Egg",
  "Rat Zygote"
)
names(CONDITIONS_NAMES) <- c(
  "mouseEgg",
  "mouseZygote",
  "ratEgg",
  "ratZygote"
)

# Read file functions


get_qorts_summary <- function(directory, file_name = "QC.summary.txt") {
  rnames <- row.names(read.csv(paste(list.dirs(directory, recursive = F)[[1]], file_name, sep = "/"), row.names = 1, sep = "\t", header = F))
  summary0 <- do.call(cbind, lapply(list.dirs(directory, recursive = F), FUN = function(x) {
    if (length(list.files(x)) > 2) {
      new <- read.csv(paste(x, file_name, sep = "/"), row.names = 1, sep = "\t", header = F)
      rnames <- intersect(row.names(new), rnames)
      new <- new[rnames, 1]
      return(new)
    } else {
      print(x)
      return()
    }
  }))
  row.names(summary0) <- rnames

  colnames(summary0) <- do.call(c, lapply(list.dirs(directory, recursive = F), FUN = function(x) {
    strsplit(x, split = "/")[[1]][length(strsplit(x, split = "/")[[1]])]
  }))

  return(summary0)
}


write_to_xlsx <- function(tables, file_name, table_names = NULL, overwrite = T){
  wb <- createWorkbook()
  
  if(is.data.table(tables)){
    tables <- as.data.frame(tables)
  }
  if(is.data.frame(tables)){
    if(is.null(table_names)){
      table_names <- 'sheet1'
    }
    addWorksheet(wb, table_names)
    
    # Write data to the sheet
    writeData(wb, table_names, tables)
    

  }
  else if(is.list(tables)){
    if(is.null(table_names) & is.null(names(tables))){
      names(tables) <- paste('sheet', 1:length(tables), sep = '')
      table_names <- names(tables)
    }
    else if(is.null(table_names)){
      table_names <- names(tables)
    }
    else if(is.null(names(tables))){
      names(tables) <- table_names
    }
    
    for(i in seq_len(length(tables))){
      
        addWorksheet(wb, table_names[i])
        
        # Write data to the sheet
        writeData(wb, table_names[i], tables[[i]])
      
    }
  }
  saveWorkbook(wb, file_name, overwrite = overwrite)
}


plot_qorts <- function(res=NULL,
                       dirc=NULL, 
                       qorts=NULL, 
                       plot_types = c('biotype.rates',
                                      'chrom.type.rates', 
                                      'clipping', 
                                      'dropped.rates', 
                                      'gene.assignment.rates', 
                                      'genebody.coverage',
                                      'genebody.coverage.UMQuartile',
                                      'mapping.rates',
                                      'insert.size'),
                       sep = 'condition',
                       plot_path = './qort_plot.pdf', colorby_colors = NULL
                       ){
  if(!is.null(res)){
    res = res
    qorts = res@decoder
  }else{
    res <- read.qc.results.data(infile.dir = dirc, decoder=qorts, autodetectMissingSamples = TRUE, debugMode = F, calc.DESeq2 = T, calc.edgeR = T)
  }
  #print(qorts$condition)
  colorby <- qorts[res@decoder$unique.ID,sep]
  names(colorby) <- res@decoder$unique.ID
  #print(colorby)
  if(!is.null(coorby_colors)){
  all_plot <- build.plotter.advanced(res, colorBy = as.character(colorby), color.title = 'cellType', plotter.params = list(contrasting.colors = colorby_colors[sort(as.character(unique(colorby)))]))
  }else{
    all_plot <- build.plotter.advanced(res, colorBy = as.character(colorby), color.title = 'cellType')
    
  }
  group_list = list()
  for(c in unique(colorby)){
    if(!is.null(colorby_colors)){
      group_list[[c]] = build.plotter.advanced(res, plotter.params = list(std.color = colorby_colors[c]), highlightBy = colorby, highlight = c , highlightTitle.singular = 'Condition', outgroup.title = 'Others')
    }else{
    group_list[[c]] = build.plotter.advanced(res, highlightBy = colorby, highlight = c , highlightTitle.singular = 'Condition', outgroup.title = 'Others')
    }
    }
  pdf(plot_path)
  plot.new()
  makePlot.legend.box(all_plot)
  for(p in plot_types){
    if(p == 'biotype.rates'){
      params = '(plotter, count.type = "unambigOnly", showTypes = c("protein_coding", "ncRNA", "rRNA", "pseudogene", "UNK"))'
    }else{
      params = '(plotter)'
    }
    par(mfrow=c(1,1))
    plotter = all_plot
    eval(parse(text=paste("makePlot.", p,params, sep="")))
    
    par(mfrow=c(2,2))
    for(c in names(group_list)){
      plotter = group_list[[c]]
      eval(parse(text=paste("makePlot.", p, params, sep="")))
    }
  }
  
  #makePlot.biotype.rates(all_plot)
  #makePlot.chrom.type.rates(all_plot)
  #makePlot.clipping(all_plot)
  #makePlot.dropped.rates(all_plot)
  #makePlot.gene.assignment.rates(all_plot)
  #makePlot.genebody.coverage(all_plot)
  #makePlot.genebody.coverage.UMQuartile(all_plot)
  #makePlot.mapping.rates(all_plot)
  dev.off()
}



read_expression <- function(dir, mode = "salmon", tx2gene = NULL, tpmType = "no", dropInfReps = T) {
  if (mode %in% c("star", "hisat")) {
    files <- list.files(dir, pattern = paste(mode, "htseq.ct", sep = "."), full.names = T)
    sample_names <- as.character(data.frame(strsplit(list.files(dir, pattern = paste(mode, "htseq.ct", sep = ".")), split = "_"))[1, ])
    cts <- do.call(cbind, lapply(files, FUN = function(x) {
      read.table(x, row.names = 1, header = F)
    }))
    colnames(cts) <- sample_names
    # cts <- cts[1:(nrow(cts)-5),]
    cts <- cts[rowSums(cts != 0) > 0, ]
    return(cts)
  } else if (mode == "salmon") {
    tx2gene <- read.csv(tx2gene, sep = "\t", header = T, col.names = c("TXNAME", "GENEID"))
    salmon_files <- list.files(dir, recursive = T, pattern = "*quant.sf", full.names = T)
    names(salmon_files) <- as.character(as.data.frame(strsplit(salmon_files, "/"))[length(strsplit(salmon_files, "/")[[1]]) - 1, ])
    gene.salmon <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = tpmType, dropInfReps = dropInfReps, importer = read.delim)
    gene.salmon$counts <- gene.salmon$counts[rowSums(gene.salmon$counts != 0) > 0, ]
    gene.salmon$abundance <- gene.salmon$abundance[rowSums(gene.salmon$abundance != 0) > 0, ]
    transcript.salmon <- tximport(salmon_files, type = "salmon", txOut = TRUE, countsFromAbundance = tpmType, dropInfReps = dropInfReps, importer = read.delim)
    transcript.salmon$counts <- transcript.salmon$counts[rowSums(transcript.salmon$counts != 0) > 0, ]
    transcript.salmon$abundance <- transcript.salmon$abundance[rowSums(transcript.salmon$abundance != 0) > 0, ]
    return(list(gene = gene.salmon, transcript = transcript.salmon))
  } else if (mode == "rsem") {
    file_names <- list.files(dir, recursive = T, pattern = "*genes.results", full.names = F)
    gene_files <- list.files(dir, recursive = T, pattern = "*genes.results", full.names = T)
    names(gene_files) <- as.character(as.data.frame(strsplit(file_names, "\\."))[1, ])
    gene.rsem <- tximport(gene_files, type = "rsem", txIn = FALSE, txOut = FALSE, importer = read.delim)
    gene.rsem$counts <- gene.rsem$counts[rowSums(gene.rsem$counts != 0) > 0, ]
    gene.rsem$abundance <- gene.rsem$abundance[rowSums(gene.rsem$abundance != 0) > 0, ]
    file_names <- list.files(dir, recursive = T, pattern = "*isoforms.results", full.names = F)
    transcript_files <- list.files(dir, recursive = T, pattern = "*isoforms.results", full.names = T)
    names(transcript_files) <- as.character(as.data.frame(strsplit(file_names, "\\."))[1, ])
    transcript.rsem <- tximport(transcript_files, type = "rsem", txIn = TRUE, txOut = TRUE, importer = read.delim)
    transcript.rsem$counts <- transcript.rsem$counts[rowSums(transcript.rsem$counts != 0) > 0, ]
    transcript.rsem$abundance <- transcript.rsem$abundance[rowSums(transcript.rsem$abundance != 0) > 0, ]
    return(list(gene = gene.rsem, transcript = transcript.rsem))
  } else {
    cat("only support salmon, rsem, hisat (htseq) and star (htseq)\n")
  }
}




# read from velocyto loom files using scvelo and reticulate
read_filter_splice_loom <- function(loomd, samples = NULL, filt_samp_by_cor = 0, f_format = c('kb_loom', 'velocyto_loom', 'csv'), splice_files = NULL) {
  f_format <- match.arg(f_format)
  if(f_format == 'velocyto_loom'){
    loom <- scv$read_loom(loomd)
    spliced <- t(as.matrix(loom$layers["spliced"]))
    colnames(spliced) <- sapply(strsplit(sapply(strsplit(loom$obs_names$values, ":"), function(x) x[2]), "_"), function(x) x[1])
    row.names(spliced) <- loom$var_names$values
    unsplic <- t(as.matrix(loom$layers["unspliced"]) + as.matrix(loom$layers["spanning"]))
    colnames(unsplic) <- colnames(spliced)
    row.names(unsplic) <- loom$var_names$values
  }else if(f_format == 'kb_loom'){
    loom <- scv$read_loom(paste(loomd, 'adata.loom', sep = '/'))
    spliced <- t(as.matrix(loom$layers["spliced"]) + as.matrix(loom$layers["ambiguous"]))
    colnames(spliced) <- row.names(read.csv(paste(loomd, 'samples.txt', sep = '/'), sep ='\t', row.names = 1, header = F))
    row.names(spliced) <- loom$var$target_name
    unsplic <- t(as.matrix(loom$layers["unspliced"]))
    colnames(unsplic) <- colnames(spliced)
    row.names(unsplic) <- loom$var$target_name
  }else{
    spliced <- do.call(cbind, lapply(splice_files['splice'], function(x) {read.csv(x, header = T, row.names = 1)}))
    unsplic <- do.call(cbind, lapply(splice_files['unspliced'], function(x) {read.csv(x, header = T, row.names = 1)}))
  }
  spliced_raw <- spliced
  unsplic_raw <- unsplic
  if (!is.null(samples)) {
    unsplic <- unsplic[, samples]
    spliced <- spliced[, samples]
  }
  

  grp1 <- which(grepl(pattern = "F", colnames(spliced)))
  grp2 <- which(grepl(pattern = "U", colnames(spliced)))

  grp1 <- grp1[colMeans(cor(unsplic[, grp1], unsplic[, grp1])) > filt_samp_by_cor]
  grp2 <- grp2[colMeans(cor(unsplic[, grp2], unsplic[, grp2]))> filt_samp_by_cor]


  grp1_filt_s <- rowMeans(spliced[, grp1]) > 5 & rowSums(spliced[, grp1] > 0) > length(grp1) * 0.2
  grp1_filt_u <- rowMeans(unsplic[, grp1]) > 1 & rowSums(unsplic[, grp1] > 0) > length(grp1) * 0.2

  grp2_filt_s <- rowMeans(spliced[, grp2]) > 5 & rowSums(spliced[, grp2] > 0) > length(grp2) * 0.2
  grp2_filt_u <- rowMeans(unsplic[, grp2]) > 1 & rowSums(unsplic[, grp2] > 0) > length(grp2) * 0.2


  spliced <- spliced[(grp1_filt_s & grp1_filt_u) | (grp2_filt_s & grp2_filt_u), ]
  unsplic <- unsplic[(grp1_filt_s & grp1_filt_u) | (grp2_filt_s & grp2_filt_u), ]

  unsplic <- unsplic[, c(grp1, grp2)]
  spliced <- spliced[, c(grp1, grp2)]


  test_cor <- do.call(rbind, lapply(1:nrow(spliced), FUN = function(x) {
    c(
      "pc" = cor(spliced[x, ], unsplic[x, ]),
      "sc" = cor(spliced[x, ], unsplic[x, ], method = "spearman"),
      "kc" = cor(spliced[x, ], unsplic[x, ], method = "kendall"),
      "rsq" = summary(lm(spliced[x, ] ~ unsplic[x, ]))$adj.r.squared,
      "err.rat" = sd(unsplic[x, ]) / sd(spliced[x, ]),
      "expr.rat" = mean(unsplic[x, ]) / mean(spliced[x, ])
    )
  }))
  row.names(test_cor) <- row.names(spliced)
  test_cor <- data.frame(test_cor)
  genes <- row.names(subset(test_cor, rsq > 0.1 & kc > 0.1 & err.rat > 0.005 & err.rat < 5))
  return(list(spliced_raw = spliced_raw, unspliced_raw = unsplic_raw, spliced = spliced, unspliced=unsplic, metrics = test_cor, genes = genes))
}





read_rnasplice_dex_dtu <- function(dexseq_ds_rds, drim_seq_filt_rds) {
  dexseq_ds <- readRDS(dexseq_ds_rds)
  drimseq_ds <- readRDS(drim_seq_filt_rds)
  dex_norm <- cbind(as.data.frame(stringr::str_split_fixed(rownames(counts(dexseq_ds)), ":", 2)), as.data.frame(counts(dexseq_ds, normalized = TRUE))[, 1:(nrow(dexseq_ds@colData) / 2)])
  colnames(dex_norm) <- c("groupID", "featureID", as.character(colData(dexseq_ds)$sample.1)[1:(nrow(dexseq_ds@colData) / 2)])
  row.names(dex_norm) <- NULL
  obj <- list()
  obj$dexseq <- dexseq_ds
  obj$drimseq <- drimseq_ds
  obj$dex_norm <- dex_norm
  return(obj)
}



prepareCount <- function(cds, dirc) {
  cells <- row.names(cds$meta)
  # tpm expression file
  # genes_table <- read.table(paste(FILES,dirc, "genes.csv", sep=""), row.names = 1, header=T, sep = ',')

  # if(ribo_filter){
  # genes_table <- genes_table[!grepl('ribosom', genes_table$Gene.description) & !grepl('rRNA', genes_table$Gene.description),]
  # }
  # genes <- row.names(subset(genes_table, Gene.featureType == 'ORF'))
  if (dir.exists(paste(FILES, dirc, "salmon", sep = "/"))) {
    salmon <- read_expression(paste(FILES, dirc, "salmon", sep = "/"), mode = "salmon", tx2gene = paste(FILES, dirc, "tx2gene.tsv", sep = "/"))
    cds[["salmon"]] <- salmon
    tpm <- salmon$gene$abundance[, cells]
    cds[["tpm"]] <- list()
    cds[["tpm"]][["bio"]] <- tpm[which(!grepl("ERCC", row.names(tpm))), ]
    cds[["tpm"]][["bio"]] <- t(t(cds[["tpm"]][["bio"]]) * 1e6 / colSums(cds[["tpm"]][["bio"]]))
    cds[["tpm"]][["spike"]] <- tpm[which(grepl("ERCC", row.names(tpm))), ]
  }
  if (dir.exists(paste(FILES, dirc, "hisat", sep = "/"))) {
    hisat <- read_expression(paste(FILES, dirc, "hisat", sep = "/"), mode = "hisat", tx2gene = paste(FILES, dirc, "tx2gene.tsv", sep = "/"))
    if (is.null(cds[["ct"]])) {
      cds[["ct"]] <- list()
    }
    cds[["ct"]][["bio"]] <- hisat[which(!grepl("ERCC", row.names(hisat))), cells]
    cds[["ct"]][["spike"]] <- hisat[which(grepl("ERCC", row.names(hisat))), cells]
  }
  if (dir.exists(paste(FILES, dirc, "star", sep = "/"))) {
    star <- read_expression(paste(FILES, dirc, "star", sep = "/"), mode = "star", tx2gene = paste(FILES, dirc, "tx2gene.tsv", sep = "/"))
    if (is.null(cds[["ct"]])) {
      cds[["ct"]] <- list()
    }
    cds[["ct"]][["bio_star"]] <- star[which(!grepl("ERCC", row.names(star))), cells]
    cds[["ct"]][["spike_star"]] <- star[which(grepl("ERCC", row.names(star))), cells]
  }
  if ("tpm" %in% names(cds) & !("ct" %in% names(cds))) {
    ct <- round(cds[["salmon"]]$gene$counts[, cells])
    cds[["ct"]] <- list()
    cds[["ct"]][["bio"]] <- round(ct[which(!grepl("ERCC", row.names(ct))), ])

    cds[["ct"]][["spike"]] <- round(ct[which(grepl("ERCC", row.names(ct))), ])
  }
  return(cds)
}


prepareGeneFeatures <- function(cds, gtf_file = NULL, dirc = "mouse_data", org_type = "mouse") {
  if (org_type == "mouse") {
    org_type <- "mmusculus_gene_ensembl"
    kegg_org <- "mmu"
  } else if (org_type == "rat") {
    org_type <- "rnorvegicus_gene_ensembl"
    kegg_org <- "rno"
  }
  gtf <- rtracklayer::import(gtf_file)
  gtf <- data.frame(gtf[gtf$type == "gene", ])
  gtf <- gtf[, colSums(is.na(gtf)) != nrow(gtf)]
  gtf <- data.frame(row.names = gtf$gene_id, gtf)
  gtf$gene_short_name <- gtf$gene
  genes <- c()
  if ("tpm" %in% names(cds)) {
    genes <- union(rownames(cds$tpm$bio), genes)
    cds$features_tpm <- gtf[row.names(cds$tpm$bio), ]
  }
  if ("ct" %in% names(cds)) {
    genes <- union(genes, rownames(cds$ct$bio))
    cds$features_ct <- gtf[row.names(cds$ct$bio), ]
  }
  gtf <- gtf[genes, ]
  cds$features <- gtf
  return(cds)
}

read_htseq_intergenic <- function(fn = "./dataset/mouse_data/mouse_intergenic_1000.hisat.ct.txt", ct, meta, cov_thresh = 5) {
  samples <- colnames(ct)

  intergenic <- read.csv(fn, row.names = 1, header = T, sep = "\t")
  intergenic1 <- intergenic[-c(1:5), ]
  intergenic1_u <- intergenic1[, grepl(".u.", colnames(intergenic1))]
  intergenic1_n <- intergenic1[, grepl("int_n", colnames(intergenic1))]
  intergenic1_s <- intergenic1[, grepl("int_s", colnames(intergenic1))]
  colnames(intergenic1_u) <- sapply(strsplit(colnames(intergenic1_n), "_"), function(x) {
    x[1]
  })
  colnames(intergenic1_n) <- sapply(strsplit(colnames(intergenic1_n), "_"), function(x) {
    x[1]
  })
  colnames(intergenic1_s) <- sapply(strsplit(colnames(intergenic1_n), "_"), function(x) {
    x[1]
  })
  intergenic1_u <- intergenic1_u[, samples]
  intergenic1_s <- intergenic1_s[, samples]
  intergenic1_n <- intergenic1_n[, samples]
  u_1k <- grepl("_u_0", row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), "_"), length) == 3

  d_1k <- grepl("_d_0", row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), "_"), length) == 3

  s_1k <- grepl("_0", row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), "_"), length) > 3

  u_2to10k <- grepl("_u", row.names(intergenic1)) & !grepl("_0", row.names(intergenic1)) & !grepl("far", row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), "_"), length) == 3
  d_2to10k <- grepl("_d", row.names(intergenic1)) & !grepl("_0", row.names(intergenic1)) & !grepl("_far", row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), "_"), length) == 3

  s_2to10k <- !grepl("_0", row.names(intergenic1)) & !grepl("far", row.names(intergenic1)) & sapply(strsplit(row.names(intergenic1), "_"), length) > 3

  far <- grepl("_far", row.names(intergenic1))

  num_genes <- colSums(ct > 0)
  print(num_genes)
  num_reads <- colSums(intergenic1_u)
  print(num_reads)

  norm_df <- data.frame(
    row.names = colnames(intergenic1_u),
    up1k = colSums(intergenic1_u[u_1k, ] > cov_thresh) / num_genes,
    down1k = colSums(intergenic1_u[d_1k, ] > cov_thresh) / num_genes,
    u_2to10k = colSums(intergenic1_u[u_2to10k, ] > cov_thresh) * 100 / num_reads,
    d_2to10 = colSums(intergenic1_u[d_2to10k, ] > cov_thresh) * 100 / num_reads,
    s_1k = colSums(intergenic1_u[s_1k, ] > cov_thresh) / num_genes,
    far = colSums(intergenic1_u[far, ] > cov_thresh) * 1000 / num_reads,
    extend_u = colSums((intergenic1_n[u_1k, ] - intergenic1_s[u_1k, ]) > cov_thresh) / num_genes,
    extend_d = colSums((intergenic1_n[d_1k, ] - intergenic1_s[d_1k, ]) > cov_thresh) / num_genes,
    cellType = meta$cellType
  )

  df <- data.frame(
    row.names = colnames(intergenic1_u),
    up1k = colSums(intergenic1_u[u_1k, ] > cov_thresh),
    down1k = colSums(intergenic1_u[d_1k, ] > cov_thresh),
    u_2to10k = colSums(intergenic1_u[u_2to10k, ] > cov_thresh),
    d_2to10 = colSums(intergenic1_u[d_2to10k, ] > cov_thresh),
    s_1k = colSums(intergenic1_u[s_1k, ] > cov_thresh),
    far = colSums(intergenic1_u[far, ] > cov_thresh),
    extend_u = colSums((intergenic1_n[u_1k, ] - intergenic1_s[u_1k, ]) > cov_thresh),
    extend_d = colSums((intergenic1_n[d_1k, ] - intergenic1_s[d_1k, ]) > cov_thresh),
    all = colSums(intergenic1_u > cov_thresh),
    cellType = meta$cellType
  )


  small_df <- data.frame(
    row.names = colnames(intergenic1_u),
    far = colSums(intergenic1_u[far, ] > cov_thresh),
    all = colSums(intergenic1_u > cov_thresh),
    cellType = meta$cellType
  )
  return(list(df = df, norm_df = norm_df, small_df = small_df, raw = intergenic, u = intergenic1_u, s = intergenic1_s, n = intergenic1_n, u_1k = u_1k, d_1k = d_1k, s_1k = s_1k, u_2to10k = u_2to10k, d_2to10k = d_2to10k, s_2to10k = s_2to10k, far = far))
}

# Exploratory plotting functions (Mainly just UMAPs)
sample_PCA <- function(cds,
                       meta,
                       color_by = "cellType",
                       umap = F,
                       labeling = FALSE,
                       point_size = 4,
                       dimension = 2,
                       reduce_noise = FALSE,
                       return_matrix = FALSE,
                       highlight = NULL,
                       main = "plot",
                       umap.config = umap_config,
                       legend.label = c("oocyte", "zygote"),
                       legend.position = c(1.05, 1.05)) {
  if (!umap) {
    test <- rsvd::rpca(t(cds), center = T, scale = T)
    percentVar <-
      c(
        100 * round(test$sdev[1]^2 / sum(test$sdev^2), 3),
        100 * round(test$sdev[2]^2 / sum(test$sdev^2), 3),
        100 * round(test$sdev[3]^2 / sum(test$sdev^2), 3)
      )
    axis_x <- paste0("PC1: ", percentVar[1], "% variance")
    axis_y <- paste0("PC2: ", percentVar[2], "% variance")
    axis_z <- paste0("PC3: ", percentVar[3], "% variance")
    test <- test$x
  } else {
    umap_test <-
      umap(
        t(cds),
        n_components = dimension,
        n_neighbors = umap.config[["n_neighbors"]],
        min_dist = umap.config[["min_dist"]],
        metric = umap.config[["metric"]],
        random_state = umap.config[["seed"]]
      )
    axis_x <- paste0("UMAP 1")
    axis_y <- paste0("UMAP 2")
    axis_z <- paste0("UMAP 3")
    test <- umap_test$layout
  }
  alpha <- 1
  if (dimension == 2) {
    theme0 <-
      theme_bw() + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.position = legend.position,
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.key.size = unit(0.1, "lines")
      )
    scale_color_manual()
    if (color_by == "cellType") {
      dot_color <- DOT_COLOR
      colors_used <- dot_color[as.character(unique(meta[[color_by]]))]
      meta$condition <- CONDITIONS_NAMES[as.character(meta[[color_by]])]
      alpha <- 1
    }
    test <- data.frame(test[, c(1, 2)], condition = as.factor(meta[["cellType"]]), org = as.factor(meta[["organism"]]))
    colnames(test) <- c("PC1", "PC2", "condition", "org")
    if (labeling == TRUE) {
      plot0 <-
        ggplot(test,
          label = T,
          aes(PC1, PC2, color = condition, label = row.names(test), shape = NULL)
        ) +
        geom_text(size = 2) +
        geom_point(size = 0) +
        xlab(axis_x) +
        ylab(axis_y) +
        coord_fixed() +
        theme0 +
        ggtitle(main)
    } else {
      plot0 <-
        ggplot(test, label = T, aes(PC1, PC2, color = condition, shape = NULL)) +
        geom_point(size = point_size, alpha = alpha) +
        xlab(axis_x) +
        ylab(axis_y) +
        theme0 +
        ggtitle(main)
    }
    if (color_by == "cellType") {
      plot0 <- plot0 + scale_color_manual(values = colors_used, labels = legend.label)
    } else {
      plot0 <- plot0
    }
  }else{ # 3D plots in plotly
    test <-
      data.frame(test[, c(1, 2, 3)], condition = as.factor(meta[[color_by]]))
    colnames(test) <- c("PC1", "PC2", "PC3", "condition")
    if (labeling == TRUE) {
      plot0 <- plot_ly(
        test,
        x = ~PC1,
        y = ~PC2,
        z = ~PC3,
        color = ~condition,
        text = rownames(test)
      ) %>%
        add_text() %>%
        layout(scene = list(
          xaxis = list(title = axis_x),
          yaxis = list(title = axis_y),
          zaxis = list(title = axis_z)
        ))
    } else {
      plot0 <- plot_ly(
        test,
        x = ~PC1,
        y = ~PC2,
        z = ~PC3,
        color = ~condition
      ) % >% layout(scene = list(
        xaxis = list(title = axis_x),
        yaxis = list(title = axis_y),
        zaxis = list(title = axis_z)
      ))
    }
  }
  umap_axis <- ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(2, "cm")
  )
  plot0 <- plot0 + guides(x = umap_axis, y = umap_axis) +
    theme(axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, "pt"))), axis.title = element_text(hjust = 0)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)
  return(plot0)
}



options(ucscChromosomeNames = FALSE)




plot_utr_coverage <- function(gene, utr_res, txdb, bw_file_list, cols, cell_names, y_margin = -0.05, y_height = 2.5, file_name_suffix = "", loci = NULL) {

  if (is.null(loci)) {
    loci <- utr_res[utr_res$gene_short_names == gene, "loci"]
  }
  # cols <- cols[cell_names]
  chrom_range <- strsplit(loci, split = ":")[[1]]
  chrom <- chrom_range[1]
  rg <- strsplit(chrom_range[2], "-")[[1]]
  strt <- as.numeric(rg[1])
  nd <- as.numeric(rg[2])
  pas <- utr_res[utr_res$gene_short_names == gene, "Predicted_Proximal_APA"]
  strd <- utr_res[utr_res$gene_short_names == gene, "strand"]

  if (strd == "+") {
    start <- strt - 100
    end <- nd + 10
  } else {
    start <- strt - 10
    end <- nd + 100
  }

  gr <- GRanges(chrom, IRanges(start, end), strand = strd)
  # grW <- parse2GRanges(utr_res[gene,'loci'])
  ids <- getGeneIDsFromTxDb(gr, txdb)
  print(ids)
  symbols <- ids
  genes <- geneTrack(ids, txdb,
    symbols,
    asList = FALSE
  )
  # temp_score <- importScore(bw_file_list[1],
  # bw_file_list[2],
  # format="BigWig", ranges = gr)
  temp_score <- sapply(bw_file_list, FUN = function(x) {
    importScore(x, format = "BigWig", ranges = gr)
  })
  print('a')
  # setTrackStyleParam(temp_score, "color", c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']))
  # strand(trackList[['Mouse']]@dat) <- '+'
  # strand(trackList[['Mouse']]@dat2) <- '-'

  # temp <- geneModelFromTxdb(mouse_gtf, org.Mm.eg.db, gr = gr)
  trackList <- trackList(c(genes, temp_score))
  # names(trackList) <- c(gene, 'Mouse')

  optSty <- optimizeStyle(trackList, theme = "bw")
  # viewerStyle <- trackViewerStyle()
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  names(trackList) <- c(paste(gene, " (", strd, ")", sep = ""), cell_names)
  for (i in 1:length(bw_file_list)) {
    setTrackStyleParam(trackList[[i + 1]], "color", cols[i])
    setTrackStyleParam(trackList[[i + 1]], "ylabpos", "bottomright")
    setTrackStyleParam(trackList[[i + 1]], "ylabgp", list(cex = 1.5, col = "black"))
    trackList[[i + 1]]@style@yaxis@main <- T
    trackList[[i + 1]]@style@yaxis@gp$cex <- 0.8
    trackList[[i + 1]]@style@xscale@draw <- F
  }
  setTrackViewerStyleParam(viewerStyle, "xaxis", T)
  # setTrackViewerStyleParam(viewerStyle, "xgp", list(cex = 1.3, col = 'black'))
  setTrackViewerStyleParam(viewerStyle, "margin", c(0.14, 0.25, y_margin, 0.05))
  setTrackXscaleParam(trackList[[1]], "draw", F)
  setTrackStyleParam(trackList[[1]], "ylabpos", "bottomright")
  # setTrackStyleParam(trackList[[2]], "ylabpos", "bottomright")
  # setTrackStyleParam(trackList[[3]], "ylabpos", "bottomright")
  # setTrackStyleParam(trackList[[4]], "ylabpos", "bottomright")
  # setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=1.5, col="black"))
  # setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=1.5, col="black"))
  # setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=1.5, col="black"))
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex = 0, col = "black"))
  setTrackStyleParam(trackList[[1]], "height", 0.1)
  trackList[[1]]@style@yaxis@main <- T
  # trackList[[2]]@style@yaxis@main <- T
  # trackList[[3]]@style@yaxis@main <- T
  # trackList[[4]]@style@yaxis@main <- T
  # trackList[[2]]@style@yaxis@gp$cex <- 0.8
  # trackList[[3]]@style@yaxis@gp$cex <- 0.8
  # trackList[[4]]@style@yaxis@gp$cex <- 0.8
  # trackList
  png(paste(paste(gene, file_name_suffix, sep = "_"), ".png", sep = ""), width = 4.2, height = y_height, units = "in", res = 300)
  vp <- viewTracks(trackList, gr = gr, viewerStyle = viewerStyle)
  addGuideLine(c(strt, pas, nd), vp = vp, col = c("black", "red", "black"), lwd = 2.3)
  grid.text(paste(gene, " (", strd, ")", sep = ""), 0.05, 0.55, rot = 90, gp = gpar(fontsize = 28, fontface = "bold"))
  dev.off()
  vp
}


plot_range_coverage <- function(txdb, bw_file_list, cols, cell_names, gene_name = NULL, range = NULL, strand = NULL, log = T, file_name_suffix = "", add_intron_line = T, y_lim = c(0, 1)) {
  # chrom_range <- strsplit(range, split = ':')[[1]]
  # chrom <- chrom_range[1]
  # rg <- strsplit(chrom_range[2], '-')[[1]]
  # strt <- as.numeric(rg[1])
  # nd <- as.numeric(rg[2])

  # if(strd == '+'){
  # start <- strt-200
  # end <- nd+10
  # }else{
  # start <- strt-10
  #  end <- nd+200
  # }

  # gr <- GRanges(chrom, IRanges(start, end), strand=strd)

  if (!is.null(range)) {
    grW <- parse2GRanges(range)
    ids <- getGeneIDsFromTxDb(grW, txdb)
    gene <- ids[1]
  } else {
    gene <- gene_name
  }
  if (!is.null(gene_name)) {
    gene <- gene_name
  }
  symbols <- gene
  # gene <- ifelse(is.null(gene_name), ids[1], gene_name)

  genes <- geneTrack(gene, txdb,
    gene,
    asList = FALSE
  )
  strand <- genes@dat@strand@values[1]
  strd <- strand
  if (!is.null(range)) {
    min_st <- max(min(genes@dat@ranges@start), min(grW@ranges@start)) + 1
    max_en <- min(max(genes@dat@ranges@start + genes@dat@ranges@width), max(grW@ranges@start + grW@ranges@width)) - 1
  } else {
    min_st <- min(genes@dat@ranges@start) + 1
    max_en <- max(genes@dat@ranges@start + genes@dat@ranges@width) - 1
  }
  gr <- GRanges(seqnames = c(genes@dat@seqnames@values), IRanges(start = c(min_st), end = c(max_en)), strand = c("*"), mcols = c("a"))

  gr_intron <- GenomicRanges::setdiff(gr, genes@dat, ignore.strand = T)
  # return(list(gr=gr, genes = genes))
  intron_pos <- sort(c(gr_intron@ranges@start + 5, gr_intron@ranges@start + gr_intron@ranges@width - 5))
  # temp_score <- importScore(bw_file_list[1],
  # bw_file_list[2],
  # format="BigWig", ranges = gr)
  temp_score <- sapply(bw_file_list, FUN = function(x) {
    importScore(x, format = "BigWig", ranges = gr)
  })
  # setTrackStyleParam(temp_score, "color", c(DOT_COLOR['mouseZygote'], DOT_COLOR['mouseEgg']))
  # strand(trackList[['Mouse']]@dat) <- '+'
  # strand(trackList[['Mouse']]@dat2) <- '-'

  # temp <- geneModelFromTxdb(mouse_gtf, org.Mm.eg.db, gr = gr)
  trackList <- trackList(c(genes, temp_score))
  if (log) {
    trackList[[2]]@dat$score <- log10(trackList[[2]]@dat$score + 1)
    trackList[[3]]@dat$score <- log10(trackList[[3]]@dat$score + 1)
  }
  # names(trackList) <- c(gene, 'Mouse')

  optSty <- optimizeStyle(trackList, theme = "bw")
  # viewerStyle <- trackViewerStyle()
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  names(trackList) <- c(paste(gene, " (", strd, ")", sep = ""), cell_names[1], cell_names[2])
  setTrackStyleParam(trackList[[2]], "color", cols[1])
  setTrackStyleParam(trackList[[3]], "color", cols[2])
  setTrackViewerStyleParam(viewerStyle, "xaxis", T)
  setTrackViewerStyleParam(viewerStyle, "xgp", list(cex = 0.8, col = "black"))
  setTrackViewerStyleParam(viewerStyle, "margin", c(0.14, 0.15, -0.25, 0.05))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackStyleParam(trackList[[1]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[2]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[3]], "ylabpos", "bottomright")
  setTrackStyleParam(trackList[[3]], "ylabgp", list(cex = 1.5, col = "black"))
  setTrackStyleParam(trackList[[2]], "ylabgp", list(cex = 1.5, col = "black"))
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex = 0, col = "black"))
  setTrackStyleParam(trackList[[1]], "height", 0.1)
  trackList[[1]]@style@yaxis@main <- T
  trackList[[2]]@style@yaxis@main <- T
  trackList[[3]]@style@yaxis@main <- T
  trackList[[2]]@style@yaxis@gp$cex <- 0.8
  trackList[[3]]@style@yaxis@gp$cex <- 0.8
  # trackList[[1]]@style@xaxis@gp$cex <- 0.8
  # return(trackList)
  png(paste(paste(gene, file_name_suffix, sep = "_"), ".intron.png", sep = ""), width = 4.3, height = 2.5, units = "in", res = 300)
  vp <- viewTracks(trackList, gr = gr, viewerStyle = viewerStyle)
  grid.text(paste(gene, " (", strd, ")", sep = ""), 0.03, 0.65, rot = 90, gp = gpar(fontsize = 24, fontface = "bold"))
  if (add_intron_line) {
    addGuideLine2(intron_pos, vp = vp, col = rep("red", length(intron_pos)), lty = rep(c("dashed", "dotted"), length(intron_pos) / 2), lwd = 2.5, y_lim = y_lim)
  }
  dev.off()
  return(list(vp = vp, gene = genes, gr_intron = gr_intron))
}







mast_diff <- function(obj = NULL, plot = F, ct = NULL, meta = NULL, FCThresh = log2(1.25), normFactor = NULL, control = "mouseEgg", tpm = T,
                      freq = 0.5, max_thres = 3, bin_by = "median", nbins = 20, min_per_bin = 30,
                      correct_wild_coef = T, corr_det = T, min_cell_grp = 3, min_cell = 6, include_filt_as_NA = F) {
  # ct <- if(tpm){obj$tpm$bio}else{t(t(obj$tpm$bio)/edgeR::calcNormFactors(obj$tpm$bio))}
  res <- list()
  if (!is.null(obj)) {
    ct <- if (tpm) {
      obj$tpm$bio
    } else {
      obj$ct$bio
    }
    meta <- obj$meta[colnames(ct), ]
  }
  if (!tpm) {
    if (!is.null(normFactor)) {
      ct <- t(t(ct) / normFactor)
    } else {
      ct <- edgeR::cpm(ct)
    }
  }
  ct <- ct[rowSums(ct > 0) > ncol(ct) * 0, ]

  target <- setdiff(unique(meta$cellType), c(control))
  ctrl_mat <- ct[, row.names(subset(meta, cellType == control))]
  ctrl_mat <- ctrl_mat[rowSums(ctrl_mat > 0) > max(ncol(ctrl_mat) * freq, min_cell_grp), ]
  trgt_mat <- ct[, row.names(subset(meta, cellType == target))]
  trgt_mat <- trgt_mat[rowSums(trgt_mat > 0) > max(ncol(trgt_mat) * freq, min_cell_grp), ]
  genes <- union(row.names(ctrl_mat), row.names(trgt_mat))
  #print(length(genes))

  genes <- intersect(row.names(ct)[which(rowSums(ct > 0) > min_cell)], genes)
  filt_genes <- setdiff(row.names(ct), genes)
  ct <- ct[genes, ]
  #print(length(filt_genes))

  #print(dim(ct))
  gene_f <- data.frame(row.names = row.names(ct), features = row.names(ct))
  freq_expressed <- freq
  FCTHRESHOLD <- FCThresh
  sca0 <- FromMatrix(as.matrix(log2(ct + 1)), meta, gene_f)
  if (nbins > 0) {
    thres <- thresholdSCRNACountMatrix(assay(sca0)[rowMedians(assay(sca0)) < max_thres, ], conditions = meta$cellType, nbins = nbins, min_per_bin = min_per_bin, bin_by = bin_by)
    thres_ct <- assay(sca0)
    thres_ct[rowMedians(assay(sca0)) < max_thres, ] <- thres$counts_threshold
    # thres_ct[thres$original_data > 2] <- thres$original_data[thres$original_data > 2]
    assays(sca0) <- list(thresh = thres_ct, tpm = assay(sca0))
    if (plot) {
      par(mfrow = c(ceiling(sqrt(nbins)), ceiling(sqrt(nbins))))
      plot(thres)
    }
  } else {
    assays(sca0) <- list(thresh = assay(sca0), tpm = assay(sca0))
  }

  # assays(sca0) <- list(thresh=assay(sca0), tpm=assay(sca0))
  # expressed_genes <- freq(sca0) > freq_expressed
  # sca0 <- sca0[expressed_genes,]
  # print(dim(sca0))
  # res['sca'] <- sca0
  comparisons <- list(comp = c(control, target))
  for (cnd in comparisons) {
    sca <- sca0[, colData(sca0)$cellType %in% cnd]
    cond <- factor(colData(sca)$cellType)
    cond <- relevel(cond, control)
    other <- paste("cellType", cnd[which(cnd != control)], sep = "")
 
    colData(sca)$cellType <- cond
    if (corr_det) { 
      zlmCond <- zlm(~ cellType + sensitivity, sca, useContinuousBayes = TRUE)
    } else {
      zlmCond <- zlm(~cellType, sca, useContinuousBayes = TRUE)
    }
    summaryCond <- summary(zlmCond, doLRT = other)
    summaryDt <- summaryCond$datatable
    # return(summaryDt)
    fcHurdle <- merge(summaryDt[contrast == other & component == "H", .(primerid, `Pr(>Chisq)`)], # hurdle P values
      summaryDt[contrast == other & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
      by = "primerid"
    ) # logFC coefficients
    # print(dim(fcHurdle))
    fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, "fdr")]
    fcHurdleSig <- merge(fcHurdle, data.table::as.data.table(mcols(sca)), by = "primerid")
    data.table::setorder(fcHurdleSig, fdr)
    row.names(fcHurdleSig) <- fcHurdleSig$primerid
    fcHurdleSig <- data.frame(fcHurdleSig, row.names = 1)
    if (correct_wild_coef) {
      con <- subset(data.frame(summaryDt), component == "C" & contrast == summaryDt$contrast[1])
      con <- data.frame(con, row.names = con$primerid)
      fc <- subset(data.frame(summaryDt), component == "logFC" & contrast == summaryDt$contrast[1])
      fc <- data.frame(fc, row.names = fc$primerid)
      sub_fc <- data.frame(row.names = fc$primerid, primerid = fc$primerid, coef_diff = abs(con$coef - fc$coef), z = abs(con$z) > abs(fc$z))
      congenes <- row.names(subset(sub_fc, coef_diff > 0.1 & z))
      fcHurdleSig[, "orig_logfc"] <- fc[row.names(fcHurdleSig), "coef"]
      fcHurdleSig[, "con_logfc"] <- con[row.names(fcHurdleSig), "coef"]
      fcHurdleSig[congenes, c("coef", "ci.hi", "ci.lo")] <- con[congenes, c("coef", "ci.hi", "ci.lo")]
    }

    fcHurdleSig$Log2FC <- fcHurdleSig$coef
    fcHurdleSig$Log2FC_uncorr<- fcHurdleSig$orig_logfc
    fcHurdle <- data.frame(fcHurdle, row.names = 1)
    fcHurdle$Log2FC <- fcHurdle$coef
    if (include_filt_as_NA) { 
      filt_gene_res <- matrix(nrow = length(filt_genes), ncol = ncol(fcHurdleSig))
      row.names(filt_gene_res) <- filt_genes
      colnames(filt_gene_res) <- colnames(fcHurdleSig)
      filt_gene_res <- data.frame(filt_gene_res)
      filt_gene_res[, 1] <- 1
      filt_gene_res[, "Log2FC"] <- 0
      fcHurdleSig <- rbind(fcHurdleSig, filt_gene_res)
      fcHurdleSig$fdr <- p.adjust(fcHurdleSig[, 1], "BH")
      fcHurdleSig[is.na(fcHurdleSig$Log2FC), "Log2FC"] <- 0
    }
    res[[paste(cnd[1], "_v_", cnd[2], sep = "")]] <- list(DESig = fcHurdleSig, model = zlmCond, DEfull = fcHurdle, samples = colData(sca0)$cellType %in% cnd, data = sca0, summaryDt = summaryCond)
  }

  return(res)
}

# get organism stuff for clusterprofiler
get_organism_items <- function(organisms){
  if(organisms == 'mouse'){
    library(org.Mm.eg.db)
    orgdb = org.Mm.eg.db
    orgabv = 'mmu'
    orgname = "Mus musculus"
  }else if(organisms == 'rat'){
    library(org.Rn.eg.db)
    orgdb = org.Rn.eg.db
    orgabv = 'rno'
    orgname = "Rattus norvegicus"
  }
  else if(organisms == 'celegans'){
    library(org.Ce.eg.db)
    orgdb = org.Ce.eg.db
    orgabv = 'cel'
    orgname = 'Caenorhabditis elegans'
  }
  else if(organisms == 'human'){
    library(org.Hs.eg.db)
    orgdb = org.Hs.eg.db
    orgabv = 'hsa'
    orgname = 'Homo sapiens'
  }
  return(list(orgdb = orgdb, orgkegg = orgabv, orgname = orgname))
}

enrich_CP <- function(ora_genes, organisms, n_type = 'ALIAS',universe = NULL, classic = T, GO_BP_only = F, enrich_all = T, Msig = NULL, alpha = 1, full_combine = T){
  
  items = get_organism_items(organisms = organisms)
  orgabv = items$orgkegg
  orgname = items$orgname
  orgdb = items$orgdb
  
  if(n_type == 'ENSEMBL'){
    ora_genes <- ENSEMBL2ALIAS[ora_genes]
    universe <- ENSEMBL2ALIAS[universe]
    n_type = 'ALIAS'
  }
  
  oraL <- tryCatch({
    egdf <- bitr(ora_genes, n_type, 'ENTREZID', orgdb) %>% distinct(eval(as.name('ALIAS')), .keep_all = T) %>% data.frame(row.names = 1)
    egdf$ENTREZID
  },error=function(cond){return(NULL)})
  universe <-  tryCatch({
    egdf <- bitr(universe, n_type, 'ENTREZID', orgdb) %>% distinct(eval(as.name('ALIAS')), .keep_all = T) %>% data.frame(row.names = 1)
    egdf$ENTREZID
  },error=function(cond){return(NULL)})
  
  if(is.null(oraL) || length(oraL) == 0)
  {return(NULL)}
  
  if(organisms == 'celegans'){
    oraL_kegg <- paste('CELE_', WORM.GENES[ora_genes,]$Sequence.Name, sep ='')
  }else{
    oraL_kegg <- oraL
  }
  #return(list(g = oraL, u=universe))
  
  GSE_results <- list()
  # GO Enrichment
  if(GO_BP_only | classic){
    GSE_results[['GO_BP_ora']] <- tryCatch({setReadable(enrichGO(gene = oraL,
                                                                 universe      = universe,
                                                                 OrgDb         = orgdb,
                                                                 ont           = "BP",
                                                                 pAdjustMethod = "BH",
                                                                 pvalueCutoff  = alpha, maxGSSize = 500,minGSSize = 10, 
                                                                 qvalueCutoff = alpha), OrgDb = orgdb)},error=function(cond){return(NULL)})
  }
  if(!GO_BP_only & classic){
    GSE_results[["WKP_ora"]] <- tryCatch({setReadable(enrichWP(oraL, organism = orgname, maxGSSize = 500, minGSSize = 10, universe = universe, pvalueCutoff  = alpha, qvalueCutoff = alpha ), OrgDb = orgdb)},error=function(cond){return(NULL)})
    
    GSE_results[['GO_CC_ora']] <- tryCatch({setReadable(enrichGO(gene = oraL,
                                                                 universe      = universe,
                                                                 OrgDb         = orgdb,
                                                                 ont           = "CC",
                                                                 pAdjustMethod = "BH",
                                                                 pvalueCutoff  = alpha,
                                                                 qvalueCutoff = alpha,maxGSSize = 500,minGSSize = 10, 
    ), OrgDb = orgdb)},error=function(cond){return(NULL)})
    
    GSE_results[['GO_MF_ora']] <- tryCatch({setReadable(enrichGO(gene = oraL,
                                                                 universe      = universe,
                                                                 OrgDb         = orgdb,
                                                                 ont           = "MF",
                                                                 pAdjustMethod = "BH",
                                                                 pvalueCutoff  = alpha, maxGSSize = 500,minGSSize = 10, 
                                                                 qvalueCutoff = alpha), OrgDb = orgdb)},error=function(cond){return(NULL)})
    
    
    #KEGG Enrichment
    
    GSE_results[['KEGG_ora']] <- tryCatch({setReadable(enrichKEGG(gene = oraL, universe = universe,
                                                                  organism     = orgabv, maxGSSize = 500,minGSSize = 10, 
                                                                  pvalueCutoff  = alpha,qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
    
    
    GSE_results[['MKEGG_ora']] <- tryCatch({setReadable(enrichMKEGG(gene = oraL, universe = universe, maxGSSize = 500,minGSSize = 10, 
                                                                    organism = orgabv, pvalueCutoff  = alpha,qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
    
    # Reactome Enrichment
    GSE_results[['REACT_ora']] <- tryCatch({setReadable(enrichPathway(gene=oraL,organism = organisms, maxGSSize = 500,minGSSize = 10,  universe = universe, pvalueCutoff  = alpha,qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
  }
  
  Msig_res <- NULL
  if(!is.null(Msig)){
    Msig_res <- lapply(Msig, function(x){
      cat = strsplit(x, '-')[[1]][1]
      sub_cat = ''
      if(length(strsplit(x, '-')[[1]]) > 1){
        sub_cat = strsplit(x, '-')[[1]][2]
      }
      df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
      tryCatch({setReadable(enricher(oraL, TERM2GENE=df, maxGSSize = 500,minGSSize = 10,  universe = universe, pvalueCutoff = alpha, qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
    })
  }
  
  if(!is.null(Msig_res)){
    names(Msig_res) <- Msig
    for(n in names(Msig_res)){
      GSE_results[[n]] <- Msig_res[[n]]
    }
  }

  if(full_combine == T){
    alpha = 1
      all_sets <- NULL
      all_sets_n <- NULL
      if(classic){
        wiki <- data.frame(clusterProfiler:::get_wp_data(orgname))
        wiki_t2g <- wiki[, c('wpid', 'gene')]
        colnames(wiki_t2g) <- c('term', 'gene')
        wiki_t2n <- unique(wiki[, c('wpid', 'name')])
        colnames(wiki_t2n) <- c('term', 'name')
        row.names(wiki_t2n) <- wiki_t2n$wpid
        react <- as.list(ReactomePA:::get_Reactome_DATA(organisms))
        go_bp <- as.list(clusterProfiler:::get_GO_data(orgdb, 'BP', "ENTREZID"))
        kegg <- as.list(clusterProfiler:::prepare_KEGG(orgabv, "KEGG", "ncbi-geneid"))
        go_kegg_react_list <- c(react$PATHID2EXTID, go_bp$PATHID2EXTID, kegg$PATHID2EXTID)
        go_kegg_react_names <- c(react$PATHID2NAME, go_bp$PATHID2NAME, kegg$PATHID2NAME)
        go_kegg_react_P2G <- data.frame(do.call(rbind, lapply(names(go_kegg_react_list), FUN = function(x){
          cbind(rep(x, length(go_kegg_react_list[[x]])), go_kegg_react_list[[x]])
        })))
        colnames(go_kegg_react_P2G) <- c('term', 'gene')
        
        go_kegg_react_P2N <- data.frame(term = names(go_kegg_react_names), name = go_kegg_react_names)
        
        all_sets <- rbind(go_kegg_react_P2G, wiki_t2g)
        all_sets_n <- rbind(go_kegg_react_P2N, wiki_t2n)
      }
      
      Msig_df <- NULL
      if(!is.null(Msig)){
        Msig_df <- do.call(rbind, lapply(Msig, function(x){
          cat = strsplit(x, '-')[[1]][1]
          sub_cat = ''
          if(length(strsplit(x, '-')[[1]]) > 1){
            sub_cat = strsplit(x, '-')[[1]][2]
          }
          df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
          colnames(df) <- c('term', 'gene')
          df
        }))
        Msig_df <- data.frame(Msig_df)
        View(Msig_df)
        colnames(Msig_df) <- c('term', 'gene')
        all_sets <- rbind(all_sets, Msig_df)
        all_sets_n <- rbind(all_sets_n, data.frame(term=unique(Msig_df[,1]), name=unique(Msig_df[,1])))
      }
      if(!is.null(all_sets)){
        GSE_results[['combined_full']] <- tryCatch({setReadable(enricher(oraL, TERM2GENE = all_sets,
                                                                         TERM2NAME = all_sets_n, 
                                                                         maxGSSize = 500,minGSSize = 10,  
                                                                         universe = universe, pvalueCutoff = alpha, 
                                                                         qvalueCutoff = alpha), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})}
  }
  
  
  return(GSE_results)
}




gse_CP <- function( organisms, logFC=NULL, n_type = 'ALIAS', classic = T, simplify_go = T, full_combine = T, Msig = NULL, alpha = 1, disease= F){
  if(full_combine){
    alpha = 1
  }
  
  items = get_organism_items(organisms = organisms)
  orgabv = items$orgkegg
  orgname = items$orgname
  orgdb = items$orgdb
  GSE_results <- list()
  gse_list0 <- NULL
  gse_list <- NULL
  if(!is.null(logFC)){ 
    gse_list <- logFC
    if(n_type == 'ENSEMBL'){
      gse_list <- tapply(gse_list, ENSEMBL2ALIAS[names(gse_list)], FUN = function(x){
        x[which(abs(x) == max(abs(x)))]
      })
      n_type ='ALIAS'
      gnames <- names(gse_list)
      gse_list <- as.vector(gse_list)
      names(gse_list) <- gnames
    }
    gse_list <- sort(gse_list, T)
    gse_list0 <- sort(gse_list0, T)
    egdf <- bitr(names(gse_list), n_type, 'ENTREZID', orgdb) %>% distinct(eval(as.name('ALIAS')), .keep_all = T) %>% data.frame(row.names = 1)
    gse_list <- gse_list[intersect(names(gse_list), row.names(egdf))]
    names(gse_list) <- egdf[names(gse_list),]$ENTREZID
    #return(gse_list)
    if(classic){
      GSE_results[["WKP_gse"]] <- setReadable(gseWP(gse_list, eps = 0, organism = orgname, minGSSize = 10, nPermSimple = 100000, maxGSSize = 500, pvalueCutoff=alpha), OrgDb = orgdb, keyType = 'ENTREZID')
      #GSE_results[["BIO_gse"]] <- setReadable0(GSEA(gse_list, TERM2GENE = BIOCYC[,c(1,2)], TERM2NAME = BIOCYC[,c(1,3)], nPerm = 1000, minGSSize = 5, maxGSSize = 500), gene2symbol = entrez2symbol, keyType = 'ENTREZ')
      GSE_results[['GO_BP_gse']] <- setReadable(gseGO(geneList = gse_list,
                                                      OrgDb        = orgdb, 
                                                      keyType = "ENTREZID", nPermSimple = 10000,
                                                      ont          = "BP",eps = 0,
                                                      minGSSize    = 10,
                                                      maxGSSize    = 500,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE, by = 'fgsea'), OrgDb = orgdb, keyType = 'ENTREZID')
      
      GSE_results[['GO_CC_gse']] <- setReadable(gseGO(geneList = gse_list,
                                                      OrgDb        = orgdb,
                                                      keyType = "ENTREZID",
                                                      ont          = "CC",nPermSimple = 10000,
                                                      minGSSize    = 10,eps = 0,
                                                      maxGSSize    = 500,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE), OrgDb = orgdb, keyType = 'ENTREZID')
      
      GSE_results[['GO_MF_gse']] <- setReadable(gseGO(geneList     = gse_list,
                                                      OrgDb        = orgdb,keyType = "ENTREZID",
                                                      ont          = "MF",
                                                      minGSSize    = 10,eps = 0,
                                                      maxGSSize    = 500,nPermSimple = 10000,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE), OrgDb = orgdb, keyType = 'ENTREZID')
      GSE_results[['KEGG_gse']]<- setReadable(gseKEGG(geneList = gse_list, 
                                                      organism     = orgabv,
                                                      minGSSize = 10, eps = 0,
                                                      maxGSSize = 500,nPermSimple = 10000,
                                                      pvalueCutoff = alpha,
                                                      verbose      = FALSE, keyType = 'ncbi-geneid'),  OrgDb = orgdb, keyType = 'ENTREZID')
      GSE_results[['MKEGG_gse']] <- setReadable(gseMKEGG(gene = gse_list,  minGSSize = 10, eps = 0, nPermSimple = 10000, maxGSSize = 500, organism = orgabv, pvalueCutoff=alpha),  OrgDb = orgdb, keyType = 'ENTREZID')
      GSE_results[['REACT_gse']] <- setReadable(gsePathway(geneList =  gse_list, organism = organisms, nPermSimple = 10000, minGSSize = 10, maxGSSize = 500, 
                                                           pvalueCutoff=alpha,eps = 0,
                                                           pAdjustMethod="BH", verbose=FALSE), OrgDb = orgdb, keyType = 'ENTREZID')
    }
    Msig_res <- NULL
    if(!is.null(Msig)){
      Msig_res <- lapply(Msig, function(x){
        cat = strsplit(x, '-')[[1]][1]
        sub_cat = ''
        if(length(strsplit(x, '-')[[1]]) > 1){
          sub_cat = strsplit(x, '-')[[1]][2]
        }
        df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
        tryCatch({setReadable(GSEA(
          gse_list,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 0,nPermSimple = 10000,
          pvalueCutoff = alpha,
          pAdjustMethod = "BH",
          TERM2GENE = df,
          verbose = TRUE,
          seed = FALSE,
          by = "fgsea"), OrgDb = orgdb, keyType = 'ENTREZID')},error=function(cond){return(NULL)})
      })
    }
    
    if(!is.null(Msig_res)){
      names(Msig_res) <- Msig
      for(n in names(Msig_res)){
        GSE_results[[n]] <- Msig_res[[n]]
      }
    }
    #for(r in names(GSE_results)){
    #GSE_results[[r]]@result <- subset(GSE_results[[r]]@result, qvalue < 0.05)
    #}

    if(full_combine == T){
        all_sets <- NULL
        all_sets_n <- NULL
        if(classic){
          wiki <- data.frame(clusterProfiler:::get_wp_data(orgname))
          wiki_t2g <- wiki[, c('wpid', 'gene')]
          colnames(wiki_t2g) <- c('term', 'gene')
          wiki_t2n <- unique(wiki[, c('wpid', 'name')])
          colnames(wiki_t2n) <- c('term', 'name')
          row.names(wiki_t2n) <- wiki_t2n$wpid
          react <- as.list(ReactomePA:::get_Reactome_DATA(organisms))
          go_bp <- as.list(clusterProfiler:::get_GO_data(orgdb, 'BP', "ENTREZID"))
          kegg <- as.list(clusterProfiler:::prepare_KEGG(orgabv, "KEGG", "ncbi-geneid"))
          go_kegg_react_list <- c(react$PATHID2EXTID, go_bp$PATHID2EXTID, kegg$PATHID2EXTID)
          go_kegg_react_names <- c(react$PATHID2NAME, go_bp$PATHID2NAME, kegg$PATHID2NAME)
          go_kegg_react_P2G <- data.frame(do.call(rbind, lapply(names(go_kegg_react_list), FUN = function(x){
            cbind(rep(x, length(go_kegg_react_list[[x]])), go_kegg_react_list[[x]])
          })))
          colnames(go_kegg_react_P2G) <- c('term', 'gene')
          
          go_kegg_react_P2N <- data.frame(term = names(go_kegg_react_names), name = go_kegg_react_names)
          
          all_sets <- rbind(go_kegg_react_P2G, wiki_t2g)
          all_sets_n <- rbind(go_kegg_react_P2N, wiki_t2n)
        }
        Msig_df <- NULL
        if(!is.null(Msig)){
          Msig_df <- do.call(rbind, lapply(Msig, function(x){
            cat = strsplit(x, '-')[[1]][1]
            sub_cat = ''
            if(length(strsplit(x, '-')[[1]]) > 1){
              sub_cat = strsplit(x, '-')[[1]][2]
            }
            df = get_msig(organisms, cat = cat, sub_cat = sub_cat, gmt_dir = './dataset/Msigdb/') 
            colnames(df) <- c('term', 'gene')
            df
          }))
          Msig_df <- data.frame(Msig_df)
          colnames(Msig_df) <- c('term', 'gene')
          all_sets <- rbind(all_sets, Msig_df)
          all_sets_n <- rbind(all_sets_n, data.frame(term=unique(Msig_df[,1]), name=unique(Msig_df[,1])))
        }
        GSE_results[['combined_full']] <- setReadable(GSEA(
          gse_list,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 0,nPermSimple = 10000,
          pvalueCutoff = alpha,
          pAdjustMethod = "BH",
          gson = NULL,
          TERM2GENE = all_sets,
          TERM2NAME = all_sets_n,
          verbose = TRUE,
          seed = FALSE,
          by = "fgsea",
        ), OrgDb = orgdb, keyType = 'ENTREZID')
      }
    }
    
  GSE_results[['gfc0']] <- gse_list0
  GSE_results[['gfc']] <- gse_list
  
  return(GSE_results)
}



deg_utr2 <- function(file, ct, compare, meta, impute = F, 
                     method = "fisher.test", filter_by_PAS_motif = TRUE,
                     filter_by_APA_dist = TRUE, 
                     alpha = 0.05, combine_p = NULL, diff_thresh = 0.05,
                     fasta_file = "../mouse_rat_proposal//dataset/igv/mouse/mouse_spike.fa",
                     expr_filt = 10, min_samp_filt = 5, na_filt = 5, 
                     utr_cov_filt = 5, old_apa_dist = T, filter_by_pas_dist = NULL) {
  
  ## read in the new dapars file that includes long. short and PDUI values
  dapars <- read.csv(file, sep = "\t", header = T, row.names = 1)
  # dapars <- dapars[,!grepl(regex('^X'), colnames(dapars))]
  # dapars <- cbind(dapars[,1:3],dapars[,grepl(regex(paste(utr_samples, collapse = '|')), colnames(dapars))])
  # colnames(dapars) <- str_remove(colnames(dapars), '^X')
  samples <- row.names(subset(meta, cellType %in% compare))
  dapars <- cbind(dapars[, c(1, 2, 3)], dapars[, grepl(regex(paste(samples, collapse = "|")), colnames(dapars))])
  ct <- ct[, samples]
  
  dapars$gene_short_names <- sapply(row.names(dapars), FUN = function(x) {
    strsplit(x, "\\|")[[1]][2]
  })
  dapars_orig <- dapars
  all_genes <- unique(dapars$gene_short_names)
  cat("Total genes assessed: ", length(unique((dapars$gene_short_names))), "genes\n")
  ## Filter first based on gene expression overall in both conditions (we don't want differential expression to affect results too much)
  genes <- row.names(ct)[rowSums(ct[,row.names(subset(meta, cellType == compare[1]))] > expr_filt) > min_samp_filt]
  genes <- intersect(genes, row.names(ct)[rowSums(ct[,row.names(subset(meta, cellType == compare[2]))] > expr_filt) > min_samp_filt])
  
  
  
  
  ## Only account for UTR coverage in genes that are assigned at have average expression of at least 5 uniquely mapped reads across all samples
  genes <- row.names(ct)[rowMeans(ct) > expr_filt]
  
  dapars <- dapars[dapars$gene_short_names %in% genes, ]
  cat("filtering based on gene expression: left with ", length(unique((dapars$gene_short_names))), "genes\n")
  # dapars <- subset(dapars, fit_value >= 10) # Maybe filter also based on regression fit value
  cat("filtering based on fit value: left with ", length(unique((dapars$gene_short_names))), "genes\n")

  # vector to store gene name and UTR region correspondence in case gene names is lost with imputation
  gene2region <- dapars$gene_short_names
  names(gene2region) <- sapply(row.names(dapars), FUN = function(x) {
    strsplit(x, "\\|")[[1]][1]
  })
  colnames(dapars) <- str_remove(colnames(dapars), pattern = "wig.")
  # split file into long, short and pdui
  d_long <- dapars[, grepl("long_exp", colnames(dapars))]
  d_short <- dapars[, grepl("short_exp", colnames(dapars))]
  d_pdui <- dapars[, grepl("PDUI", colnames(dapars))]

  # change column names
  colnames(d_long) <- sapply(strsplit(colnames(d_long), "_"), FUN = function(x) {
    strsplit(x[1], "\\.")[[1]][1]
  })
  colnames(d_short) <- sapply(strsplit(colnames(d_short), "_"), FUN = function(x) {
    strsplit(x[1], "\\.")[[1]][1]
  })
  colnames(d_pdui) <- sapply(strsplit(colnames(d_pdui), "_"), FUN = function(x) {
    strsplit(x[1], "\\.")[[1]][1]
  })
  d_long <- d_long[,row.names(meta)]
  d_short <- d_short[,row.names(meta)]
  d_pdui<- d_pdui[,row.names(meta)]
  # select filtering genes based on number of passes (Non NAs) and overall coverage (average > 2 either in long or short UTR in both conditions)
  grp <- meta[colnames(d_pdui), "cellType"]
  btch <- meta[colnames(d_pdui), "experiment"]
  cond1_ind <- which(grp == compare[1])
  cond2_ind <- which(grp == compare[2])
  # print(d_pdui[dapars$gene_short_name == 'Cdk1',])
  ## NA FILTER
  # onyly keep genes that have at least than na_filt non-NAs in terms of coverage in both conditions
  na.filt.genes <- rowSums(!is.na(d_pdui[, cond1_ind])) >= na_filt & rowSums(!is.na(d_pdui[, cond2_ind])) >= na_filt

  dapars <- dapars[na.filt.genes, ]

  cat("filtering based on number of non NAs in each condition ( > ", na_filt," in each condition): left with ", length(unique((dapars$gene_short_names))), "genes\n")


  d_long <- d_long[na.filt.genes, ]
  d_short <- d_short[na.filt.genes, ]
  d_pdui <- d_pdui[na.filt.genes, ]

  # Filter based on coverage of UTR regions
  c1.filt.genes <- rowMeans(d_long[, cond1_ind], na.rm = T) > utr_cov_filt | rowMeans(d_short[, cond1_ind], na.rm = T) > utr_cov_filt
  c2.filt.genes <- rowMeans(d_long[, cond2_ind], na.rm = T) > utr_cov_filt | rowMeans(d_short[, cond2_ind], na.rm = T) > utr_cov_filt

  final.filt.genes <- c1.filt.genes & c2.filt.genes

  # filtering matrices with genes selected prior
  dapars <- dapars[final.filt.genes, ]
  d_long <- d_long[final.filt.genes, ]
  d_short <- d_short[final.filt.genes, ]
  #d_pdui <- d_pdui[final.filt.genes, ]
  #d_long1 <- d_long[final.filt.genes, ] %>% mutate(across(everything(), replace_na, 0))
  #d_short1 <- d_short[final.filt.genes, ] %>% mutate(across(everything(), replace_na, 0))
  d_pdui <- d_long / (d_long + d_short)
  cat("filtering based on number of mean coverage of long/short UTR in both condition: left with ", length(unique((dapars$gene_short_names))), "genes\n")
  if (impute) {
    dapars_out <- data.frame(Gene = row.names(dapars), dapars[, 1:3], d_pdui)
    write.table(dapars_out, file = "./temp.dp.tsv", sep = "\t", quote = F, row.names = F)
    d_pdui <- scDaPars(
      raw_PDUI_file = "./temp.dp.tsv",
      out_dir = "apa/scDaPars_result",
      filter_gene_thre = 0.2,
      filter_cell_thre = 0.1
    )
    return(d_pdui)
    method <- "ks.test"
  }
  cat("performing tests\n")
  if (method == "ks.test") {
    test <- apply(d_pdui, 1, FUN = function(x) {
      if (sum(x[!is.na(x)]) == 0 | sum(!is.na(x[cond1_ind])) <= 3 | sum(!is.na(x[cond2_ind])) <= 3) {
        c(1, 0)
      } else {
        c(ks.test(x[cond1_ind][!is.na(x[cond1_ind])], x[cond2_ind][!is.na(x[cond2_ind])])$p.value, mean(x[cond2_ind][!is.na(x[cond2_ind])]) - mean(x[cond1_ind][!is.na(x[cond1_ind])]))
      }
    })
    test <- t(test)
  } else if (method == "binom") {
    test <- do.call(rbind, pblapply(row.names(d_pdui), FUN = function(n) {
      x <- d_pdui[n, ]
      w <- d_long[n, ] + d_short[n, ]
      not_na <- !is.na(x)
      x0 <- c(x[not_na])
      w0 <- c(w[not_na])
      data <- data.frame(pdui = x0, cellType = grp[not_na])

      test_0 <- tryCatch(
        {
          mylogit <- glm(pdui ~ cellType, data = data, family = "quasibinomial", weights = w0)
          mylogit0 <- glm(pdui ~ 1, data = data, family = "quasibinomial", weights = w0)
          pval <- anova(mylogit, mylogit0, test = "F")[, "Pr(>F)"][2]
          diff <- predict(mylogit, data.frame(cellType = compare[2], batch = "A"), type = "response") - predict(mylogit, data.frame(cellType = compare[1], batch = "A"), type = "response")
          return(c(pval, diff))
        },
        error = function(cond) {
          return(c(1, 0))
        }
      )
      return(test_0)
    }))
    test[is.na(test[, 1]), ] <- c(1, 0)
  } else if (method == "betab") {
    test <- do.call(rbind, pblapply(row.names(d_pdui), FUN = function(n) {
      xx <- d_pdui[n, ]
      x <- d_long[n, ]
      w <- d_long[n, ] + d_short[n, ]
      not_na <- !is.na(xx)
      x0 <- round(c(x[not_na]))
      w0 <- round(c(w[not_na]))
      data <- data.frame(y = x0, n = w0, cellType = grp[not_na])
      if (sum(x0 == w0) == length(x0)) {
        return(c(1, 0))
      }
      test_0 <- tryCatch(
        {
          mylogit <- betabin(formula = cbind(y, n - y) ~ cellType, random = ~1, data = data)
          mylogit0 <- betabin(formula = cbind(y, n - y) ~ 1, random = ~1, data = data)
          pval <- anova(mylogit, mylogit0)@anova.table[, "P(> Chi2)"][2]
          preds <- predict(mylogit, data.frame(cellType = c(compare[2], compare[1])))
          diff <- preds[2] - preds[1]
          return(c(pval, diff))
        },
        error = function(cond) {
          print(cond)
          print(data)
          return(c(1, 0))
        }
      )
      return(test_0)
    }))
    test[is.na(test[, 1]), ] <- c(1, 0)
  } else if (method == "betareg") {
    test <- do.call(rbind, pblapply(row.names(d_pdui), FUN = function(n) {
      x <- d_pdui[n, ]
      not_na <- !is.na(x)
      x0 <- c(x[not_na])
      if (sum(x0 == 0) > 0 | sum(x0 == 1) > 0) {
        x0 <- (x0 * (length(x0) - 1) + 0.5) / length(x0)
      }
      if (length(unique(x0)) == 1) {
        return(c(1, 0))
      }
      data <- data.frame(pdui = x0, cellType = grp[not_na])
      data$cellType <- as.factor(data$cellType)
      test_0 <- tryCatch(
        {
          mylogit <- betareg(pdui ~ cellType, data = data, link = "log")
          mylogit0 <- betareg(pdui ~ 1, data = data, link = "log")
          pval <- lrtest(mylogit, mylogit0)[, "Pr(>Chisq)"][2]
          diff <- predict(mylogit, data.frame(cellType = compare[2])) - predict(mylogit, data.frame(cellType = compare[1]))
          return(c(pval, diff))
        },
        error = function(cond) {
          print(cond)
          return(c(1, 0))
        }
      )
      return(test_0)
    }))
  } else {
    l1 <- length(cond1_ind)
    l2 <- length(cond2_ind)
    d_long1 <- d_long
    #d_long1[is.na(d_long1)] <- 0
    d_short1 <- d_short
    #d_short1[is.na(d_short1)] <- 0
    d_pdui1 <- d_pdui
    utrl1_mean <- round(rowMeans(d_long1[, cond1_ind], na.rm = T))
    utrl2_mean <- round(rowMeans(d_long1[, cond2_ind], na.rm = T))
    utrs1_mean <- round(rowMeans(d_short1[, cond1_ind], na.rm = T))
    utrs2_mean <- round(rowMeans(d_short1[, cond2_ind], na.rm = T))
    pdui1_mean <- rowMeans(d_pdui1[, cond1_ind], na.rm = T)
    pdui2_mean <- rowMeans(d_pdui1[, cond2_ind], na.rm = T)
    #pdui1_mean <- utrl1_mean/(utrl1_mean + utrs1_mean )
    #pdui2_mean <- utrl2_mean/(utrl2_mean + utrs2_mean )
    test <- do.call(rbind, pblapply(row.names(d_long1), FUN = function(x) {
      # utr_l1 <- d_long[x,cond1_ind][!is.na(d_long[x,cond1_ind])] # long utr coverage in condition 1
      # utr_l2 <- d_long[x,cond2_ind][!is.na(d_long[x,cond2_ind])] # short utr coverage in condition 2
      # utr_s1 <- d_short[x,cond1_ind][!is.na(d_short[x,cond1_ind])] # long utr coverage in condition 1
      # utr_s2 <- d_short[x,cond2_ind][!is.na(d_short[x,cond2_ind])] # short utr coverage in condition 2
      pdui_1 <- pdui1_mean[x]
      pdui_2 <- pdui2_mean[x]
    
      
      # c(fisher.test(x = rbind(c(mean(utr_l1), mean(utr_s1)), c(mean(utr_l2), mean(utr_s2))))$p.value, pdui_2 - pdui_1)
      twobytwo <- rbind(c(utrl1_mean[x], utrs1_mean[x]), c(utrl2_mean[x], utrs2_mean[x]))
      c(fisher.test(x = twobytwo)$p.value, pdui_2 - pdui_1)


      # if(x == 'XM_039113041.1|Pou2f2|NC_051336.1|-'){
      # print(rbind(c(sum(utr_l1)/l1, sum(utr_s1)/l1), c(sum(utr_l2)/l1, sum(utr_s2)/l2)))
      # }
    }))
  }
  cat("finished\n")
  # print(dim(test))
  # print(dim(d_long))
  row.names(test) <- row.names(d_long)
  colnames(test) <- c("pval", "mean.diff")
  test <- data.frame(test)
  test$pval[test$pval > 1] <- 1
  # test$padj <- p.adjust(test$pval, method = 'BH')
  # test$fdr <- qvalue(test$pval)$qvalue

  if (impute) {
    test$gene_short_names <- gene2region[row.names(test)]
  } else {
    test$gene_short_names <- sapply(row.names(test), FUN = function(x) {
      strsplit(x, "\\|")[[1]][2]
    })
  }
  # test$diff <- abs(test$mean.diff) > diff_thresh & test$fdr < alpha
  test$fit_value <- dapars[row.names(test), ]$fit_value
  test$predicted_p_APA <- dapars_orig[row.names(test), ]$Predicted_Proximal_APA
  test$loci <- dapars_orig[row.names(test), ]$Loci
  test$strand <- sapply(strsplit(row.names(test), "\\|"), FUN = function(x) {x[4]})
  test$APA_dist <- 0
  if(!old_apa_dist){
    test$utr_length <- sapply(row.names(test), function(x){as.numeric(strsplit(x, "\\|")[[1]][3])})
    test$APA_dist <- sapply(seq_len(nrow(test)), function(i){
      loci = as.numeric(strsplit(strsplit(test[i, 'loci'], ":")[[1]][2], '-')[[1]])
      if(test[i, 'strand'] == '+'){
        test[i, 'predicted_p_APA']-(loci[2]-test[i, 'utr_length'])-1
      }else{
        loci[1]+test[i, 'utr_length']- test[i, 'predicted_p_APA']-1
     }
    })
    test$offset <- sapply(seq_len(nrow(test)), function(i){
      loci = as.numeric(strsplit(strsplit(test[i, 'loci'], ":")[[1]][2], '-')[[1]])
      if(test[i, 'strand'] == '+'){
        loci[2]- test[i, 'utr_length'] - loci[1]
      }else{
        loci[2] - (loci[1]+test[i, 'utr_length'])
      }
    })
  }else{
    test[test$strand == "+", ]$APA_dist <- abs(sapply(strsplit(test[test$strand == "+", ]$loci, "-"),
      FUN = function(x) {
        as.numeric(strsplit(x[1], ":")[[1]][2])
      }
    ) - test[test$strand == "+", ]$predicted_p_APA) - 1
    test[test$strand == "-", ]$APA_dist <- abs(sapply(strsplit(test[test$strand == "-", ]$loci, "-"),
      FUN = function(x) {
        as.numeric(x[2])
      }
    ) - test[test$strand == "-", ]$predicted_p_APA) - 1
  }

  if (filter_by_PAS_motif) {
    test <- post_dapars_pas_filter(test, fasta_file, up_range = 80, down_range = 120, offset = 0)$PAS_motif
    test <- subset(test, num_motif > 0)
  }
  if (!is.null(filter_by_pas_dist)){
    test <- subset(test, APA_dist > filter_by_pas_dist)
  }
  test$padj <- p.adjust(test$pval, method = "BH")
  test$fdr <- qvalue(test$pval)$qvalue
  test$diff <- abs(test$mean.diff) > diff_thresh & test$fdr < alpha
  cat("adjusting p values and combining p values\n")
  # test$APA_dist <- abs(sapply(strsplit(test$loci, '-'), FUN = function(x){as.numeric(x[2])}) - test$predicted_p_APA)-1

  gene_res <- data.frame(do.call(rbind, tapply(row.names(test), test$gene_short_names, function(x) {
    df <- test[x, ]
    min_pval <- min(df[, "pval"])
    min_pval_ind <- which(df[, "pval"] == min_pval)
    min_pval_ind <- min_pval_ind[which(abs(df[min_pval_ind, "mean.diff"]) == max(abs(df[min_pval_ind, "mean.diff"])))][1]
    if (!is.null(combine_p)) {
      df[min_pval_ind, "pval"] <- metapod::combineParallelPValues(as.list(df[, "pval"]), method = combine_p)$p.value
    }
    return(cbind(df[min_pval_ind, ], dapars[x[min_pval_ind], c(1, 2, 3)]))
  })))
  gene_res$padj <- p.adjust(gene_res$pval)
  gene_res$fdr <- qvalue(gene_res$pval)$qvalue
  gene_res$diff <- abs(gene_res$mean.diff) > diff_thresh & gene_res$fdr < alpha
  pdui <- dapars_orig[, grepl("PDUI", colnames(dapars_orig))]
  colnames(pdui) <- colnames(d_long)
  pdui <- pdui[row.names(d_long), ]
  pdui_impute <- t(apply(pdui, 1, FUN = function(x) {
    x[is.na(x)] <- mean(x, na.rm = T)
    x
  }))

  return(list( deg = test, long = d_long, gene_res = gene_res, short = d_short, df = dapars_orig, pdui = pdui, pdui_imp = pdui_impute, gene_universe = all_genes))
}



fisher_proportion_test <- function(nascent, mature, groups, control = "mouseEgg") {
  grp1 <- which(groups == control)
  grp2 <- which(groups == setdiff(groups, control))
  nascent_c <- nascent[, grp1]
  nascent_t <- nascent[, grp2]
  mature_c <- mature[, grp1]
  mature_t <- mature[, grp2]
  nascent_p1 <- rowMeans(nascent_c / (mature_c + nascent_c), na.rm = T)
  nascent_p01 <- rowSums(nascent_c) / rowSums(mature_c)

  nascent_p2 <- rowMeans(nascent_t / (mature_t + nascent_t), na.rm = T)
  nascent_p02 <- rowSums(nascent_t) / rowSums(mature_t)
  res <- data.frame(do.call(rbind, lapply(1:nrow(nascent), FUN = function(i) {
    # diff <- c(nascent_p2[i]- nascent_p1[i])/mean(nascent[i,]/c(mature[i,]+nascent[i,]), na.rm=T)
    diff <- log2(nascent_p02[i] / nascent_p01[i]) # /(sum(nascent[i,])/sum(mature[i,]))
    c(fisher.test(x = round(rbind(c(mean(nascent_c[i, ]), mean(mature_c[i, ])), c(mean(nascent_t[i, ]), mean(mature_t[i, ])))))$p.value, diff)
  })))

  row.names(res) <- row.names(nascent)
  colnames(res) <- c("fisher.p", "prop_diff")
  res[res$fisher.p >= 1, 1] <- 1
  res$qvalue <- qvalue(res[, 1])$qvalue
  res$features <- row.names(res)
  res
}


stageR_dexseqRes <- function(dex) {
  # dex = DEXSeqResults(dex_obj, independentFiltering = FALSE)
  qval <- perGeneQValue(dex)
  gene_pvals <- data.frame(groupID = names(qval), padj = qval)
  # gene_pvals <- read.csv(per_gene, header = T,  sep = '\t')
  t_pvals <- dex$pvalue
  t_pvals[is.na(t_pvals)] <- 1
  # res.t = DRIMSeq::results(drim, level = "feature")
  # res.t$pvalue <- no.na(res.t$pvalue)
  pScreen <- gene_pvals$padj
  names(pScreen) <- gene_pvals$groupID
  pConfirmation <- matrix(t_pvals, ncol = 1)
  dimnames(pConfirmation) <- list(dex$featureID, "transcript")
  # View(pConfirmation)
  tx2gene <- data.frame(dex[, c("featureID", "groupID")], dex[, c("featureID", "groupID")])
  # View(tx2gene)
  stageRObj <- stageRTx(
    pScreen = pScreen,
    pConfirmation = pConfirmation,
    pScreenAdjusted = T,
    tx2gene = tx2gene[, 1:2]
  )

  stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
  drim.padj <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = F)
}

prop_diff <- function(props, group) {
  lvl <- unique(group)
  diff <- abs(rowMeans(props[, group == lvl[1]], na.rm = T) - rowMeans(props[, group == lvl[2]], na.rm = T))
}


get_dexseq_prop <- function(dex_obj) {
  drim.prop <- reshape2::melt(dex_obj$dex_norm, id = c("groupID", "featureID"))
  drim.prop <- drim.prop[order(drim.prop$groupID, drim.prop$variable, drim.prop$featureID), ]

  # Calculate proportions from counts
  system.time({
    drim.prop <- drim.prop %>%
      group_by(groupID, variable) %>%
      mutate(total = sum(value)) %>%
      group_by(variable, add = TRUE) %>%
      mutate(prop = value / total)
  })

  # Convert the data.frame to wide-format data with reshape2::dcast
  drim.prop <- reshape2::dcast(drim.prop[, c(1, 2, 3, 6)], groupID + featureID ~ variable)
}



plotDEXSeqDTU <- function(expData = NULL, geneID = NULL, samps = NULL, isProportion = FALSE) {
  colnames(expData)[1:2] <- c("gid", "tid")
  sub <- subset(expData, gid == geneID)
  colnames(samps) <- c("sample_id", "group")
  sub <- reshape2::melt(sub, id = c("gid", "tid"))
  sub <- merge(samps, sub, by.x = "sample_id", by.y = "variable")
  if (!isProportion) {
    sub$value <- log2(sub$value + 1)
  }

  # clrs = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
  # "deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
  # "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")
  clrs <- DOT_COLOR
  p <- ggplot(sub, aes(tid, value, color = group, fill = group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.8, lwd = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", shape = 5, size = 3, position = position_dodge(width = 0.8)) +
    scale_color_manual(values = clrs) +
    scale_fill_manual(values = clrs) +
    geom_quasirandom(size = 1, dodge.width = 0.8, alpha = 0.8) +
    theme_bw() +
    ggtitle(geneID) +
    xlab("Transcripts") +
    theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust = 0.5, size = 10), 
          axis.title = element_text(size = 12), title = element_text(size = 12),
          legend.position = 'none')

  if (!isProportion) {
    p <- p + ylab("log(Expression)")
  } else {
    p <- p + ylab("Proportions")+ylim(c(0,1))
  }
  p
}


cp_tree_ridge_plot <- function(res, n_cat = 50, nclust = 8, alpha = 0.05, geneSet = NULL) {

  res@result <- subset(res@result, qvalue < alpha)
  if (sum(grepl("GO", res@result$ID))) {
    res@setType <- "BP"
    # res <- clusterProfiler::simplify(res)
  }
  gcolors <- paletteer::paletteer_d("ggsci::springfield_simpsons")[1:nclust]
  sets <- geneSet
  if (is.null(geneSet)) {
    sets <- 1:20
  }

  # res@result$Description[nchar(res@result$Description) > 60] <- res@result$ID[nchar(res@result$Description) > 60]
  res@result$Description <- stringr::str_wrap(res@result$Description, 30)
  res@result <- res@result[sets, ]
  res@result <- res@result[tapply(res@result$ID, res@result$Description, function(x) x[1]), ]
  res <- enrichplot::pairwise_termsim(res)
  # res@result <- res@result[tapply(res@result$ID, res@result$Description, function(x) x[order(res@result[x, 'qvalue'])][1]),]
  res_tree <- addSmallLegend(enrichplot::treeplot(res,
    showCategory = min(n_cat, nrow(res@result)),
    geneClusterPanel = "pie",
    cluster.params = list(method = "ward.D2", n = nclust, color = gcolors, label_words_n = 4, label_format = 25),
    color = NULL, offset_tiplab = 0.8, fontsize = 4
  ) +
    geom_tiplab(
      offset = 0.8, hjust = 0,
      show.legend = FALSE,
      align = TRUE, size = 3.5, lineheight = 0.75
    ) + xlim(c(0, 20)) + scale_size(
      name = "number of genes",
      limits = c(10,110),
      breaks = c(25,50,75,100),
      range = c(1, 5)
    ), pointSize = 2, textSize = 10)
  tree <- res_tree
  res_tree$layers[c(7, 8)] <- NULL
  res_tree$layers[c(3, 4)] <- NULL
  res@result$qvalues <- res@result$qvalue
  res_ridge <- addSmallLegend(ridgeplot(res, showCategory = min(n_cat, nrow(res@result)), fill = 'qvalues') + 
                                scale_fill_viridis_c(name = 'FDR', limits= c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04,0.05))+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + xlim(c(-4, 4)) + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
      xlab(expression("log"[2] * "FC")), pointSize = 3, textSize = 10, spaceLegend = 1.5)

  res_ridge$data$label <- res_ridge$data$category
  # return(list(ridge = res_ridge, tree=res_tree))
  ridge_tree <- res_ridge %>% insert_left(res_tree, width = 3)
  #print(ridge_tree)
  return(list(tree = res_tree, ridge = res_ridge, both = ridge_tree))
  # res_ridge  %>% insert_left(res_tree, width = 2)
}


volcano_plot <- function(res, pval_col = 'fdr', top_genes = NULL, alpha = 0.05, fc = log2(2), main = NULL){
  res <- subset(res, !is.na(res$Log2FC))
  res$padjust <- res[[pval_col]]
  if(!'gene_short_name' %in% colnames(res)){
    res$gene_short_name <- row.names(res)
  }
  if(is.null(top_genes)){
    top_up_genes <- row.names(res[order(res$Log2FC, decreasing = T),])[1:20]
    top_down_genes <- row.names(res[order(res$Log2FC),])[1:20]
    top_genes <- c(top_up_genes, top_down_genes)
  }
  res$deg <- "NoDE"
  res[res$fdr < alpha & res$Log2FC > fc,]$deg <- "up-regulated"
  res[res$fdr < alpha & res$Log2FC < -fc,]$deg <- "down-regulated"
  res[!row.names(res) %in% top_genes,]$gene_short_name <- NA
  res[res$deg == 'NoDE',]$gene_short_name <- NA
  res[grepl('LOC',row.names(res)) | grepl('Rik', row.names(res)) | grepl('Gm', row.names(res)),]$gene_short_name <- NA
  p <- ggplot(data=res, aes(x=Log2FC, y=-log10(padjust), color = deg, label=gene_short_name)) + geom_point(size = 0.2)+geom_label_repel(size=3, max.overlaps = 300) 
  theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 16), 
    axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 16), 
    legend.position = "none",
  )
  p <- p+theme+
    xlim(min(res$Log2FC),max(res$Log2FC))+ 
    ylim(0, 13)+ggtitle(main)+
    ylab(expression('-Log'['10']*'Pvalue'))+
    xlab(expression('-Log'['2']*'FC'))+
    scale_color_manual(values = c( "#6f95e6", "grey","#a4302a"))+
    geom_vline(xintercept = -1, linetype="dashed",color = "red", size=1)+
    geom_vline(xintercept = 1, linetype="dashed",color = "red", size=1)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",color = "red", size=1)
  return(p)
}



custom_cnet_plot <- function(cp_res, top_n_cat = 10, seed=12345, category = NULL, gene_color = NULL, gene_color2 = NULL, layout = 'fr', color_cat_pval = F){
  cp_res1 <- cp_res
  if(!is.null(category)){
    category = intersect(category, row.names(cp_res1@result))
    cp_res1@result <- cp_res1@result[category,]
  }
  top_n_cat <- min(max(length(category), top_n_cat),nrow(cp_res1))
  cp_res1@pvalueCutoff <- 1
  cp_res1@qvalueCutoff <- 1
  pal <- c(rev(c(colorRampPalette(brewer.pal(9, 'Blues'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Blues'))(128)[64:128])), 
           c(colorRampPalette(brewer.pal(9, 'Reds'))(32)[1:16],colorRampPalette(brewer.pal(9, 'Reds'))(128)[64:128]))
  set.seed(seed)
  plot1 <- cnetplot( cp_res1, showCategory = top_n_cat, foldChange = gene_color, layout = layout, node_label = 'none')+
    scale_colour_gradientn(name = expression('Log'['2']*'FC'), limits= c(-4, 4), 
                           colours = pal, 
                           na.value = 'black')
  cat_df <- plot1$data[1:top_n_cat,]
  cat_df$qval <- cp_res@result[match(cat_df$name, cp_res@result$Description),]$qvalue
  cat_df$qval[cat_df$qval > 0.1] <- NA
  if(!color_cat_pval){
    plot1 <- plot1+ggnewscale::new_scale_color() +
      ggraph::geom_node_point(aes_(size=~size), data = cat_df, color = 'black') + 
      scale_size(limits = c(1,150), breaks = c(10,20,40,80), range = c(1,5))
  }
  plot1$data$name <- stringr::str_wrap(plot1$data$name, 25)
  plot1 <- plot1+ggraph::geom_node_text(aes_(label=~name), data = plot1$data[1:top_n_cat,], size = 3, bg.color = "white", repel=TRUE)
  
  
  if(!is.null(gene_color2)){
    plot2 <- plot1
    plot2$data$color[-c(1:top_n_cat)] <- gene_color2[plot1$data$name[-c(1:top_n_cat)]]
    
  }
  if(color_cat_pval){
    plot1 <- plot1+ggnewscale::new_scale_color() +
      ggraph::geom_node_point(aes_(color=~qval, size=~size), data = cat_df) + 
      scale_size(limits = c(1,150), breaks = c(10,20,40,80), range = c(1,5)) +
      scale_colour_gradientn(name = "FDR", na.value = 'black', colours = colorRampPalette(rev(brewer.pal(9,  'Purples')))(255)[0:200], 
                             limits= c(0, 0.1), breaks = c(0,2.5e-2,  5e-2, 7.5e-2, 1e-1),  oob = scales::oob_squish)#+
    #theme(legend.box = "horizontal", legend.position="bottom")
    #theme(plot.margin=unit(c(0,0,0,0),"mm"), aspect.ratio = 1)
    if(!is.null(gene_color2)){
      plot2 <- plot2+ggnewscale::new_scale_color() +
        ggraph::geom_node_point(aes_(color=~qval, size=~size), data = cat_df) + 
        scale_size(limits = c(1,150), breaks = c(10,20,40,80), range = c(1,5),) +
        scale_colour_gradientn(name = "FDR",na.value = 'black', 
                               colours = colorRampPalette(rev(brewer.pal(9,  'Purples')))(255)[0:200], limits= c(0, 0.1), breaks = c(0, 2.5e-2,  5e-2, 7.5e-2, 1e-1), 
                               oob = scales::oob_squish)#+
      # theme(legend.box = "vertical", legend.position="bottom")
    }
  }
  return(list(plot1 = plot1, plot2 = plot2))
}




make_very_custom_DAP_vs_DEG_composite_plot <- function(dap_res, deg_res){
  # dap_res and deg_res are the data frames containing final results of either dap and deg respectively, row.names of both are gene_names
  ## Compare all genes that have DAP analysis D+PDUI change vs Log2FC, no correlation but association
  genes_shared <- intersect(row.names(dap_res), row.names(deg_res))

  utr_v_deg2 <- data.frame(row.names = genes_shared, utr=dap_res[genes_shared, 'mean.diff'], 
                           deg = deg_res[genes_shared, 'Log2FC']
  )
  
  utr_v_deg2$deg_sig <- 'Non'
  utr_v_deg2[intersect(genes_shared, row.names(subset(deg_res, fdr < 0.05 & Log2FC > log2(2)))), 'deg_sig'] <- 'DEG up'
  utr_v_deg2[intersect(genes_shared, row.names(subset(deg_res, fdr < 0.05 & Log2FC < -log2(2)))), 'deg_sig'] <- 'DEG down'
  
  utr_v_deg2$utr_sig <- 'Non'
  utr_v_deg2[intersect(genes_shared, row.names(subset(dap_res, mean.diff < -0.2 & fdr < 0.05))), 'utr_sig'] <- 'Shortened UTR'
  utr_v_deg2[intersect(genes_shared, row.names(subset(dap_res, mean.diff > 0.2 & fdr < 0.05))), 'utr_sig'] <- 'Lengthened UTR'
  
  
  utr_v_deg2$utr_sig2 <- 'Non'
  utr_v_deg2[intersect(genes_shared, row.names(subset(dap_res, abs(mean.diff) > 0.2 & fdr < 0.05))), 'utr_sig2'] <- 'Sig DAP'
  
  
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
    scale_size_manual(values = c('Sig DAP'=1.5, 'Non'=0.5))+
    annotate('text', x=0.5, y=2.4, size = 6, label = paste('PCC: ', round(cor(utr_v_deg2$utr, utr_v_deg2$deg), 3), sep = ''))+
    geom_hline(yintercept = c(1,-1), color = c('red', 'blue'), linetype = 'dashed')
  
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
  return(utr_deg_mix_plot)
}

#'graphopt' 'fr' 'kk' 'drl'  'lgl'