# Load packages
suppressPackageStartupMessages({
  library(betareg)
  library(tidyr)
  library(qvalue)
  library(pbapply)
  library(lmtest)
  library(biomaRt)
  library(BASiCS)
  library(scDD)
  library(BPSC)
  library(ggrepel)
  library(apeglm)
  library(heatmap3)
  library(ggplot2)
  library(ggfortify)
  library(stringr)
  library(RColorBrewer)
  library(MKmisc)
  library(DESeq2)
  library(Rtsne)
  library(MAST)
  library(reticulate)
  library(edgeR)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
  library(pcaMethods)
  library(segmented)
  library(robust)
  library(MASS)
  library(umap)
  library(QoRTs)
  library(tximport)
  library(Seurat)
  require(ReactomePA)
  require(clusterProfiler)
  require(org.Sc.sgd.db)
  require(meshes)
  library(msigdbr)
  library(Seurat)
  library(doParallel)
  library(grid)
  library(gridExtra)
  library(pathview)
  library(igraph)
  library(ggraph)
  library(qvalue)
  library(ggrepel)
  library(DRIMSeq)
  library(ggbeeswarm)
  library(stageR)
  library(DEXSeq)
  library(trackViewer)
  library(GenomicFeatures)
})
reticulate::use_condaenv("scvelo")
scv <- import("scvelo")

FILES <- "/Users/ysu13/My Drive/mouse_rat/star_based/dataset/"
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


read_filter_splice <- function(splice_f, unsplic_f) {
  spliced <- t(read.csv(splice_f, header = T, row.names = 1))
  spliced_raw <- t(read.csv(splice_f, header = T, row.names = 1))
  unsplic <- t(read.csv(unsplic_f, header = T, row.names = 1))
  unsplic_raw <- t(read.csv(unsplic_f, header = T, row.names = 1))
  grp1 <- grepl(pattern = "F", colnames(spliced))
  grp2 <- grepl(pattern = "U", colnames(spliced))
  grp1_filt_s <- rowMeans(spliced[, grp1]) > 5 & rowSums(spliced[, grp1] > 0) > sum(grp1) * 0.1
  grp1_filt_u <- rowMeans(unsplic[, grp1]) > 1 & rowSums(unsplic[, grp1] > 0) > sum(grp1) * 0.1

  grp2_filt_s <- rowMeans(spliced[, grp2]) > 5 & rowSums(spliced[, grp2] > 0) > sum(grp2) * 0.1
  grp2_filt_u <- rowMeans(unsplic[, grp2]) > 1 & rowSums(unsplic[, grp2] > 0) > sum(grp2) * 0.1

  spliced <- spliced[(grp1_filt_s & grp1_filt_u) | (grp2_filt_s & grp2_filt_u), ]
  unsplic <- unsplic[(grp1_filt_s & grp1_filt_u) | (grp2_filt_s & grp2_filt_u), ]


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
  return(list(spliced = spliced_raw, unspliced = unsplic_raw, metrics = test_cor, genes = genes))
}


read_filter_splice_loom <- function(loomf, samples = NULL) {
  loom <- scv$read_loom(loomf)
  spliced <- t(as.matrix(loom$layers["spliced"]))
  colnames(spliced) <- sapply(strsplit(sapply(strsplit(loom$obs_names$values, ":"), function(x) x[2]), "_"), function(x) x[1])
  row.names(spliced) <- loom$var_names$values
  unsplic <- t(as.matrix(loom$layers["unspliced"]) + as.matrix(loom$layers["spanning"]))
  colnames(unsplic) <- colnames(spliced)
  row.names(unsplic) <- loom$var_names$values
  if (!is.null(samples)) {
    unsplic <- unsplic[, samples]
    spliced <- spliced[, samples]
  }
  spliced_raw <- spliced
  unsplic_raw <- unsplic

  grp1 <- which(grepl(pattern = "F", colnames(spliced)))
  grp2 <- which(grepl(pattern = "U", colnames(spliced)))

  grp1 <- grp1[colMeans(cor(unsplic_raw[, grp1], unsplic_raw[, grp1])) > 0.7]
  grp2 <- grp2[colMeans(cor(unsplic_raw[, grp2], unsplic_raw[, grp2])) > 0.7]
  spliced_raw <- spliced_raw[, c(grp1, grp2)]
  unsplic_raw <- unsplic_raw[, c(grp1, grp2)]


  grp1_filt_s <- rowMeans(spliced[, grp1]) > 5 & rowSums(spliced[, grp1] > 0) > length(grp1) * 0.5
  grp1_filt_u <- rowMeans(unsplic[, grp1]) > 1 & rowSums(unsplic[, grp1] > 0) > length(grp1) * 0.5

  grp2_filt_s <- rowMeans(spliced[, grp2]) > 5 & rowSums(spliced[, grp2] > 0) > length(grp2) * 0.5
  grp2_filt_u <- rowMeans(unsplic[, grp2]) > 1 & rowSums(unsplic[, grp2] > 0) > length(grp2) * 0.5

  print(dim(unsplic))
  print(dim(spliced))
  spliced <- spliced[(grp1_filt_s & grp1_filt_u) | (grp2_filt_s & grp2_filt_u), ]
  unsplic <- unsplic[(grp1_filt_s & grp1_filt_u) | (grp2_filt_s & grp2_filt_u), ]

  unsplic <- unsplic[, c(grp1, grp2)]
  spliced <- spliced[, c(grp1, grp2)]
  print(dim(unsplic))
  print(dim(spliced))

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
  return(list(spliced = spliced_raw, unspliced = unsplic_raw, metrics = test_cor, genes = genes))
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



prepareMeta <- function(cds, dirc, filter_conditions = vector(), filter_samples = c(), QC = F) {
  dot_color <- DOT_COLOR
  stat <- read.table(paste(FILES, dirc, "meta.tsv", sep = "/"), header = TRUE, row.names = 1)
  cells <- row.names(stat)
  ct <- get_qorts_summary(paste(FILES, dirc, "/QCData/hisat/", sep = "/"), file_name = "QC.geneCounts.formatted.for.DESeq.txt.gz")
  ct <- ct[1:(nrow(ct) - 5), ]
  spike.reads <- ct[which(grepl("ERCC", row.names(ct))), ]
  bio.reads <- ct[which(!grepl("ERCC", row.names(ct))), ]
  spike.rate <- round(colSums(spike.reads) / colSums(ct), 5)
  names(spike.rate) <- colnames(ct)
  stat$spike.rate <- spike.rate[row.names(stat)]
  eff.reads.ratio <- (colSums(ct) - colSums(spike.reads)) / (stat[colnames(ct), ]$input.read.pair.count * (1 - stat[colnames(ct), ]$spike.rate))
  names(eff.reads.ratio) <- colnames(ct)
  stat$eff.reads.ratio <- eff.reads.ratio[row.names(stat)]
  stat$num_bio_reads <- colSums(bio.reads)[row.names(stat)]
  stat$sensitivity <- colSums(bio.reads > 0)[row.names(stat)]
  meta <- data.frame(stat[cells, ])
  remove_samples <- union(filter_samples, row.names(meta[meta$condition %in% filter_conditions, ]))
  meta <- meta[!row.names(meta) %in% remove_samples, ]
  meta$color <- dot_color[meta$cellType]
  cds$bio.reads <- bio.reads
  cds$meta <- meta
  if (QC) {
    cds$meta <- subset(cds$meta, sensitivity > 500 & complexity > 0.01 & gap < 1.5 & eff.reads.ratio > 0.05)
    cds$meta <- subset(cds$meta, intraDiffinterRank > 1)
  }
  return(cds)
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
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.margin = margin(6, 10, 6, 6),
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
        coord_fixed() +
        theme0 +
        ggtitle(main)
    }
    if (color_by == "cellType") {
      plot0 <- plot0 + scale_color_manual(values = colors_used, labels = legend.label)
    } else {
      plot0 <- plot0
    }
  } else {
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

makePCA_UMAP <- function(cds, cellTypes, counts = F, output_dir = "./", color_by = "cellType") {
  umap_config <- list(n_neighbors = 15, min_dist = 0.5, metric = "pearson")
  samples <- row.names(cds$meta[cds$meta$cellType %in% cellTypes, ])
  meta <- cds$meta[samples, ]
  mat <- cds$tpm$bio
  # mat <- mat[rowSums(mat > 0) > ncol(mat)*0.1,]
  # mat <- mat[rowMeans(mat) >= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.025) & rowMeans(mat) <= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.975), ]
  mat <- log(mat + 1)
  if (counts) {
    ct <- cds$ct$bio
    # ct <- cds$ct$bio[setdiff(row.names(cds$ct$bio), c('Gm26917', 'Lars2')),]
    mat <- DESeq2::varianceStabilizingTransformation(as.matrix(ct))
    # mat <- log(t(t(mat[,samples])/edgeR::calcNormFactors(mat[,samples]), )+1)
    # mat <- log(edgeR::cpm(mat)+1)
    # mat <- mat[rowSums(mat > 0) > ncol(mat)*0.1,]
  }
  plts <- list()
  print(dim(mat))
  # mat <- mat[rowMeans(mat) >= quantile(rowMeans(mat)[rowMeans(mat) > 0], 0.025), ]
  sample_PCA(mat, meta, umap = T, dimension = 2, main = "UMAP of Log Expression", labeling = T, point_size = 2, color_by = color_by, umap.config = umap_config)
  ggsave(paste(output_dir, "_uMAP_label.tiff", sep = ""), units = "in", width = 5, height = 4, dpi = 300, compression = "lzw")
  plts[["umap"]] <- sample_PCA(mat, meta, umap = T, dimension = 2, main = "UMAP of Log Expression", labeling = T, point_size = 4, color_by = color_by, umap.config = umap_config)
  ggsave(paste(output_dir, "_uMAP.tiff", sep = ""), units = "in", width = 5, height = 4, dpi = 300, compression = "lzw")
  sample_PCA(mat, meta, dimension = 2, main = "PCA of Log Expression", labeling = T, point_size = 2, color_by = color_by)
  ggsave(paste(output_dir, "_PCA_label.tiff", sep = ""), units = "in", width = 5, height = 4, dpi = 300, compression = "lzw")
  plts[["pca"]] <- sample_PCA(mat, meta, dimension = 2, main = "PCA of Log Expression", point_size = 4, color_by = color_by)
  ggsave(paste(output_dir, "_PCA.tiff", sep = ""), units = "in", width = 5, height = 4, dpi = 300, compression = "lzw")
  return(plts)
}


options(ucscChromosomeNames = FALSE)




plot_utr_coverage <- function(gene, utr_res, txdb, bw_file_list, cols, cell_names, y_margin = -0.05, y_height = 2.5, file_name_suffix = "", loci = NULL) {
  require(GenomicRanges)
  require(trackViewer)
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



# function for importing bigwig files
import_wig <- function(file, file2, format = c(
                         "BED", "bedGraph", "WIG",
                         "BigWig"
                       ), ranges = GRanges(), ignore.strand = TRUE) {
  if (missing(file)) {
    stop("file is required.")
  }
  format <- match.arg(format)
  if (!is(ranges, "GRanges")) {
    stop("ranges must be an object of GRanges.")
  }
  gr <- trackViewer:::orderedGR(ranges)
  seqn <- unique(as.character(seqnames(gr)))
  filterByRange <- function(r) {
    if (length(gr) > 0) {
      r <- r[r[, 1] %in% seqn, , drop = FALSE]
      nr <- nrow(r)
      if (nr > 0) {
        idx <- rep(FALSE, nr)
        l <- floor(nr / 1000)
        for (i in 0:l) {
          f <- min(i * 1000 + 1, nr)
          t <- min((i + 1) * 1000, nr)
          x <- r[f:t, , drop = FALSE]
          xgr <- GRanges(x[, 1], IRanges(start = as.numeric(x[
            ,
            2
          ]), end = as.numeric(x[, 3])))
          suppressWarnings(ol <- findOverlaps(xgr, gr,
            ignore.strand = ignore.strand
          ))
          if (length(ol) > 0) {
            idx[queryHits(ol) + i * 1000] <- TRUE
          }
        }
        r <- r[idx, , drop = FALSE]
      }
    }
    r
  }
  getWigInfo <- function(firstline) {
    firstline <- unlist(strsplit(firstline, "\\s"))
    firstline <- firstline[firstline != ""]
    structure <- firstline[1]
    firstline <- firstline[-1]
    firstline <- do.call(rbind, strsplit(firstline, "=",
      fixed = TRUE
    ))
    firstline <- firstline[match(c(
      "chrom", "span", "start",
      "step"
    ), firstline[, 1]), ]
    info <- c(structure, firstline[, 2])
    names(info) <- c(
      "structure", "chrom", "span", "start",
      "step"
    )
    return(info)
  }
  readWIG <- function(buf, lastWigInfo = NULL) {
    buf <- gsub("^\\s+", "", buf)
    buf <- gsub("\\s+$", "", buf)
    buf <- buf[grepl(
      "^(variableStep|fixedStep|([0-9]+))",
      buf
    )]
    infoLine <- grep("Step", buf)
    if (length(infoLine) > 0) {
      if (infoLine[1] != 1) {
        if (is.null(lastWigInfo[1])) {
          stop("WIG file must contain track definition line, \n                     which should start by variableStep or fixedStep.")
        } else {
          buf <- c(lastWigInfo, buf)
          infoLine <- grep("Step", buf)
        }
      }
    } else {
      if (is.null(lastWigInfo[1])) {
        stop("WIG file must contain track definition line, \n                     which should start by variableStep or fixedStep.")
      } else {
        buf <- c(lastWigInfo, buf)
        infoLine <- grep("Step", buf)
      }
    }
    lastWigInfo <- buf[infoLine[length(infoLine)]]
    while (infoLine[length(infoLine)] == length(buf)) {
      buf <- buf[-length(buf)]
    }
    block <- c(infoLine, length(buf) + 1)
    dif <- diff(block)
    block <- rep(infoLine, dif)
    buf <- split(buf, block)
    r <- lapply(buf, function(.ele) {
      wiginfo <- getWigInfo(.ele[1])
      span <- wiginfo["span"]
      step <- as.numeric(wiginfo["step"])
      if (wiginfo["structure"] == "variableStep") {
        start <- as.numeric(strsplit(.ele[2], "\\s+")[[1]][1])
        if (is.na(span)) {
          span <- 1
        }
        lastrow <- strsplit(.ele[length(.ele)], "\\s+")[[1]]
        end <- as.numeric(lastrow)[1] + as.numeric(span)
      } else {
        start <- as.numeric(wiginfo["start"])
        if (!is.na(span)) {
          end <- start + (length(.ele) - 1) * step +
            as.numeric(span) - 1
        } else {
          end <- start + length(.ele) * step - 1
        }
      }
      c(wiginfo["chrom"], start, end, span, step, wiginfo["structure"])
    })
    r <- do.call(rbind, r)
    wiginfo <- getWigInfo(lastWigInfo)
    if (wiginfo["structure"] == "fixedStep") {
      lastWigInfo <- gsub(
        "start=\\d+(\\s)", paste("start=",
          as.numeric(r[nrow(r), 3]) + 1, "\\1",
          sep = ""
        ),
        lastWigInfo
      )
    }
    buf <- lapply(buf, "[", -1)
    buf <- CharacterList(buf, compress = TRUE)
    if (length(gr) > 0) {
      r <- cbind(r, rid = 1:nrow(r))
      r <- filterByRange(r)
      buf <- buf[as.numeric(r[, "rid"])]
    }
    list(gr = GRanges(seqnames = r[, 1], ranges = IRanges(start = as.numeric(r[
      ,
      2
    ]), end = as.numeric(r[, 3])), score = buf, span = as.numeric(r[
      ,
      4
    ]), step = as.numeric(r[, 5]), structure = r[
      ,
      6
    ]), lastWigInfo = lastWigInfo)
  }
  readBED <- function(buf) {
    buf <- strsplit(buf, "\t", fixed = TRUE)
    len <- sapply(buf, length)
    buf <- buf[len > 2]
    len <- len[len > 2]
    if (length(buf) < 1) {
      return(GRanges(score = numeric(0)))
    }
    maxLen <- max(len)
    if (all(len == maxLen)) {
      buf <- do.call(rbind, buf)
    } else {
      NAs <- rep("", maxLen)
      buf <- do.call(rbind, lapply(buf, function(.ele) {
        c(
          .ele,
          NAs
        )[1:maxLen]
      }))
    }
    if (ncol(buf) == 3) {
      buf <- cbind(buf, ".")
    }
    if (ncol(buf) == 4) {
      if (all(grepl("^[\\d\\.]+$", buf[, 4])) && length(unique(nchar(buf[
        ,
        4
      ]))) > 1) {
        buf <- cbind(buf, buf[, 4])
      } else {
        buf <- cbind(buf, 1)
      }
    }
    if (ncol(buf) == 5) {
      buf <- cbind(buf, "*")
    }
    buf[!buf[, 6] %in% c("+", "-"), 6] <- "*"
    buf <- filterByRange(buf)
    if (nrow(buf) > 0) {
      GRanges(seqnames = buf[, 1], ranges = IRanges(start = as.numeric(buf[
        ,
        2
      ]) + 1, end = as.numeric(buf[, 3])), strand = buf[
        ,
        6
      ], score = as.numeric(buf[, 5]))
    } else {
      GRanges(seqnames = buf[, 1], ranges = IRanges(start = as.numeric(buf[
        ,
        2
      ]), end = as.numeric(buf[, 3])), strand = buf[
        ,
        6
      ], score = as.numeric(buf[, 5]))
    }
  }
  readFourCols <- function(buf) {
    buf <- gsub("^\\s+", "", buf)
    buf <- gsub("\\s+$", "", buf)
    buf <- buf[!grepl("^(browser|track|#)", buf)]
    buf <- strsplit(buf, "\t", fixed = TRUE)
    len <- sapply(buf, length)
    buf <- buf[len == 4]
    if (length(buf) < 1) {
      return(GRanges(score = numeric(0)))
    }
    buf <- do.call(rbind, buf)
    buf <- filterByRange(buf)
    if (nrow(buf) > 0) {
      GRanges(seqnames = buf[, 1], ranges = IRanges(start = as.numeric(buf[
        ,
        2
      ]) + 1, end = as.numeric(buf[, 3])), score = as.numeric(buf[
        ,
        4
      ]))
    } else {
      GRanges(seqnames = buf[, 1], ranges = IRanges(start = as.numeric(buf[
        ,
        2
      ]), end = as.numeric(buf[, 3])), score = as.numeric(buf[
        ,
        4
      ]))
    }
  }
  readbedGraph <- function(buf) {
    readFourCols(buf)
  }
  readBigWig <- function(file) {
    if (length(gr) > 0) {
      import(con = file, format = "BigWig", which = gr)
    } else {
      import(con = file, format = "BigWig")
    }
  }
  readFile <- function(file, format, FUN) {
    if (format == "WIG") {
      res <- NULL
      con <- file(file, open = "r")
      on.exit(close(con))
      lastWigInfo <- NULL
      while (length(buf <- readLines(con, n = 1e+06, warn = FALSE)) >
        0) {
        buf <- FUN(buf, lastWigInfo)
        lastWigInfo <- buf$lastWigInfo
        if (length(res) < 1) {
          res <- buf$gr
        } else {
          suppressWarnings(res <- c(res, buf$gr))
        }
      }
    } else {
      s <- file.info(file)$size
      if (s < 1e+08) {
        buf <- readChar(file, s, useBytes = TRUE)
        buf <- strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]]
        res <- FUN(buf)
      } else {
        message("file is too huge. Please consider to use bedtools or bedops to subset the data.")
        res <- NULL
        con <- file(file, open = "r")
        on.exit(close(con))
        while (length(buf <- readLines(con,
          n = 1e+06,
          warn = FALSE
        )) > 0) {
          buf <- FUN(buf)
          if (length(res) < 1) {
            res <- buf
          } else {
            suppressWarnings(res <- c(res, buf))
          }
        }
      }
    }
    res <- unique(res)
    return(res)
  }
  readFiles <- function(file, format) {
    FUN <- get(paste("read", format, sep = ""))
    if (format == "BigWig") {
      res <- unique(FUN(file))
    } else {
      res <- readFile(file, format, FUN)
    }
    return(res)
  }
  res <- readFiles(file, format)
  if (!missing(file2)) {
    res2 <- readFiles(file2, format)
    return(new("track",
      dat = trackViewer:::orderedGR(res), dat2 = trackViewer:::orderedGR(res2),
      type = "data", format = format
    ))
  } else {
    return(new("track",
      dat = trackViewer:::orderedGR(res), type = "data",
      format = format
    ))
  }
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


  genes <- intersect(row.names(ct)[which(rowSums(ct > 0) > min_cell)], genes)
  filt_genes <- setdiff(row.names(ct), genes)
  ct <- ct[genes, ]
  print(length(filt_genes))

  print(dim(ct))
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

      congenes <- row.names(subset(sub_fc, coef_diff > 0.01 & z))
      fcHurdleSig[, "orig_logfc"] <- fc[row.names(fcHurdleSig), "coef"]
      fcHurdleSig[, "con_logfc"] <- con[row.names(fcHurdleSig), "coef"]
      fcHurdleSig[congenes, c("coef", "ci.hi", "ci.lo")] <- con[congenes, c("coef", "ci.hi", "ci.lo")]
    }

    fcHurdleSig$Log2FC <- fcHurdleSig$coef
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



enrich_CP <- function(ora_genes, organisms, universe = NULL, logFC = NULL, GSE = F, GO_BP_only = F, enrich_all = T) {
  require(ReactomePA)
  require(clusterProfiler)
  require(org.Sc.sgd.db)
  require(meshes)
  require(dplyr)
  if (organisms == "mouse") {
    require(org.Mm.eg.db)
    orgdb <- org.Mm.eg.db
    orgabv <- "mmu"
    orgname <- "Mus musculus"
    n_type <- "ALIAS"
  } else if (organisms == "rat") {
    require(org.Rn.eg.db)
    orgdb <- org.Rn.eg.db
    orgabv <- "rno"
    orgname <- "Rattus norvegicus"
    n_type <- "ALIAS"
  } else if (organisms == "celegans") {
    require(org.Ce.eg.db)
    orgdb <- org.Ce.eg.db
    orgabv <- "cel"
    orgname <- "Caenorhabditis elegans"
    n_type <- "ENSEMBL"
  }

  # oraL <- unique(bitr(ora_genes, 'ALIAS', 'ENTREZID', orgdb)$ENTREZID)
  # universe <-  unique(bitr(universe, 'ALIAS', 'ENTREZID', orgdb)$ENTREZID)
  oraL <- tryCatch(
    {
      egdf <- bitr(ora_genes, n_type, "ENTREZID", orgdb) %>%
        distinct(eval(as.name("ALIAS")), .keep_all = T) %>%
        data.frame(row.names = 1)
      egdf$ENTREZID
    },
    error = function(cond) {
      return(NULL)
    }
  )
  universe <- tryCatch(
    {
      egdf <- bitr(universe, n_type, "ENTREZID", orgdb) %>%
        distinct(eval(as.name("ALIAS")), .keep_all = T) %>%
        data.frame(row.names = 1)
      egdf$ENTREZID
    },
    error = function(cond) {
      return(NULL)
    }
  )

  if (is.null(oraL) || length(oraL) == 0) {
    return(NULL)
  }

  if (organisms == "celegans") {
    oraL_kegg <- paste("CELE_", WORM.GENES[ora_genes, ]$Sequence.Name, sep = "")
  } else {
    oraL_kegg <- oraL
  }

  GSE_results <- list()
  # GO Enrichment
  if (GO_BP_only | enrich_all) {
    GSE_results[["GO_BP_ora"]] <- tryCatch(
      {
        setReadable(enrichGO(
          gene = oraL,
          universe = universe,
          OrgDb = orgdb,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.1, maxGSSize = 500, minGSSize = 10,
          qvalueCutoff = 0.05
        ), OrgDb = orgdb)
      },
      error = function(cond) {
        return(NULL)
      }
    )
  }
  if (!GO_BP_only) {
    GSE_results[["WKP_ora"]] <- tryCatch(
      {
        setReadable(enrichWP(oraL, organism = orgname, maxGSSize = 500, minGSSize = 10, universe = universe, pvalueCutoff = 0.1, qvalueCutoff = 0.05), OrgDb = orgdb)
      },
      error = function(cond) {
        return(NULL)
      }
    )

    GSE_results[["GO_CC_ora"]] <- tryCatch(
      {
        setReadable(enrichGO(
          gene = oraL,
          universe = universe,
          OrgDb = orgdb,
          ont = "CC",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.1,
          qvalueCutoff = 0.05, maxGSSize = 500, minGSSize = 10,
        ), OrgDb = orgdb)
      },
      error = function(cond) {
        return(NULL)
      }
    )

    GSE_results[["GO_MF_ora"]] <- tryCatch(
      {
        setReadable(enrichGO(
          gene = oraL,
          universe = universe,
          OrgDb = orgdb,
          ont = "MF",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.1, maxGSSize = 500, minGSSize = 10,
          qvalueCutoff = 0.05
        ), OrgDb = orgdb)
      },
      error = function(cond) {
        return(NULL)
      }
    )


    # KEGG Enrichment

    GSE_results[["KEGG_ora"]] <- tryCatch(
      {
        setReadable(enrichKEGG(
          gene = oraL, universe = universe,
          organism = orgabv, maxGSSize = 500, minGSSize = 10,
          pvalueCutoff = 0.1, qvalueCutoff = 0.05
        ), OrgDb = orgdb, keyType = "ENTREZID")
      },
      error = function(cond) {
        return(NULL)
      }
    )


    GSE_results[["MKEGG_ora"]] <- tryCatch(
      {
        setReadable(enrichMKEGG(
          gene = oraL, universe = universe, maxGSSize = 500, minGSSize = 10,
          organism = orgabv, pvalueCutoff = 0.1, qvalueCutoff = 0.05
        ), OrgDb = orgdb, keyType = "ENTREZID")
      },
      error = function(cond) {
        return(NULL)
      }
    )

    # Reactome Enrichment
    GSE_results[["REACT_ora"]] <- tryCatch(
      {
        setReadable(enrichPathway(gene = oraL, organism = organisms, maxGSSize = 500, minGSSize = 10, universe = universe, pvalueCutoff = 0.1, qvalueCutoff = 0.05), OrgDb = orgdb, keyType = "ENTREZID")
      },
      error = function(cond) {
        return(NULL)
      }
    )
  }

  # Mesh Enrichment
  # GSE_results[['MESH_ora']] <- setReadable0(enrichMeSH(oraL, MeSHDb = "MeSH.Sce.S288c.eg.db", database="gene2pubmed", category = 'G'), gene2symbol = entrez2symbol, keyType = 'ENTREZ')


  gse_list0 <- NULL
  gse_list <- NULL
  if (!is.null(logFC)) {
    gse_list <- logFC
    gse_list0 <- logFC
    gse_list <- sort(gse_list, T)
    gse_list0 <- sort(gse_list0, T)
    egdf <- bitr(names(gse_list), "ALIAS", "ENTREZID", orgdb) %>%
      distinct(eval(as.name("ALIAS")), .keep_all = T) %>%
      data.frame(row.names = 1)
    gse_list <- gse_list[intersect(names(gse_list), row.names(egdf))]
    names(gse_list) <- egdf[names(gse_list), ]$ENTREZID
  }
  GSE_results[["gfc0"]] <- gse_list0
  GSE_results[["gfc"]] <- gse_list
  return(GSE_results)
}



gse_CP <- function(ora_genes, organisms, logFC = NULL, GSE = F, simplify_go = T, combine = T, simple_combine = T) {
  require(ReactomePA)
  require(clusterProfiler)
  require(org.Sc.sgd.db)
  require(meshes)
  if (organisms == "mouse") {
    require(org.Mm.eg.db)
    orgdb <- org.Mm.eg.db

    orgabv <- "mmu"
    orgname <- "Mus musculus"
    n_type <- "ALIAS"
  } else if (organisms == "rat") {
    require(org.Rn.eg.db)
    orgdb <- org.Rn.eg.db
    orgabv <- "rno"
    orgname <- "Rattus norvegicus"
    n_type <- "ALIAS"
  } else if (organisms == "celegans") {
    require(org.Ce.eg.db)
    orgdb <- org.Ce.eg.db
    orgabv <- "cel"
    orgname <- "Caenorhabditis elegans"
    n_type <- "ENSEMBL"
  }
  require(dplyr)

  GSE_results <- list()
  gse_list0 <- NULL
  gse_list <- NULL
  if (!is.null(logFC)) {
    gse_list <- logFC
    print(length(gse_list))
    gse_list0 <- logFC
    gse_list <- sort(gse_list, T)
    gse_list0 <- sort(gse_list0, T)
    egdf <- bitr(names(gse_list), n_type, "ENTREZID", orgdb) %>%
      distinct(eval(as.name("ALIAS")), .keep_all = T) %>%
      data.frame(row.names = 1)
    gse_list <- gse_list[intersect(names(gse_list), row.names(egdf))]
    names(gse_list) <- egdf[names(gse_list), ]$ENTREZID
    # names(gse_list) <- sapply(names(gse_list), function(x){ifelse(x %in% genes_table$Gene.symbol, unique(row.names(genes_table[genes_table$Gene.symbol == x,])), x)})
    # names(gse_list0) <-  paste('CELE_', WORM.GENES[names(gse_list0),]$Sequence.Name, sep ='')

    if (GSE) {
      print(length(gse_list))
      GSE_results[["WKP_gse"]] <- setReadable(gseWP(gse_list, eps = 0, organism = orgname, minGSSize = 10, nPermSimple = 100000, maxGSSize = 250, pvalueCutoff = 0.5), OrgDb = orgdb, keyType = "ENTREZID")
      # GSE_results[["BIO_gse"]] <- setReadable0(GSEA(gse_list, TERM2GENE = BIOCYC[,c(1,2)], TERM2NAME = BIOCYC[,c(1,3)], nPerm = 1000, minGSSize = 5, maxGSSize = 500), gene2symbol = entrez2symbol, keyType = 'ENTREZ')
      GSE_results[["GO_BP_gse"]] <- setReadable(gseGO(
        geneList = gse_list,
        OrgDb = orgdb,
        keyType = "ENTREZID", nPermSimple = 100000,
        ont = "BP", eps = 0,
        minGSSize = 10,
        maxGSSize = 250,
        pvalueCutoff = 0.5,
        verbose = FALSE, by = "fgsea"
      ), OrgDb = orgdb, keyType = "ENTREZID")

      GSE_results[["GO_CC_gse"]] <- setReadable(gseGO(
        geneList = gse_list,
        OrgDb = orgdb,
        keyType = "ENTREZID",
        ont = "CC", nPermSimple = 100000,
        minGSSize = 10, eps = 0,
        maxGSSize = 250,
        pvalueCutoff = 0.5,
        verbose = FALSE
      ), OrgDb = orgdb, keyType = "ENTREZID")

      GSE_results[["GO_MF_gse"]] <- setReadable(gseGO(
        geneList = gse_list,
        OrgDb = orgdb, keyType = "ENTREZID",
        ont = "MF",
        minGSSize = 10, eps = 0,
        maxGSSize = 250, nPermSimple = 100000,
        pvalueCutoff = 0.5,
        verbose = FALSE
      ), OrgDb = orgdb, keyType = "ENTREZID")
      GSE_results[["KEGG_gse"]] <- setReadable(gseKEGG(
        geneList = gse_list,
        organism = orgabv,
        minGSSize = 10, eps = 0,
        maxGSSize = 250, nPermSimple = 100000,
        pvalueCutoff = 0.5,
        verbose = FALSE, keyType = "ncbi-geneid"
      ), OrgDb = orgdb, keyType = "ENTREZID")
      GSE_results[["MKEGG_gse"]] <- setReadable(gseMKEGG(gene = gse_list, minGSSize = 10, eps = 0, nPermSimple = 100000, maxGSSize = 250, organism = orgabv, pvalueCutoff = 0.5), OrgDb = orgdb, keyType = "ENTREZID")
      GSE_results[["REACT_gse"]] <- setReadable(gsePathway(
        geneList = gse_list, organism = organisms, nPermSimple = 100000, minGSSize = 10, maxGSSize = 250,
        pvalueCutoff = 0.5, eps = 0,
        pAdjustMethod = "BH", verbose = FALSE
      ), OrgDb = orgdb, keyType = "ENTREZID")
      # GSE_results[['MESH_gse']] <- setReadable0(gseMeSH(gse_list, MeSHDb = "MeSH.Sce.S288c.eg.db", database = "gene2pubmed", category = "G", nPerm = 1000, minGSSize = 5, maxGSSize = 500), gene2symbol = entrez2symbol, keyType = 'ENTREZ')
    }
    for (r in names(GSE_results)) {
      GSE_results[[r]]@result <- subset(GSE_results[[r]]@result, qvalue < 0.05)
    }
    if (combine) {
      if (simple_combine == T) {
        combined <- GSE_results[["GO_BP_gse"]]
        res_ <- c("WKP_gse", "GO_BP_gse", "KEGG_gse", "MKEGG_gse", "REACT_gse", "GO_CC_gse", "GO_MF_gse")
        combined@result <- do.call(rbind, lapply(res_, FUN = function(x) {
          GSE_results[[x]]@result
        }))
        row.names(combined@result) <- combined@result$ID
        combined@geneSets <- do.call(c, lapply(res_, FUN = function(x) {
          GSE_results[[x]]@geneSets
        }))
        names(combined@geneSets) <- do.call(c, lapply(res_, FUN = function(x) {
          names(GSE_results[[x]]@geneSets)
        }))

        GSE_results[["combined"]] <- combined
        GSE_results[["combined_up"]] <- combined
        GSE_results[["combined_up"]]@result <- subset(GSE_results[["combined_up"]]@result, NES > 0)
        GSE_results[["combined_down"]] <- combined
        GSE_results[["combined_down"]]@result <- subset(GSE_results[["combined_down"]]@result, NES < 0)
      } else {
        react <- as.list(ReactomePA:::get_Reactome_DATA(organisms))
        # print(length(react$PATHID2EXTID))
        go_bp <- as.list(clusterProfiler:::get_GO_data(orgdb, "BP", "ENTREZID"))
        # print(length(go_bp$PATHID2EXTID))
        # print(length(react$PATHID2EXTID))
        kegg <- as.list(clusterProfiler:::prepare_KEGG(orgabv, "KEGG", "ncbi-geneid"))
        # print(length(kegg$PATHID2EXTID))
        # print(length(go_bp$PATHID2EXTID))
        # print(length(react$PATHID2EXTID))
        go_kegg_react_list <- c(react$PATHID2EXTID, go_bp$PATHID2EXTID, kegg$PATHID2EXTID)
        # print(length(go_kegg_react_list))
        go_kegg_react_names <- c(react$PATHID2NAME, go_bp$PATHID2NAME, kegg$PATHID2NAME)
        # print(print(length(go_kegg_react_names)))
        go_kegg_react_P2G <- data.frame(do.call(rbind, lapply(names(go_kegg_react_list), FUN = function(x) {
          cbind(rep(x, length(go_kegg_react_list[[x]])), go_kegg_react_list[[x]])
        })))
        # print(dim(go_kegg_react_P2G))
        colnames(go_kegg_react_P2G) <- c("term", "gene")
        print(length(unique(go_kegg_react_P2G$term)))
        go_kegg_react_P2N <- data.frame(term = names(go_kegg_react_names), name = go_kegg_react_names)
        View(go_kegg_react_P2N)
        GSE_results[["combined"]] <- setReadable(GSEA(
          gse_list,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 250,
          eps = 0, nPermSimple = 100000,
          pvalueCutoff = 0.5,
          pAdjustMethod = "BH",
          gson = NULL,
          TERM2GENE = go_kegg_react_P2G,
          TERM2NAME = go_kegg_react_P2N,
          verbose = TRUE,
          seed = FALSE,
          by = "fgsea",
        ), OrgDb = orgdb, keyType = "ENTREZID")
      }
    }
  }



  GSE_results[["gfc0"]] <- gse_list0
  GSE_results[["gfc"]] <- gse_list

  return(GSE_results)
}

deg_utr2 <- function(file, ct, compare, meta, impute = F, method = "fisher.test", filter_by_PAS_motif = TRUE, alpha = 0.05, combine_p = NULL, diff_thresh = 0.05, fasta_file = "../mouse_rat_proposal//dataset/igv/mouse/mouse_spike.fa") {
  ## read in the new dapars file that includes long. short and PDUI values
  dapars <- read.csv(file, sep = "\t", header = T, row.names = 1)
  # dapars <- dapars[,!grepl(regex('^X'), colnames(dapars))]
  # dapars <- cbind(dapars[,1:3],dapars[,grepl(regex(paste(utr_samples, collapse = '|')), colnames(dapars))])
  # colnames(dapars) <- str_remove(colnames(dapars), '^X')
  samples <- row.names(subset(meta, cellType %in% compare))
  dapars <- cbind(dapars[, c(1, 2, 3)], dapars[, grepl(regex(paste(samples, collapse = "|")), colnames(dapars))])

  ct <- ct[, samples]
  ## Filter first based on coverage of the each gene's entire body
  ## Only account for UTR coverage in genes that are assigned at have average expression of at least 5 uniquely mapped reads across all samples
  genes <- row.names(ct)[rowMeans(ct) > 5]
  # genes <- intersect(genes, row.names(ct)[rowSums(ct[,row.names(subset(meta, cellType == compare[2]))] > 10) > 0])
  dapars$gene_short_names <- sapply(row.names(dapars), FUN = function(x) {
    strsplit(x, "\\|")[[1]][2]
  })
  dapars_orig <- dapars
  all_genes <- unique(dapars$gene_short_names)
  cat("Total genes assessed: ", length(unique((dapars$gene_short_names))), "genes\n")
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
  # select filtering genes based on number of passes (Non NAs) and overall coverage (average > 2 either in long or short UTR in both conditions)
  grp <- meta[colnames(d_pdui), "cellType"]
  btch <- meta[colnames(d_pdui), "experiment"]
  cond1_ind <- which(grp == compare[1])
  cond2_ind <- which(grp == compare[2])
  # print(d_pdui[dapars$gene_short_name == 'Cdk1',])
  ## NA FILTER
  # onyly keep genes that have at least than 3 non NAs in terms of coverage in both conditions
  na.filt.genes <- rowSums(!is.na(d_pdui[, cond1_ind])) >= 3 & rowSums(!is.na(d_pdui[, cond2_ind])) >= 3

  dapars <- dapars[na.filt.genes, ]

  cat("filtering based on number of non NAs in each condition ( > 3 in each condition): left with ", length(unique((dapars$gene_short_names))), "genes\n")


  d_long <- d_long[na.filt.genes, ]
  d_short <- d_short[na.filt.genes, ]
  d_pdui <- d_pdui[na.filt.genes, ]

  # Filter based on coverage of UTR regions
  c1.filt.genes <- rowMeans(d_long[, cond1_ind], na.rm = T) > 1 | rowMeans(d_short[, cond1_ind], na.rm = T) > 1
  c2.filt.genes <- rowMeans(d_long[, cond2_ind], na.rm = T) > 1 | rowMeans(d_short[, cond2_ind], na.rm = T) > 1

  final.filt.genes <- c1.filt.genes & c2.filt.genes

  # filtering matrices with genes selected prior
  dapars <- dapars[final.filt.genes, ]
  d_long <- d_long[final.filt.genes, ] %>% mutate(across(everything(), replace_na, 0))
  d_short <- d_short[final.filt.genes, ] %>% mutate(across(everything(), replace_na, 0))
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
    d_long1[is.na(d_long1)] <- 0
    d_short1 <- d_short
    d_short1[is.na(d_short1)] <- 0
    utrl1_mean <- round(rowMeans(d_long1[, cond1_ind], na.rm = T))
    utrl2_mean <- round(rowMeans(d_long1[, cond2_ind], na.rm = T))
    utrs1_mean <- round(rowMeans(d_short1[, cond1_ind], na.rm = T))
    utrs2_mean <- round(rowMeans(d_short1[, cond2_ind], na.rm = T))
    test <- do.call(rbind, pblapply(row.names(d_long), FUN = function(x) {
      # utr_l1 <- d_long[x,cond1_ind][!is.na(d_long[x,cond1_ind])] # long utr coverage in condition 1
      # utr_l2 <- d_long[x,cond2_ind][!is.na(d_long[x,cond2_ind])] # short utr coverage in condition 2
      # utr_s1 <- d_short[x,cond1_ind][!is.na(d_short[x,cond1_ind])] # long utr coverage in condition 1
      # utr_s2 <- d_short[x,cond2_ind][!is.na(d_short[x,cond2_ind])] # short utr coverage in condition 2
      pdui_1 <- mean(d_pdui[x, cond1_ind][!is.na(d_pdui[x, cond1_ind])])
      pdui_2 <- mean(d_pdui[x, cond2_ind][!is.na(d_pdui[x, cond2_ind])])
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
  test$strand <- sapply(strsplit(row.names(test), "\\|"), FUN = function(x) {
    x[4]
  })
  test$APA_dist <- 0
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
  # test <- subset(test, APA_dist >= 75)
  if (filter_by_PAS_motif) {
    test <- post_dapars_pas_filter(test, fasta_file, up_range = 80, down_range = 120, offset = 0)$PAS_motif
    test <- subset(test, num_motif > 0)
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

  return(list(deg = test, long = d_long, gene_res = gene_res, short = d_short, df = dapars_orig, pdui = pdui, pdui_imp = pdui_impute, gene_universe = all_genes))
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
    c(fisher.test(x = rbind(c(mean(nascent_c[i, ]), mean(mature_c[i, ])), c(mean(nascent_t[i, ]), mean(mature_t[i, ]))))$p.value, diff)
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
    theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust = 0.5, size = 8), axis.title = element_text(size = 12), title = element_text(size = 12))

  if (!isProportion) {
    p <- p + ylab("log(Expression)")
  } else {
    p <- p + ylab("Proportions")
  }
  p
}

find_key_gs <- function(res, keys = NULL, key_length = 5, alpha = 0.1) {
  if (is.null(keys)) {
    return(NULL)
  }
  if (length(key_length) != length(keys)) {
    key_length <- rep(key_length[1], length(keys))
  }
  res <- subset(res, qvalue < alpha)
  gs <- unique(do.call(c, lapply(1:length(keys), function(x) {
    ids <- res$ID[grepl(regex(keys[x]), res$Description, ignore.case = T)]
    ids <- ids[1:min(length(ids), key_length[x])]
  })))
  return(gs)
}


addSmallLegend <- function(myPlot, pointSize = 1.5, textSize = 7, spaceLegend = 0.3) {
  myPlot +
    guides(
      shape = guide_legend(override.aes = list(size = pointSize)),
      color = guide_legend(override.aes = list(size = pointSize))
    ) +
    theme(
      legend.title = element_text(size = textSize),
      legend.text = element_text(size = textSize),
      legend.key.size = unit(spaceLegend, "lines")
    )
}

cp_tree_ridge_plot <- function(res, n_cat = 50, nclust = 8, alpha = 0.05, geneSet = NULL) {
  require(ggtree)
  require(aplot)
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
    color = "p.adjust", offset_tiplab = 0.8, fontsize = 4
  ) +
    geom_tiplab(
      offset = 0.8, hjust = 0,
      show.legend = FALSE,
      align = TRUE, size = 3.5, lineheight = 0.75
    ) + xlim(c(0, 20)) + scale_size(
      name = "number of genes",
      range = c(1, 4)
    ))
  tree <- res_tree
  res_tree$layers[c(7, 8)] <- NULL
  res_tree$layers[c(3, 4)] <- NULL
  res_ridge <- addSmallLegend(ridgeplot(res, showCategory = min(n_cat, nrow(res@result))) + scale_fill_viridis_c() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + xlim(c(-4, 4)) + geom_vline(xintercept = 0, linetype = "dashed", color = "red") + xlab(expression("log"[2] * "FC")), pointSize = 3, textSize = 10, spaceLegend = 0.9)

  res_ridge$data$label <- res_ridge$data$category
  # return(list(ridge = res_ridge, tree=res_tree))
  ridge_tree <- res_ridge %>% insert_left(res_tree, width = 3)
  print(ridge_tree)
  return(list(tree = res_tree, ridge = res_ridge, both = ridge_tree))
  # res_ridge  %>% insert_left(res_tree, width = 2)
}
