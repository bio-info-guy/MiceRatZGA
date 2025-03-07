# corrected dotplot function for clusterProfiler
test_dotplot <- function (object, x = "geneRatio", color = "p.adjust", showCategory = 10, 
                       size = NULL, split = NULL, font.size = 12, title = "", orderBy = "x", 
                       label_format = 30, decreasing = TRUE) 
{
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size)) 
      size <- "Count"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size)) 
      size <- "GeneRatio"
  }
  else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size)) 
      size <- "Count"
  }
  else {
    if (is.null(size)) 
      size <- "Count"
  }
  if (inherits(object, c("enrichResultList", "gseaResultList"))) {
    ldf <- lapply(object, test_fortify, showCategory = showCategory, 
                  split = split)
    df <- dplyr::bind_rows(ldf, .id = "category")
    df$category <- factor(df$category, levels = names(object))
  }
  else {
    df <- test_fortify(object, showCategory = showCategory, split = split)
  }
  
  if (orderBy != "x" && !orderBy %in% colnames(df)) {
    message("wrong orderBy parameter; set to default `orderBy = \"x\"`")
    orderBy <- "x"
  }
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text = x)))
  }
  label_func <- enrichplot:::default_labeller(label_format)
  if (is.function(label_format)) {
    label_func <- label_format
  }
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))
  p <- ggplot(df, aes_string(x = x, y = "Description", size = size, 
                             fill = colorBy)) + geom_point() + aes(shape = I(enrichplot:::enrichplot_point_shape)) + 
    enrichplot:::set_enrichplot_color(type = "fill", name = color) + 
    scale_y_discrete(labels = label_func) + ylab(NULL) + 
    ggtitle(title) + enrichplot:::theme_dose(font.size) + scale_size(range = c(3, 
                                                                  8))
  class(p) <- c("enrichplotDot", class(p))
  return(p)
}




# Corrected version of fortify function from clusterProfiler
test_fortify <- function(model, data, showCategory=5, by = "Count",
                             order=FALSE, drop=FALSE, split=NULL, ...) {
  res <- as.data.frame(model@result)

  res <- res[!is.na(res$Description), ]
  if (inherits(model, "gseaResult")) {
    res$Count <- enrichplot:::str_count(res$core_enrichment, "/")
    res$.sign <- "activated"
    res$.sign[res$NES < 0] <- "suppressed"
  }
  if (drop) {
    res <- res[res$Count != 0, ]
  }
  if (inherits(model, "gseaResult")) {
    res$GeneRatio <- res$Count / res$setSize
  } else if (inherits(model, "enrichResult")) {
    res$GeneRatio <- enrichplot:::parse_ratio(res$GeneRatio)
    if ("BgRatio" %in% colnames(res)) {
      ## groupGO output doesn't have this column
      res$BgRatio <- enrichplot:::parse_ratio(res$BgRatio)
    }
  }
  
  if (order) {
    if (by == "Count") {
      idx <- order(res$Count, decreasing=TRUE)
    } else {
      idx <- order(res$GeneRatio, decreasing=TRUE)
    }
    res <- res[idx,]
  }
  
  topN <- function(res, showCategory) {
    if ( is.numeric(showCategory) ) {
      if ( showCategory <= nrow(res) ) {
        res <- res[1:showCategory,]
      }
    } else { ## selected categories
      res <- res[res$Description %in% showCategory,]
    }
    return(res)
  }
  
  if (is.null(split)) {
    res <- topN(res, showCategory)
  } else {
    lres <- split(res, as.character(res[, split]))
    lres <- lapply(lres, topN, showCategory = showCategory)
    res <- do.call('rbind', lres)
  }
  
  res$Description <- factor(res$Description,
                            levels=rev(unique(res$Description)))
  
  return(res)
}

#Modified importScore function from trackViewer for importing bigwig files on Windows Machines??
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


# helper function for changing ggplot legend element sizes 
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



# modified version of addGuideLine function from trackViewer
addGuideLine2 <- function(guideLine, col = "gray", lty = "dashed", lwd = 1, 
                          vp = NULL, y_lim =c(0,1)) 
{
  if (missing(guideLine) | !(inherits(guideLine, c("numeric", 
                                                   "integer"))) | length(guideLine) < 1) 
    stop("guideLine is required as a numeric vector of coordinates of genome")
  len <- length(guideLine)
  trimLen <- function(obj, len) {
    if (length(obj) < len) 
      obj <- rep(obj, len)[1:len]
    obj
  }
  vpmultiple <- FALSE
  if (length(vp) > 0) {
    stopifnot(is(vp, "viewport"))
    if (is(vp, "vpTree")) {
      vpmultiple <- TRUE
    }
  }
  selectVP <- function(x, tree) {
    xscales <- sapply(tree$children, function(.ele) {
      seekViewport(names(.ele$children))
      current.viewport()$xscale
    }, simplify = FALSE)
    xscales <- do.call(cbind, xscales)
    i <- which(x >= xscales[1, ] & x <= xscales[2, ])
    if (length(i) < 1) {
      message(x, " out of the range.")
      return(NULL)
    }
    seekViewport(paste0("panel.", i))
    current.viewport()
  }
  col <- trimLen(col, len)
  lty <- trimLen(lty, len)
  lwd <- trimLen(lwd, len)
  for (i in seq_along(guideLine)) {
    if (vpmultiple) {
      currentVP <- selectVP(guideLine[i], vp)
      grid.lines(x = guideLine[i], y =  y_lim , gp = gpar(col = col[i], 
                                                          lty = 'dashed', lwd = lwd[i], alpha = 0.8), default.units = "native")
    }
    else {
      currentVP <- vp
      grid.lines(x = guideLine[i], y =  y_lim , gp = gpar(col = col[i], 
                                                          lty = 'dashed', lwd = lwd[i], alpha = 0.8), default.units = "native", 
                 vp = currentVP)
    }
  }
  i <- 1
  while (i < length(guideLine)) {
    if (vpmultiple) {
      currentVP <- selectVP(guideLine[i], vp)
      grid.lines(x = c(guideLine[i], guideLine[i+1]), y =  y_lim[1] , gp = gpar(col = col[i], 
                                                                                lty = lty[i], lwd = lwd[i], alpha = 0.8), default.units = "native", vp = currentVP)
      grid.lines(x = c(guideLine[i], guideLine[i+1]), y =  y_lim[2] , gp = gpar(col = col[i], 
                                                                                lty = lty[i], lwd = lwd[i], alpha = 0.8), default.units = "native", vp = currentVP)
    }
    else {
      currentVP <- vp
      grid.lines(x = c(guideLine[i], guideLine[i+1]), y =  y_lim[1] , gp = gpar(col = col[i], 
                                                                                lty = 'dashed', lwd = lwd[i], alpha = 0.8), default.units = "native", vp = currentVP)
      grid.lines(x = c(guideLine[i], guideLine[i+1]), y =  y_lim[2], gp = gpar(col = col[i], 
                                                                               lty = 'dashed', lwd = lwd[i], alpha = 0.8), default.units = "native", vp = currentVP)
      #grid.lines(x=c(guideLine[1], guideLine[2]), y = 0.03, gp = gpar(col = 'black', 
      #lty = 'solid', lwd = 2, alpha = 1), default.units = "native", arrow = arrow(type = 'closed', length = unit(0.1, "inches")), vp = currentVP)
    }
    i <- i+2
  }
  
  return(invisible())
  
}



extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  
  if (inherits(x, 'list')) {
    geneSets <- x
  } else {
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x@result)
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description        
  }
  #print(geneSets)
  if (is.numeric(n)) {
    print(n)
    return(geneSets)
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    if (inherits(x, 'list')) {
      showCategory <- showCategory[showCategory %in% names(x)]
    }
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (inherits(x, 'list')) {
    nn <- length(x)
  } else {
    nn <- nrow(x@result)
  }
  if (nn < n) {
    n <- nn
  }
  
  return(n)
}

setReadable0 <- function (x, gene2symbol, keyType = "auto") {
  if(is.null(x)){return(x)}
  if (!(is(x, "enrichResult") || is(x, "groupGOResult") || 
        is(x, "gseaResult"))) 
    stop("input should be an 'enrichResult' or 'gseaResult' object...")
  isGSEA <- FALSE
  if (is(x, "gseaResult")) 
    isGSEA <- TRUE
  if (keyType == "auto") {
    keyType <- x@keytype
    if (keyType == "UNKNOWN") {
      stop("can't determine keyType automatically; need to set 'keyType' explicitly...")
    }
  }
  if (x@readable) 
    return(x)
  gc <- geneInCategory(x)
  if (isGSEA) {
    genes <- names(x@geneList)
  }
  else {
    genes <- x@gene
  }
  gn <- gene2symbol
  gc <- lapply(gc, function(i) gn[i])
  res <- x@result
  gc <- gc[as.character(res$ID)]
  geneID <- sapply(gc, paste0, collapse = "/")
  if (isGSEA) {
    res$core_enrichment <- unlist(geneID)
  }
  else {
    res$geneID <- unlist(geneID)
  }
  x@gene2Symbol <- gn
  x@result <- res
  x@keytype <- keyType
  x@readable <- TRUE
  return(x)
}

list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
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


# number of ortho genes
# Function to convert genes from organism1 to orthology_class to organism2, ie. Rat to OrthoClass to Mice; results are interpreted as mice genes
# requires named list of orthology_class to organism 2 genes, currently only allow one gene per orthology class for organism 2
# requires named list of organism 1 genes to orthology_class, multiple organism 1 genes can be in a orthology class
# The current limitation accounts for the following two cases of orthology:
## org_1_genes (one) --> org2_genes (many)
## org_1_genes (one) --> org2_genes (one)
# But does not account for:
## org_1_genes (many) --> org2_genes (many)
## org_1_genes (many) --> org2_genes (one)

ortho_convert <- function(genes, org_to_ortho, ortho_2_org2){
  #no_ortho <- setdiff(genes, names(org_to_ortho))
  #ortho <- intersect(genes, names(org_to_ortho))
  #org2_gene <- ortho_2_org2[as.character(org_to_ortho[ortho])]
  
  #return(list(no_ortho = no_ortho, ortho = ortho, conv_genes = org2_gene))
  res <- ortho_2_org2[as.character(org_to_ortho[genes])]
  res[is.na(res)] <- paste(genes[is.na(res)], '_no_ortho', sep = '')
  res
}

# convert rat gene names to mice ortholog and intersect with list of mice genes
rat_2_mice_intersect <- function(rat_genes, mice_genes){
  conv <- ortho_convert(rat_genes, Rat2HomClass, HomClass2Mouse)
  intsct_genes <- intersect(conv, mice_genes)
}


rat_2_mice_wrapper <- function(rat_genes){
  conv <- ortho_convert(rat_genes, Rat2HomClass, HomClass2Mouse)
}



## Deprecated differential UTR function, now using deg_utr2
# deg_utr <- function(file, ct, compare, meta, impute = F, method = 'fisher.test', combine_p = 'fisher'){
#   ##read in the new dapars file that includes long. short and PDUI values
#   dapars <- read.csv(file, sep = '\t', header = T, row.names = 1)
#   dapars_orig <- dapars
#   dapars$strand = sapply(strsplit(row.names(dapars), '\\|'), FUN = function(x){x[4]})
#   dapars$APA_dist = 0
#   dapars[dapars$strand == '+',]$APA_dist <- abs(sapply(strsplit(dapars[dapars$strand == '+',]$Loci, '-'), 
#                                                        FUN = function(x){as.numeric(strsplit(x[1], ':')[[1]][2])}) - dapars[dapars$strand == '+',]$Predicted_Proximal_APA)-1
#   dapars[dapars$strand == '-',]$APA_dist <- abs(sapply(strsplit(dapars[dapars$strand == '-',]$Loci, '-'), 
#                                                        FUN = function(x){as.numeric(x[2])}) - dapars[dapars$strand == '-',]$Predicted_Proximal_APA)-1
#   
#   ## Filter first based on coverage of the each gene's entire body
#   ## Only account for UTR coverage in genes that are assigned at least 10 uniquely mapped reads
#   
#   genes <- row.names(ct)[rowSums(ct[,row.names(subset(meta, cellType == compare[1]))] > 10) > 5]
#   genes <- intersect(genes, row.names(ct)[rowSums(ct[,row.names(subset(meta, cellType == compare[2]))] > 10) > 5])
#   dapars$gene_short_names <- sapply(row.names(dapars), FUN = function(x){strsplit(x,"\\|")[[1]][2]})
#   print(length(unique((dapars$gene_short_names))))
#   
#   
#   all_genes <- unique(dapars$gene_short_names)
#   dapars <- dapars[dapars$gene_short_names %in% genes,]
#   #dapars <- subset(dapars, fit_value >= 2) # Maybe filter also based on regression fit value
#   print(length(unique((dapars$gene_short_names))))
#   # vector to store gene name and UTR region correspondence in case gene names is lost with imputation
#   gene2region <- dapars$gene_short_names
#   names(gene2region) <- sapply(row.names(dapars), FUN = function(x){strsplit(x,"\\|")[[1]][1]})
#   
#   #split file into long, short and pdui
#   d_long <- dapars[,grepl('long_exp', colnames(dapars))]
#   d_short <- dapars[,grepl('short_exp', colnames(dapars))]
#   d_pdui <- dapars[,grepl('PDUI', colnames(dapars))]
#   
#   # change column names
#   colnames(d_long) <- sapply(strsplit(colnames(d_long), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
#   colnames(d_short) <- sapply(strsplit(colnames(d_short), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
#   colnames(d_pdui) <- sapply(strsplit(colnames(d_pdui), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
#   
#   # select filtering genes based on number of passes (Non NAs) and overall coverage (average > 2 either in long or short UTR in both conditions)
#   grp <- meta[colnames(d_pdui), "cellType"]
#   btch <- meta[colnames(d_pdui), "experiment"]
#   cond1_ind <- which(grp == compare[1])
#   cond2_ind <- which(grp == compare[2])
#   
#   #print(d_pdui[dapars$gene_short_name == 'Cdk1',])
#   ## NA FILTER
#   # First filter out all genes that have at least 5 non-NAs in terms of coverage in both conditions
#   na.filt.genes <- rowSums(!is.na(d_pdui[,cond1_ind])) >= 5 & rowSums(!is.na(d_pdui[,cond2_ind])) >= 5
#   dapars <- dapars[na.filt.genes, ]
#   
#   print(length(unique((dapars$gene_short_names))))
#   
#   d_long <- d_long[na.filt.genes, ]
#   d_short <- d_short[na.filt.genes, ]
#   d_pdui <- d_pdui[na.filt.genes, ]
#   # Filter based on coverage of UTR regions
#   c1.filt.genes <- rowMeans(d_long[, cond1_ind], na.rm = T) > 1 | rowMeans(d_short[, cond1_ind], na.rm = T) > 1
#   c2.filt.genes <- rowMeans(d_long[, cond2_ind], na.rm = T) > 1 | rowMeans(d_short[, cond2_ind], na.rm = T) > 1
#   
#   #print(sum(c2.filt.genes))
#   final.filt.genes <- c1.filt.genes & c2.filt.genes 
#   #print(sum(final.filt.genes))
#   # filtering matrices with genes selected prior
#   dapars <- dapars[final.filt.genes, ]
#   d_long <- d_long[final.filt.genes, ]
#   d_short <- d_short[final.filt.genes, ]
#   d_pdui <- d_pdui[final.filt.genes, ]
#   print(length(unique((dapars$gene_short_names))))
#   #print('Cdk1' %in% dapars$gene_short_name)
#   if(impute){
#     dapars_out <- data.frame(Gene = row.names(dapars), dapars[1:3,], d_pdui)
#     write.table(dapars_out, file = './temp.dp.tsv', sep = '\t', quote = F, row.names = F)
#     d_pdui = scDaPars(raw_PDUI_file = './temp.dp.tsv',
#                       out_dir = "apa/scDaPars_result",
#                       filter_gene_thre = 0.2,
#                       filter_cell_thre = 0.1, k= 8)
#     method = 'ks.test'
#   }
#   #dapars_ratio <- dapars_ratio[rowSums(!is.na(dapars_ratio)) > 20 ,]
#   
#   #dapars_ratio[is.na(dapars_ratio)] <- 0
#   
#   #print(sum(rowSums(!is.na(dapars_ratio[,cond1_ind])) >= 10 & rowSums(!is.na(dapars_ratio[,cond2_ind])) >=10))
#   #print(sum(rowSums(!is.na(dapars_ratio)) > 20))
#   #dapars_ratio <- dapars_ratio[rowSums(!is.na(dapars_ratio[,cond1_ind])) >= 7 & rowSums(!is.na(dapars_ratio[,cond2_ind])) >=7,]
#   if(method == 'ks.test'){
#     test <- apply(d_pdui, 1, FUN = function(x){
#       if(sum(x[!is.na(x)]) == 0 | sum(!is.na(x[cond1_ind])) <= 3 | sum(!is.na(x[cond2_ind])) <= 3){
#         c(1, 0)
#       }else{
#         c(ks.test(x[cond1_ind][!is.na(x[cond1_ind])], x[cond2_ind][!is.na(x[cond2_ind])])$p.value, mean(x[cond2_ind][!is.na(x[cond2_ind])]) - mean(x[cond1_ind][!is.na(x[cond1_ind])]))
#       }
#     })
#     test <- t(test)
#   }else{
#     
#     l1 <- length(cond1_ind)
#     l2 <- length(cond2_ind)
#     d_long1 <- d_long
#     #d_long1[is.na(d_long1)] <- 0
#     d_short1 <- d_short
#     #d_short1[is.na(d_short1)] <- 0
#     utrl1_mean <- round(rowMeans(d_long1[, cond1_ind], na.rm = T))
#     utrl2_mean <- round(rowMeans(d_long1[, cond2_ind], na.rm = T))
#     utrs1_mean <- round(rowMeans(d_short1[, cond1_ind], na.rm = T))
#     utrs2_mean <- round(rowMeans(d_short1[, cond2_ind], na.rm = T))
#     pdui1_mean <- rowMeans(d_pdui[, cond1_ind], na.rm = T)
#     pdui2_mean <- rowMeans(d_pdui[, cond2_ind], na.rm = T)
#     test <- do.call(rbind, pblapply(row.names(d_long), FUN = function(x) {
#       # utr_l1 <- d_long[x,cond1_ind][!is.na(d_long[x,cond1_ind])] # long utr coverage in condition 1
#       # utr_l2 <- d_long[x,cond2_ind][!is.na(d_long[x,cond2_ind])] # short utr coverage in condition 2
#       # utr_s1 <- d_short[x,cond1_ind][!is.na(d_short[x,cond1_ind])] # long utr coverage in condition 1
#       # utr_s2 <- d_short[x,cond2_ind][!is.na(d_short[x,cond2_ind])] # short utr coverage in condition 2
#       pdui_1 <- pdui1_mean[x]#mean(d_pdui[x, cond1_ind], na.rm = T)
#       pdui_2 <- pdui2_mean[x]#mean(d_pdui[x, cond2_ind], na.rm = T)
#       # c(fisher.test(x = rbind(c(mean(utr_l1), mean(utr_s1)), c(mean(utr_l2), mean(utr_s2))))$p.value, pdui_2 - pdui_1)
#       twobytwo <- rbind(c(utrl1_mean[x], utrs1_mean[x]), c(utrl2_mean[x], utrs2_mean[x]))
#       c(fisher.test(x = twobytwo)$p.value, pdui_2 - pdui_1)
#       
#       
#       # if(x == 'XM_039113041.1|Pou2f2|NC_051336.1|-'){
#       # print(rbind(c(sum(utr_l1)/l1, sum(utr_s1)/l1), c(sum(utr_l2)/l1, sum(utr_s2)/l2)))
#       # }
#     }))
#     cat('done')
#   }
#   #print(dim(test))
#   #print(dim(d_long))
#   row.names(test) <- row.names(d_long)
#   colnames(test) <- c('pval', 'mean.diff')
#   test <- data.frame(test)
#   test$pval[test$pval > 1] = 1
#   test$padj <- p.adjust(test$pval)
#   test$fdr <- qvalue(test$pval)$qvalue
#   
#   if(impute){
#     test$gene_short_names <- gene2region[row.names(test)]
#   }else{
#     test$gene_short_names <- sapply(row.names(test), FUN = function(x){strsplit(x,"\\|")[[1]][2]})
#   }
#   test$diff <- abs(test$mean.diff) > 0.2 & test$fdr < 0.05
#   test$fit_value <- dapars[row.names(test),]$fit_value
#   test$predicted_p_APA <- dapars_orig[row.names(test),]$Predicted_Proximal_APA
#   test$loci <- dapars_orig[row.names(test),]$Loci
#   test$strand = sapply(strsplit(row.names(test), '\\|'), FUN = function(x){x[4]})
#   test$APA_dist = 0
#   test[test$strand == '+',]$APA_dist <- abs(sapply(strsplit(test[test$strand == '+',]$loci, '-'), 
#                                                    FUN = function(x){as.numeric(strsplit(x[1], ':')[[1]][2])}) - test[test$strand == '+',]$predicted_p_APA)-1
#   test[test$strand == '-',]$APA_dist <- abs(sapply(strsplit(test[test$strand == '-',]$loci, '-'), 
#                                                    FUN = function(x){as.numeric(x[2])}) - test[test$strand == '-',]$predicted_p_APA)-1
#   
#   #test$APA_dist <- abs(sapply(strsplit(test$loci, '-'), FUN = function(x){as.numeric(x[2])}) - test$predicted_p_APA)-1
#   
#   gene_res <- data.frame(do.call(rbind, tapply(row.names(test), test$gene_short_names, function(x){
#     df <- test[x,]
#     min_pval <- min(df[,'pval'])
#     min_pval_ind = which(df[,'pval'] == min_pval)
#     min_pval_ind <- min_pval_ind[which(abs(df[min_pval_ind, 'mean.diff']) == max(abs(df[min_pval_ind, 'mean.diff'])))][1]
#     if(!is.null(combine_p)){
#       df[min_pval_ind,'pval'] <- metapod::combineParallelPValues(as.list(df[,'pval']), method = combine_p)$p.value
#     }
#     return(cbind(df[min_pval_ind,], dapars[x[min_pval_ind],c(1,2,3)]))
#   })))
#   gene_res$padj <- p.adjust(gene_res$pval)
#   gene_res$fdr <- qvalue(gene_res$pval)$qvalue
#   gene_res$diff <- abs(gene_res$mean.diff) > 0.2 & gene_res$fdr < 0.05
#   print(length(unique((test$gene_short_names))))
#   pdui <- dapars_orig[,grepl('PDUI', colnames(dapars_orig))]
#   colnames(pdui) <- colnames(d_long)
#   pdui <- pdui[row.names(d_long),]
#   pdui_impute <- t(apply(pdui, 1, FUN = function(x){x[is.na(x)] = mean(x, na.rm =T); x}))
#   
#   return(list(deg= test, long = d_long, gene_res = gene_res, short = d_short, df = dapars_orig, pdui = pdui, pdui_imp = pdui_impute, gene_universe = all_genes))
# }
