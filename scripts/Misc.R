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

