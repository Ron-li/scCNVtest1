#' Create offtarget and ontarget reads
#'
#' This function uses the terminal to find overlap between fragments file and
#' peaks. Bedtools and Tabix need to be installed in order for this
#' preprocessing to occur
#'
#' @param fragments_file path to the fragments file
#' @param peak_path path to the peaks BED file
#' @param blacklist_path path to a blacklist BED file, if it exists
#' @return A list containing the paths to the offtarget and ontarget files
#' @export
getOverlap <- function(fragments_file, peak_path, blacklist_path=NULL){
  runcmd <- T
  if(!file.exists("/usr/bin/bedtools")){
    runcmd <- F
    print(paste("WARNING: bedtools not installed in system.",
                "https://bedtools.readthedocs.io/en/latest/content/installation.html"))
  }
  if(!file.exists("/usr/bin/tabix")){
    runcmd <- F
    print(paste("WARNING: tabix not installed in system.",
                "http://wiki.wubrowse.org/How_to_install_tabix"))
  }

  #removed blacklisted sections
  new_path = substr(fragments_file,1,nchar(fragments_file)-3)
  if(runcmd){
    cmd <- paste("gunzip -k", fragments_file)

    if(file.exists(blacklist_path)){
      cmd <-paste(cmd, "&& bedtools intersect -a",
                  substr(fragments_file,1,nchar(fragments_file)-3),
                  "-b", blacklist_path, "-v -wa > fragments.blacklisted.tsv")
      system(cmd)
      new_path <- "fragments.blacklisted.tsv"
    }

    #create offtarget file
    cmd <- paste("bedtools intersect -a", new_path, peak_path, "-v -wa > fragments.offtarget.tsv",
                 "&& bgzip fragments.offtarget.tsv",
                 "&& tabix -p bed fragments.offtarget.tsv.gz")
    system(cmd)

    #create ontarget file
    cmd <- paste("bedtools intersect -a", new_path, peak_path, "-wa > fragments.ontarget.tsv",
                 "&& bgzip fragments.ontarget.tsv",
                 "&& tabix -p bed fragments.ontarget.tsv.gz")
  }
}

#' Remove unmappable regions from genomic tiles
#'
#' This function removes any unmapped region defined by a .csv
#' from a GenomicRanges tile object
#'
#' @param tiles GRanges tiles object
#' @param unmapped_region path to a .csv defining the unmappable regions
#' @return GRanges object with the unmappable regions removed
#' @export
removeUnmapped <- function(tiles,unmapped_region){
  gr_b <- read.csv(unmapped_region,header = F)
  colnames(gr_b) <- c("chr", "start", "end")
  pairs <- subsetByOverlaps(tiles, makeGRangesFromDataFrame(gr_b), invert = TRUE, ignore.strand = TRUE)
  return(pairs)
}


#' Convert Genomic Tiles to a dataframe
#'
#' This function converts a GRanges object to a dataframe
#' with appropriate row names
#'
#' @param tiles GRanges Genomic Tiles Object
#' @return A dataframe with genomic features as row names
#' @export
tileTodf <- function(tiles){
  df <- as.data.frame(tiles)[,c(1:4)]
  row_names <- do.call(paste, c(df[,c(1:3)], sep="-"))
  rownames(df) <- row_names
  return(df)
}


#' Calculate GC content and append to dataframe
#'
#' This function the GC content of GenomicRanges Tiles and appends it to
#' an existing dataframe. This function uses the Hsapiens genome from
#' BSgenome.Hsapiens.UCSC.hg19.
#'
#' @param df existing dataframe
#' @param tiles GRanges tiles object
#' @param GC_cutoff Maximum GC content to be included
#' @param ref.genome Choose the reference, 'hg19' or 'hg38'
#' @return dataframe containing GC content, count, and frequency.
#' @export
calcGC <- function(df, tiles, GC_cutoff, ref.genome = "hg38"){
  GC_Content <- c()
  N_Count <- c()
  N_Freq <- c()

  if (ref.genome == "hg38") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38 ::Hsapiens
    for (i in 1:length(tiles))
    {
      GC_Content <- append(GC_Content,biovizBase::GCcontent(Hsapiens, tiles[i]))
      N_Count <- append(N_Count,Biostrings::alphabetFrequency(Biostrings::getSeq(Hsapiens,tiles[i])[[1]],as.prob=F)["N"])
      N_Freq <- append(N_Freq,Biostrings::alphabetFrequency(Biostrings::getSeq(Hsapiens,tiles[i])[[1]],as.prob=T)["N"])
    }
  } else if (ref.genome == "hg19") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19 ::Hsapiens
    for (i in 1:length(tiles))
    {
      GC_Content <- append(GC_Content,biovizBase::GCcontent(Hsapiens, tiles[i]))
      N_Count <- append(N_Count,Biostrings::alphabetFrequency(Biostrings::getSeq(Hsapiens,tiles[i])[[1]],as.prob=F)["N"])
      N_Freq <- append(N_Freq,Biostrings::alphabetFrequency(Biostrings::getSeq(Hsapiens,tiles[i])[[1]],as.prob=T)["N"])
    }
  }

  df$GC_Content <- GC_Content
  df$GC_rm_idex <- (df$GC_Content <= GC_cutoff)
  df$N_Count <- N_Count
  df$N_Freq <- N_Freq
  return(df)
}

#' Calculate overlap between feature bin and peaks
#'
#' This function calculates the length of each Genomic feature not occupied by
#' peaks, and checks if the percentage of peaks is above a certain threshold.
#' The lengths and boolean values are appended as columns to the provided dataframe
#'
#' @param df dataframe to contain the lengths without peaks
#' @param tiles GRanges tiles object
#' @param ref.genome Choose the reference, 'hg19' or 'hg38'
#' @param peak_off_cutoff Maximum percentage of peaks in a given feature
#' @return dataframe containing the lengths of features without peaks and column of features above the peak cutoff
#' @export
calcOverlap <- function(df, tiles, peaks, ref.genome = "hg38", peak_occ_cutoff){
  lengths = c()
  starts <-  df[,2]
  ends <- df[,3]
  if (ref.genome == "hg38"){
    genome <- Seqinfo(genome = "hg38")
    peak <- rtracklayer::import(peaks, genome = genome)
  } else if (ref.genome == "hg19") {
    genome <- Seqinfo(genome = "hg19")
    peak <- rtracklayer::import(peaks, genome = genome)
  }

  for (i in 1:length(tiles)){
    pairs <- IRanges::findOverlapPairs(tiles[i], peak, ignore.strand = TRUE)
    peak_width = sum(BiocGenerics::width(IRanges::pintersect(pairs, ignore.strand = TRUE)))
    window_length = strtoi(ends[i]) - strtoi(starts[i]) + 1
    lengths = append(lengths, window_length - peak_width)
  }

  df$length <- lengths
  df$peak_occ_abn <- (((df$width-df$length)/df$width) >= peak_occ_cutoff)
  return(df)
}


#' Calculate read counts in each bin
#'
#' This function calculates the read count of specific cells in a fragments file.
#' These counts are returned as a matrix
#'
#' @param cells .tsv file of cells to include in the calculation
#' @param df dataframe containing GCcontent, peak content and width information
#' @param fragments path to a .tsv.gz fragments file with corresponding .tbi file
#' @param binsize the size of the bin
#' @param tiles GRanges tiles object
#' @return A matrix containing the read counts. Row names are features, column names are cells
#' @export
getBinMatrix <- function(cells, df, fragments, binsize, tiles){
  cell <- utils::read.table(cells, header = TRUE, sep = ",")
  #cell <- cell[cell$cell_id == "Normal", ]

  df$tile_keep <- (!(df$GC_rm_idex | df$peak_occ_abn) & (df$width == binsize))
  tiles_clean <- tiles[df$tile_keep]

  bin_matrix <- Signac::FeatureMatrix(
    fragments = CreateFragmentObject(fragments),
    features = tiles_clean,
    cells = cell[, 1],
    #chunk = 10,
    sep = c('-', '-'),
    process_n = 200,
    verbose = T

  )
  data <- as.matrix(bin_matrix)
  return(data)
}

# Don't export this function
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

# Don't export this function
.wrap <- function(...) {
  file.sep <- .Platform$file.sep
  splitted <- paste(unlist(strsplit(paste(...), split = file.sep)),
                    collapse = paste0(file.sep, " "))
  splitted.pasted <- paste(strwrap(paste(splitted),
                                   exdent = 2), collapse = "\n")
  gsub(paste0(file.sep, " "), file.sep, splitted.pasted)
}

#' Calculate mappability and append to dataframe
#'
#' This function calculate mappability based on the bin size and appends it to
#' an existing dataframe.
#'
#' @param df dataframe containing a feature lengths column
#' @param binsize The size of each bin
#' @param GC.mappa.grange mappability reference file
#' @return dataframe containing mappability.
#' @export
addMappa <- function(df, binsize, GC.mappa.grange){
  if (!.is.wholenumber(binsize/1000)) {
    stop(.wrap("Please provide a bin size which is a multiple of 1000."))
  }
  # GC.mappa.grange <- NULL
  # black.GC.mappa.folder <- CopyhelpeR::getPathHelperFiles(ref.genome)
  # load(file.path(black.GC.mappa.folder, "GC_mappability.rda"))

  MERGEBINNUMBER <- binsize/1000

  custom.bin <- data.frame()
  for (chr in seqlevels(GC.mappa.grange)) {
    selection <- as(seqnames(GC.mappa.grange) == chr, "vector")
    num <- sum(selection)
    ncol <- num %/% MERGEBINNUMBER + 1
    repn <- ncol*MERGEBINNUMBER - num

    start <- c(start(GC.mappa.grange)[selection], rep(NA, repn))
    start.m <- matrix(start, nrow = MERGEBINNUMBER)
    start.bin <- colMins(start.m, na.rm = TRUE)

    end <- c(end(GC.mappa.grange)[selection], rep(NA, repn))
    end.m <- matrix(end, nrow = MERGEBINNUMBER)
    end.bin <- colMaxs(end.m, na.rm = TRUE)

    mappa <- c(GC.mappa.grange$mappability[selection], rep(NA, repn))
    mappa.m <- matrix(mappa, nrow = MERGEBINNUMBER)
    mappability.bin <- colMeans(mappa.m, na.rm = TRUE)

    GC <- c(GC.mappa.grange$GCcontent[selection], rep(NA, repn))
    GC.m <- matrix(GC, nrow = MERGEBINNUMBER)
    GCcontent.bin <- colMeans(GC.m, na.rm = TRUE)

    chr.bin <- cbind(seqnames = paste0("chr", chr),
                     start = start.bin,
                     end = end.bin,
                     mappability = mappability.bin,
                     GC = GCcontent.bin)
    custom.bin <- rbind(custom.bin, chr.bin)
  }
  custom.bin[, 2:3] <-
    apply(custom.bin[, 2:3], c(1, 2), as.integer)
  custom.bin[, 4:5] <-
    apply(custom.bin[, 4:5], c(1, 2), as.numeric)
  custom.bin <- na.omit(custom.bin)

  a <- left_join(df, custom.bin, by = c('seqnames', 'start', "end"))
  rownames(a) <- do.call(paste, c(df[, 1:3], sep = "-"))

  return(a)
}



#' Correct the read counts according to the GC contents and mappability of
#' each bin.
#'
#' @param df dataframe to contain the bin information
#' @param data A matrix containing the read counts. Row names are features, column names are cells
#' @return A matrix containing the corrected read counts. Row names are features, column names are cells.
#' @export
gcmappaCorrection <- function(df, data){
  dtavg <- data.frame(count.avg = apply(X = data, MARGIN = 1, FUN = mean))
  dtavg$GC_Content <- df[rownames(dtavg), "GC_Content"]
  dtavg$GC_Content %<>% round(4)
  dtavg$mappa <- df[rownames(dtavg), "mappability"]
  dtavg$mappa %<>% round(4)

  rough <- loess(count.avg ~ GC_Content, data = dtavg, span = 0.03)
  i <- seq(0, 1, by = 0.0001)
  final <- loess(predict(rough, i) ~ i, span = 0.3)
  normv <- predict(final, dtavg$GC_Content)

  if(any(is.na(normv))){
    missmin <- which(dtavg$GC_Content == min(dtavg$GC_Content))
    missmax <- which(dtavg$GC_Content == max(dtavg$GC_Content))
    normv[missmin] <- min(final$fitted)
    normv[missmax] <- max(final$fitted)
  }

  data_gccorrect <- data/(normv/median(normv, na.rm = TRUE))

  dtavg$gccount.avg <- apply(X = data_gccorrect, MARGIN = 1, FUN = mean)
  rough <- loess(gccount.avg ~ mappa, data = dtavg, span = 0.03)
  i <- seq(0, 1, by = 0.0001)
  final <- loess(predict(rough, i) ~ i, span = 0.3)
  normv <- predict(final, dtavg$mappa)

  if(any(is.na(normv))){
    missmin <- which(dtavg$mappa == min(dtavg$mappa))
    missmax <- which(dtavg$mappa == max(dtavg$mappa))
    normv[missmin] <- min(final$fitted)
    normv[missmax] <- max(final$fitted)
  }

  data_mappacorrect <- data_gccorrect/(normv/median(normv, na.rm = TRUE))

  return(data_mappacorrect)
}



#' Plot reads per cell per bin distributons
#'
#' This function creates a distribution plot for reads per cell
#'
#' @param data A matrix containing the read counts. Row names are features, column names are cells
#' @param project name for the plot
#' @return Creates a distribution plot.
#' @export
plotReadsDistributon <- function(data, project){
  reads_cell_bin <- unlist(data, use.names = FALSE)
  print(paste("Save reads per cell per bin distributons to ",project,"_reads_distributon.pdf",sep=""))
  grDevices::pdf(paste(project,"_reads_distributon.pdf",sep=""))
  par(mar=c(5,5,7,5))
  hist(reads_cell_bin,breaks = 3000,xlim = c(0,50),main=paste("Distribution of reads per cell per bin:",project,"\nMean:",round(mean(reads_cell_bin),1),"\nMedian:",median(reads_cell_bin),"\nMin:",min(reads_cell_bin),"\nMax:",max(reads_cell_bin)))
  abline(v=c(mean(reads_cell_bin),median(reads_cell_bin)), col=c("blue", "red"), lty=c(1,2), lwd=c(1, 3))
  grDevices::dev.off()
}

#' Plot Off-target Reads Coverage per cell
#'
#' This function creates a distribution plot for off-target Reads Coverage per cell
#'
#' @param data A matrix containing the read counts. Row names are features, column names are cells
#' @param project name for the plot
#' @return Creates a distribution plot and return an array of cells
#' @export
plotCoverageDistributon <- function(data, project,cov_cutoff){
  print(paste("Save reads per cell distributons to ",project,"_reads_distributon.pdf",sep=""))
  grDevices::pdf(paste(project,"_coverage_distributon.pdf",sep=""))
  hist(colSums(data),breaks=600,main = "Off-target Reads Coverage per cell",xlab = "Off-target Reads Coverage per cell")
  abline(v=cov_cutoff, col="red", lty=1, lwd=3 )
  grDevices::dev.off()
  return(colnames(data[,colSums(data)>=cov_cutoff]))
}



# plotCoverageDistributon2 <- function(data, project,cov_cutoff){
#   print(paste("Save reads per cell distributons to ",project,"_reads_distributon.pdf",sep=""))
#   grDevices::pdf(paste(project,"_coverage_distributon.pdf",sep=""))
#   hist(colSums(data),breaks=600,main = "Off-target Reads Coverage per cell",xlab = "Off-target Reads Coverage per cell")
#   abline(v=cov_cutoff, col="red", lty=1, lwd=3 )
#   grDevices::dev.off()
#   return(colnames(data[,colSums(data)>=cov_cutoff]))
# }



#' Calculate coverage
#'
#' This function calculates coverage by multiplying read counts by the average
#' fragment length and dividing by the length of the feature
#'
#' @param data Matrix of copy numbers
#' @param df dataframe containing a feature lengths column
#' @param avg_frag_length Average length of a genomic framgnet
#' @return coverage matrix
#' @export
calcCoverage <- function(data, df, avg_frag_length){
  coverage <- data
  # Get acutal lengh (peaks length removed) for each window
  lengths <- df[rownames(coverage),"length"]
  for( j in 1:ncol(data)){
    coverage[,j] = (data[,j] * avg_frag_length) / lengths
  }
  return(coverage)
}


#' Select GC matched background
#'
#' This function matches each feature to a specified number of its closest
#' GC matched neighbors. A certain number of these neighbors are then
#' sampled to create a GC matched background
#'
#' @param coverage coverage matrix (as produced by calcCoverage())
#' @param df dataframe containing the GC content of genomic features
#' @param GC_radius half the number of neighbors that will be matched to each feature
#' @param GC_sample_size number of features that will be sampled to produce the background
#' @return A matrix containing the names of GC matched features
#' @export
getGCBackground <- function(coverage, df, GC_radius, GC_sample_size){
  row_names <- rownames(coverage)
  GC_Content <- df[row_names,"GC_Content"]
  names_with_GC <- cbind(row_names,GC_Content)
  names_with_GC <- names_with_GC[order(names_with_GC[,2]),]

  name_table <- c()
  GC_table <- c()

  #sort by GC content
  for(i in 1:length(row_names)){
    target_GC = GC_Content[i]
    target_index = match(target_GC,  names_with_GC[,2])
    RADIUS = GC_radius
    SAMPLE_SIZE = GC_sample_size
    DIAMETER = RADIUS * 2

    if(target_index == 1){
      index_range = 2:(DIAMETER+1)
    }else if(target_index <= (RADIUS+1)){
      index_range = append(1:(target_index - 1),(target_index+1):(DIAMETER+1))
    }else if(target_index == length(row_names)){
      index_range = (length(row_names)-DIAMETER):(length(row_names)-1)
    }else if(target_index >= length(row_names)-RADIUS){
      index_range = append((length(row_names)-DIAMETER):(target_index - 1),(target_index+1):length(row_names))
    }else{
      index_range = append((target_index-RADIUS):(target_index-1), (target_index+1):(target_index+RADIUS))
    }
    name_row <- c()
    GC_row <- c()

    ## remove all background that have "chrX" or "chrY".
    index_range_chr <- do.call('rbind', strsplit(names_with_GC[index_range,1],'-',fixed=TRUE))[,1]
    index_range <-index_range[!(index_range_chr %in% c("chrX","chrY"))]

    sample_indexes = sample(index_range, SAMPLE_SIZE, replace=FALSE)
    for(index in sample_indexes){
      name_row <- append(name_row, names_with_GC[index,1])
      GC_row <- append(GC_row, names_with_GC[index,2])
    }

    name_table <- rbind(name_table,name_row)
    GC_table <- rbind(GC_table, GC_row)
  }
  rownames(name_table) <- row_names
  GC_table <- cbind(target=GC_Content, GC_table)
  return(name_table)
}

#' Calculate random selected background
#'
#' Calculates the average coverage of each row of feature names in name_table
#' in order to produce a background average. Any averages of zero are replaced
#' with the smallest value for that cell.
#'
#' @param coverage Coverage matrix (feature x cell)
#' @param name_table Table of feature names to include in background
#' @return average coverage table for each feature and cell
#' @export
calcAvgTable <- function(coverage, name_table){
  avg_table <- c()
  for(window in 1:nrow(coverage)){
    names <- name_table[rownames(coverage)[window],]
    avg_col <- rowMeans(coverage[names,])
    qu <- quantile(avg_col, c(0:4/4))
    names_median <- names[((avg_col < qu[3]) & (avg_col > qu[1]))]
    avg_row <- colMeans(coverage[names_median,])
    avg_table <- rbind(avg_table, avg_row)
  }

  ### Replace all 0 in background to the smallest number in the cell #######
  avg_table_no_zero <- avg_table
  for(cell in 1:ncol(avg_table_no_zero)){
    col <- avg_table_no_zero[,cell]
    col_min <- min(col[col > 0])
    avg_table_no_zero[col==0,cell] <- col_min
  }
  row.names(avg_table_no_zero) <- row.names(name_table)
  return(avg_table_no_zero)
}



# GCplot <- function(fc,df,project,cormethod){
#   GC_Content <- df[rownames(fc),"GC_Content"]
#   grDevices::pdf(paste(project,"_CNV_vs_GC_spearman.pdf",sep=""))
#   print(paste("Save GC-CNV corraltion to ",project,"CNV_vs_GC.pdf",sep=""))
#   co <- stats::cor.test(GC_Content,rowMeans(fc),method = cormethod)
#   graphics::plot(GC_Content,rowMeans(fc),xlab="GC",ylab="CNV fold change",pch = 20,main=paste("Correlation",round(co$estimate,2)))
#   grDevices::dev.off()
#   return(GC_Content)
# }
#
#
#
# ### Plot Corrlation plot for potential bias factors ######
# TADplot <- function(fc,tiles,df,project,TAD_cluster){
#   gr_a <- tiles
#   gr_b <- utils::read.table(TAD_cluster)
#   colnames(gr_b) <- c("chr", "start", "end", "name","width")
#   pairs <- findOverlapPairs(gr_a, makeGRangesFromDataFrame(gr_b,keep.extra.columns=T), ignore.strand = TRUE)
#   mcols(pairs)$overlap_width <- BiocGenerics::width(IRanges::pintersect(pairs, ignore.strand = TRUE))
#   tiles_tad <- as.data.frame(pairs)
#   tiles_tad$name <- paste(tiles_tad[,1],tiles_tad[,2],tiles_tad[,3],sep="-")
#   isIDmax <- with(tiles_tad, stats::ave(overlap_width, name, FUN=function(x) seq_along(x)==which.max(x)))==1
#   tiles_tad_uniq <- tiles_tad[isIDmax, ]
#   rownames(tiles_tad_uniq) <- tiles_tad_uniq$name
#   tiles_tad_uniq <- tiles_tad_uniq[,c("name","second.name")]
#   tiles_tad_uniq$tad <- "2_Intermediate"
#   tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_1","cluster_2"),]$tad <- "1_Cold"
#   tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_9","cluster_10"),]$tad <- "3_Hot"
#   tiles_tad_uniq$tad_num <- 2
#   tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_1","cluster_2"),]$tad_num <- 1
#   tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_9","cluster_10"),]$tad_num <- 3
#   tad <- as.character(tiles_tad_uniq[rownames(fc),]$tad)
#   tad_num <- as.numeric(tiles_tad_uniq[rownames(fc),]$tad_num)
#   pdata <- as.data.frame(cbind(tad,cnv=rowMeans(fc)))
#   pdata <- pdata[!is.na(pdata$tad),]
#   pdata$cnv <- as.numeric(levels(pdata$cnv))[pdata$cnv]
#   print(paste("Save TAD-CNV plot to ",project,"CNV_vs_TAD.pdf",sep=""))
#   grDevices::pdf(paste(project,"_CNV_vs_TAD.pdf",sep=""))
#   print(ggplot2::ggplot(pdata, ggplot2::aes_string(x='tad',y='cnv'))+
#           ggplot2::geom_boxplot() + ggplot2::geom_jitter(shape=16))
#   grDevices::dev.off()
#   return(tad_num)
# }
#
#
# ### Plot Corrlation plot for potential bias factors ######
# ATACplot <- function(cells, target_fragments, df, fc, project, cormethod){
#   target <- getBinMatrix(cells, df, target_fragments)
#   target <- t(target)
#   target_c <- colSums(target)
#   target_t <- target_c[rownames(fc)]
#   print(paste("Save ATAC-CNV corraltion to ",project,"CNV_vs_ATAC.pdf",sep=""))
#   grDevices::pdf(paste(project,"_CNV_vs_ATAC.pdf",sep=""))
#   co <- stats::cor.test(target_t,rowMeans(fc),method = cormethod)
#   plot(target_t,rowMeans(fc),xlab="ATAC",ylab="CNV fold change",pch = 20,main=paste("Correlation",round(co$estimate,2)))
#   grDevices::dev.off()
#   return(target_t)
# }



#' Clustering cells
#'
#' This function clusters cells using Louvain tree cluster.
#'
#' @param fc Fold change matrix
#' @param cov_cutoff_bc An array of cells
#' @param k Used for tree clustering. Smaller 'k' usually yields finer clusters.
#' @param d An integer scalar specifying the number of dimensions to use for a PCA on the expression matrix prior to the nearest neighbor search.
#' @param cluster_cell_cut_off Remove clusters that have number of cells less than this parameter
#' @param project names for the file to be saved
#' @return creates a dataframe with cell names and groups
#' @export
snnTreeCluster <- function(fc, cov_cutoff_bc, k = 10, d = 50, cluster_cell_cut_off, project){
  fc <- fc[, cov_cutoff_bc]

  print("Start Louvain tree clustering...")
  # Smaller 'k' usually yields finer clusters:
  # d: An integer scalar specifying the number of dimensions to use for a PCA on
  # the expression matrix prior to the nearest neighbor search.
  graph <- buildSNNGraph(fc, k = k, type = "number", d = d)
  lc <- cluster_louvain(graph)
  hc.clus <- membership(lc)
  print("Tree clustering finished.")

  SNN_results <- as.data.frame(cbind(colnames(fc),hc.clus))
  colnames(SNN_results) <- c("barcode","SNNTREEcluster")

  trclust <- SNN_results$SNNTREEcluster
  names(trclust) <- SNN_results$barcode

  trclustsort <- sort(trclust) # order by 1,2,3....

  print(paste("Remove clusters that have number of cells less than ", cluster_cell_cut_off,".",sep=""))
  cluster <- summary(as.factor(trclustsort))
  trclustsort <- trclustsort[trclustsort %in% names(cluster[cluster > cluster_cell_cut_off])]
  trclustsort_order <- trclustsort
  j=1
  for(i in names(sort(summary(as.factor(trclustsort)),decreasing=T))){
    trclustsort_order[names(trclustsort[trclustsort==i])] <- j
    j <- j+1
  }
  trclustsort_order <- factor(trclustsort_order,levels=names(summary(as.factor(trclustsort_order))))
  trclustsort_order <- sort(trclustsort_order)
  output <- data.frame(barcode=names(trclustsort_order),SNNTREEcluster=as.numeric(trclustsort_order))

  return(output)
}

#' Fit Gaussian mixture model
#'
#' This function calculates some statistics for each of the cell to help identify the tumor and normal group.
#'
#' @param fc Fold change matrix
#' @param SNN_results A dataframe with group information for every cell
#' @param GMMpara Pretrained parameters for Gaussian mixture model
#' @param null_dist Pre-calculated statistics distribution of a normal cohort
#' @param project names for the file to be saved
#' @return creates a dataframe with cell names, neutral proportion, neutral probability and p-value
#' @export
fitGMM <- function(fc, SNN_results, GMMpara, null_dist, project){
  chr_group_name = factor(sapply(rownames(fc), function(x) strsplit(x, '-')[[1]][1]),
                          levels = c(paste("chr", 1:22, sep=""), "chrX", "chrY")) ## only for plot

  df <- read.csv(GMMpara)
  cellnames = SNN_results$barcode
  dat = fc[, cellnames]

  pi1 <- df$pi1
  pi2 <- df$pi2
  mu1 <- df$mu1
  mu2 <- df$mu2
  sigma1 <- df$sigma1
  sigma2 <- df$sigma2

  compute_post_prob <- function(dat)
  {
    total_prob <- pi1*dnorm(as.matrix(dat),mu1,sigma1)/(pi1*dnorm(as.matrix(dat),mu1,sigma1) +
                                                          pi2*dnorm(as.matrix(dat),mu2,sigma2))
    #p2.post <- 1-p1.post
    total_prob[is.na(total_prob)]=0  # some NA when the read count is too large
    return(total_prob)
  }

  total_prob = compute_post_prob(dat)

  cutoff=0.5
  results=data.frame("neutral_prop" = colMeans(total_prob > cutoff), "neutral_prob" = colMeans(total_prob),
                     "cluster" = as.factor(SNN_results$SNNTREEcluster))
  print(paste("Save Gaussian mixture model boxplot to ", project, "_GMM_boxplot.pdf", sep = ""))
  pdf(paste(project, "_GMM_boxplot.pdf", sep = ""))
  print(ggplot(data=results,aes(x=cluster,y=neutral_prob,colour=cluster))+geom_boxplot())
  dev.off()

  # norm.cutoff= quantile(results$neutral_prop[results$cluster==norm_clust],0.05)
  # results %>% group_by(cluster) %>% summarise(prop=mean(neutral_prop>norm.cutoff))

  # compare neutral_prop with normal cluster, and get pvalue for each cell
  # null_dist = results$neutral_prob[results$cluster==norm_clust]
  null_dist = read.csv(null_dist)
  null_dist <- null_dist$null_dist
  results$neutral_prob_pval = sapply(results$neutral_prob, function(x) mean(null_dist < x))
  results %>% group_by(cluster) %>% summarise(prop=mean(neutral_prob_pval<0.05))

  return(results)
}

#' Plot CNV heatmap
#'
#' This function creates a heatmap using CNV data, where each row is a cell
#' and reach column is a genomic feature. Cells can be clustered together
#' hierarchically based on the tree clustering result
#'
#' @param fc Fold change matrix
#' @param data Read counts before correction
#' @param ginicutoff Used as a cutoff for the gini index
#' @param cov_cutoff_bc An array of cells
#' @param SNN_results a dataframe with cell names and groups
#' @param GMMresults a dataframe with cell names, neutral proportion, neutral probability and p-value
#' @param gz chromosome length data, produced by seqlengths()
#' @param hg19.genes GRanger object that shows the genes
#' @param ref.genome The genome reference, 'hg19' or 'hg38'
#' @param k Used in the plot file name to show the number of groups
#' @param project Used in the plot file name
#' @param a change the cutoff of the color bar
#' @param b change the cutoff of the color bar
#' @return creates heatmap and returns corresponding dataframe
#' @export
plotCNV <- function(fc, data, ginicutoff = 0.75, cov_cutoff_bc, SNN_results,
                    GMMresults, ref.genome = "hg38", k, project, a = 1, b = 2){
  if (ref.genome == "hg38"){
    gz<- seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)[1:24]
    hg19.genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  } else if (ref.genome == "hg19"){
    gz<- seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)[1:24]
    hg19.genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  }

  data <- data[, SNN_results$barcode]
  SNN_results$gini <- apply(data, 2, gini.wtd)
  SNN_results$seq_dep <- apply(data, 2, sum)
  GMMresults <- GMMresults[SNN_results$barcode, ]


  # filter out 0.75 gini
  SNN_results <- SNN_results[SNN_results$gini <= ginicutoff, ]
  # write.csv(SNN_results, paste(project, "_", k, "_rmGini_SNNTREEcluster.csv", sep = ""), row.names = FALSE)

  cov_cutoff_bc <- SNN_results$barcode
  GMMresults <- GMMresults[SNN_results$barcode, ]

  fc <- fc[,cov_cutoff_bc]
  mat2 <- t(fc)

  trclustsort_order <- SNN_results$SNNTREEcluster
  names(trclustsort_order) <- SNN_results$barcode

  hmdata <- mat2[SNN_results$barcode, ]
  rowclusterparam = FALSE
  # right_anno = ComplexHeatmap::rowAnnotation(Treecluster = as.factor(SNN_results$SNNTREEcluster),width = grid::unit(1, "mm"))

  ######## Define color list for each chromsome
  chrs <- do.call('rbind', strsplit(rownames(fc),'-',fixed=TRUE))[,1]
  ann <- factor(chrs,levels = names(gz)) ###### ann: factor: chr1, chr2, ..., chr24
  ann_color <- rep(c("gray97","gray30"),length(gz)/2)
  names(ann_color) <- names(gz)
  col1 <- list(chr = ann_color)

  ## Get average fc for each cluster
  cluster_group <- names(summary(as.factor(trclustsort_order)))
  cluster_cnv <- c()
  for(i in cluster_group){
    cluster_cells <- names(trclustsort_order[trclustsort_order==i])
    cluster_cnv <- rbind(cluster_cnv,colMeans(hmdata[cluster_cells,]))
  }

  ## Percent of cell numbers in each cluster
  percent <- paste(round(summary(as.factor(trclustsort_order))/sum(summary(as.factor(trclustsort_order))),2)*100,"%",sep="")

  chr <- do.call('rbind', strsplit(rownames(fc),'-',fixed=TRUE))[,1]
  start <- as.numeric(do.call('rbind', strsplit(rownames(fc),'-',fixed=TRUE))[,2])
  end <- as.numeric(do.call('rbind', strsplit(rownames(fc),'-',fixed=TRUE))[,3])
  loc <- as.data.frame(cbind(chr,start,end))

  ## Perform circular binary segmentation (CBS) on each cluster & create bottom annotation
  seg_cnv <- c()
  ha = HeatmapAnnotation(df = as.data.frame(t(cluster_cnv)))

  for(i in 1:nrow(cluster_cnv)){
    cluster_seg <- suppressMessages(segment(CNA(cluster_cnv[i,], chr, start)))$output
    gr_b <- cluster_seg[,c("chrom","loc.start","loc.end","seg.mean")]
    colnames(gr_b) <- c("chr", "start", "end", "seg.mean")
    pairs <- findOverlapPairs(makeGRangesFromDataFrame(loc), makeGRangesFromDataFrame(gr_b,keep.extra.columns=T), ignore.strand = TRUE)
    cluster_seg_window <- as.data.frame(pairs)[,c("first.seqnames","first.start","first.end","second.seg.mean")]
    cluster_seg_cnv <- as.numeric(cluster_seg_window$second.seg.mean)
    seg_cnv <- rbind(seg_cnv,cluster_seg_cnv)
    ha@anno_list[[i]] <- ComplexHeatmap::SingleAnnotation(fun = anno_points(cluster_seg_cnv, add_points = TRUE,
                                                                            height = grid::unit(2, "cm"),ylim=c(0,10)))
    ## Get MAD socre for each cluster
    MAD <- madDiff(cluster_cnv[i,])
    ha@anno_list[[i]]@name <- paste("cluster",i," / ",percent[i]," \nMAD ",round(MAD,2),sep = "")
    ha@anno_list[[i]]@label <- paste("anno",i," / ",percent[i]," \nMAD ",round(MAD,2),sep = "")
  }

  ha@anno_size <- unit(rep("20",nrow(cluster_cnv)), "mm", data=NULL)
  ha@height <- unit(22*nrow(cluster_cnv), "mm", data=NULL)

  # if (hg19.genes == "hg38") {
  #   transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene")
  # } else if (hg19.genes == "hg19") {
  #   transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
  # }

  gene_symbol<- AnnotationDbi::select(org.Hs.eg.db, keys=hg19.genes$gene_id,
                                      columns="SYMBOL", keytype="ENTREZID")
  all.equal(hg19.genes$gene_id, gene_symbol$ENTREZID)
  hg19.genes$gene_id<- gene_symbol$SYMBOL


  rownames(seg_cnv) <- rownames(cluster_cnv)
  colnames(seg_cnv) <- colnames(cluster_cnv)

  print(paste("Save segmentation CNV to ",project,"_seg_cnv_clusters.rds",sep=""))

  saveRDS(seg_cnv, file = paste(project,"_seg_cnv_clusters.rds",sep=""))


  cnv_cluster_gene <- c()
  for(i in (1:nrow(seg_cnv))){
    gr_b <- cbind(loc,cnv = seg_cnv[i,])
    pairs <- findOverlapPairs(hg19.genes, makeGRangesFromDataFrame(gr_b,keep.extra.columns=T), ignore.strand = TRUE)
    cluster_gene <- data.frame(chr = seqnames(pairs@first) %>% as.character,
                               start = start(pairs@first) %>% as.numeric,
                               end = end(pairs@first) %>% as.numeric,
                               gene_id = elementMetadata(pairs@first)$gene_id,
                               cnv = elementMetadata(pairs@second)$cnv,
                               cluster = i)
    cnv_cluster_gene <- rbind(cnv_cluster_gene,cluster_gene)
  }
  print(paste("Save Gene CNV to ",project,"_Gene_Copy_Number_per_Cluster.csv",sep=""))
  write.table(cnv_cluster_gene,paste(project,"_Gene_Copy_Number_per_Cluster.csv",sep=""),sep = "," ,quote = F,row.names = F)

  pdf_height <- 12 + 1.2* nrow(seg_cnv)

  cell_order <- rownames(hmdata)

  pval <- GMMresults[cell_order,]

  Aneuploid_score <- -log10(pval$neutral_prob_pval)
  Aneuploid_score[Aneuploid_score > 3] <- 3
  summary(Aneuploid_score)

  right_anno = rowAnnotation(Ane = Aneuploid_score, col = list(Ane = circlize::colorRamp2(c(0, 1.3, 3), c("white","orange","red"))))


  print(paste("Save CNV heatmap to ",project,"_treecluster_",k,"_heatmap.pdf",sep=""))
  grDevices::pdf(paste(project,"_treecluster_",k,"_heatmap.pdf",sep=""),  width=24,height=pdf_height )

  print(ComplexHeatmap::Heatmap(hmdata, col = circlize::colorRamp2(c(0, 0.5, a, b), c("cornflowerblue","white","white","red")),
                                name = "CNV score", column_title = paste(project," - CNV from offtarget scATAC reads",sep=""),
                                show_column_names = FALSE, show_row_names = FALSE,
                                heatmap_legend_param = list(title = "CNV-fc                 ", direction = c("horizontal")),
                                cluster_columns =FALSE, cluster_rows = rowclusterparam,
                                top_annotation = ComplexHeatmap::HeatmapAnnotation(pt = ComplexHeatmap::anno_empty(border = FALSE,height = grid::unit(4, "mm")),
                                                                                   chr = ann,show_legend = F,col = col1,border = T),
                                right_annotation = right_anno,
                                row_split = trclustsort_order,row_gap = unit(0.5, "mm"),
                                row_dend_width = unit(60, "mm"),
                                bottom_annotation = ha,
                                column_split = ann,column_gap = unit(0, "mm")))

  for(i in 1:length(unique(sub("\\-.*", "", colnames(hmdata))))) { # length(names(gz))
    ComplexHeatmap::decorate_annotation("pt", slice = i, {
      grid::grid.text(gsub("chr","",paste(names(gz)[i])), just = "centre")
    })
  }

  grDevices::dev.off()
}


# ### Plot Corrlation plot for potential bias factors ######
# chrCNVplot <- function(fc, gz, project){
#   CNV <- rowMeans(fc)
#   chrs <- do.call('rbind', strsplit(rownames(fc),'-',fixed=TRUE))[,1]
#   pdata <- as.data.frame(cbind(chrs,CNV))
#   pdata$CNV <- as.numeric(levels(pdata$CNV))[pdata$CNV]
#   pdata$chrs <- factor(chrs,levels = names(gz))
#
#
#
#   print(paste("Save CNV boxplot to ",project,"_chr_CNV_boxplot.pdf",sep=""))
#   grDevices::pdf(paste(project,"_chr_CNV_boxplot.pdf",sep=""),  width=14,height=10 )
#   print(ggplot2::ggplot(pdata, ggplot2::aes_string(x='chrs',y='CNV'))+
#     ggplot2::geom_boxplot() + ggplot2::geom_jitter(shape=16) +
#     ggplot2::scale_y_continuous(breaks = seq(0, 10),limits = c(0,10))+
#     ggplot2::geom_hline(yintercept = 1,size=1.5) +
#     ggplot2::ggtitle(paste("CNV fold change distributon / MAD score:",round(madDiff(CNV),2)))+
#     theme(plot.title = element_text(hjust = 0.5)))
#   grDevices::dev.off()
#   return(pdata)
# }


#' Supercell - PCA & hierarchical clustering
#' Clusters cells into super cells using Principle Component Analysis
#'
#' @param fc Fold change matrix
#' @param cellsPerGroup number of cells per supercell
#' @return dataframe containing the supercell groups
#' @export
superCells <- function(fc,cellsPergroup){
  print("Start PCA...")
  fc.pca <- stats::prcomp(t(fc),retx = TRUE, center = TRUE, scale. = TRUE)
  cordata <- stats::cor(t(fc.pca$x[,1:50]), method="spearman")
  distance <- (1-cordata)/2
  print("Calculating distance matrix...")
  rowdistance = stats::dist(as.matrix(distance), method = "euclidean")
  print("Start hierarchical clustering...")
  rowcluster = stats::hclust(rowdistance, method = "ward.D2")
  #plot(rowcluster)
  print("Create supercell groups...")
  supercells <- stats::cutree(rowcluster, k = (ncol(fc)/cellsPergroup), h = NULL)
  print("Done")
  output <- data.frame(Barcode=names(supercells),Supercells=paste("supercell_",supercells,sep = ""))
  return(output)
}


#' Supercell - merge cells
#'
#' This function merges the cnv data of a group of cells according
#' to their supercell grouping
#' @param data matrix containing cnv data
#' @param supercells matrix containing supercell grouping assignemnts
#' @return supercells with calculated cnv data
#' @export
superCells_merge <- function(data,supercells){
  data_super <- c()

  for(i in unique(supercells$Supercells)){
    data_super <- cbind(data_super,rowSums(data[,as.character(supercells[supercells$Supercells==i,]$Barcode)]))
  }

  colnames(data_super) <- as.character(unique(supercells$Supercells))
  return(as.matrix(data_super))
}


