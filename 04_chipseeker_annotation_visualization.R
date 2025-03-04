# load in required packages
library(ChIPseeker)
# library(clusterProfiler)
# library(ReactomePA)
# library(tidyverse)
library(openxlsx)
library(optparse)


option_list <- list(
  make_option(c("--beddir"),
              type = "character",
              help = "Directory containing input BED files",
              metavar = "DIR"),

  make_option(c("--bedpattern"),
              type = "character",
              help = "Suffix pattern of BED files (e.g., .rearrangedCols.bed)",
              metavar = "PATTERN"),

  make_option(c("--out"),
              type = "character",
              default = "./04.peak.anno.plot/",
              help = "Output directory for results",
              metavar = "OUTDIR"),

  make_option(c("--hg"),
              type = "character",
              default = "hg19",
              help = "Genome assembly (hg19/hg38) [default %default]"))

opt <- parse_args(OptionParser(option_list = option_list),
                  args = commandArgs(trailingOnly = TRUE))
beddir <- normalizePath(opt$beddir)
bedpattern <- opt$bedpattern
outdir <- normalizePath(opt$out)
hg <- tolower(opt$hg)


# which human genome version to use
if(hg == "hg19"){
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  selchr <- c(paste0("chr", c(1:22, "X", "Y")))
  # target.pos <- GRanges("chr2", IRanges(start = 134877506 - 20000,
  #                                       end = 135212192 + 20000))
  # target.pos$SYMBOL <- "MGAT5"
}else if(hg == "hg38"){
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  selchr <- c(paste0("chr", c(1:22, "X", "Y")))
}else if(hg == "mm10"){
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  selchr <- c(paste0("chr", c(1:19, "X", "Y")))
  # target.pos <- GRanges("chr1",
  #                       IRanges(start = 127205018 - 20000,
  #                               end = 127488336 + 20000))
  # target.pos$SYMBOL <- "Mgat5"
}


# set the output directory
setwd(outdir)

# locations of all input files
files <- list.files(path = beddir,
                    pattern = bedpattern , # "stringent.auc.threshold.merge.bed",
                    full.names = TRUE,
                    include.dirs = FALSE)

fn <- files
fn <- gsub("^.*/", "", fn)
fn <- gsub(".rearrangedCols.bed", "", fn)
fn <- gsub("_hgAligned_stampyFiltered.mappedSorted.normalized.bed", "", fn)
fn <- gsub(".rearrangedCols.bed", "", fn)

# Make a list of all the names because.... thats what it wants
files <- as.list(files)
names(files) <- fn

# iterate through and make individual plots for each peaks file
# current_num <- 1
for( current_num in 1:length(files) ) {
    # what is the name for the current sample
    current_name <- names(files)[current_num]

    # name of the current peak file
    current_peak_file_name <- current_name
    message(current_peak_file_name,": ",basename( files[[current_num]] ) )
    #  read in a file
    peak <- readPeakFile( files[[current_num]] )
    ## only keep peaks with chr1:22 and chrXY for human
    ## only keep peaks with chr1:19 and chrXY for mouse
    index <- seqnames(peak) %in% selchr
    peak <- peak[index,]
    seqlevels(peak) <- seqlevels(peak)[seqlevels(peak) %in% selchr]


    # TSS peak binding
    promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
    tagMatrix <- getTagMatrix(peak, windows=promoter)

    # average profile straight from bed
    print("TSS average peaks profile")
    pdf(file = paste0(current_peak_file_name, "_AveragePeaks.pdf"),
        width = 8, height = 4.5, pointsize = 6)
    p <- plotAvgProf(tagMatrix, conf = 0.95, xlim = c(-3000, 3000),
                     xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
    print(p)
    dev.off()


    print("peaks profile")
    pdf(file = paste0(current_peak_file_name, "_PeaksProf.pdf"),
        width = 8, height = 4.5, pointsize = 6)
    p <- plotPeakProf2(peak, upstream = rel(0.2), downstream = rel(0.2),
                       conf = 0.95, by = "gene", type = "body", nbin = 800,
                       TxDb = txdb, weightCol = "V5",ignore_strand = F)
    print(p)
    dev.off()


    # peaks coverage plot
    print("coverage plot - overall")
    pdf(file = paste0(current_peak_file_name, "_PeaksOverall.pdf"),
        width = 8, height = 10, pointsize = 6)
    plot( covplot(peak, weightCol="V5") )
    dev.off()


    # peak annotation
    peakAnno <- annotatePeak(peak, # files[[current_num]],
                             tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
    paste("saving annotation texts")
    write.table(peakAnno, file = paste0(current_peak_file_name, "_AnnotatedPeaks.txt"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    write.table(peakAnno@annoStat, file = paste0(current_peak_file_name, "_AnnotatedPeaksDistribution.txt"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    # visualise
    print("pie chart")
    pdf(file = paste0(current_peak_file_name, "_PeakBinding_PieChart.pdf"),
        width = 6, height = 4.5, pointsize = 6)
    plotAnnoPie(peakAnno)
    dev.off()

    print("bar chart")
    pdf(file = paste0(current_peak_file_name, "_PeakBinding_BarChart.pdf"),
        width = 8, height = 4.5, pointsize = 6)
    plot(plotAnnoBar(peakAnno))
    dev.off()

    print("venn diagram")
    pdf(file = paste0(current_peak_file_name, "_PeakBinding_VennDiagram.pdf"),
        width = 8, height = 8, pointsize = 8)
    # plot(vennpie(peakAnno))
    vennpie(peakAnno)
    dev.off()

    print("upset plot")
    pdf(file = paste0(current_peak_file_name, "_PeakBinding_UpsetPlot.pdf"),
        width = 11, height = 8, pointsize = 8)
    print(upsetplot(peakAnno))
    dev.off()


    # TF binding relative to TSS
    print("tss binding")
    pdf(file = paste0(current_peak_file_name, "_PeakBinding_TSS.pdf"),
        width = 8, height = 2.5, pointsize = 6)
    plot( plotDistToTSS(peakAnno,
                  title="Distribution of transcription factor-binding loci\nrelative to TSS")
    )
    dev.off()

}
