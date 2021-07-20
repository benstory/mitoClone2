#'Create a list object from a list of single-cell BAM files where each
#' contains a matrix of the of AGCT nt counts at chosen sites
#'
#' Uses the \code{deepSNV} package to count nucleotide frequencies at
#' every position in the mitochondrial genome for every cell.
#'@param bamfiles A character vector specifying the bam file paths
#'@param sites String specifying genomic regions, defaults to the
#'entire mitochondrial genome
#'@param ncores Number of threads to use for the computation. Default
#'1
#'@param ignore_nonstandard Ignore basecalls that are not AGCTN
#'@return A list of base count matrices which can serve as an input to
#'\code{\link{mutationCallsFromExclusionlist}} or
#'\code{\link{mutationCallsFromCohort}}
#'@examples bamCounts <- baseCountsFromBamList(bamfiles =
#'list(system.file("extdata", "mm10_10x.bam",
#'package="mitoClone2")),sites="chrM:1-15000", ncores=1)
#'@export baseCountsFromBamList
baseCountsFromBamList <- function(bamfiles,
                                  sites = "chrM:1-16569",
                                  ncores = 1,
                                  ignore_nonstandard = FALSE) {
    if (!length(sites) == 1) {
        stop('Your sites parameter must be a character/GRanges object of length 1')
    }
    mito.chr <- GenomicRanges::GRanges(sites)
    mc.out <- parallel::mclapply(bamfiles, function(bampath) {
        bam.file <-
            deepSNV::bam2R(
                         bampath,
                         chr = GenomicRanges::seqnames(mito.chr),
                         start = GenomicRanges::start(mito.chr),
                         stop = GenomicRanges::end(mito.chr)
                     )
        bam.file.sub <-
            bam.file[, c("a", "t", "c", "g", "_", "n", "ins", "del")]
        bam.file <-
            bam.file[, c("A", "T", "C", "G", "-", "N", "INS", "DEL")]
        bam.file <- bam.file + bam.file.sub
        if (ignore_nonstandard) {
            bam.file <- bam.file[, c('A', 'T', 'C', 'G', 'N')]
        }
        return(bam.file)
    }, mc.cores = ncores)
    names(mc.out) <- names(bamfiles)
    return(mc.out)
}

#' Read nucleotide counts from a 10x Genomics .bam file
#'
#' This function uses a C interface to read the nucleotide counts on
#' each position of a .bam alignment. The counts are individually
#' tabulated for each cell barcode as specified by the user. The
#' counts of both strands are reported separately and nucleotides
#' below a quality cutoff are masked.
#'
#' This code is an adaption of code that was originally written by
#' Moritz Gerstung for the deepSNV package
#'
#' @param file The file location of the BAM file as a string.
#' @param sites The chromosome locations of interest in BED format as
#'a string. Alternatively a single GRanges object will also work.
#' @param q An optional cutoff for the nucleotide Phred
#'quality. Default q = 25. Nucleotides with Q < q will be masked
#'by 'N'.
#' @param mq An optional cutoff for the read mapping quality. Default
#'mq = 0 (no filter). reads with MQ < mq will be discarded.
#' @param s Optional choice of the strand. Defaults to s = 2 (both).
#' @param head.clip Should n nucleotides from the head of reads be
#'clipped? Default 0.
#' @param max.depth The maximal depth for the pileup command. Default
#'1,000,000.
#' @param verbose Boolean. Set to TRUE if you want to get additional
#'output.
#' @param mask Integer indicating which flags to filter. Default 0 (no
#'mask). Try 1796 (BAM_DEF_MASK).
#' @param keepflag Integer indicating which flags to keep. Default 0
#'(no mask). Try 3 (PAIRED|PROPERLY_PAIRED).
#' @param max.mismatches Integer indicating maximum MN value to allow
#'in a read. Default NULL (no filter).
#' @param ncores Integer indicating the number of threads to use for
#'the parallel function call that summarize the results for each
#'bam file. Default 1.
#' @param ignore_nonstandard Boolean indicating whether or not gapped
#'alignments, insertions, or deletions should be included in the
#'final output. Default FALSE. If you have an inflation of
#'spliced mitochondrial reads it is recommended to set this to
#'TRUE.
#' @return A named \code{\link{list}} of \code{\link{matrix}} with
#'rows corresponding to genomic positions and columns for the
#'nucleotide counts (A, T, C, G, -), masked nucleotides (N),
#'(INS)ertions, (DEL)etions that count how often a read begins
#'and ends at the given position, respectively. Each member of
#'the list corresponds to an invididual cells or entity based on
#'the cell barcode of interest. The names of the elements of the
#'list correspond to the respective cell barcodes.  For the
#'intents and purposes of the mitoClone2 package this object is
#'equivalent to the output from the
#'\code{\link{baseCountsFromBamList}} function.  The returned
#'list has a variable length depending on the ignore_nonstandard
#'parameter and each element contains a matrix has 8 columns and
#'(stop - start + 1) rows. The two strands have their counts
#'merged. If no counts are present in the provided sites
#'parameter nothing will be returned.  IMPORTANT: The names of
#'the list will NOT reflect the source filename and will
#'exclusively be named based on the respective the barcodes
#'extracted from said file. If merging multiple datasets, it is
#'important to change the list's names once imported to avoid
#'naming collisions.
#' @examples bamCounts <- bam2R_10x(file = system.file("extdata",
#'"mm10_10x.bam", package="mitoClone2"), sites="chrM:1-15000")
#' @author Benjamin Story (adapted from original code with permission
#'from Moritz Gerstung)
#' @export bam2R_10x
bam2R_10x <- function(file,
                      sites = "MT:1-16569",
                      q = 25,
                      mq = 0,
                      s = 2,
                      head.clip = 0,
                      max.depth = 1000000,
                      verbose = FALSE,
                      mask = 0,
                      keepflag = 0,
                      max.mismatches = NULL,
                      ncores = 1,
                      ignore_nonstandard = FALSE) {
    if (!length(sites) == 1) {
        stop('Your sites parameter must be a character/GRanges object of length 1')
    }
    mito.chr <- GenomicRanges::GRanges(sites)
    chr = GenomicRanges::seqnames(mito.chr)
    start = GenomicRanges::start(mito.chr)
    stop = GenomicRanges::end(mito.chr)
    if (is.null(max.mismatches))
        max.mismatches <- -1
    region = paste(chr, ":", start, "-", stop, sep = "")
    result = .Call(
        "bam2R_10x",
        as.character(file),
        as.character(chr),
        as.integer(start),
        as.integer(stop),
        vector("integer", (stop - start + 1) * 11 * 2),
        as.integer(q),
        as.integer(mq),
        as.integer(s),
        as.integer(head.clip),
        as.integer(max.depth),
        as.integer(verbose),
        as.integer(mask),
        as.integer(keepflag),
        as.integer(max.mismatches)
    )
    barcode.n <- names(result)
    result <- parallel::mclapply(result, function(mat) {
        bam.file <-
            matrix(mat,
                   nrow = stop - start + 1,
                   dimnames = list(
                       NULL,
                       c('A','T','C','G','-','N',
                         'INS','DEL','HEAD','TAIL','QUAL',
                         'a','t','c','g','_','n',
                         'ins','del','head','tail','qual')
                   ))
        bam.file.sub <-
            bam.file[, c("a", "t", "c", "g", "_", "n", "ins", "del")]
        bam.file <-
            bam.file[, c("A", "T", "C", "G", "-", "N", "INS", "DEL")]
        bam.file <- bam.file + bam.file.sub
        if (ignore_nonstandard) {
            bam.file <- bam.file[, c('A', 'T', 'C', 'G', 'N')]
        }
        return(bam.file)
    }, mc.cores = ncores)
    names(result) <- barcode.n
    return(result)
}
