## Version 2 of the 10x nucleotide counting (src/bam2R_v2.cpp).
##
## Independent of version 1 (bam2R_10x() in bam_functions.R and
## src/bam2R_10x.cpp), both can be used side by side.
##
## Layers:
##   bam2R_v2()           core reader: sparse counts, any regions, one pass
##     |-- mutationCallsFromExclusionlist()  \
##     |-- mutationCallsFromCohort()          > accept the sparse counts directly
##     |-- pullcountsVars()                  /  (see "INTEGRATION" at the end of this file)
##     |-- bam2R_10x_v2()       dense per cell matrices, same output as bam2R_10x()
##     |-- pullcountsVars10x()  M / N count tables of variants
##     |-- strandBias10x()      strand bias of variants
##
## Only bam2R_v2() talks to C (via .Call: bam2R_v2_barcodes, bam2R_v2_count),
## everything else works on its output.

## Internal: mclapply that fails loudly. mclapply returns NULL for a
## worker that died (e.g. killed when out of memory) and a try-error
## for a worker that failed, both would silently drop data
## otherwise. Windows does not support more than one core.
mclapply10x <- function(X, FUN, ncores) {
    if (.Platform$OS.type == 'windows') {
        ncores <- 1
    }
    res <- parallel::mclapply(X, FUN, mc.cores = ncores)
    failed <- vapply(res, function(r) is.null(r) || inherits(r, 'try-error'), logical(1))
    if (any(failed)) {
        err <- res[failed][[1]]
        if (is.null(err)) {
            stop('A parallel worker did not return a result (out of memory?). ',
                 'Try again with a lower ncores.', call. = FALSE)
        }
        stop(conditionMessage(attr(err, 'condition')), call. = FALSE)
    }
    return(res)
}

#' Read nucleotide counts from a 10x Genomics .bam file (version 2, sparse)
#'
#' This function uses a C interface to read the nucleotide counts on
#' each position of a .bam alignment. The counts are individually
#' tabulated for each cell barcode. In contrast to
#' \code{\link{bam2R_10x}} any number of regions on any chromosome can
#' be provided and all of them are read in one single pass over the
#' file. Only the site/barcode combinations that were actually
#' observed are stored, which keeps the memory footprint small. The
#' counts of both strands are reported separately and nucleotides
#' below a quality cutoff are masked.
#'
#' The result can be passed directly to the mutation calling of the
#' package: \code{\link{mutationCallsFromExclusionlist}},
#' \code{\link{mutationCallsFromCohort}} and
#' \code{\link{pullcountsVars}} accept it in place of the list with
#' one matrix per cell and return the same result. Variants are
#' called on the mitochondrial chromosome of the counts.
#' \code{\link{strandBias10x}} and \code{\link{bam2R_10x_v2}} work
#' on the output as well.
#'
#' This code is an adaption of code that was originally written by
#' Moritz Gerstung for the deepSNV package
#'
#' @param file The file location of the BAM file as a string.
#' @param sites The regions of interest. A GRanges or GRangesList
#'object with any number of ranges on any chromosome, alternatively
#'a character vector that can be converted to one
#'(e.g. c("chrM:1-16569","chr7:140753336")). Overlapping ranges
#'are merged.
#' @param barcodes Character vector of the cell barcodes to tabulate
#'(e.g. the filtered barcodes from Cell Ranger). If provided the
#'file is only read once and min_reads_per_barcode is ignored.
#'Default NULL (find barcodes using min_reads_per_barcode).
#' @param q An optional cutoff for the nucleotide Phred
#'quality. Default q = 25. Nucleotides with Q <= q will be masked
#'by 'N' (only nucleotides with Q > q are counted).
#' @param mq An optional cutoff for the read mapping quality. Default
#'mq = 0 (no filter). reads with MQ < mq will be discarded.
#' @param head.clip Number of nucleotides at the head (5' end) of
#'each read that are masked by 'N'. Default 0.
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
#' @param ncores Integer indicating the number of threads to use. The
#'sites are split into ncores chunks of similar size that are read
#'in parallel. Default 1.
#' @param min_reads_per_barcode Int defining how many reads a barcode
#'must for it to be considered in the pileup tabulations. Default 50
#' @param barcode.tag The BAM tag that contains the cell barcode.
#'Default CB.
#' @param umi.tag The BAM tag that contains the UMI (e.g. UB for Cell
#'Ranger). If provided UMIs are counted instead of reads: all reads
#'of a cell that share a UMI and show the same nucleotide on the
#'same strand at a site are counted once (removes PCR duplicates).
#'Reads without this tag are discarded. Default NULL (count reads).
#' @return A named \code{\link{list}} of 16 sparse matrices (class
#'dgCMatrix) with rows corresponding to covered genomic sites
#'(named chromosome:position) and columns corresponding to cell
#'barcodes. There is one matrix for each of the nucleotide counts
#'(A, T, C, G, -), masked nucleotides (N), (INS)ertions and
#'(DEL)etions on the forward strand and one for each of these on
#'the reverse strand (a, t, c, g, _, n, ins, del).  Sites without
#'any coverage are not included. Use
#'\code{\link{mutationCallsFromExclusionlist}} or
#'\code{\link{mutationCallsFromCohort}} to call variants,
#'\code{\link{pullcountsVars}} to extract the M and N count
#'tables of known variants and \code{\link{strandBias10x}} to
#'check variants for strand bias.
#' @examples sparseCounts <- bam2R_v2(file = system.file("extdata",
#'"mm10_10x.bam", package="mitoClone2"), sites="chrM:1-15000")
#'## call variants, with thresholds for this tiny example of 5 cells
#'mutCalls <- mutationCallsFromExclusionlist(sparseCounts, lim.cov = 5,
#'min.af = 0.1, min.num.samples = 1, universal.var.cells = 5,
#'max.var.na = 1, max.cell.na = 1.1, genome = 'mm10')
#' @author Benjamin Story (adapted from original code with permission
#'from Moritz Gerstung)
#' @export bam2R_v2
bam2R_v2 <- function(file,
                     sites = "MT:1-16569",
                     barcodes = NULL,
                     q = 25,
                     mq = 0,
                     head.clip = 0,
                     max.depth = 1000000,
                     verbose = FALSE,
                     mask = 0,
                     keepflag = 0,
                     max.mismatches = NULL,
                     ncores = 1,
                     min_reads_per_barcode = 50,
                     barcode.tag = 'CB',
                     umi.tag = NULL) {
    if (methods::is(sites, 'GRangesList')) {
        sites <- unlist(sites)
    }
    ## merge overlapping ranges so no site is counted twice
    sites <-
        GenomicRanges::reduce(GenomicRanges::GRanges(sites), ignore.strand = TRUE)
    if (is.null(max.mismatches))
        max.mismatches <- -1
    ## BAM tags are always two characters, anything else would be cut or garbage in C
    if (!is.character(barcode.tag) || length(barcode.tag) != 1 || nchar(barcode.tag) != 2) {
        stop("barcode.tag must be a two character BAM tag, e.g. 'CB'")
    }
    if (!is.null(umi.tag) &&
        (!is.character(umi.tag) || length(umi.tag) != 1 || nchar(umi.tag) != 2)) {
        stop("umi.tag must be NULL or a two character BAM tag, e.g. 'UB'")
    }
    ## an empty tag turns the UMI counting off
    if (is.null(umi.tag))
        umi.tag <- ''
    ## first pass - only needed if the user does not know the barcodes
    if (is.null(barcodes)) {
        barcode.n <- .Call(
            "bam2R_v2_barcodes",
            as.character(file),
            as.character(GenomicRanges::seqnames(sites)),
            as.integer(GenomicRanges::start(sites)),
            as.integer(GenomicRanges::end(sites)),
            as.character(barcode.tag)
        )
        ## sorted, so the order does not depend on the hash table (same order as bam2R_10x)
        barcodes <- sort(names(barcode.n)[barcode.n >= min_reads_per_barcode], method = 'radix')
        if (verbose) {
            message('Unique barcodes found: ', length(barcode.n))
            message('Barcodes with >= ', min_reads_per_barcode, ' reads: ', length(barcodes))
        }
        if (length(barcode.n) == 0) {
            warning('No read in the sites has the barcode tag ', barcode.tag,
                    '. Is barcode.tag correct (Cell Ranger: CB)?', call. = FALSE)
        } else if (length(barcodes) == 0) {
            warning('None of the ', length(barcode.n), ' barcodes has at least',
                    ' min_reads_per_barcode = ', min_reads_per_barcode,
                    ' reads in the sites (the most reads of a barcode: ', max(barcode.n), ').',
                    call. = FALSE)
        }
    }
    barcodes <- unique(as.character(barcodes))
    tiled <- tileSites10x(sites, ncores)
    result <- mclapply10x(tiled$chunks, function(chunk) {
        countRegions10x(file, chunk, barcodes, q, mq, head.clip, max.depth, verbose,
                        mask, keepflag, max.mismatches, barcode.tag, umi.tag)
    }, ncores)
    result <- assembleSparse10x(result, tiled$tiles, barcodes, ncores)
    if (length(barcodes) > 0 && nrow(result[[1]]) == 0) {
        warning(noCoverage10x(file, sites, barcodes, barcode.tag, umi.tag), call. = FALSE)
    }
    return(result)
}

## Internal: why did the selected cells not cover any site? Only run if
## that happened, it reads the barcodes of the file once more.
noCoverage10x <- function(file, sites, barcodes, barcode.tag, umi.tag) {
    msg <- 'No read of the selected cells covers the sites. '
    barcode.n <- .Call(
        "bam2R_v2_barcodes",
        as.character(file),
        as.character(GenomicRanges::seqnames(sites)),
        as.integer(GenomicRanges::start(sites)),
        as.integer(GenomicRanges::end(sites)),
        as.character(barcode.tag)
    )
    if (length(barcode.n) == 0) {
        return(paste0(msg, 'No read in the sites has the barcode tag ', barcode.tag,
                      '. Is barcode.tag correct (Cell Ranger: CB)?'))
    }
    if (!any(barcodes %in% names(barcode.n))) {
        return(paste0(msg, 'None of the barcodes occurs in the file, which contains e.g. ',
                      paste(utils::head(names(barcode.n), 2), collapse = ', '),
                      '. Are the barcodes written the same way (e.g. the -1 suffix)?'))
    }
    if (umi.tag != '') {
        return(paste0(msg, 'The cells have reads, but none of them has the UMI tag ', umi.tag,
                      ' (reads without it are skipped). Is umi.tag correct (Cell Ranger: UB)?'))
    }
    paste0(msg, 'The cells have reads, but none of them passed the filters',
           ' (q, mq, mask, keepflag, max.mismatches).')
}

## Internal: cut the (reduced) sites into ncores chunks of similar
## width. Every requested position gets a site number: 1 for the first
## position of the first tile up to the total number of positions.
tileSites10x <- function(sites, ncores) {
    chunk.width <- ceiling(sum(GenomicRanges::width(sites)) / ncores)
    tiles <- unlist(GenomicRanges::tile(sites, width = chunk.width))
    tiles$offset <-
        as.integer(cumsum(c(0, GenomicRanges::width(tiles)))[seq_along(tiles)])
    chunks <-
        split(tiles, ceiling(cumsum(GenomicRanges::width(tiles)) / chunk.width))
    list(tiles = tiles, chunks = chunks)
}

## Internal: count the regions of a chunk of tiles in one file (C code)
countRegions10x <- function(file, chunk, barcodes, q, mq, head.clip, max.depth, verbose,
                            mask, keepflag, max.mismatches, barcode.tag, umi.tag) {
    .Call(
        "bam2R_v2_count",
        as.character(file),
        as.character(GenomicRanges::seqnames(chunk)),
        as.integer(GenomicRanges::start(chunk)),
        as.integer(GenomicRanges::end(chunk)),
        chunk$offset,
        barcodes,
        as.integer(q),
        as.integer(mq),
        as.integer(head.clip),
        as.integer(max.depth),
        as.integer(verbose),
        as.integer(mask),
        as.integer(keepflag),
        as.integer(max.mismatches),
        as.character(barcode.tag),
        as.character(umi.tag)
    )
}

## Internal: the results of the C code (site number, barcode index and
## count per column) -> 16 sparse matrices (covered sites x cells)
assembleSparse10x <- function(result, tiles, cells, ncores) {
    ## only covered sites become rows: translate the site number back to chromosome:position
    sites.covered <- sort(unique(unlist(lapply(result, '[[', 'sites'), use.names = FALSE)))
    tile.n <- findInterval(sites.covered, tiles$offset + 1)
    site.name <-
        sprintf(
            '%s:%d',
            as.character(GenomicRanges::seqnames(tiles))[tile.n],
            GenomicRanges::start(tiles)[tile.n] + sites.covered - tiles$offset[tile.n] - 1L
        )
    ## one sparse matrix per column, capital letters = forward strand
    nt.names <- c('A','T','C','G','-','N','INS','DEL',
                  'a','t','c','g','_','n','ins','del')
    out <- mclapply10x(seq_along(nt.names), function(n) {
        ## merge the chunks
        site <- unlist(lapply(result, function(res) res$site[[n]]), use.names = FALSE)
        barcode <- unlist(lapply(result, function(res) res$barcode[[n]]), use.names = FALSE)
        count <- unlist(lapply(result, function(res) res$count[[n]]), use.names = FALSE)
        Matrix::sparseMatrix(
                    i = match(site, sites.covered),
                    j = barcode,
                    x = as.numeric(count),
                    dims = c(length(site.name), length(cells)),
                    dimnames = list(site.name, cells)
                )
    }, ncores)
    names(out) <- nt.names
    return(out)
}

#' Read nucleotide counts from one BAM file per cell (version 2, sparse)
#'
#' Sparse counterpart of \code{\link{baseCountsFromBamList}} for
#' experiments with one BAM file per cell (e.g. Smart-seq2). Every
#' file becomes one cell. The result has the same format as the
#' output of \code{\link{bam2R_v2}} and can be passed directly to
#' \code{\link{mutationCallsFromExclusionlist}},
#' \code{\link{mutationCallsFromCohort}},
#' \code{\link{pullcountsVars}} and \code{\link{strandBias10x}}.
#'
#' In contrast to \code{\link{baseCountsFromBamList}} any number of
#' regions on any chromosome can be read, the counts of both strands
#' are kept separately, and the deepSNV package is not needed. With
#' the default parameters the counts are identical to those of
#' \code{\link{baseCountsFromBamList}} (after merging the strands).
#'
#' @param bamfiles A character vector (or list) with the paths of the
#'BAM files, one per cell. The names of the vector are used as cell
#'names, unnamed files are named after the file (without .bam).
#' @param sites The regions of interest. A GRanges or GRangesList
#'object with any number of ranges on any chromosome, alternatively
#'a character vector that can be converted to one
#'(e.g. c("chrM:1-16569","chr7:140753336")). Overlapping ranges
#'are merged. Default: the entire mitochondrial genome.
#' @param q An optional cutoff for the nucleotide Phred
#'quality. Default q = 25. Nucleotides with Q <= q will be masked
#'by 'N' (only nucleotides with Q > q are counted).
#' @param mq An optional cutoff for the read mapping quality. Default
#'mq = 0 (no filter). reads with MQ < mq will be discarded.
#' @param head.clip Number of nucleotides at the head (5' end) of
#'each read that are masked by 'N'. Default 0.
#' @param max.depth The maximal depth for the pileup command. Default
#'1,000,000.
#' @param mask Integer indicating which flags to filter. Default 0 (no
#'mask). Try 1796 (BAM_DEF_MASK).
#' @param keepflag Integer indicating which flags to keep. Default 0
#'(no mask). Try 3 (PAIRED|PROPERLY_PAIRED).
#' @param max.mismatches Integer indicating maximum MN value to allow
#'in a read. Default NULL (no filter).
#' @param umi.tag The BAM tag that contains the UMI. If provided UMIs
#'are counted instead of reads, see \code{\link{bam2R_v2}}. Default
#'NULL (count reads).
#' @param ncores Number of files that are read in parallel. Default 1.
#' @return A named \code{\link{list}} of 16 sparse matrices (class
#'dgCMatrix) with rows corresponding to covered genomic sites
#'(named chromosome:position) and one column per BAM file, see
#'\code{\link{bam2R_v2}}.
#' @examples
#'## one BAM file per cell (here: the same file twice)
#'bam <- system.file("extdata", "mm10_10x.bam", package = "mitoClone2")
#'counts <- baseCountsFromBamList_v2(c(cell1 = bam, cell2 = bam),
#'sites = "chrM:1-15000")
#'dim(counts$A)
#' @author Benjamin Story
#' @export baseCountsFromBamList_v2
baseCountsFromBamList_v2 <- function(bamfiles,
                                     sites = "chrM:1-16569",
                                     q = 25,
                                     mq = 0,
                                     head.clip = 0,
                                     max.depth = 1000000,
                                     mask = 0,
                                     keepflag = 0,
                                     max.mismatches = NULL,
                                     umi.tag = NULL,
                                     ncores = 1) {
    cells <- names(bamfiles)
    bamfiles <- as.character(unlist(bamfiles, use.names = FALSE))
    if (length(bamfiles) == 0) {
        stop('No BAM files provided.')
    }
    if (is.null(cells)) {
        cells <- make.unique(sub("\\.(bam|cram)$", "", basename(bamfiles), ignore.case = TRUE))
    } else if (anyDuplicated(cells) || any(is.na(cells) | cells == '')) {
        stop('The names of bamfiles are used as cell names and have to be unique and not empty.')
    }
    if (methods::is(sites, 'GRangesList')) {
        sites <- unlist(sites)
    }
    sites <-
        GenomicRanges::reduce(GenomicRanges::GRanges(sites), ignore.strand = TRUE)
    if (is.null(max.mismatches))
        max.mismatches <- -1
    if (!is.null(umi.tag) &&
        (!is.character(umi.tag) || length(umi.tag) != 1 || nchar(umi.tag) != 2)) {
        stop("umi.tag must be NULL or a two character BAM tag, e.g. 'UB'")
    }
    if (is.null(umi.tag))
        umi.tag <- ''
    ## all regions in one pass per file, the files are read in parallel
    tiled <- tileSites10x(sites, 1)
    result <- mclapply10x(seq_along(bamfiles), function(i) {
        ## no barcode tag: the whole file is one cell
        res <- countRegions10x(bamfiles[i], tiled$tiles, 'cell', q, mq, head.clip, max.depth,
                               FALSE, mask, keepflag, max.mismatches, '', umi.tag)
        res$barcode <- lapply(res$barcode, function(b) rep(i, length(b)))
        res
    }, ncores)
    result <- assembleSparse10x(result, tiled$tiles, cells, ncores)
    reads <- Reduce('+', lapply(result[c('A', 'T', 'C', 'G', 'N', 'a', 't', 'c', 'g', 'n')],
                                Matrix::colSums))
    if (any(reads == 0)) {
        warning(sum(reads == 0), ' of ', length(reads), ' file(s) have no reads in the sites: ',
                paste(utils::head(cells[reads == 0], 5), collapse = ', '),
                if (sum(reads == 0) > 5) ', ...' else '',
                '. Check the sites (chromosome names)',
                if (umi.tag != '') paste0(' and the umi.tag ', umi.tag) else '', '.',
                call. = FALSE)
    }
    return(result)
}

#' Read nucleotide counts from a 10x Genomics .bam file (version 2, dense)
#'
#' Drop-in replacement for \code{\link{bam2R_10x}} that is built on
#' \code{\link{bam2R_v2}}: same arguments, same output (one dense
#' matrix per cell, strands merged), but the file is read with the
#' version 2 reader.
#'
#' The dense output needs (stop - start + 1) * 8 integers for every
#' single cell, if you have many cells or are only interested in a
#' set of variants use \code{\link{bam2R_v2}} together with
#' \code{\link{pullcountsVars10x}} instead.
#'
#' @inheritParams bam2R_v2
#' @param sites The chromosome locations of interest in BED format as
#'a string. Alternatively a single GRanges object will also work.
#' @param s Optional choice of the strand. Defaults to s = 2 (both).
#'Kept for compatibility with \code{\link{bam2R_10x}}, the strands
#'are always merged.
#' @param ignore_nonstandard Boolean indicating whether or not gapped
#'alignments, insertions, or deletions should be included in the
#'final output. Default FALSE. If you have an inflation of
#'spliced mitochondrial reads it is recommended to set this to
#'TRUE.
#' @return A named \code{\link{list}} of \code{\link{matrix}}, one
#'per cell barcode, with (stop - start + 1) rows and the columns
#'A, T, C, G, -, N, INS, DEL (A, T, C, G, N if ignore_nonstandard
#'= TRUE). Identical in format to the output of
#'\code{\link{bam2R_10x}} and \code{\link{baseCountsFromBamList}}.
#'Only barcodes with at least one count are returned.
#' @examples bamCounts <- bam2R_10x_v2(file = system.file("extdata",
#'"mm10_10x.bam", package="mitoClone2"), sites="chrM:1-15000")
#' @author Benjamin Story (adapted from original code with permission
#'from Moritz Gerstung)
#' @export bam2R_10x_v2
bam2R_10x_v2 <- function(file,
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
                         ignore_nonstandard = FALSE,
                         min_reads_per_barcode = 50,
                         barcodes = NULL,
                         barcode.tag = 'CB',
                         umi.tag = NULL) {
    if (!length(sites) == 1) {
        stop('Your sites parameter must be a character/GRanges object of length 1')
    }
    mito.chr <- GenomicRanges::GRanges(sites)
    ## not called 'stop': that would hide the function stop()
    start.pos <- GenomicRanges::start(mito.chr)
    end.pos <- GenomicRanges::end(mito.chr)
    result <- bam2R_v2(
        file = file,
        sites = mito.chr,
        barcodes = barcodes,
        q = q,
        mq = mq,
        head.clip = head.clip,
        max.depth = max.depth,
        verbose = verbose,
        mask = mask,
        keepflag = keepflag,
        max.mismatches = max.mismatches,
        ncores = ncores,
        min_reads_per_barcode = min_reads_per_barcode,
        barcode.tag = barcode.tag,
        umi.tag = umi.tag
    )
    ## merge the two strands
    nt.names <- c("A", "T", "C", "G", "-", "N", "INS", "DEL")
    result <-
        Map('+', result[nt.names], result[c("a", "t", "c", "g", "_", "n", "ins", "del")])
    ## only barcodes with at least one count are returned
    barcode.n <-
        colnames(result[[1]])[Reduce('+', lapply(result, Matrix::colSums)) > 0]
    if (length(barcode.n) == 0) {
        return(list())
    }
    if (ignore_nonstandard) {
        nt.names <- c('A', 'T', 'C', 'G', 'N')
    }
    ## rows of the dense matrix that are covered
    pos <- as.integer(gsub("^.*:", "", rownames(result[[1]]))) - start.pos + 1L
    result <- result[nt.names]
    ## the dense matrices are created for blocks of 100 cells at a time
    blocks <- split(barcode.n, ceiling(seq_along(barcode.n) / 100))
    result <- mclapply10x(blocks, function(cells) {
        ## array of covered positions x cells x nucleotides
        block <-
            sapply(result, function(x) as.matrix(x[, cells, drop = FALSE]), simplify = 'array')
        block <- array(block, dim = c(length(pos), length(cells), length(nt.names)))
        lapply(seq_along(cells), function(cell) {
            bam.file <-
                matrix(0L,
                       nrow = end.pos - start.pos + 1,
                       ncol = length(nt.names),
                       dimnames = list(NULL, nt.names))
            bam.file[pos, ] <- as.integer(block[, cell, ])
            return(bam.file)
        })
    }, ncores)
    result <- unlist(result, recursive = FALSE, use.names = FALSE)
    names(result) <- barcode.n
    return(result)
}

## Internal: split variants into the site (chromosome:position) and
## the reference / alternate allele. Accepts strings such as "1434
## G>A", "chr7:140753336 A>T" or a GRanges object with the metadata
## columns ref and alt (see mut2GR).
vars2sites10x <- function(SparseCounts, vars) {
    if (is.character(vars)) {
        bad <- !grepl("^(.+:)?[0-9]+ *[A-Za-z]+>[A-Za-z]+$", vars)
        if (any(bad)) {
            stop('Variants must be given as position reference>alternate, e.g. 1434 G>A or',
                 ' chr7:140753336 A>T. Not understood: ',
                 paste(utils::head(vars[bad], 3), collapse = ', '))
        }
        ## variants without a chromosome are on the one chromosome that was
        ## counted, or on the mitochondrial chromosome if there are several
        nochr <- !grepl(":", vars)
        chr <- NA
        if (any(nochr)) {
            chr <- mitoChr10x(SparseCounts)
            if (is.na(chr)) {
                stop(
                    paste0(
                        'Your counts cover several chromosomes and none of them is',
                        ' named like a mitochondrial chromosome.',
                        ' Please provide variants in the following format: chrM:1337 G>A'
                    )
                )
            }
        }
        var.df <- data.frame(name = vars, stringsAsFactors = FALSE)
        vars <- ifelse(nochr, paste0(chr, ":", vars), vars)
        alleles <- strsplit(gsub("^.*[0-9] *", "", vars), ">")
        var.df$site <- gsub(" *[A-Za-z]+>[A-Za-z]+$", "", vars)
        var.df$ref <- sapply(alleles, '[', 1)
        var.df$alt <- sapply(alleles, '[', 2)
    } else {
        var.df <- data.frame(
            site = paste0(GenomicRanges::seqnames(vars), ":", GenomicRanges::start(vars)),
            ref = as.character(vars$ref),
            alt = as.character(vars$alt),
            stringsAsFactors = FALSE
        )
        var.df$name <- paste0(var.df$site, " ", var.df$ref, ">", var.df$alt)
    }
    ## lower case letters are reverse strand columns in the counts, alleles are always upper case
    var.df$ref <- toupper(var.df$ref)
    var.df$alt <- toupper(var.df$alt)
    if (!all(c(var.df$ref, var.df$alt) %in% c('A', 'T', 'C', 'G'))) {
        stop('Only single nucleotide variants (A/T/C/G) are supported.')
    }
    return(var.df)
}

## Internal: subset the cells of the sparse counts, with a clear error
## for cells that are not in the counts
subsetcells10x <- function(SparseCounts, cells) {
    if (is.null(cells)) {
        return(SparseCounts)
    }
    missing <- !(cells %in% colnames(SparseCounts[[1]]))
    if (any(missing)) {
        stop(
            paste0(
                'You are attempting to subset for ',
                sum(missing),
                ' cell(s) that do not appear in your input SparseCounts!'
            )
        )
    }
    lapply(SparseCounts, function(x) x[, cells, drop = FALSE])
}

## Internal: pull the counts of one allele per variant (variants x
## cells). No looping over variants: a selection matrix that links
## each variant to its site is multiplied with the count matrix of
## each nucleotide.
pullallele10x <- function(SparseCounts, site, allele, strand = c('both', 'forward', 'reverse')) {
    strand <- match.arg(strand)
    hit <- match(site, rownames(SparseCounts[[1]]))
    out <- Reduce('+', lapply(c('A', 'T', 'C', 'G'), function(nt) {
        ## variants without coverage have no site and remain 0
        use <- which(allele == nt & !is.na(hit))
        sel <-
            Matrix::sparseMatrix(
                        i = use,
                        j = hit[use],
                        x = 1,
                        dims = c(length(site), nrow(SparseCounts[[1]]))
                    )
        nt.counts <- switch(
            strand,
            "both" = SparseCounts[[nt]] + SparseCounts[[tolower(nt)]],
            "forward" = SparseCounts[[nt]],
            "reverse" = SparseCounts[[tolower(nt)]]
        )
        sel %*% nt.counts
    }))
    return(out)
}

#'Pull variant counts (sparse 10x counts)
#'
#'Equivalent of \code{\link{pullcountsVars}} for the output of
#' \code{\link{bam2R_v2}}. All variants are pulled at once,
#' they can be located on any chromosome.
#'@param SparseCounts A list of sparse count matrices as produced by
#'\code{\link{bam2R_v2}}
#'@param vars Character vector of variants to pull, in format 5643G>T
#'or chr7:140753336 A>T. The chromosome may only be omitted if
#'all counts are from the same chromosome. Alternatively a GRanges
#'object with the metadata columns ref and alt.
#'@param cells Character vector for cells to select, or NULL if all
#'cells from the input are to be used
#'@return A list with two entries, M (count table on the variant
#'allele) and N (count table on the reference allele). The result
#'can be passed to \code{\link{mutationCallsFromMatrix}} as
#'\code{mutationCallsFromMatrix(t(MN$M), t(MN$N))}.
#'@examples sparseCounts <- bam2R_v2(file = system.file("extdata",
#'"mm10_10x.bam", package="mitoClone2"), sites="chrM:1-15000")
#'counts.known.vars <- pullcountsVars10x(sparseCounts, vars=c("13575 T>C","13778 A>T"))
#'@export
pullcountsVars10x <- function(SparseCounts, vars, cells = NULL) {
    var.df <- vars2sites10x(SparseCounts, vars)
    ## subset cells of interest if appropriate
    if (!is.null(cells)) {
        message('Subsetting for specific cells...')
    }
    SparseCounts <- subsetcells10x(SparseCounts, cells)
    ## pull counts alt
    M <- pullallele10x(SparseCounts, var.df$site, var.df$alt)
    ## pull total counts per position
    N <- Reduce('+', lapply(c('A', 'T', 'C', 'G'), function(nt) {
        pullallele10x(SparseCounts, var.df$site, rep(nt, nrow(var.df)))
    }))
    ## maintain consistency with previous code for N by excluding the mutant allele calls
    N <- N - M
    M <- as.matrix(M)
    N <- as.matrix(N)
    rownames(M) <- row.names(N) <- var.df$name
    colnames(M) <- colnames(N) <- colnames(SparseCounts[[1]])
    return(list(M = M, N = N))
}

#'Strand bias of variants (sparse 10x counts)
#'
#'Tabulates how many reads support the reference and the alternate
#' allele of each variant on the forward and on the reverse strand
#' and checks in two ways whether a variant is supported by both
#' strands. POOLED: the reads of all cells are summed and it is tested
#' whether the alternate allele is distributed differently across the
#' two strands than the reference allele. Variants that are almost
#' exclusively seen on one strand, while the reference allele is seen
#' on both, are often sequencing or alignment artifacts. PER CELL:
#' for a true variant the cells with many alternate reads on the
#' forward strand should also be the cells with many alternate reads
#' on the reverse strand (strand.cor), and there should be cells in
#' which the variant is seen on both strands (n.cells.conf). NOTE: 10x
#' libraries are strand specific, therefore many sites are covered by
#' one strand only. For these sites both alleles are on the same
#' strand, the variant will not be flagged and strand.cor is NA.
#'@param SparseCounts A list of sparse count matrices as produced by
#'\code{\link{bam2R_v2}}
#'@param vars Character vector of variants, in format 5643G>T or
#'chr7:140753336 A>T. Alternatively a GRanges object with the
#'metadata columns ref and alt.
#'@param cells Character vector for cells to select, or NULL if all
#'cells from the input are to be used
#'@param max.p Variants with a Fisher's exact test p-value below this
#'value are flagged. Default 0.001.
#'@param min.sor Variants additionally need a strand odds ratio of at
#'least this value to be flagged, i.e. how many times more skewed
#'towards one strand is the alternate allele in comparison to the
#'reference allele. Default 4.
#'@param min.reads.conf Minimum number of alternate reads a cell needs
#'on the forward AND on the reverse strand to be counted as a cell
#'in which the variant was confidently detected. Default 2.
#'@return A data.frame with one row per variant containing the read
#'counts of the reference and alternate allele on the forward and
#'reverse strand, the fraction of reads on the forward strand for
#'each allele, the strand odds ratio (sor, always >= 1), the
#'p-value of a Fisher's exact test and the strand.flag (TRUE =
#'variant shows strand bias in the pooled reads). The per cell
#'columns are: n.cells.alt (cells with at least one alternate
#'read), n.cells.conf (cells with at least min.reads.conf
#'alternate reads on both strands), strand.cor (Pearson
#'correlation between the forward and the reverse alternate
#'counts across the cells with at least one alternate read),
#'mean.af (mean allele frequency of the covered cells) and vmr
#'(variance of the allele frequency across the covered cells
#'divided by its mean, high for variants that are restricted to a
#'clone and low for variants/artifacts that are present at a
#'similar level in all cells).
#'@examples sparseCounts <- bam2R_v2(file = system.file("extdata",
#'"mm10_10x.bam", package="mitoClone2"), sites="chrM:1-15000")
#'strand.bias <- strandBias10x(sparseCounts, vars=c("13575 T>C","13778 A>T"))
#'@export
strandBias10x <- function(SparseCounts,
                          vars,
                          cells = NULL,
                          max.p = 0.001,
                          min.sor = 4,
                          min.reads.conf = 2) {
    var.df <- vars2sites10x(SparseCounts, vars)
    ## the variants are the row names of the result, so they have to be unique
    if (anyDuplicated(var.df$name)) {
        warning('Duplicated variants are only reported once.')
        var.df <- var.df[!duplicated(var.df$name), , drop = FALSE]
    }
    SparseCounts <- subsetcells10x(SparseCounts, cells)
    ## reads per allele and strand: variants x cells
    ref.fwd <- pullallele10x(SparseCounts, var.df$site, var.df$ref, 'forward')
    ref.rev <- pullallele10x(SparseCounts, var.df$site, var.df$ref, 'reverse')
    alt.fwd <- pullallele10x(SparseCounts, var.df$site, var.df$alt, 'forward')
    alt.rev <- pullallele10x(SparseCounts, var.df$site, var.df$alt, 'reverse')
    ## POOLED: summed over all cells
    out <- data.frame(
        row.names = var.df$name,
        ref.fwd = Matrix::rowSums(ref.fwd),
        ref.rev = Matrix::rowSums(ref.rev),
        alt.fwd = Matrix::rowSums(alt.fwd),
        alt.rev = Matrix::rowSums(alt.rev)
    )
    out$ref.frac.fwd <- out$ref.fwd / (out$ref.fwd + out$ref.rev)
    out$alt.frac.fwd <- out$alt.fwd / (out$alt.fwd + out$alt.rev)
    ## strand odds ratio with a pseudocount, symmetric so that it is always >= 1
    sor <-
        ((out$ref.fwd + 1) * (out$alt.rev + 1)) / ((out$ref.rev + 1) * (out$alt.fwd + 1))
    out$sor <- pmax(sor, 1 / sor)
    out$p.value <- apply(out[, c('ref.fwd', 'ref.rev', 'alt.fwd', 'alt.rev')], 1, function(x) {
        stats::fisher.test(matrix(x, nrow = 2))$p.value
    })
    out$strand.flag <- out$p.value < max.p & out$sor >= min.sor
    ## PER CELL: in how many cells is the variant seen (on both strands)
    out$n.cells.alt <- Matrix::rowSums((alt.fwd + alt.rev) > 0)
    out$n.cells.conf <-
        Matrix::rowSums((alt.fwd >= min.reads.conf) & (alt.rev >= min.reads.conf))
    ## correlation of the forward and reverse alternate counts across the
    ## cells with an alternate read. Cells without one are 0 on both
    ## strands and do not change the sums, so only n has to be adjusted
    n <- out$n.cells.alt
    cov.fr <- n * Matrix::rowSums(alt.fwd * alt.rev) - out$alt.fwd * out$alt.rev
    var.f <- n * Matrix::rowSums(alt.fwd ^ 2) - out$alt.fwd ^ 2
    var.r <- n * Matrix::rowSums(alt.rev ^ 2) - out$alt.rev ^ 2
    out$strand.cor <-
        ifelse(n >= 3 & var.f > 0 & var.r > 0, cov.fr / sqrt(var.f * var.r), NA)
    ## allele frequency per cell (NA if the cell is not covered)
    alt <- as.matrix(alt.fwd + alt.rev)
    af <- alt / (alt + as.matrix(ref.fwd + ref.rev))
    out$mean.af <- rowMeans(af, na.rm = TRUE)
    out$vmr <- apply(af, 1, stats::var, na.rm = TRUE) / out$mean.af
    return(out)
}

## ---------------------------------------------------------------------
## INTEGRATION with the mutation calling of the package
##
## mutationCallsFromExclusionlist(), mutationCallsFromCohort() and
## pullcountsVars() check their input with isSparseCounts10x() and hand
## the sparse counts of bam2R_v2() to the functions below. They apply
## exactly the same filters as the code for the dense counts, but never
## create a dense table of all sites x all cells: the cells are
## processed in blocks.
## ---------------------------------------------------------------------

## Internal: is this the output of bam2R_v2 (and not a list with one
## matrix per cell)?
isSparseCounts10x <- function(x) {
    is.list(x) && length(x) > 0 &&
        all(c('A', 'T', 'C', 'G', 'N', 'a', 't', 'c', 'g', 'n') %in% names(x)) &&
        all(vapply(x, function(m) methods::is(m, 'Matrix'), logical(1)))
}

## Internal: the chromosome that variants without a chromosome
## ("1337 G>A") refer to. The only chromosome of the counts, or the one
## that is named like a mitochondrial chromosome. NA if there is none.
mitoChr10x <- function(SparseCounts) {
    chr <- unique(sub(":[0-9]+$", "", rownames(SparseCounts[[1]])))
    if (length(chr) > 1) {
        chr <- chr[toupper(sub("^chr", "", chr, ignore.case = TRUE)) %in% c('M', 'MT')]
    }
    if (length(chr) != 1) {
        return(NA_character_)
    }
    return(chr)
}

## Internal: counts of the mitochondrial chromosome with both strands
## merged, one sparse matrix (sites x cells) per nucleotide. The sites
## are sorted by position. If a reference sequence is provided, sites
## outside of it are dropped.
mitoCounts10x <- function(SparseCounts, nts, mito.dna = NULL) {
    site <- rownames(SparseCounts[[1]])
    chr <- mitoChr10x(SparseCounts)
    if (length(site) > 0 && is.na(chr)) {
        stop(
            paste0(
                'Your counts cover several chromosomes and none of them is named like',
                ' a mitochondrial chromosome (chrM, MT). Please subset the counts first.'
            )
        )
    }
    keep <- which(sub(":[0-9]+$", "", site) == chr)
    if (length(keep) < length(site)) {
        message('Using the counts on ', chr, ' (', length(keep), ' of ', length(site), ' sites).')
    }
    pos <- as.integer(sub("^.*:", "", site[keep]))
    if (!is.null(mito.dna)) {
        inside <- pos <= length(mito.dna)
        if (!all(inside)) {
            warning(sum(!inside), ' sites are located behind the end of the reference',
                    ' sequence (', length(mito.dna), ' bp) and are ignored.',
                    ' Is the genome parameter correct?')
        }
        keep <- keep[inside]
        pos <- pos[inside]
    }
    keep <- keep[order(pos)]
    counts <- lapply(nts, function(nt) {
        SparseCounts[[nt]][keep, , drop = FALSE] + SparseCounts[[tolower(nt)]][keep, , drop = FALSE]
    })
    names(counts) <- nts
    return(list(chr = chr, pos = sort(pos), counts = counts))
}

## Internal: stop with a clear message if the counts have no cells or no sites
checkNotEmpty10x <- function(n.cells, n.sites) {
    if (n.cells == 0) {
        stop('The counts do not contain any cell.', call. = FALSE)
    }
    if (n.sites == 0) {
        stop('The counts do not contain any covered site on the mitochondrial chromosome',
             ' (or on any chromosome with all.chromosomes = TRUE).', call. = FALSE)
    }
}

## Internal: mutationCallsFromCohort (dense and sparse), no site has an
## alternate allele in at least MINCELL cells
noCandidates10x <- function(MINCELL, MINREADS) {
    stop('No variant is supported by at least MINCELL = ', MINCELL, ' cells with at least',
         ' MINREADS = ', MINREADS, ' reads on the variant. Try a lower MINCELL or MINREADS.',
         call. = FALSE)
}

## Internal: mutationCallsFromCohort (dense and sparse), warn if no
## candidate variant is specific for a single patient
noPatientVariants10x <- function(n.single, n.candidates, n.shared) {
    if (n.single == 0) {
        warning('None of the ', n.candidates, ' candidate variants is specific for a single',
                ' patient (', n.shared, ' are shared between patients and are on the',
                ' exclusionlist). Try a lower MINCELLS.PATIENT or MINFRAC.PATIENT.',
                call. = FALSE)
    }
}

## Internal: split the cells into blocks, so that a dense table of
## sites x cells of a block has at most 5e6 entries (40 Mb)
cellBlocks10x <- function(n.cells, n.sites) {
    per.block <- max(1, floor(5e6 / max(1, n.sites)))
    split(seq_len(n.cells), ceiling(seq_len(n.cells) / per.block))
}

## Internal: mutationCallsFromExclusionlist for sparse counts. See
## there for the parameters, mito.dna is the reference sequence.
exclusionlist10x <- function(SparseCounts,
                             lim.cov,
                             min.af,
                             min.num.samples,
                             min.af.universal,
                             universal.var.cells,
                             exclusionlists.use,
                             max.var.na,
                             max.cell.na,
                             mito.dna,
                             ncores,
                             ...) {
    nts <- c('A', 'T', 'C', 'G')
    mito <- mitoCounts10x(SparseCounts, nts, mito.dna)
    cells <- colnames(SparseCounts[[1]])
    n.sites <- length(mito$pos)
    checkNotEmpty10x(length(cells), n.sites)
    blocks <- cellBlocks10x(length(cells), n.sites)
    ## a cell is covered at a site if it has at least lim.cov reads
    is.covered <- function(x) {
        cov <- Reduce('+', x)
        cov >= lim.cov & cov > 0
    }
    ## PASS 1: per site and nucleotide in how many cells is the allele
    ## frequency above the cutoffs (4 columns per statistic)
    stats <- mclapply10x(blocks, function(block) {
        x <- lapply(mito$counts, function(m) as.matrix(m[, block, drop = FALSE]))
        cov <- Reduce('+', x)
        covered <- is.covered(x)
        af <- lapply(x, function(m) {
            m <- m / cov
            m[!covered] <- NA
            m
        })
        list(
            n.mut = vapply(af, function(m) rowSums(m > min.af, na.rm = TRUE), numeric(n.sites)),
            n.any = vapply(af, function(m) rowSums(m > 0, na.rm = TRUE), numeric(n.sites)),
            n.universal = vapply(af, function(m) rowSums(m >= min.af.universal, na.rm = TRUE), numeric(n.sites)),
            n.covered = rowSums(covered)
        )
    }, ncores)
    total <- function(what) {
        matrix(Reduce('+', lapply(stats, '[[', what)), nrow = n.sites)
    }
    ## one row per possible variant, same order as for the dense counts:
    ## by alternate allele (A, T, C, G), then by position
    ref <- mito.dna[mito$pos]
    vars <- data.frame(
        site = rep(seq_len(n.sites), length(nts)),
        alt = rep(nts, each = n.sites),
        n.mut = as.vector(total('n.mut')),
        n.any = as.vector(total('n.any')),
        n.universal = as.vector(total('n.universal')),
        n.covered = rep(as.vector(total('n.covered')), length(nts)),
        stringsAsFactors = FALSE
    )
    vars$ref <- ref[vars$site]
    vars$name <- paste0(mito$pos[vars$site], ' ', vars$ref, '>', vars$alt)
    ## no reference allele, no sites where the reference is unknown
    vars <- vars[!is.na(vars$ref) & vars$ref != 'N' & vars$ref != vars$alt, , drop = FALSE]
    ## number of variants left after each filter, reported if none is left
    left <- c(candidates = nrow(vars))
    ## part 1 - enough cells with the variant
    vars <- vars[vars$n.mut >= min.num.samples, , drop = FALSE]
    left['min.af/min.num.samples'] <- nrow(vars)
    ## part 2 - filter based on the exclusionlist
    is.names <- vapply(exclusionlists.use, function(x) typeof(x) == "character", logical(1))
    if (sum(is.names) > 0) {
        removal.names.list <- unique(unlist(exclusionlists.use[is.names]))
        vars <- vars[!vars$name %in% removal.names.list, , drop = FALSE]
    }
    if (sum(!is.names) > 0 && nrow(vars) > 0) {
        removal.ranges.list <-
            unique(unlist(GenomicRanges::GRangesList(exclusionlists.use[!is.names])))
        ## the exclusionlists are located on chrM
        var.gr <- GenomicRanges::GRanges(paste0('chrM:', mito$pos[vars$site]))
        hits <- S4Vectors::queryHits(GenomicRanges::findOverlaps(var.gr, removal.ranges.list))
        vars <- vars[!seq_len(nrow(vars)) %in% hits, , drop = FALSE]
    }
    left['exclusionlists'] <- nrow(vars)
    ## remove variants that were never seen in a covered cell
    vars <- vars[vars$n.any > 0, , drop = FALSE]
    ## remove variants that are at a certain universal af in a certain number of cells
    vars <- vars[!vars$n.universal >= universal.var.cells, , drop = FALSE]
    left['universal.var.cells'] <- nrow(vars)
    ## vars must have less than X % NA's
    vars <- vars[(length(cells) - vars$n.covered) < max.var.na * length(cells), , drop = FALSE]
    left['max.var.na'] <- nrow(vars)
    if (nrow(vars) == 0) {
        stop('No variant passed the filters. Variants left after each filter: ',
             paste0(names(left), ' = ', left, collapse = ', '),
             '. Relax the filter that removes them (a site only counts as covered',
             ' in a cell with at least lim.cov = ', lim.cov, ' reads).',
             call. = FALSE)
    }
    ## PASS 2: cells must have less than X % NA's in the remaining variants
    var.site <- unique(vars$site)
    n.na <- unlist(mclapply10x(blocks, function(block) {
        covered <- is.covered(lapply(mito$counts, function(m) {
            as.matrix(m[var.site, block, drop = FALSE])
        }))
        colSums(!covered[match(vars$site, var.site), , drop = FALSE])
    }, ncores), use.names = FALSE)
    cells <- cells[n.na < max.cell.na * nrow(vars)]
    if (length(cells) == 0) {
        stop('No cell passed the filters. Try a higher max.cell.na.', call. = FALSE)
    }
    MN <- suppressMessages(pullcountsVars10x(SparseCounts, paste0(mito$chr, ':', vars$name), cells))
    rownames(MN$M) <- rownames(MN$N) <- vars$name
    mutationCallsFromMatrix(t(MN$M), t(MN$N), ...)
}

## Internal: counts of ALL sites with both strands merged, one sparse
## matrix (sites x cells) per nucleotide. Same format as mitoCounts10x,
## but chr has one entry per site.
allCounts10x <- function(SparseCounts, nts) {
    site <- rownames(SparseCounts[[1]])
    counts <- lapply(nts, function(nt) {
        SparseCounts[[nt]] + SparseCounts[[tolower(nt)]]
    })
    names(counts) <- nts
    list(chr = sub(":[0-9]+$", "", site),
         pos = as.integer(sub("^.*:", "", site)),
         counts = counts)
}

## Internal: reference nucleotide of sites (chr, pos) from a FASTA file
## (path) or from an object with a getSeq method (BSgenome, FaFile).
## Returns ref (the nucleotide, NA if unknown) and known (is the
## sequence of the site part of the reference?).
fetchReference10x <- function(reference, chr, pos) {
    if (length(chr) == 0) {
        return(list(ref = character(0), known = logical(0)))
    }
    site.gr <- GenomicRanges::GRanges(seqnames = chr, ranges = paste0(pos, '-', pos))
    if (is.character(reference)) {
        if (length(reference) != 1 || !file.exists(reference)) {
            stop('reference must be the path of ONE existing FASTA file or a BSgenome object.')
        }
        ## nearby sites are fetched together: one stretch of sequence per region
        region <- GenomicRanges::reduce(site.gr, min.gapwidth = 1000)
        seqs <- .Call(
            "bam2R_v2_fasta",
            path.expand(reference),
            as.character(GenomicRanges::seqnames(region)),
            as.integer(GenomicRanges::start(region)),
            as.integer(GenomicRanges::end(region))
        )
        hit <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(site.gr, region))
        offset <- pos - GenomicRanges::start(region)[hit] + 1L
        ref <- substr(seqs[hit], offset, offset)
        known <- !is.na(seqs[hit])
        ## behind the end of the sequence
        ref[!is.na(ref) & ref == ''] <- NA
    } else {
        if (!requireNamespace('Biostrings', quietly = TRUE)) {
            stop('The package Biostrings is required to use a BSgenome object as reference.')
        }
        known <- tryCatch(GenomicRanges::seqnames(GenomicRanges::seqinfo(reference)),
                          error = function(e) NULL)
        if (is.null(known)) {
            stop('reference must be the path of a FASTA file or a BSgenome object.')
        }
        ref <- rep(NA_character_, length(chr))
        known <- chr %in% known
        if (any(known)) {
            ref[known] <- as.character(Biostrings::getSeq(reference, site.gr[known]))
        }
    }
    ## soft masked (lower case) sequence, everything but A/C/G/T is unknown
    ref <- toupper(ref)
    ref[!is.na(ref) & !ref %in% c('A', 'C', 'G', 'T')] <- 'N'
    return(list(ref = ref, known = known))
}

## Internal: reference nucleotide for the sites of ALL chromosomes
## (mutationCallsFromCohort with all.chromosomes = TRUE). sites: output
## of allCounts10x, mito.dna: built-in mitochondrial reference or NULL
## (USE.REFERENCE = FALSE), reference: FASTA / BSgenome or NULL.
reference10x <- function(sites, mito.dna, reference) {
    nts <- names(sites$counts)
    n.sites <- length(sites$pos)
    ## the most abundant nucleotide of each site
    abundant <- function() {
        totalntCounts <- vapply(sites$counts, Matrix::rowSums, numeric(n.sites))
        nts[apply(matrix(totalntCounts, nrow = n.sites), 1, which.max)]
    }
    chrs <- unique(sites$chr)
    mito.chr <- chrs[toupper(sub("^chr", "", chrs, ignore.case = TRUE)) %in% c('M', 'MT')]
    is.mito <- sites$chr %in% mito.chr
    ## the built-in mitochondrial reference (NA behind its end)
    builtin <- function(idx) {
        ref <- mito.dna[sites$pos[idx]]
        if (any(is.na(ref))) {
            warning(sum(is.na(ref)), ' mitochondrial sites are located behind the end of the',
                    ' reference sequence and are ignored. Is the genome parameter correct?',
                    call. = FALSE)
        }
        ref
    }
    if (!is.null(reference)) {
        fetched <- fetchReference10x(reference, sites$chr, sites$pos)
        ref <- fetched$ref
        ## the mitochondrial chromosome may be missing from the reference:
        ## use the built-in sequence
        fill <- !fetched$known & is.mito
        if (any(fill) && !is.null(mito.dna)) {
            ref[fill] <- builtin(which(fill))
            fetched$known[fill] <- TRUE
        }
        missing <- unique(sites$chr[!fetched$known])
        if (length(missing) > 0) {
            stop('The reference does not contain the following sequence(s) of your counts: ',
                 paste(utils::head(missing, 5), collapse = ', '),
                 '. The names have to be identical to those in the BAM file (chr7 is not 7).',
                 call. = FALSE)
        }
        return(ref)
    }
    ref <- abundant()
    if (is.null(mito.dna)) {
        message('USE.REFERENCE = FALSE: the most abundant nucleotide of every site is used as',
                ' the reference allele.')
        return(ref)
    }
    if (any(is.mito)) {
        ref[is.mito] <- builtin(which(is.mito))
    }
    if (any(!is.mito)) {
        warning(
            '\n',
            '#######################################################################\n',
            '## NO REFERENCE GENOME PROVIDED FOR THE NON-MITOCHONDRIAL SITES !!!\n',
            '## Only the mitochondrial reference is part of this package. For the\n',
            '## ', sum(!is.mito), ' site(s) on ', paste(utils::head(setdiff(chrs, mito.chr), 3), collapse = ', '),
            if (length(setdiff(chrs, mito.chr)) > 3) ', ...' else '', '\n',
            '## the MOST ABUNDANT nucleotide is used as the reference allele:\n',
            '##  - homozygous variants and variants carried by most cells are\n',
            '##    INVISIBLE (they look like the reference)\n',
            '##  - REF and ALT are SWAPPED where the variant is the majority\n',
            '## Provide the genome your BAM file was aligned to:\n',
            '##    reference = "/path/to/genome.fa"   (FASTA file)\n',
            '##    reference = BSgenome.Hsapiens.UCSC.hg38::Hsapiens   (BSgenome)\n',
            '#######################################################################',
            call. = FALSE
        )
    }
    return(ref)
}

## Internal: mutationCallsFromCohort for sparse counts. See there for
## the parameters, mito.dna is the mitochondrial reference sequence or
## NULL if the most abundant nucleotide of a site is the reference.
## all.chromosomes = FALSE: only the mitochondrial chromosome, as for
## the dense counts. TRUE: every site of the counts, the reference of
## the other chromosomes comes from reference (FASTA / BSgenome).
cohort10x <- function(SparseCounts,
                      patient,
                      MINREADS,
                      MINCELL,
                      MINFRAC,
                      MINCELLS.PATIENT,
                      MINFRAC.PATIENT,
                      MINFRAC.OTHER,
                      mito.dna,
                      all.chromosomes = FALSE,
                      reference = NULL,
                      ncores = 1) {
    cells <- colnames(SparseCounts[[1]])
    if (length(patient) != length(cells)) {
        stop('patient must have one entry per cell (column) of the counts: ',
             length(cells), ' cells, ', length(patient), ' patients.')
    }
    if (MINREADS < 1) {
        stop('MINREADS must be at least 1.')
    }
    nts <- c('A', 'T', 'C', 'G', 'N')
    if (all.chromosomes) {
        sites <- allCounts10x(SparseCounts, nts)
        n.sites <- length(sites$pos)
        reference <- reference10x(sites, mito.dna, reference)
        ## mitochondrial variants keep their name without a chromosome
        chrs <- unique(sites$chr)
        mito.chr <- chrs[toupper(sub("^chr", "", chrs, ignore.case = TRUE)) %in% c('M', 'MT')]
        site.label <- ifelse(sites$chr %in% mito.chr, sites$pos, paste0(sites$chr, ':', sites$pos))
    } else {
        sites <- mitoCounts10x(SparseCounts, nts, mito.dna)
        n.sites <- length(sites$pos)
        sites$chr <- rep(sites$chr, n.sites)
        if (is.null(mito.dna)) {
            ## the most abundant nucleotide is the reference
            totalntCounts <- vapply(sites$counts, Matrix::rowSums, numeric(n.sites))
            totalntCounts <- matrix(totalntCounts, nrow = n.sites)
            reference <- nts[apply(totalntCounts, 1, which.max)]
        } else {
            reference <- mito.dna[sites$pos]
        }
        site.label <- as.character(sites$pos)
    }
    checkNotEmpty10x(length(cells), n.sites)
    patients <- unique(patient)
    ## cells x patients
    patient.cells <- vapply(patients, function(pa) as.numeric(patient == pa), numeric(length(cells)))
    patient.cells <- matrix(patient.cells, nrow = length(cells))
    blocks <- cellBlocks10x(length(cells), n.sites)
    stats <- mclapply10x(blocks, function(block) {
        x <- lapply(sites$counts, function(m) as.matrix(m[, block, drop = FALSE]))
        total_cov_pos <- Reduce('+', x)
        x <- x[c('A', 'T', 'C', 'G')]
        list(
            ## cells with at least MINREADS reads on the allele
            support = vapply(x, function(m) rowSums(m >= MINREADS), numeric(n.sites)),
            ## mutant cells per patient: at least MINREADS reads and
            ## MINFRAC % of the total coverage at this position
            mutant = lapply(x, function(m) {
                (m >= MINREADS & m >= (MINFRAC * total_cov_pos)) %*%
                    patient.cells[block, , drop = FALSE]
            })
        )
    }, ncores)
    support <- matrix(Reduce('+', lapply(stats, '[[', 'support')), nrow = n.sites)
    ## one row per candidate variant: by site, then by allele (A, T, C, G)
    vars <- data.frame(
        site = rep(seq_len(n.sites), 4),
        alt = rep(c('A', 'T', 'C', 'G'), each = n.sites),
        support = as.vector(support),
        stringsAsFactors = FALSE
    )
    vars$ref <- reference[vars$site]
    varcount.bypatient <- do.call(rbind, lapply(c('A', 'T', 'C', 'G'), function(nt) {
        Reduce('+', lapply(stats, function(st) st$mutant[[nt]]))
    }))
    keep <- which(!is.na(vars$ref) & vars$ref != vars$alt & vars$support >= MINCELL)
    keep <- keep[order(vars$site[keep], match(vars$alt[keep], c('A', 'T', 'C', 'G')))]
    if (length(keep) == 0) {
        noCandidates10x(MINCELL, MINREADS)
    }
    vars <- vars[keep, , drop = FALSE]
    varcount.bypatient <- varcount.bypatient[keep, , drop = FALSE]
    alleles <- paste0(vars$ref, '>', vars$alt)
    vars$name <- paste0(site.label[vars$site], alleles)
    ## the name with a space is used on the exclusionlist, the full name to pull the counts
    vars$name.space <- paste0(site.label[vars$site], ' ', alleles)
    vars$name.full <- paste0(sites$chr[vars$site], ':', sites$pos[vars$site], ' ', alleles)
    dimnames(varcount.bypatient) <- list(vars$name, patients)
    ## total cells per patient
    patient.totalcells <- colSums(patient.cells)
    ## fraction of all cells containg said variant per patient
    varcount.bypatient.fraction <- t(t(varcount.bypatient) / patient.totalcells)
    ## find variants with at least MINCELLS.PATIENT in any patient AND at least MINFRAC.PATIENT % of the cells of that patient
    filter <- apply(varcount.bypatient, 1, function(ptid) {
        max(ptid) >= MINCELLS.PATIENT &
            max(ptid) / patient.totalcells[which.max(ptid)] >= MINFRAC.PATIENT
    })
    n.enriched <- apply(varcount.bypatient.fraction, 1, function(x) sum(x >= MINFRAC.OTHER * max(x)))
    n.patients <- rowSums(varcount.bypatient >= MINCELLS.PATIENT)
    if (nrow(vars) == 0) {
        filter <- n.enriched <- n.patients <- logical(0)
    }
    ## variants in a single patient / in more than one patient
    singlepatient <- filter & n.enriched == 1 & n.patients == 1
    multipatient <- filter & n.enriched > 1 & n.patients > 1
    noPatientVariants10x(sum(singlepatient), length(singlepatient), sum(multipatient))
    mutation.bypatient <-
        patients[apply(varcount.bypatient[singlepatient, , drop = FALSE], 1, which.max)]
    selected <- vars[singlepatient, , drop = FALSE]
    out <- lapply(patients, function(ptid) {
        if (sum(mutation.bypatient == ptid) == 0) {
            return(NULL)
        }
        use <- selected[mutation.bypatient == ptid, , drop = FALSE]
        MN <- suppressMessages(pullcountsVars10x(SparseCounts, use$name.full, cells[patient == ptid]))
        rownames(MN$M) <- rownames(MN$N) <- use$name
        mutationCallsFromMatrix(t(MN$M), t(MN$N))
    })
    names(out) <- patients
    ## mutations that are shared between patients go to the exclusionlist
    out$exclusionlist <- vars$name.space[multipatient]
    return(out)
}
