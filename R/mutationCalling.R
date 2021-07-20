#'Create a mutationCalls objects from nucleotide base calls and
#' defines a exclusionlist (cohort)
#'
#'Identifies relevant mitochondrial somatic variants from raw counts
#' of nucleotide frequencies measured in single cells from several
#' individuals. Applies two sets of filters: In the first step,
#' filters on coverage to include potentially noisy variants; in the
#' second step, compares allele frequencies between patients to remove
#' variants that were observed in several individuals and that
#' therefore are unlikely to represent true somatic variants (e.g. RNA
#' editing events). The exclusionlist derived from the original Velten
#' et al. 2021 dataset is available internal and can be used on single
#' individuals using \code{\link{mutationCallsFromExclusionlist}}
#'@param BaseCounts A list of base call matrices (one matrix per cell)
#'as produced by \code{\link{baseCountsFromBamList}} or
#'\code{\link{bam2R_10x}}.
#'@param patient A character vector associating each cell / entry in
#'the \code{BaseCount} list with a patient
#'@param sites Vector specifying genomic regions, defaults to the
#'entire mitochondrial genome. Excepts a string but may be
#'included as a GRanges object.
#'@param MINREADS Minimum number of reads on a site in a single cell
#'to qualify the site as covered
#'@param MINCELL Minimum number of cells across the whole data set to
#'cover a site
#'@param MINFRAC Fraction of reads on the mutant allele to
#'provisionally classify a cell as mutant
#'@param MINCELLS.PATIENT Minimum number of mutant cells per patient
#'to classify the mutation as relevant in that patient, AND
#'@param MINFRAC.PATIENT Minimum fraction of mutant cells per patient
#'to classify the mutation as relevant in that patient
#'@param MINFRAC.OTHER Minimum fraction of mutant cells identified in
#'a second patient for the mutation to be excluded. Fraction
#'relative to the fraction of of cells from the patient where a
#'variant is enriched.
#'@param USE.REFERENCE Boolean. The variant calls will be of the
#'format REF>ALT where REF is decided based on the selected
#'\code{genome} annotation. If set to FALSE, the reference allele
#'will be the most abundant.
#'@param genome The mitochondrial genome of the sample being
#'investigated. Please note that this is the UCSC standard
#'chromosome sequence. Default: hg38.
#'@return A list of \code{\link{mutationCalls}} objects (one for each
#'\code{patient}) and an entry named \code{exclusionlist}
#'containing a exclusionlist of sites with variants in several
#'individuals
#'@examples sites.gr <- GenomicRanges::GRanges("chrM:1-15000")
#'BaseCounts <- bam2R_10x(file = system.file("extdata",
#'"mm10_10x.bam", package="mitoClone2"), sites=sites.gr)
#'mutCalls <- mutationCallsFromCohort(BaseCounts,
#'patient=c('sample2','sample1','sample2','sample2','sample1','sample2'),
#'MINCELL=1, MINFRAC=0, MINCELLS.PATIENT=1, genome='mm10',
#'sites=sites.gr)
#'@export
mutationCallsFromCohort <- function(BaseCounts,
                                    sites,
                                    patient,
                                    MINREADS = 5,
                                    MINCELL = 20,
                                    MINFRAC = 0.1,
                                    MINCELLS.PATIENT = 10,
                                    MINFRAC.PATIENT = 0.01,
                                    MINFRAC.OTHER = 0.1,
                                    USE.REFERENCE = TRUE,
                                    genome = 'hg38') {
    message("Making sure 'sites' parameter is set correctly...")
    if (!length(sites) == 1) {
        stop('Your sites parameter must be a character vector or GRanges object of length 1')
    }
    sites <- GenomicRanges::GRanges(sites)
    ## read in the
    ntcountsArray <- simplify2array(BaseCounts)
    ntcountsArray <- aperm(ntcountsArray, c(1, 3, 2))
    ntcountsArray <- ntcountsArray[, , c('A', 'T', 'C', 'G', 'N')]
    if (USE.REFERENCE) {
        reference <-
            switch(genome,
                   "hg38" = hg38.dna,
                   "hg19" = hg19.dna,
                   "mm10" = mm10.dna)
        reference <-
            reference[GenomicRanges::start(sites):GenomicRanges::end(sites)]
        message(paste0(
            'Looks good. Using the UCSC ',
            genome,
            ' genome as a reference for variants.',
            ' Be wary, only mitochondrial annotations are available!'
        ))
    } else{
        message(paste0(
            'Looks good. However, the mutation names are run specific and may include N.'
        ))
        totalntCounts <- apply(ntcountsArray, c(1, 3), sum)
        reference <-
            colnames(totalntCounts)[apply(totalntCounts, 1, which.max)]
    }
    ## get total reads per cell
    total_cov_per_cell <- rowSums(colSums(ntcountsArray))
    
    ## make the mutation calls per var
    variant_calls <- lapply(seq_along(reference), function(pos) {
        ## check which position have at least MINREADS total reads in each cell
        support <- apply(ntcountsArray[pos, , ] >= MINREADS, 2, sum)
        ## check which position have at least MINREADS total reads in each cell
        ## remove reference calls and sites where unknown N calls dominate - perhapse modify later
        support <-
            support[!names(support) %in% c(reference[pos], "N", "INS", "DEL", "-")]
        ## only keep variants showing up in at least MINCELL
        candidates <- names(support)[support >= MINCELL]
        ## return null if no variants
        if (length(candidates) == 0) {
            NULL
        } else {
            total_cov_pos <- rowSums(ntcountsArray[pos, , ])
            out_matrix <- sapply(candidates, function(mutid) {
                ## mutations have at least MINREADS reads and exist as MINFRAC % of the total coverage at this nt position
                is_mut <-
                    ntcountsArray[pos, , mutid] >= MINREADS &
                    ntcountsArray[pos, , mutid] >= (MINFRAC * total_cov_pos)
                ## reference have at least MINREADS reads and exist as MINFRAC % of the total coverage at this nt position
                is_ref <-
                    ntcountsArray[pos, , reference[pos]] >= MINREADS &
                    ntcountsArray[pos, , reference[pos]] >= (MINFRAC * total_cov_pos)
                ## label cell depending on mutation and reference allele detection
                ifelse(is_mut,
                ifelse(is_ref, "BOTH", "MUT"),
                ifelse(is_ref, "WT", "DROP"))
            })
            colnames(out_matrix) <-
                paste0(pos, reference[pos], ">", candidates)
            return(out_matrix)
        }
    })
    ## merge the variant call matricies
    variant_calls <- do.call(cbind, variant_calls)
    ## count number of cell that are mutant per patient
    varcount.bypatient <- sapply(unique(patient), function(pa) {
        as.integer(colSums(variant_calls[patient == pa, ] == "BOTH" |
                           variant_calls[patient == pa, ] == "MUT"))
    })
    row.names(varcount.bypatient) <- colnames(variant_calls)
    ## total cells per patient
    patient.totalcells <-
        as.integer(table(patient)[colnames(varcount.bypatient)])
    names(patient.totalcells) <- colnames(varcount.bypatient)
    ## fraction of all cells containg said variant per patient
    varcount.bypatient.fraction <-
        t(t(varcount.bypatient) / patient.totalcells)
    ## find variants with at least MINCELLS.PATIENT in any patient AND at least MINFRAC.PATIENT % of all reads at a given position for that patient
    filter <-  apply(varcount.bypatient, 1, function(ptid) {
        max(ptid) >= MINCELLS.PATIENT &
            max(ptid) / patient.totalcells[which.max(ptid)] >= MINFRAC.PATIENT
    })
    ## identify enriched variants that appear in at least MINCELLS.PATIENT and at a fraction of over MINFRAC.OTHER % of the maximum observed  frequency of this variant in only a single patient
    ## CONSIDER CHANGING OR ADDING A SPECIFIC CUTOFF FRAC filter
    singlepatient <-
        filter &
        apply(varcount.bypatient.fraction, 1, function(x)
            sum(x >= MINFRAC.OTHER * max(x))) == 1 &
        rowSums(varcount.bypatient >= MINCELLS.PATIENT) == 1
    ## identify enriched variants that appear in at least MINCELLS.PATIENT and at a fraction of over MINFRAC.OTHER % of the maximum observed  frequency of this variant in more than one patient
    multipatient <-
        filter &
        apply(varcount.bypatient.fraction, 1, function(x)
            sum(x >= MINFRAC.OTHER * max(x))) > 1 &
        rowSums(varcount.bypatient >= MINCELLS.PATIENT) > 1
    ## extract mutations that are unique to individual patients
    mutation.bypatient <-
        colnames(varcount.bypatient)[apply(varcount.bypatient[singlepatient, ], 1, which.max)]
    ## pull the variant cells
    variant_calls_selected <- variant_calls[, singlepatient]
    ## extract the baseCounts for variants of interest
    out <- lapply(unique(patient), function(ptid) {
        ##retrieve matrices of allele counts for patient specific variants
        if (sum(mutation.bypatient == ptid) == 0) {
            return(NULL)
        }
        MN <-
            pullcountsVars(BaseCounts[patient == ptid], colnames(variant_calls_selected)[mutation.bypatient == ptid])
        ##create mutationCalls object
        o <- mutationCallsFromMatrix(t(MN$M), t(MN$N))
    })
    names(out) <- unique(patient)
    ## extract mutations that are shared in one or more patient for exclusionlist
    out$exclusionlist <- rownames(varcount.bypatient[multipatient, ])
    out$exclusionlist <-
        gsub("(\\d+)(.+)", "\\1 \\2", out$exclusionlist)
    return(out)
}

#'Create a mutationCalls object from nucleotide base calls using a
#' exclusionlist (single individual)
#'
#'Identifies relevant mitochondrial somatic variants from raw counts
#' of nucleotide frequencies. Applies two sets of filters: In the
#' first step, filters on coverage and minimum allele frequency to
#' exclude potentially noisy variants; in the second step, filters
#' against a exclusionlist of variants that were observed in several
#' individuals and that therefore are unlikely to represent true
#' somatic variants (e.g. RNA editing events). These exclusionlists
#' are created using \code{\link{mutationCallsFromCohort}}
#'@param BaseCounts A list of base call matrices (one matrix per cell)
#'as produced by \code{\link{baseCountsFromBamList}}
#'@param lim.cov Minimal coverage required per cell for a cell to be
#'classified as covered
#'@param min.af Minimal allele frequency for a cell to be classified
#'as mutant
#'@param min.num.samples Minimal number of cells required to be
#'classified as covered and mutant according to the thresholds
#'set in \code{lim.cov} and \code{min.af}. Usually specified as a
#'fraction of the total number of cells.
#'@param universal.var.cells Maximum number of cells required to be
#'classified as mutant according to the threshold set in
#'\code{min.af.universal}.  Usually specified as a fraction of
#'the total number of cells; serves to avoid e.g. germline
#'variants.
#'@param min.af.universal Minimal allele frequency for a cell to be
#'classified as mutant, in the context of removing universal
#'variants. Defaults to \code{min.af}, but can be set to lower
#'values.
#'@param exclusionlists.use List of sites to exclude for variants
#'calling. The default exclusionlists object included with this
#'package contains exclude or hardmask in GRanges format. The
#'four exclusionlists included in this case are: "three" (hg38
#'sites that are part of homopolymer(e.g. AAA) of at least 3 bp
#'in length), "mutaseq" (sites discovered to be overrepresented
#'in AML SmartSeq2 data analysis from Velten et al 2021),
#'"masked" (sites that are softmasked in either the UCSC or
#'Refseq genome annotations), and "rnaEDIT" which are sites that
#'are subjected to RNA-editing according to the REDIportal. These
#'lists can also be input manually by a researcher and provided
#'as either coordinates (as a string) or as a GRanges objects.
#'@param max.var.na Final filtering step: Remove all mutations with no
#'coverage in more than this fraction of cells
#'@param max.cell.na Final filtering step: Remove all cells with no
#'coverage in more than this fraction of mutations
#'@param genome The mitochondrial genome of the sample being
#'investigated. Please note that this is the UCSC standard
#'chromosome sequence. Default: hg38.
#'@param ncores number of cores to use for tabulating potential
#'variants (defaults to 2)
#'@param ... Parameters passed to
#'\code{\link{mutationCallsFromMatrix}}
#'@return An object of class \code{\link{mutationCalls}}
#'@examples load(system.file("extdata/example_counts.Rda",package = "mitoClone2"))
#'Example <- mutationCallsFromExclusionlist(example.counts,
#' min.af=0.05, min.num.samples=5,
#' universal.var.cells = 0.5 * length(example.counts),
#' binarize = 0.1)
#'@export
mutationCallsFromExclusionlist <- function(BaseCounts,
                                           lim.cov = 20,
                                           min.af = 0.2,
                                           min.num.samples = 0.01 * length(BaseCounts),
                                           min.af.universal = min.af,
                                           universal.var.cells = 0.95 * length(BaseCounts),
                                           exclusionlists.use = exclusionlists,
                                           max.var.na = 0.5,
                                           max.cell.na = 0.95,
                                           genome = 'hg38',
                                           ncores = 1,
                                           ...) {
    mito.dna <- switch(genome,
                       "hg38" = hg38.dna,
                       "hg19" = hg19.dna,
                       "mm10" = mm10.dna)
    varaf <- parallel::mclapply(BaseCounts, function(x) {
        ## focus on A,G,C,T
        x <- x[, c('A', 'T', 'C', 'G')]
        ## find cell that have less than 100 cov over agct at a given pos
        zeroes <- rowSums(x) < lim.cov
        ## af calc
        x.af <- x / rowSums(x)
        ## alleleRatio
        ##x.af <- x / (x+apply(x,1,max))
        x.af <- reshape2::melt(x.af)
        colnames(x.af) <- c('pos', 'nt', 'af')
        ## remove reference af's
        x.af <- x.af[!(mito.dna[x.af$pos] == x.af$nt), ]
        ## remove N site
        x.af <- x.af[!(mito.dna[x.af$pos] == 'N'), ]
        x.af$name <- paste0(x.af$pos, ' ', mito.dna[x.af$pos], '>', x.af$nt)
        ## find dominant NT
        x.af$af[x.af$pos %in% which(zeroes)] <- NA
        x <- x.af$af
        names(x) <- x.af$name
        return(x)
    }, mc.cores = ncores) #remove parallelism here
    varaf <- do.call(cbind, varaf)
    ## you could allow for only sites with coverage!
    ## currently you filter at a rate of 10% cells dropping out max
    ##varaf <- varaf[rowSums(is.na(varaf))/length(BaseCounts) < max.fraction.na,]
    varaf <-
        varaf[rowSums(varaf > min.af, na.rm = TRUE) >= min.num.samples, ]
    
    is.names <-
        sapply(exclusionlists.use, function(x)
            typeof(x) == "character")
                                        #part 2 - filter based on the exclusionlist
    if (sum(is.names) > 0) {
        removal.names.list <- unique(unlist(exclusionlists.use[is.names]))
        varaf <- varaf[!row.names(varaf) %in% removal.names.list, ]
    }
    if (sum(!is.names) > 0) {
        removal.ranges.list <-
            unique(unlist(GenomicRanges::GRangesList(exclusionlists.use[!is.names])))
        varaf <-
            varaf[-c(S4Vectors::queryHits(
                                    GenomicRanges::findOverlaps(mut2GR(row.names(varaf)), removal.ranges.list)
                                )), ]
    }
    ##if(drop.empty){
    ##colSums(varaf,na.rm=TRUE) > 0
    ##}
    ## remove empty rows
    varaf <- varaf[rowSums(varaf, na.rm = TRUE) > 0, , drop = FALSE]
    ## remove variants that are at a certain universal af in a certain number of cells
    varaf <-
        varaf[!rowSums(varaf >= min.af.universal, na.rm = TRUE) >= universal.var.cells, ]
    ## vars must have less than X % NA's
    varaf <- varaf[rowSums(is.na(varaf)) < max.var.na * NCOL(varaf), ]
    ## cells must have less than X % NA's
    varaf <-
        varaf[, colSums(is.na(varaf)) < max.cell.na * NROW(varaf)]
    MN <- pullcountsVars(BaseCounts, rownames(varaf), colnames(varaf))
    mutationCallsFromMatrix(t(MN$M), t(MN$N), ...)
}

#'Convert mutation string to GRanges
#'
#'@param mut The mutation to convert to a GRanges in the format of
#'"position reference>alternate".
#'@return Returns a GRanges object containg the site of the variant
#'along with reference/alternate allele data in the metacolumns
#'@examples mutation.as.granges <- mut2GR('1434 G>A')
#'@examples mutation.as.granges.no.space <- mut2GR('1434G>A')
#'@export
mut2GR <- function(mut) {
    gr <-
        GenomicRanges::GRanges(paste0('chrM:', as.numeric(gsub(
                                                   " *[A-Z].*", "", mut
                                               ))))
    gr$ref <-
        sapply(mut, function(x) {
            unlist(strsplit(gsub("\\d+ *", "", x), ">"))[[1]]
        })
    gr$alt <-
        sapply(mut, function(x) {
            unlist(strsplit(gsub("\\d+ *", "", x), ">"))[[2]]
        })
    return(gr)
}

#'Pull variant counts
#'
#'@param BaseCounts A list of base call matrices (one matrix per cell)
#'as produced by \code{\link{baseCountsFromBamList}}
#'@param vars Character vector of variants to pull, in format 5643G>T
#'@param cells Character vector for cells to select, or NULL if all
#'cells from the input are to be used
#'@return A list with two entries, M (count table on the variant
#'allele) and N (count table on the reference allele)
#'@examples load(system.file("extdata/example_counts.Rda",package = "mitoClone2"))
#'known.variants <- c("9 T>C","12 G>A","13 G>A")
#'counts.known.vars <- pullcountsVars(example.counts, vars=known.variants)
#'@export
pullcountsVars <- function(BaseCounts, vars, cells = NULL) {
    var.gr <- mut2GR(vars)
    ## subset cells of interest if appropriate
    if (!is.null(cells)) {
        message('Subsetting for specific cells...')
        if (!all(cells %in% names(BaseCounts))) {
            stop(
                paste0(
                    'You are attempting to subset for ',
                    sum(!(cells %in% names(BaseCounts))),
                    'cell(s) that do not appear in your input BaseCounts!'
                )
            )
        }
        BaseCounts <- BaseCounts[names(BaseCounts) %in% cells]
    }
    ## pull counts alt
    M <- sapply(BaseCounts, function(cell) {
        mapply(function(p, x)
            cell[p, x],
            GenomicRanges::start(var.gr),
            S4Vectors::mcols(var.gr)$alt)
    })
    ## pull total counts per position
    N <- sapply(BaseCounts, function(cell) {
        rowSums(cell[GenomicRanges::start(var.gr), c('A', 'G', 'C', 'T'), drop =
                                                                              FALSE])
    })
    if (!is.matrix(M) | !is.matrix(N)) {
        M <- matrix(M, ncol = length(M), dimnames = list(vars, names(M)))
        N <- matrix(N, ncol = length(N), dimnames = list(vars, names(N)))
    }
    ## maintain consistency with previous code for N by excluding the mutant allele calls
    N <- N - M
    rownames(M) <- row.names(N) <- vars
    return(list(M = M, N = N))
}
