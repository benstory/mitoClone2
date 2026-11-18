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
#'\code{\link{bam2R_10x}}. Alternatively the sparse counts
#'produced by \code{\link{bam2R_v2}}.
#'@param patient A character vector associating each cell / entry in
#'the \code{BaseCount} list with a patient
#'@param sites Vector specifying genomic regions, defaults to the
#'entire mitochondrial genome. Excepts a string but may be
#'included as a GRanges object. Not needed for the sparse counts
#'of \code{\link{bam2R_v2}}.
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
#'@param customGenome A GRanges object containing a custom annotation. If
#'provided, this genome will be used instead of the predefined options
#'specified by the `genome` parameter. Default is NULL.
#'@param all.chromosomes Boolean, only for the sparse counts of
#'\code{\link{bam2R_v2}}. FALSE (default): variants are called on
#'the mitochondrial chromosome only. TRUE: variants are called on
#'every site of the counts, on any chromosome. Variants outside of
#'the mitochondrial chromosome are named with their chromosome,
#'e.g. chr7:140753336A>T. The thresholds are the same for all
#'sites, please note that the defaults were chosen for the high
#'coverage of the mitochondrial genome.
#'@param reference Only used if \code{all.chromosomes = TRUE}. The
#'genome that the BAM file was aligned to, either the path of a
#'FASTA file (an index is created if there is none) or a BSgenome
#'object. The sequence names have to be identical to the BAM file.
#'Default NULL: the mitochondrial sites use the \code{genome}
#'annotation and for all other sites the most abundant nucleotide
#'is used as reference allele (with a warning: homozygous
#'variants can not be detected that way and REF/ALT are swapped
#'where the variant is the majority).
#'@return A list of \code{\link{mutationCalls}} objects (one for each
#'\code{patient}) and an entry named \code{exclusionlist}
#'containing a exclusionlist of sites with variants in several
#'individuals
#'@examples sites.gr <- GenomicRanges::GRanges("chrM:1-15000")
#'BaseCounts <- bam2R_10x(file = system.file("extdata",
#'"mm10_10x.bam", package="mitoClone2"), sites=sites.gr)
#'mutCalls <- mutationCallsFromCohort(BaseCounts,
#'patient=c('sample2','sample1','sample2','sample2','sample1'),
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
                                    genome = 'hg38',
                                    customGenome = NULL,
                                    all.chromosomes = FALSE,
                                    reference = NULL) {
    if (!isSparseCounts10x(BaseCounts) && (all.chromosomes || !is.null(reference))) {
        stop('all.chromosomes and reference are only available for the sparse counts',
             ' of bam2R_v2, the list with one matrix per cell does not know the chromosome.')
    }
    if (!all.chromosomes && !is.null(reference)) {
        stop('reference is only used together with all.chromosomes = TRUE.')
    }
    ## sparse counts from bam2R_v2: the sites are part of the counts
    if (isSparseCounts10x(BaseCounts)) {
        mito.dna <- NULL
        if (USE.REFERENCE) {
            if (!is.null(customGenome)) {
                mito.dna <- customGenome
            } else {
                mito.dna <- switch(genome,
                                   "hg38" = hg38.dna,
                                   "hg19" = hg19.dna,
                                   "mm10" = mm10.dna)
            }
        }
        return(
            cohort10x(
                BaseCounts,
                patient = patient,
                MINREADS = MINREADS,
                MINCELL = MINCELL,
                MINFRAC = MINFRAC,
                MINCELLS.PATIENT = MINCELLS.PATIENT,
                MINFRAC.PATIENT = MINFRAC.PATIENT,
                MINFRAC.OTHER = MINFRAC.OTHER,
                mito.dna = mito.dna,
                all.chromosomes = all.chromosomes,
                reference = reference
            )
        )
    }
    message("Making sure 'sites' parameter is set correctly...")
    if (!length(sites) == 1) {
        stop('Your sites parameter must be a character vector or GRanges object of length 1')
    }
    sites <- GenomicRanges::GRanges(sites)
    if (length(patient) != length(BaseCounts)) {
        stop('patient must have one entry per cell (element of BaseCounts): ',
             length(BaseCounts), ' cells, ', length(patient), ' patients.')
    }
    if (GenomicRanges::width(sites) > NROW(BaseCounts[[1]])) {
        stop('sites covers ', GenomicRanges::width(sites), ' positions, but the counts only have ',
             NROW(BaseCounts[[1]]), ' rows. Please use the sites that the counts were read with.')
    }
    ## read in the
    ntcountsArray <- simplify2array(BaseCounts)
    ntcountsArray <- aperm(ntcountsArray, c(1, 3, 2))
    ntcountsArray <- ntcountsArray[, , c('A', 'T', 'C', 'G', 'N')]
    if (USE.REFERENCE) {
      if (!is.null(customGenome)) {
        # Use custom genome if provided
        reference <- customGenome
      } else {
        reference <- switch(genome,
                            "hg38" = hg38.dna,
                            "hg19" = hg19.dna,
                            "mm10" = mm10.dna)
      }
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
    if (is.null(variant_calls)) {
        noCandidates10x(MINCELL, MINREADS)
    }
    ## count number of cell that are mutant per patient
    varcount.bypatient <- sapply(unique(patient), function(pa) {
        ## drop = FALSE: a patient with a single cell crashed before
        as.integer(colSums(variant_calls[patient == pa, , drop = FALSE] == "BOTH" |
                           variant_calls[patient == pa, , drop = FALSE] == "MUT"))
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
    noPatientVariants10x(sum(singlepatient), length(singlepatient), sum(multipatient))
    ## extract mutations that are unique to individual patients
    mutation.bypatient <-
        colnames(varcount.bypatient)[apply(varcount.bypatient[singlepatient, , drop = FALSE], 1, which.max)]
    ## pull the variant cells
    variant_calls_selected <- variant_calls[, singlepatient, drop = FALSE]
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
    ## drop = FALSE: a single shared variant was lost before
    out$exclusionlist <- rownames(varcount.bypatient[multipatient, , drop = FALSE])
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
#'as produced by \code{\link{baseCountsFromBamList}} or
#'\code{\link{bam2R_10x}}. Alternatively the sparse counts
#'produced by \code{\link{bam2R_v2}}.
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
#'@param customDNA A character vector containing a custom DNA sequence. If
#'provided, this sequence will be used instead of the predefined options
#'specified by the `genome` parameter. Default is NULL.
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
                                           customDNA = NULL,
                                           ncores = 1,
                                           ...) {
    ## sparse counts from bam2R_v2: the defaults refer to the number of
    ## cells, which is the number of columns and not the length of the list
    sparse <- isSparseCounts10x(BaseCounts)
    if (sparse) {
        if (missing(min.num.samples)) {
            min.num.samples <- 0.01 * ncol(BaseCounts[[1]])
        }
        if (missing(universal.var.cells)) {
            universal.var.cells <- 0.95 * ncol(BaseCounts[[1]])
        }
    }
    if (!is.null(customDNA)) {
    # Use custom genome sequence if provided
    mito.dna <- customDNA
    } else {
    mito.dna <- switch(genome,
                       "hg38" = hg38.dna,
                       "hg19" = hg19.dna,
                       "mm10" = mm10.dna)
    
    }
    if (!is.null(customDNA) | genome != 'hg38'){
      print('Be aware! The provided exclusionlists are defined based on HUMAN GENOME hg38. Provide custom lists if this is not your genome.')
    }
    if (sparse) {
        return(
            exclusionlist10x(
                BaseCounts,
                lim.cov = lim.cov,
                min.af = min.af,
                min.num.samples = min.num.samples,
                min.af.universal = min.af.universal,
                universal.var.cells = universal.var.cells,
                exclusionlists.use = exclusionlists.use,
                max.var.na = max.var.na,
                max.cell.na = max.cell.na,
                mito.dna = mito.dna,
                ncores = ncores,
                ...
            )
        )
    }
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
        varaf[rowSums(varaf > min.af, na.rm = TRUE) >= min.num.samples, , drop = FALSE]
    
    ## vapply: an empty list of exclusionlists is allowed
    is.names <-
        vapply(exclusionlists.use, function(x)
            typeof(x) == "character", logical(1))
                                        #part 2 - filter based on the exclusionlist
    if (sum(is.names) > 0) {
        removal.names.list <- unique(unlist(exclusionlists.use[is.names]))
        varaf <- varaf[!row.names(varaf) %in% removal.names.list, , drop = FALSE]
    }
    if (sum(!is.names) > 0 && nrow(varaf) > 0) {
        removal.ranges.list <-
            unique(unlist(GenomicRanges::GRangesList(exclusionlists.use[!is.names])))
        ## (a negative index would drop ALL variants if none of them overlaps)
        on.list <-
            S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2GR(row.names(varaf)), removal.ranges.list))
        varaf <- varaf[!seq_len(nrow(varaf)) %in% on.list, , drop = FALSE]
    }
    ##if(drop.empty){
    ##colSums(varaf,na.rm=TRUE) > 0
    ##}
    ## remove empty rows
    varaf <- varaf[rowSums(varaf, na.rm = TRUE) > 0, , drop = FALSE]
    ## remove variants that are at a certain universal af in a certain number of cells
    varaf <-
        varaf[!rowSums(varaf >= min.af.universal, na.rm = TRUE) >= universal.var.cells, , drop = FALSE]
    ## vars must have less than X % NA's
    varaf <- varaf[rowSums(is.na(varaf)) < max.var.na * NCOL(varaf), , drop = FALSE]
    ## cells must have less than X % NA's
    varaf <-
        varaf[, colSums(is.na(varaf)) < max.cell.na * NROW(varaf), drop = FALSE]
    if (nrow(varaf) == 0 || ncol(varaf) == 0) {
        stop('No variant or no cell passed the filters. Try a lower lim.cov, min.af or',
             ' min.num.samples, or a higher max.var.na or max.cell.na.', call. = FALSE)
    }
    MN <- pullcountsVars(BaseCounts, rownames(varaf), colnames(varaf))
    mutationCallsFromMatrix(t(MN$M), t(MN$N), ...)
}

#'Convert mutation string to GRanges
#'
#'@param mut The mutation to convert to a GRanges in the format of
#'"position reference>alternate". The chromosome may be given in
#'front, e.g. "chr7:140753336 A>T", variants without a chromosome
#'are located on chrM.
#'@return Returns a GRanges object containg the site of the variant
#'along with reference/alternate allele data in the metacolumns
#'@examples mutation.as.granges <- mut2GR('1434 G>A')
#'@examples mutation.as.granges.no.space <- mut2GR('1434G>A')
#'@examples mutation.on.chr7 <- mut2GR('chr7:140753336 A>T')
#'@export
mut2GR <- function(mut) {
    ## an optional chromosome in front of the position (the chromosome
    ## name may contain a colon itself, e.g. HLA-A*01:01:25 C>T)
    pattern <- "^(.*):([0-9]+ *[A-Za-z]+>[A-Za-z]+)$"
    has.chr <- grepl(pattern, mut)
    chr <- ifelse(has.chr, sub(pattern, "\\1", mut), 'chrM')
    site <- ifelse(has.chr, sub(pattern, "\\2", mut), mut)
    pos <- as.numeric(gsub(" *[A-Z].*", "", site))
    if (all(chr == 'chrM')) {
        gr <- GenomicRanges::GRanges(paste0('chrM:', pos))
    } else {
        gr <- GenomicRanges::GRanges(seqnames = chr, ranges = paste0(pos, '-', pos))
    }
    names(site) <- mut
    gr$ref <-
        sapply(site, function(x) {
            unlist(strsplit(gsub("\\d+ *", "", x), ">"))[[1]]
        })
    gr$alt <-
        sapply(site, function(x) {
            unlist(strsplit(gsub("\\d+ *", "", x), ">"))[[2]]
        })
    return(gr)
}

#'Pull variant counts
#'
#'@param BaseCounts A list of base call matrices (one matrix per cell)
#'as produced by \code{\link{baseCountsFromBamList}}. Alternatively
#'the sparse counts produced by \code{\link{bam2R_v2}}, see
#'\code{\link{pullcountsVars10x}}.
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
    ## sparse counts from bam2R_v2
    if (isSparseCounts10x(BaseCounts)) {
        return(pullcountsVars10x(BaseCounts, vars, cells))
    }
    bad <- !grepl("^(.+:)?[0-9]+ *[A-Z]>[A-Z]$", vars)
    if (any(bad)) {
        stop('Variants must be given as position reference>alternate, e.g. 1434 G>A or',
             ' 1434G>A. Not understood: ', paste(utils::head(vars[bad], 3), collapse = ', '))
    }
    var.gr <- mut2GR(vars)
    if (length(BaseCounts) > 0 && max(GenomicRanges::start(var.gr)) > NROW(BaseCounts[[1]])) {
        stop('The variant position ', max(GenomicRanges::start(var.gr)), ' is behind the end of',
             ' the counts (', NROW(BaseCounts[[1]]), ' positions).')
    }
    ## the list with one matrix per cell only contains the mitochondrial genome
    if (any(as.character(GenomicRanges::seqnames(var.gr)) != 'chrM')) {
        stop('The base counts (one matrix per cell) only contain the mitochondrial genome.',
             ' Use the sparse counts of bam2R_v2 or baseCountsFromBamList_v2 for variants',
             ' on other chromosomes.')
    }
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
    ## always variants x cells, also for a single variant (the cell
    ## names were changed in that case before)
    M <- matrix(M, nrow = length(vars), dimnames = list(vars, names(BaseCounts)))
    N <- matrix(N, nrow = length(vars), dimnames = list(vars, names(BaseCounts)))
    ## maintain consistency with previous code for N by excluding the mutant allele calls
    N <- N - M
    rownames(M) <- row.names(N) <- vars
    return(list(M = M, N = N))
}
