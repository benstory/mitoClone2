library(mitoClone2)

context("10x file reader, version 2")

bam <- system.file("extdata", "mm10_10x.bam", package = "mitoClone2")
cells <- c('AAATCGAAAATCCG',
           'ACTGGACGTCTTCG',
           'GGAAGACGCTCTTT',
           'GTGACCCCCTCTTG',
           'TGGACGCATTATAC')

test_that("bam2R_v2 returns 16 sparse matrices", {
    sparse <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    expect_length(sparse, 16)
    expect_named(sparse, c('A','T','C','G','-','N','INS','DEL',
                           'a','t','c','g','_','n','ins','del'))
    for (m in sparse) {
        expect_s4_class(m, "dgCMatrix")
        expect_equal(colnames(m), cells)
    }
    expect_true(all(grepl("^chrM:[0-9]+$", rownames(sparse[[1]]))))
})

test_that("bam2R_10x_v2 gives the same counts as bam2R_10x", {
    v1 <- bam2R_10x(file = bam, sites = "chrM:1-15000")
    v2 <- bam2R_10x_v2(file = bam, sites = "chrM:1-15000")
    expect_named(v2, cells)
    for (cell in cells) {
        expect_equal(dim(v2[[cell]]), c(15000, 8))
        expect_equal(v2[[cell]], v1[[cell]])
    }
})

test_that("several regions in one pass give the same counts", {
    one <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    ## a gap at 13301, so the two regions stay separate (adjacent ranges are merged)
    two <- bam2R_v2(file = bam, sites = c("chrM:1-13300", "chrM:13302-15000"),
                    barcodes = cells)
    expect_false("chrM:13301" %in% rownames(two[[1]]))
    common <- intersect(rownames(one[[1]]), rownames(two[[1]]))
    expect_equal(length(common), nrow(two[[1]]))
    for (n in names(one)) {
        expect_equal(as.matrix(two[[n]][common, ]), as.matrix(one[[n]][common, ]))
    }
})

test_that("several cores give the same result", {
    ## mclapply can only use one core on Windows
    skip_on_os("windows")
    one <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    two <- bam2R_v2(file = bam, sites = "chrM:1-15000", barcodes = cells, ncores = 2)
    for (n in names(one)) {
        expect_equal(as.matrix(two[[n]]), as.matrix(one[[n]]))
    }
})

test_that("invalid tags are rejected", {
    expect_error(bam2R_v2(file = bam, sites = "chrM:1-100", barcode.tag = "CBX"))
    expect_error(bam2R_v2(file = bam, sites = "chrM:1-100", umi.tag = ""))
})

test_that("variant counts and strand bias", {
    sparse <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    vars <- c("13575 T>C", "13778 A>T")
    MN <- pullcountsVars10x(sparse, vars)
    expect_equal(dim(MN$M), c(2, 5))
    expect_equal(rownames(MN$M), vars)
    expect_gt(sum(MN$M), 0)
    sb <- strandBias10x(sparse, vars)
    expect_equal(rownames(sb), vars)
    expect_equal(sb$alt.fwd + sb$alt.rev, unname(rowSums(MN$M)))
    expect_true(all(sb$sor >= 1))
    ## lower case alleles, duplicated variants, unknown cells
    counts <- c('ref.fwd', 'ref.rev', 'alt.fwd', 'alt.rev')
    expect_equal(unlist(strandBias10x(sparse, "13575 t>c")[, counts]),
                 unlist(sb[1, counts]))
    expect_warning(dup <- strandBias10x(sparse, vars[c(1, 1)]))
    expect_equal(nrow(dup), 1)
    expect_error(strandBias10x(sparse, vars, cells = "NOT_A_CELL"), "do not appear")
    expect_error(pullcountsVars10x(sparse, vars, cells = "NOT_A_CELL"), "do not appear")
})

context("mutation calling with the sparse counts of version 2")

## compare two mutationCalls objects
expect_same_calls <- function(a, b) {
    expect_equal(a@M, b@M)
    expect_equal(a@N, b@N)
    expect_equal(a@ternary, b@ternary)
    expect_equal(a@cluster, b@cluster)
}

test_that("mutationCallsFromExclusionlist: sparse counts give the same result as dense counts", {
    dense <- bam2R_10x_v2(file = bam, sites = "chrM:1-15000")
    sparse <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    for (lists in list(list(), exclusionlists)) {
        call.variants <- function(counts) {
            mutationCallsFromExclusionlist(counts, lim.cov = 5, min.af = 0.1,
                                           min.num.samples = 1, universal.var.cells = 5,
                                           max.var.na = 1.1, max.cell.na = 1.1,
                                           genome = 'mm10', exclusionlists.use = lists)
        }
        from.sparse <- call.variants(sparse)
        expect_s4_class(from.sparse, "mutationCalls")
        expect_gt(ncol(from.sparse@M), 0)
        expect_same_calls(from.sparse, call.variants(dense))
    }
    ## the defaults refer to the number of cells, not to the 16 matrices
    expect_same_calls(mutationCallsFromExclusionlist(sparse, genome = 'mm10'),
                      mutationCallsFromExclusionlist(dense, genome = 'mm10'))
    expect_error(mutationCallsFromExclusionlist(sparse, genome = 'mm10', lim.cov = 1e6),
                 "No variant passed the filters")
})

test_that("mutationCallsFromCohort: sparse counts give the same result as dense counts", {
    dense <- bam2R_10x_v2(file = bam, sites = "chrM:1-15000")
    sparse <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    patient <- c('sample2', 'sample1', 'sample2', 'sample2', 'sample1')
    from.dense <- mutationCallsFromCohort(dense, patient = patient, MINCELL = 1, MINFRAC = 0,
                                          MINCELLS.PATIENT = 1, genome = 'mm10',
                                          sites = "chrM:1-15000")
    ## no sites needed
    from.sparse <- mutationCallsFromCohort(sparse, patient = patient, MINCELL = 1, MINFRAC = 0,
                                           MINCELLS.PATIENT = 1, genome = 'mm10')
    expect_named(from.sparse, c('sample2', 'sample1', 'exclusionlist'))
    expect_same_calls(from.sparse$sample2, from.dense$sample2)
    expect_null(from.sparse$sample1)
    expect_equal(from.sparse$exclusionlist, from.dense$exclusionlist)
    expect_error(mutationCallsFromCohort(sparse, patient = patient[1:3]), "one entry per cell")
})

test_that("pullcountsVars accepts sparse counts", {
    dense <- bam2R_10x_v2(file = bam, sites = "chrM:1-15000")
    sparse <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    vars <- c("13575 T>C", "13778 A>T")
    expect_equal(pullcountsVars(sparse, vars), pullcountsVars(dense, vars))
    expect_equal(pullcountsVars(sparse, vars, cells = cells[2:3]),
                 pullcountsVars(dense, vars, cells = cells[2:3]))
})

test_that("mutationCallsFromCohort: all.chromosomes calls variants outside of chrM", {
    sparse <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    patient <- c('sample2', 'sample1', 'sample2', 'sample2', 'sample1')
    call.variants <- function(counts, ...) {
        mutationCallsFromCohort(counts, patient = patient, MINCELL = 1, MINFRAC = 0,
                                MINCELLS.PATIENT = 1, genome = 'mm10', ...)
    }
    mito.only <- call.variants(sparse)
    ## the same counts a second time on a made up chromosome, shifted by 1000 bp
    two.chr <- lapply(sparse, function(m) {
        copy <- m
        rownames(copy) <-
            paste0("chr7:", as.integer(sub("^.*:", "", rownames(m))) + 1000L)
        rbind(m, copy)
    })
    ## default: only the mitochondrial chromosome, nothing changes
    expect_same_calls(call.variants(two.chr)$sample2, mito.only$sample2)
    ## a FASTA file with the made up chromosome only (soft masked, N in front)
    mm10.dna <- get("mm10.dna", envir = asNamespace("mitoClone2"))
    fasta <- tempfile(fileext = ".fa")
    writeLines(c(">chr7", paste(c(rep("N", 1000), tolower(mm10.dna)), collapse = "")), fasta)
    all.chr <- call.variants(two.chr, all.chromosomes = TRUE, reference = fasta)$sample2
    vars <- colnames(all.chr@M)
    on.chr7 <- grepl("^chr7", vars)
    expect_equal(sum(on.chr7), ncol(mito.only$sample2@M))
    ## the mitochondrial calls are unchanged, the chr7 calls are a copy of them
    expect_equal(all.chr@M[, !on.chr7], mito.only$sample2@M)
    expect_equal(unname(all.chr@M[, on.chr7]), unname(mito.only$sample2@M))
    expect_equal(unname(all.chr@N[, on.chr7]), unname(mito.only$sample2@N))
    ## without a reference: a warning, the most abundant allele is the reference
    expect_warning(no.ref <- call.variants(two.chr, all.chromosomes = TRUE),
                   "NO REFERENCE GENOME PROVIDED")
    expect_s4_class(no.ref$sample2, "mutationCalls")
    ## wrong sequence names, wrong input
    wrong <- tempfile(fileext = ".fa")
    writeLines(c(">7", paste(mm10.dna, collapse = "")), wrong)
    expect_error(call.variants(two.chr, all.chromosomes = TRUE, reference = wrong),
                 "does not contain")
    expect_error(call.variants(two.chr, reference = fasta), "all.chromosomes = TRUE")
    dense <- bam2R_10x_v2(file = bam, sites = "chrM:1-15000")
    expect_error(call.variants(dense, sites = "chrM:1-15000", all.chromosomes = TRUE),
                 "only available for the sparse counts")
})

context("one BAM file per cell, version 2")

test_that("baseCountsFromBamList_v2 gives the same counts as baseCountsFromBamList", {
    sites <- "chrM:1-15000"
    old <- baseCountsFromBamList(list(cell1 = bam), sites = sites)
    new <- baseCountsFromBamList_v2(c(cell1 = bam), sites = sites)
    expect_length(new, 16)
    expect_equal(colnames(new$A), "cell1")
    ## dense, strands merged, rows = positions 1..15000
    pos <- as.integer(sub("^.*:", "", rownames(new$A)))
    merged <- matrix(0, 15000, 8, dimnames = list(NULL, colnames(old$cell1)))
    for (nt in colnames(merged)) {
        rev.nt <- c(A = 'a', T = 't', C = 'c', G = 'g', '-' = '_', N = 'n',
                    INS = 'ins', DEL = 'del')[[nt]]
        merged[pos, nt] <- as.vector(new[[nt]][, 1] + new[[rev.nt]][, 1])
    }
    expect_equal(merged, old$cell1, check.attributes = FALSE)
})

test_that("baseCountsFromBamList_v2: one column per file, any number of regions", {
    ## all reads of the example file have a cell barcode: a file read as one
    ## cell is the sum of its barcodes
    per.barcode <- bam2R_v2(file = bam, sites = c("chrM:1-5000", "chrM:6000-15000"),
                            min_reads_per_barcode = 1)
    files <- baseCountsFromBamList_v2(c(a = bam, b = bam),
                                      sites = c("chrM:1-5000", "chrM:6000-15000"))
    expect_equal(colnames(files$A), c("a", "b"))
    for (n in names(files)) {
        expect_equal(rownames(files[[n]]), rownames(per.barcode[[n]]))
        expect_equal(as.vector(files[[n]][, "a"]), as.vector(Matrix::rowSums(per.barcode[[n]])))
        expect_equal(files[[n]][, "a"], files[[n]][, "b"])
    }
    ## unnamed files are named after the file
    expect_equal(colnames(baseCountsFromBamList_v2(bam, sites = "chrM:13000-13100")$A), "mm10_10x")
    ## a region without reads
    expect_warning(baseCountsFromBamList_v2(bam, sites = "chrM:1-100"), "no reads in the sites")
    expect_error(baseCountsFromBamList_v2(c(a = bam, a = bam)), "unique")
    ## the result goes straight into the mutation calling (two identical
    ## files: every variant is shared between the patients)
    expect_warning(calls <- mutationCallsFromCohort(files, patient = c("P1", "P2"), MINCELL = 1,
                                                    MINFRAC = 0, MINCELLS.PATIENT = 1,
                                                    genome = 'mm10'),
                   "specific for a single patient")
    expect_named(calls, c("P1", "P2", "exclusionlist"))
})

test_that("mut2GR accepts variants with a chromosome", {
    gr <- mut2GR(c("1434 G>A", "chr7:140753336 A>T", "HLA-A*01:01:25C>T"))
    expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chrM", "chr7", "HLA-A*01:01"))
    expect_equal(GenomicRanges::start(gr), c(1434, 140753336, 25))
    expect_equal(unname(gr$ref), c("G", "A", "C"))
    expect_equal(unname(gr$alt), c("A", "T", "T"))
    ## variants of different chromosomes are not merged by removeWindow
    expect_equal(removeWindow(c("100 A>G", "chr1:100 C>T", "chr1:101 G>A"), window = 3),
                 "100 A>G")
    ## the dense counts and mitoPlot are mitochondrial only
    load(system.file("extdata/example_counts.Rda", package = "mitoClone2"))
    expect_error(pullcountsVars(example.counts, "chr7:9 T>C"), "only contain the mitochondrial")
    expect_error(mitoPlot("chr7:9 T>C"), "only shows mitochondrial")
})

context("helpful messages")

test_that("common mistakes give clear errors and warnings", {
    ## wrong chromosome name: the message names the one of the file
    expect_error(bam2R_v2(file = bam, sites = "MT:1-100"), "called chrM in this file")
    expect_error(bam2R_10x(file = bam, sites = "MT:1-100"), "called chrM in this file")
    ## nothing counted
    expect_warning(bam2R_v2(file = bam, sites = "chrM:1-100", barcode.tag = "XB"),
                   "barcode tag XB")
    expect_warning(bam2R_v2(file = bam, sites = "chrM:1-15000", min_reads_per_barcode = 1e6),
                   "min_reads_per_barcode")
    expect_warning(bam2R_v2(file = bam, sites = "chrM:1-15000", barcodes = "AAAA-1"),
                   "None of the barcodes occurs in the file")
    expect_warning(bam2R_v2(file = bam, sites = "chrM:1-15000", umi.tag = "UB"),
                   "UMI tag UB")
    ## calling
    sparse <- bam2R_v2(file = bam, sites = "chrM:1-15000")
    dense <- bam2R_10x_v2(file = bam, sites = "chrM:1-15000")
    patient <- c('sample2', 'sample1', 'sample2', 'sample2', 'sample1')
    expect_error(mutationCallsFromCohort(sparse, patient = patient, genome = 'mm10'),
                 "No variant is supported by at least MINCELL")
    expect_error(mutationCallsFromCohort(dense, sites = "chrM:1-15000", patient = patient,
                                         genome = 'mm10'),
                 "No variant is supported by at least MINCELL")
    expect_error(mutationCallsFromCohort(dense, sites = "chrM:1-15000", patient = patient[1:3]),
                 "one entry per cell")
    expect_error(mutationCallsFromCohort(dense, sites = "chrM:1-16000", patient = patient),
                 "only have 15000 rows")
    ## a patient with a single cell
    load(system.file("extdata/example_counts.Rda", package = "mitoClone2"))
    one.cell <- suppressWarnings(
        mutationCallsFromCohort(example.counts, sites = "chrM:1-24",
                                patient = c("A", rep("B", length(example.counts) - 1)),
                                MINCELL = 2, MINCELLS.PATIENT = 2))
    expect_named(one.cell, c("A", "B", "exclusionlist"))
    expect_warning(mutationCallsFromCohort(sparse, patient = patient, MINCELL = 1, MINFRAC = 0,
                                           MINCELLS.PATIENT = 5, genome = 'mm10'),
                   "specific for a single patient")
    ## variants and clones
    expect_error(pullcountsVars(sparse, "1434GA"), "Not understood: 1434GA")
    expect_error(pullcountsVars(dense, "1434GA"), "Not understood: 1434GA")
    expect_error(pullcountsVars(dense, "20000 A>G"), "behind the end of the counts")
    expect_error(mutationCallsFromMatrix(matrix(1, 3, 2), matrix(1, 3, 3)), "same dimensions")
    one <- mutationCallsFromMatrix(matrix(1, 3, 1, dimnames = list(NULL, "X1A.G")),
                                   matrix(1, 3, 1, dimnames = list(NULL, "X1A.G")))
    expect_error(varCluster(one, tempfolder = tempfile()), "at least 2 variants")
})
