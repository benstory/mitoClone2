library(mitoClone2)

context("10x file reader")

test_that("correct input read from BAM", {
    input <- bam2R_10x(file = system.file("extdata", "mm10_10x.bam", package="mitoClone2"), sites="chrM:1-15000")
    expect_length(input,6)
    expect_named(expected=c('ZAAATCGAAAATCCG',
                            'ZACTGGACGTCTTCG',
                            'ZGGAAGACGCTCTTT',
                            'ZGTGACCCCCTCTTG',
                            'ZTCACCCCAAGCTCT',
                            'ZTGGACGCATTATAC'),
                 input)
    for (i in seq_along(input)) {
        expect_equal(NROW(input[[i]]), expected=15000)
        expect_equal(NCOL(input[[i]]), expected=8)
    }
})
