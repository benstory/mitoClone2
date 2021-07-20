#'mutationCalls class
#'
#'To create this class from a list of bam files (where each bam file corresponds
#'to a single cell), use \code{\link{mutationCallsFromCohort}} or
#'\code{\link{mutationCallsFromExclusionlist}}. To create this class if you
#'already have the matrices of mutation counts, use its contstructor, i.e.
#'\code{mutationCallsFromMatrix(M = data1, N = data2)}.
#'
#'@slot M A matrix of read counts mapping to the \emph{mutant} allele. Columns
#'  are genomic sites and rows and single cells.
#'@slot N A matrix of read counts mapping to the \emph{nonmutant} alleles.
#'  Columns are genomic sites and rows and single cells.
#'@slot ternary Discretized version describing the mutational status of each
#'  gene in each cell, where 1 signfiies mutant, 0 signifies reference, and ?
#'  signifies dropout
#'@slot cluster Boolean vector of length \code{ncol(M)} specifying if the given
#'  mutation should be included for clustering (\code{TRUE}) or only used for
#'  annotation.
#'@slot metadata Metadata frame for annotation of single cells (used for
#'  plotting). Row names should be the same as in \code{M}
#'@slot tree Inferred mutation tree
#'@slot cell2clone Probability matrix of single cells and their assignment to
#'  clones.
#'@slot mut2clone Maps mutations to main clones
#'@slot mainClone Probability matrix of single cells and their assignment to
#'  main clones
#'@slot treeLikelihoods Likelihood matrix underlying the inference of main
#'  clones, see \code{\link{clusterMetaclones}}
#'@export
mutationCalls <- setClass(
    "mutationCalls",
    slots = c(
        M = "matrix",
        N = "matrix",
        metadata = "data.frame",
        ternary = "matrix",
        cluster = "logical",
        tree = "list",
        cell2clone = "matrix",
        mut2clone = "integer",
        mainClone = "matrix",
        treeLikelihoods = "matrix"
        
    ),
    validity = function(object) {
        if (!identical(dim(object@M), dim(object@N))) {
            return("Matrices M and N must have identical dimensions")
        }
        return(TRUE)
    }
)

#'mutationCalls constructor
#'
#'To be used when allele-specific count matrices are available.
#'@param M A matrix of read counts mapping to the \emph{mutant}
#'allele. Columns are genomic sites and rows and single cells.
#'@param N A matrix of read counts mapping to the \emph{referece}
#'allele. Columns are genomic sites and rows and single cells.
#'@param cluster If \code{NULL}, only mutations with coverage in 20
#'percent of the cells or more will be used for the clustering,
#'and all other mutations will be used for cluster annotation
#'only. Alternatively, a boolean vector of length \code{ncol(M)}
#'that specifies the desired behavior for each genomic site.
#'@param metadata A data.frame of metadata that will be transfered to
#'the final output where the \code{row.names(metadata)}
#'correspond to the the \code{row.names(M)}.
#'@param binarize Allele frequency threshold to define a site as
#'mutant (required for some clustering methods)
#'@return An object of class \code{\link{mutationCalls}}.
#'@examples load(system.file("extdata/example_counts.Rda",package = "mitoClone2"))
#' ## we have loaded the example.counts object
#' known.variants <- c("8 T>C","4 G>A","11 G>A","7 A>G","5 G>A","15 G>A","14 G>A")
#' known.subset <- pullcountsVars(example.counts, known.variants)
#' known.subset <- mutationCallsFromMatrix(t(known.subset$M), t(known.subset$N),
#' cluster = rep(TRUE, length(known.variants)))
#'@export
mutationCallsFromMatrix <- function(M,
                                    N,
                                    cluster = NULL,
                                    metadata = data.frame(row.names = rownames(M)),
                                    binarize = 0.05) {
    colnames(M) <- make.names(colnames(M))
    colnames(N) <- make.names(colnames(N))
    binfun <- function(M, N) {
        alleleRatio <-
            M / (M + N)
        apply(alleleRatio, 2, function(x)
            ifelse(is.na(x), "?", ifelse(x > binarize, "1", "0")))
    }
    ## if (!is.null(cluster)){
    ##     ##out@cluster <- cluster
    ## }else {
    ##   out@cluster <-
    ## apply(out@ternary!="?", 2, mean) > 0.2
    ## & apply(out@ternary=="1", 2, mean) > 0.04
    ## the last filter was not used when I made the figure
    ## there was a filter on the allele freq. in RNA.
    ## Should maybe include this in the other routines? 
    ## }
    ternary <- binfun(M, N)
    if (is.null(cluster)) {
        cluster <- apply(ternary != "?", 2, mean) > 0.2
    }
    out <-
        methods::new(
                     "mutationCalls",
                     M = M,
                     N = N,
                     metadata = metadata,
                     ternary = ternary,
                     cluster = cluster
                 )
    return(out)
}

#'Plot clonal assignment of single cells
#'
#'Creates a heatmap of single cell mutation calls, clustered using
#' PhISCS.
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param what One of the following: \emph{alleleFreq}: The fraction of
#'reads mapping to the mutant allele or \emph{ternary}:
#'Ternarized mutation status
#'@param show boolean vector specifying for each mutation if it should
#'be plotted on top of the heatmap as metadata; defaults to
#'mutations not used for the clustering \code{!mutcalls@cluster}
#'@param ... any arguments passed to \code{\link[pheatmap]{pheatmap}}
#'@examples P1 <-
#'readRDS(system.file("extdata/sample_example1.RDS",package =
#'"mitoClone2"))
#'plotClones(P1)
#'@return Returns TRUE only used for generating a PostScript tree
#'image of the putative mutation tree
#'@export
plotClones <- function(mutcalls,
                       what = c("alleleFreq", "ternary"),
                       show = c(),
             ...) {
    what <- match.arg(what)
    if (what == "alleleFreq")
        plotData <- mutcalls@M / (mutcalls@M + mutcalls@N)
    if (what == "ternary")
        plotData <-
            apply(mutcalls@ternary, 2, function(x)
                ifelse(x == "1", 1, ifelse(x == "?", 0,-1)))
    plotData <-
        t(plotData[, getNodes(mutcalls@tree)[-1]]) #how to order rows?
    if (length(show) > 1)
        annos <-
            data.frame(row.names = rownames(mutcalls@M),
                       mutcalls@ternary[, show],
                       mutcalls@metadata)
    if (length(show) == 1) {
        annos <-
            data.frame(row.names = rownames(mutcalls@M),
                       ann = mutcalls@ternary[, show],
                       mutcalls@metadata)
        colnames(annos)[2] <- show
    }
    if (length(show) == 0)
        annos <-
            data.frame(row.names = rownames(mutcalls@M), mutcalls@metadata)
    
    if (length(mutcalls@mut2clone) > 0) {
        annos$mainClone <-
            as.factor(apply(mutcalls@mainClone, 1, which.max))
        annos$confidence <- apply(mutcalls@mainClone, 1, max)
        plotData <- plotData[, order(annos$mainClone)]
    }
    
    pheatmap::pheatmap(
                  plotData,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE,
                  show_colnames = FALSE,
                  color = colorRampPalette(rev(
                      c("#9B0000", "#FFD72E", "#FFD72E", "#00009B")
                  ))(100),
                  annotation_col = annos,
                  ...
              )
    
}

#'mutationCalls accessors
#'
#'Retrieves the full matrix of likelihoods associating single cells
#' with clones
#'@param mutcall object of class \code{\link{mutationCalls}}.
#'@param mainClones Retrieve likelihoods associated with the main
#'Clones. Defaults to \code{TRUE} if
#'\code{\link{clusterMetaclones}} has been run.
#'@return Return \code{TRUE} if \code{\link{clusterMetaclones}} has
#'been run otherwise returns the cell by clone matrix of
#'likelihood associating each cell to a given clone.
#'@examples load(system.file("extdata/LudwigFig7.Rda",package =
#'"mitoClone2"))
#'likelihood_matrix <- getCloneLikelihood(LudwigFig7)
#'@export
getCloneLikelihood <- function(mutcall,
                               mainClones = length(mutcall@mut2clone) > 0)
    mutcall@cell2clone

#' @describeIn getCloneLikelihood Retrieve the most likely clone
#'associate with each cell.
getMainClone <-
    function(mutcall,
             mainClones = length(mutcall@mut2clone) > 0)
        as.factor(apply(
            getCloneLikelihood(mutcall, mainClones = mainClones),
            1,
            which.max
        ))

#' @describeIn getCloneLikelihood Retrieve the likelihood of the most
#'likely clone for each cell.
getConfidence <-
    function(mutcall,
             mainClones = length(mutcall@mut2clone) > 0)
        as.factor(apply(getCloneLikelihood(mutcall,mainClones = mainClones),
                        1,
                        max))

#' @describeIn getCloneLikelihood Retrieve the assignment of mutations
#'to clones, once \code{\link{clusterMetaclones}} has been run.
getMut2Clone <- function(mutcall)
    mutcall@mut2clone

#'mutationCalls cluster accessor
#'
#'Extracts all the putative variants that we want to use for
#' clustering
#'@param mutcall object of class \code{\link{mutationCalls}}.
#'@examples load(system.file("extdata/LudwigFig7.Rda",package =
#'"mitoClone2"))
#'mutations_to_cluster <- getVarsCandidate(LudwigFig7)
#'@return Returns a character vector including all the variants to be
#'used for clustering
#'@export
getVarsCandidate <- function(mutcall)
    mutcall@cluster

#'mutationCalls cluster setter
#'
#'Sets the putative variants that we want to use for clustering
#'@param mutcall object of class \code{\link{mutationCalls}}.
#'@param varlist vector of booleans with the names set to the variants
#'to use for clustering
#'@examples load(system.file("extdata/LudwigFig7.Rda",package =
#'"mitoClone2"))
#'mutations_to_cluster <- getVarsCandidate(LudwigFig7)
#'mutations_to_cluster[] <- rep(c(TRUE,FALSE),each=19)
#'LudwigFig7 <- setVarsCandidate(LudwigFig7,mutations_to_cluster)
#'@return Sets the cluster slot on a mutationCalls object
#'@export
setVarsCandidate <- function(mutcall, varlist) {
    methods::slot(mutcall, 'cluster') <- varlist
    return(mutcall)
}

#'mutationCalls counts accessor
#'
#'Extracts the counts of allele for either the mutant or all the
#' non-mutant alleles
#'@param mutcall object of class \code{\link{mutationCalls}}.
#'@param type character that is either `mutant` or `nonmutant`
#'depending on which allele count the user wants to access
#'@examples load(system.file("extdata/LudwigFig7.Rda",package = "mitoClone2"))
#'mutantAllele_count <- getAlleleCount(LudwigFig7,type='mutant')
#'@return Returns matrix of either mutant or non-mutant allele counts
#'@export
getAlleleCount <- function(mutcall, type = c('mutant', 'nonmutant')) {
    message('Extracting sum of all ', type, ' alleles')
    pulledslot <- switch(type, "mutant" = "M", "nonmutant" = "N")
    return(methods::slot(mutcall, pulledslot))
}
