#'Inference of mutational trees by of single cell mutational status
#'
#'From data on the observed mutational status of single cells at a
#' number of genomic sites, computes a likely phylogenetic tree using
#' PhISCS (https://github.com/sfu-compbio/PhISCS) and associates
#' single cells with leaves of the tree.  The function
#' \code{\link{clusterMetaclones}} should be called on the output in
#' order to group mutations into clones using a likelihood-based
#' approach.
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param fn false negative rate, i.e. the probability of only
#'observing the reference allele if there is a mutation. #add
#'gene-wise
#'@param fp false positive, i.e. the probability of observing the
#'mutant allele if there is no mutation.
#'@param cores number of cores to use for PhISCS (defaults to 1)
#'@param time maximum time to be used for PhISCS optimization, in
#'seconds (defaults to 10000)
#'@param tempfolder temporary folder to use for PhISCS output
#'@param python_env Any shell commands to execute in order to make the
#'gurobi python package available. The easiest solution is
#'running R from an environment where the gurobi python package
#'is avaiable. In some settings (e.g. RStudio Server), this
#'parameter can be used instead. \code{muta_clone} executes
#'PhISCS using a \code{system} call to python. The value of this
#'parameter is prepended to the call. If you have a conda
#'environment \code{myenv} that contains gurobipy, \code{source
#'activate myenv} can work. Occassionally RStudio Server modifies
#'your PATH so that that the conda and source commands are not
#'available. In that case you can for example use \code{export
#'PATH=/path/to/conda/:$PATH; source activate myenv}. easybuild
#'users can \code{module load anaconda/v3; source activate myenv}
#'@param force_recalc Rerun PhISCS even if the \code{tempfolder}
#'contains valid PhISCS output
#'@param method A string variable of either PhISCS or SCITE depending
#'on the tree-inferring software the user wants to use. Default:
#'PhISCS
#'@examples load(system.file("extdata/LudwigFig7.Rda",package =
#'"mitoClone2"))
#'LudwigFig7 <- varCluster(LudwigFig7,
#'python_env = "",method='SCITE')
#'@return an object of class \code{\link{mutationCalls}}, with an
#'inferred tree structure and cell to clone assignment added.
#'@export

varCluster <- function(mutcalls,
                       fn = 0.1,
                       fp = 0.02,
                       cores = 1,
                       time = 10000,
                       tempfolder = tempdir(),
                       python_env = '',
                       force_recalc = FALSE,
                       method = 'SCITE') {
    ##prepare data and run PhISCS
    suppressWarnings(dir.create(tempfolder))
    
    usedata <- mutcalls@ternary[, mutcalls@cluster]
    if (!method %in% c('PhISCS', 'SCITE')) {
        message('No method selected, defaulting to SCITE')
        method <- 'SCITE'
    }
    if (!file.exists(file.path(tempfolder, "out"))) {
        message(paste0(
            'Creating temporary dir: ',
            file.path(tempfolder, 'out'),
            ' for tree-building output.'
        ))
        dir.create(file.path(tempfolder, 'out'))
    }
    if (method == 'PhISCS') {
        write.table(
            usedata,
            file = file.path(tempfolder, "in.txt"),
            quote = FALSE,
            sep = "\t",
            row.names = gsub("[><_]", "", rownames(usedata)),
            col.names = gsub("[><_]", "", colnames(usedata))
        )
        ## if(!system("which PhISCS", intern = FALSE,
        ##ignore.stderr = TRUE, ignore.stdout = TRUE) == "0"){
        ##     stop("PhISCS not detected on your system!")
        ## }else{
        ##     ##base <-
        ##system.file("extdata/python/PhISCS-I",package = "mitoClone")
        ## }
        base <- 'PhISCS-I'
        command <-
            sprintf(
                "%spython %s -SCFile %s -fn %.2f -fp %.2f -o %s -threads %d -time %d --drawTree",
                ifelse(python_env == "", "", paste0(python_env, "; ")),
                base,
                file.path(tempfolder, "in.txt"),
                fn,
                fp,
                file.path(tempfolder, "out"),
                cores,
                time
            )
        if (!file.exists(file.path(tempfolder, "out", "in.CFMatrix")) |
            force_recalc) {
            ##message("Now running the following command:", command)
            message("Please run the following command with PhISCS-I:",
                    command)
            ##tryCatch(system(command), error = function(e){
            ## stop("PhISCS error: ",e,
            ## "Make sure that the gurobi python package is available
            ## and consider specifying python_env."))
            ##}
        } else {
            message("Results found for PhISCS run")
        }
        ##read in the result and create tree data structure
        physics <-
            read.table(
                file.path(tempfolder, "out", "in.CFMatrix"),
                header = TRUE,
                row.names = 1,
                sep = "\t"
            )
        cell.mutations <- physics
        txtCon <-
            textConnection(unique(apply(physics, 1, paste, collapse = ",")))
        clones <-
            read.csv(txtCon,
                     header = FALSE ,
                     col.names = colnames(physics))
        mutcalls@tree <- clones2tree(clones)
        nodes.order <- getNodes(mutcalls@tree)[-1]
        clone.names <-
            apply(clones[, nodes.order], 1, function(x)
                max(which(x == 1)))
        clone.names <- nodes.order[clone.names]
        clone.names[is.na(clone.names)] <- "root"
        clones <- as.matrix(clones)
        clone.names -> rownames(clones)
    }
    ## if (method == 'SACS'){
    ##     base <- 'sacs'
    ##       command <-
    ##         sprintf(
    ##             "%s %s -i %s -e %s -E %s -m %d -n %d -r 2 -a %.2f -b %.2f -k 0 > %s",
    ##             ifelse(python_env == "", "", paste0(python_env, "; ")),
    ##             base,
    ##             file.path(tempfolder, "input_scite.txt"),
    ##             file.path(tempfolder, "mut_sacs_names.txt"),
    ##             file.path(tempfolder, "cell_sacs_names.txt"),
    ##             NCOL(usedata),
    ##             NROW(usedata),
    ##             fn,
    ##             fp,
    ##             file.path(tempfolder, "sacs_out")
    ##         )
    ##     mut.raw <- usedata
    ##     mut.raw[mut.raw == '?'] <- '2'
    ##     cell.names <- row.names(mut.raw)
    ##     mutation.names <- colnames(mut.raw)
    ##     mut.raw <- apply(mut.raw, 2, as.numeric)
    ##     write.table(
    ##         x = mut.raw,
    ##         file.path(tempfolder, "input_sacs.txt"),
    ##         sep = ' ',
    ##         row.names = FALSE,
    ##         col.names = FALSE,
    ##         quote = FALSE
    ##     )##,sep=' ')
    ##     write.table(
    ##         mutation.names,
    ##         file.path(tempfolder, "mut_sacs_names.txt"),
    ##         row.names = FALSE,
    ##         col.names = FALSE,
    ##         quote = FALSE
    ##     )
    ##     write.table(
    ##         cell.names,
    ##         file.path(tempfolder, "cell_sacs_names.txt"),
    ##         row.names = FALSE,
    ##         col.names = FALSE,
    ##         quote = FALSE
    ##     )
    ##     if (!file.exists(file.path(tempfolder, "sacs_out")))
    ##         {
    ##             message("Please run the following command with SACS:",
    ##                     command)
    ##         } else {
    ##             message("Results found for SACS run")
    ##     }
    ## }
    if (method == 'SCITE') {
        if (file.exists(system.file("SCITE/scite", package = "mitoClone2"))) {
            base <- system.file("SCITE/scite", package = "mitoClone2")
        } else{
            base <- system.file("SCITE/scite.exe", package = "mitoClone2")
        }
        command <-
            sprintf(
                "%s %s -i %s -names %s -n %d -max_treelist_size 1 -a -m %d -r 1 -l 200000 -ad %.2f -fd %.2f -o %s",
                ifelse(python_env == "", "", paste0(python_env, "; ")),
                base,
                file.path(tempfolder, "input_scite.txt"),
                file.path(tempfolder, "mut_scite_names.txt"),
                NCOL(usedata),
                NROW(usedata),
                fn,
                fp,
                file.path(tempfolder, "scite_out")
            )
        mut.raw <- usedata
        mut.raw[mut.raw == '?'] <- '3'
        cell.names <- row.names(mut.raw)
        mutation.names <- colnames(mut.raw)
        mut.raw <- t(apply(mut.raw, 2, as.numeric))
        write.table(
            x = mut.raw,
            file.path(tempfolder, "input_scite.txt"),
            sep = ' ',
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )##,sep=' ')
        write.table(
            mutation.names,
            file.path(tempfolder, "mut_scite_names.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
        if (!file.exists(file.path(tempfolder, "scite_out_ml0.gv")) |
            force_recalc) {
            message("If you use this method for publicaiton, please make sure to cite :")
            message(paste0("Jahn, K., Kuipers, J. & Beerenwinkel, N. Tree inference for ",
                           "single-cell data. Genome Biol 17, 86 (2016). ",
                           "https://doi.org/10.1186/s13059-016-0936-x"))
            message("Now running the following command:", command)
            tryCatch(
                system(command),
                error = function(e)
                    stop(
                        "SCITE error: ",
                        e,
                        "Make sure that the SCITE package is installed."
                    )
            )
        } else{
            message(
                paste0("Results found, skipping SCITE run. Warning! As SCITE will not overwrite previous results",
                       " if you are attempting to run this function again you may potentially need to delete",
                       " your previous run or provide a new tempfolder.")
            )
        }
        scite.in <-
            readLines(file.path(tempfolder, "scite_out_ml0.gv"))
        begin.si <- grep("node \\[color", scite.in)
        scite.tree <- scite.in[seq(from = 3, to = begin.si[2] - 1)]
        scite.tree <- data.frame(t(sapply(scite.tree, function(s) {
            gsub("\\;$", "", unlist(strsplit(s, split = " -> ")))
        })))
        colnames(scite.tree) <- c('start', 'end')
        
        scite.cells <-
            scite.in[seq(from = begin.si[2] + 1, to = length(scite.in) - 1)]
        scite.cells <- data.frame(t(sapply(scite.cells, function(s) {
            gsub("\\;$", "", unlist(strsplit(s, split = " -> ")))
        })))
        ## PLEASE MAKE SURE TO DOUBLE CHECK THAT THIS IS ACCURATE!
        scite.cells$real <-
            cell.names[as.numeric(gsub("^s", "", scite.cells$X2)) + 1]
        ##  states = number of mutations possible
        mutcalls@tree <- list(
            states = rep(0, nrow(scite.tree)),
            children = list(),
            mutation = "root"
        )
        names(mutcalls@tree$states) <- scite.tree$end
        ## consider merging clones into the scite.tree
        ##clones <- data.frame(scite.tree)
        mutcalls@tree <-
            addnodesSCITE(mutcalls@tree, scite.tree = scite.tree)
        ##mutcalls@tree <- clones2tree(clones)
        nodes.order <- getNodes(mutcalls@tree)[-1]
        flat.tree <-
            lapply(rapply(mutcalls@tree, enquote, how = "unlist"), eval)
        flat.states.df <-
            unique(do.call(rbind, flat.tree[grep('states', names(flat.tree))]))
        row.names(flat.states.df) <- NULL
        clones <- as.matrix(flat.states.df)
        clone.names <-
            apply(flat.states.df[, nodes.order], 1, function(x)
                max(which(x == 1)))
        ##
        clone.names <- nodes.order[clone.names]
        clone.names[is.na(clone.names)] <- "root"
        row.names(clones) <- clone.names
        ##clones <- data.frame(scite.tree)
        scite.cells.multi <- split(scite.cells, scite.cells$real)
        scite.cells.multi.pre <-
            lapply(scite.cells.multi, function(x) {
                x$X1 <- gsub("^Root$", "root", x$X1)
                if (NROW(x) == 1) {
                    x.ret <- clones[x[1, ]$X1, , drop = FALSE]
                } else{
                    clone.reference <- clones[row.names(clones) %in% x$X1, ]
                    oldest.parent <-
                        which.min(apply(clone.reference, 1, function(z) {
                            ## associating multi-associated cells with the oldest parent cell possible
                            ##which.max(apply(clones,1,function(z){
                            ##raw.calls <- suppressWarnings(as.numeric(usedata[unique(x$real),]))
                            ##raw.calls[is.na(raw.calls)] <- 0
                            ## fixed dropped mutations in clustering
                            raw.calls <-
                                as.numeric(gsub("^\\?$", "3", usedata[unique(x$real), colnames(clone.reference)]))
                            z <- z[!is.na(raw.calls)]
                            ##raw.calls <- raw.calls[!is.na(raw.calls)]
                            ## use -1* for which.max
                            ##print(sum(raw.calls != z))
                            sum(raw.calls != z)
                        }))
                    x.ret <- clone.reference[oldest.parent, , drop = FALSE]
                }
                row.names(x.ret) <- unique(x$real)
                return(x.ret)
            })
        scite <- data.frame(do.call(rbind, scite.cells.multi.pre))
        cell.mutations <- scite
    }
    ##retrieve the assignment for every cell
    ##compute likelihood of assignments
    
    cell2clone.prob <-
        apply(usedata[, nodes.order], 1, function(cell) {
            apply(clones[, nodes.order], 1, function(clone) {
                                        #likelihood
                sum(log10(ifelse(
                    cell == "0" & clone == 1,
                    fn,
                          ifelse(
                              cell == "1" & clone == 0,
                              fp,
                          ifelse(
                              cell == "1" & clone == 1,
                              1 - fn,
                          ifelse(
                              cell == "0" & clone == 0,
                              1 - fp,
                          ifelse(cell == "?", 1, (1 -
                                                  fn) * (1 - fp))
                          )
                          )
                          )
                )))
            })
        })
    mutcalls@cell2clone <- t(apply(cell2clone.prob, 2, function(x) {
        10 ^ x / sum(10 ^ x)
    }))
    ##finally, determine which mutations to group:
    evaluate_likelihood <- function(data, idealized) {
        sum(log10(ifelse(
            data == "0" & idealized == 1,
            fn,
                  ifelse(
                      data == "1" & idealized == 0,
                      fp,
                  ifelse(
                      data == "1" & idealized == 1,
                      1 - fn,
                  ifelse(
                      data == "0" & idealized == 0,
                      1 - fp,
                  ifelse(data == "?", 1, (1 -
                                          fn) * (1 - fp))
                  )
                  )
                  )
        )))
        
    }
    ref <-
        evaluate_likelihood(usedata[, colnames(cell.mutations)], cell.mutations)
    mutcalls@treeLikelihoods <-
        sapply(colnames(cell.mutations), function(node1) {
            sapply(c(colnames(cell.mutations), "root"), function(node2) {
                newmodel <- cell.mutations
                if (node2 == "root")
                    newmodel[, node1] <-
                        rep(0, nrow(newmodel))
                else
                    newmodel[, node1] <- newmodel[, node2]
                evaluate_likelihood(usedata[, colnames(cell.mutations)], newmodel) - ref
            })
        })
    return(mutcalls)
}

#'Remove mutations that occuring at the same site
#'
#'Mutations co-occuring at the same genomic position may often be the
#' result of sequencing artifacts or technical biases. In cases where
#' the user which to drop these from a result this function may be
#' used. ONLY WORKS FOR MITOCHONDRIAL MUTATIONS.
#'@param x A list of strings that comprise sites that will be filtered
#'@param window Integer of how close mutations must be to one another
#'(in bp) to be removed
#'@return Returns the same list of mutations excluding those, if any,
#'that fall within the same window =
#'@examples P1.muts <- rep(TRUE,3)
#'names(P1.muts) <- c("X2537GA","X3351TC","X3350TC")
#'names(P1.muts) <- gsub("^X","",
#'gsub("(\\d+)([AGCT])([AGCT])","\\1 \\2>\\3",names(P1.muts)))
#'P1.muts <- P1.muts[removeWindow(names(P1.muts))]
#'@export
removeWindow <- function(x, window = 1) {
    if (any(sapply(x, grepl, '^[A-Z]'))) {
        stop(paste0(
            'Non-standard variant name detected.',
            ' Please only provide variant coordinates',
            ' in the following format: 1337 G>A')
        )
    }
    firstround <- sort(mut2GR(x))
    removeidentical <-
        unique(firstround[which(GenomicRanges::countOverlaps(firstround) > 1)])
    firstround <- GenomicRanges::resize(firstround, width = window)
    bigreg <- GenomicRanges::reduce(firstround)
    bigreg.w <- as.numeric(GenomicRanges::width(bigreg))
    bigreg <- bigreg[bigreg.w > window]
    bigreg <- GenomicRanges::reduce(c(bigreg, removeidentical))
    if (length(c(S4Vectors::queryHits(
                                GenomicRanges::findOverlaps(mut2GR(x), bigreg)
                            ))) > 0) {
        x <-
            x[-c(S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2GR(x), bigreg)))]
    }
    return(x)
}


#'Quick clustering of mutations
#'
#'Performs a quick hierarchical clustering on a object of class
#' \code{\link{mutationCalls}}. See \code{\link{varCluster}} for an
#' alternative that infers mutational trees and uses sound models of
#' dropout.
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param binarize If \code{FALSE}, will use raw allele frequencies for
#'the clustering. If \code{TRUE}, will use binarized
#'mutation/reference/dropout calls.
#'@param drop_empty Remove all rows in the provided mutcalls object
#'where no cells exhibit a mutation.
#'@param ... Parameters passed to \code{\link[pheatmap]{pheatmap}}
#'@return The result of running \code{\link[pheatmap]{pheatmap}}
#'@examples load(system.file("extdata/LudwigFig7.Rda",package = "mitoClone2"))
#' quickCluster <- quick_cluster(LudwigFig7)
#'@export
quick_cluster <- function(mutcalls,
                          binarize = FALSE,
                          drop_empty = TRUE,
                          ...) {
    if (drop_empty)
        mutcalls@ternary <-
            mutcalls@ternary[apply(mutcalls@ternary, 1, function(x)
                any(x == "1")), ]
    if (binarize)
        converted <-
            t(apply(mutcalls@ternary, 2, function(x)
                ifelse(x == "1", 1, ifelse(x == "0", -1, 0))))
    if (!binarize)
        converted <- t(mutcalls@M / (mutcalls@M + mutcalls@N))
    clustered <-
        pheatmap::pheatmap(converted[mutcalls@cluster, ], ...)
}

#'Cluster mutations into clones - following the tree structure
#'
#'PhISCS orders all mutations into a hierarchical mutational tree; in
#' many cases, the exact order of the acquisition of individual
#' mutations in not unanimously determined from the data. This
#' function computes the change in likelihood of the infered clonal
#' assignment if two mutations are merged into a clone. Hierarchical
#' clustering is then used to determine the clonal structure. The
#' result is visualized and should be fine-tuned using the
#' \code{min.lik} parameter.
#'@param mutcalls mutcalls object of class \code{\link{mutationCalls}}
#'for which \code{\link{varCluster}} has been run
#'@param min.lik specifies the minimum difference in likelihood
#'required. This parameter is set arbitrarily, see the vignette
#'"Computation of clonal hierarchies and clustering of mutations"
#'for more information.
#'@param plot whether dendrograms should be plotted.
#'@return Returns the provided \code{\link{mutationCalls}} class
#'object with an additional 'mainClone' metadata which allows for
#'further refinement of clonal population and association of
#'cells with a cluster of mutations (in this case clones).
#'@examples P1 <- readRDS(system.file("extdata/sample_example1.RDS",package = "mitoClone2"))
#' P1 <- clusterMetaclones(P1)
#' ## access via mainClone metadata
#'@export
clusterMetaclones <- function(mutcalls,
                              min.lik = 1,
                              plot = TRUE) {
                                        #split the tree into branches with no further splits
    branches <- getBranches(mutcalls@tree)
    mutcalls@mut2clone <-
        as.integer(rep(0, nrow(mutcalls@treeLikelihoods)))
    names(mutcalls@mut2clone) <- rownames(mutcalls@treeLikelihoods)
    
    par(mfrow = c(ceiling(sqrt(
            length(branches)
        )), ceiling(sqrt(
                length(branches)
            ))))
    
    for (i in 1:length(branches)) {
        if (sum(branches[[i]] != "root") <= 1) {
            mutcalls@mut2clone[branches[[i]]] <-
                as.integer(max(mutcalls@mut2clone) + 1)
        } else {
            ##d <- as.dist(1-cor(t(mutcalls@treeLikelihoods[branches[[i]],])))
            ub <- branches[[i]][branches[[i]] != "root"]
            d <- dist(t(mutcalls@treeLikelihoods[ub, ub]))
            ##pheatmap::pheatmap(mutcalls@treeLikelihoods[branches[[i]],branches[[i]]],
            ##clustering_distance_cols =
            ## as.dist(1-cor(mutcalls@treeLikelihoods[branches[[i]],branches[[i]]])),
            ##clustering_distance_rows =
            ##as.dist(1-cor(t(mutcalls@treeLikelihoods[branches[[i]],branches[[i]]]))))
            cl <- hclust(d)
            tryCatch(
                plot(cl),
                error = function(e)
                    "Failed to plot"
            )
            mutcalls@mut2clone[branches[[i]]] <-
                as.integer(max(mutcalls@mut2clone) + cutree(cl, h = nrow(mutcalls@M) * min.lik))
        }
        
    }
    use <- mutcalls@mut2clone[colnames(mutcalls@cell2clone)]
    mutcalls@mainClone <- sapply(unique(use), function(mainClone) {
        if (sum(use == mainClone) == 1)
            mutcalls@cell2clone[, use == mainClone]
        else
            apply(mutcalls@cell2clone[, use == mainClone], 1, sum)
    })
    
    return(mutcalls)
}


#'Plot clone-specific variants in circular plots
#'
#'@param variants Character vector of variants to plot in format
#'5643G>T or 5643 G>T.
#'@param patient Characet vector identifying which variant belongs to
#'what clone. The order should match that of the 'vars' parameter
#'and shoul dbe of identical length. If none is provided, the
#'function assumes all variants are from one single sample which
#'will be named "Main Clone". Default: NULL.
#'@param genome The mitochondrial genome of the sample being
#'investigated. Please note that this is the UCSC standard
#'chromosome sequence. Default: hg38.
#'@param customGenome A GRanges object containing a custom annotation.
#'If provided, this genome will be used instead of the predefined options
#'specified by the `genome` parameter. Default is NULL.
#'@param showLegend Boolean for whether or not the gene legend should
#'be present in the final output plot. Default: TRUE.
#'@param showLabel Boolean for whether or not the name of the variant
#'should be shown as a label in the final output plot. Default:
#'TRUE.
#'@return A ggplot object illustrating the clone specific mutations.
#'@examples known.variants <- c("9001 T>C","12345 G>A","1337 G>A")
#'mitoPlot(known.variants)
#'@export
mitoPlot <- function(variants,
                     patient = NULL,
                     genome = 'hg38',
                     customGenome = NULL,
                     showLegend = TRUE,
                     showLabel = TRUE) {
    mito.var <- mut2GR(variants)
    if (!is.null(customGenome)) {
      # Use custom genome if provided
      mito.gr <- customGenome
    } else {
      mito.gr <- switch(genome,
                        "hg38" = hg38.mito,
                        "hg19" = hg19.mito,
                        "mm10" = mm10.mito)
    }
    mito.gr <- mito.gr[mito.gr$gene_biotype != 'Mt_tRNA']
    S4Vectors::mcols(mito.gr) <-
        S4Vectors::mcols(mito.gr)[, c('external_gene_name'), drop = FALSE]
    mito.genes.df <-
        data.frame(c(
            GenomicRanges::resize(mito.gr, width = 1, fix = 'start'),
            c(GenomicRanges::resize(
                                 mito.gr, width = 1, fix = 'end'
                             ))
        ))
    mito.genes.df$gene <-
        gsub("mt-|MT-", "", mito.genes.df$external_gene_name)
    mito.genes.df$gene <-
        factor(mito.genes.df$gene, levels = rev(sort(unique(
                                       mito.genes.df$gene
                                   ))))
    mito.genes.df$external_gene_name <- NULL
    mito.genes.df$type <- 'mito'
    mito.gene.color <-
        c(grDevices::colorRampPalette(
                         c(
                             "springgreen4",
                             "chartreuse4",
                             "chartreuse",
                             "darkolivegreen2",
                             "darkolivegreen"
                         )
                     )(length(unique(
                                         subset(mito.genes.df, mito.genes.df$type == 'mito')$gene
                                     ))), rev(c(
                                              "#F9B90AFF", "#E6352FFF", "#3D79F3FF"
                                          )))
    mito.genes.df$sample <- 'Main Clone'
    ## setup the color scheme
    var.df <- data.frame(mito.var)
    var.df$ref <-  var.df$alt <- NULL
    var.df$gene <- variants
    var.df$type <- 'mutation'
    if (!is.null(patient)) {
        var.df$sample <- patient
        patcol <- grDevices::rainbow(length(unique(patient)))
        names(patcol) <- unique(patient)
        plot.df.list <- lapply(unique(patient), function(z) {
            plot.df.sub <-
                rbind(mito.genes.df,
                      subset(var.df[, colnames(mito.genes.df)], sample == z, drop = FALSE))
            plot.df.sub$sample <- z
            plot.df.sub$strand <- patcol[z]
            return(plot.df.sub)
        })
        plot.df <- do.call(rbind, plot.df.list)
    } else{
        var.df$sample <- 'Main Clone'
        var.df$strand <- 'red'
        plot.df <-
            rbind(mito.genes.df, var.df[, colnames(mito.genes.df)])
    }
    ## prepare to plot
    legendPos <- ifelse(showLegend, "top", "none")
    mito.plot.df <- subset(plot.df,plot.df$type == 'mito')
    p <- ggplot2::ggplot(data=mito.plot.df, ggplot2::aes(x = start, y=12, color=gene)) +
      ggplot2::geom_hline(yintercept=12, color = "black",alpha=1) +
      ggplot2::geom_line(size=4) +
      ggplot2::theme_void(base_size=24) + ggplot2::xlab('') +
      ggplot2::ylab('') +
      ggplot2::theme(legend.position = legendPos, axis.text.x = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values=mito.gene.color) +
      ggplot2::geom_point(data=subset(plot.df,plot.df$type == 'mutation'),size=5, ggplot2::aes(x = start, y = 12), color=subset(plot.df,plot.df$type == 'mutation')$strand) +
      ggplot2::coord_polar() +
      ggplot2::facet_wrap(~sample) +
      ggplot2::ylim(0,13)
  if(showLabel){
      p <- p + ggplot2::geom_text(data=subset(plot.df,plot.df$type == 'mutation'),ggplot2::aes(x = start, y = 12, label = subset(plot.df,plot.df$type == 'mutation')$gene),color='black',nudge_y = -3)
  }
  p
}


#'Manually overwrite clustering of mutations into clones
#'
#'The function \code{\link{clusterMetaclones}} provides an automated
#' way to group mutations into clones for subsequent analyses (such as
#' differential expression analyses). In practice, it may make sense
#' to overwrite these results manually. See the vignette 'Computation
#' of clonal hierarchies and clustering of mutations' for an example.
#'@param mutcalls mutcalls object of class \code{\link{mutationCalls}}
#'for which \code{\link{clusterMetaclones}} has been run
#'@param mutation2clones Named integer vector that assigns mutations
#'to clones. See the vignette 'Computation of clonal hierarchies
#'and clustering of mutations' for an example.
#'@return Returns the provided \code{\link{mutationCalls}} class
#'object with the 'mainClone' metadata overwritten with the
#'manual values provided by the user.
#'@examples P1 <- readRDS(system.file("extdata/sample_example1.RDS",package = "mitoClone2"))
#' new.n <- seq(17)
#' names(new.n) <- names(getMut2Clone(P1))
#' P1.newid <- overwriteMetaclones(P1,new.n)
#'@export
overwriteMetaclones <- function(mutcalls, mutation2clones) {
                                        #split the tree into branches with no further splits
    mutcalls@mut2clone <- mutation2clones
    use <- mutcalls@mut2clone[colnames(mutcalls@cell2clone)]
    mutcalls@mainClone <- sapply(unique(use), function(mainClone) {
        if (sum(use == mainClone) == 1)
            mutcalls@cell2clone[, use == mainClone]
        else
            apply(mutcalls@cell2clone[, use == mainClone], 1, sum)
    })
    
    return(mutcalls)
}


getBranches <- function(tree) {
    current <- c()
    while (length(tree$children) == 1) {
        current <- c(current, tree$mutation)
        tree <- tree$children[[1]]
    }
    if (length(tree$children) == 0)
        return(c(current, tree$mutation))
    result <- removeDepth(lapply(tree$children, getBranches))
    result[[length(result) + 1]] <- c(current, tree$mutation)
    return(result)
}

removeDepth <- function(list) {
    out <- list()
    for (i in 1:length(list)) {
        if (is.list(list[[i]])) {
            for (j in 1:length(list[[i]])) {
                out[[length(out) + 1]] <- list[[i]][[j]]
            }
        } else
            out[[length(out) + 1]] <- list[[i]]
    }
    out
}
