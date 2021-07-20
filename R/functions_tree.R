attach2tree <- function(tree, clone, mutation = NULL) {
    if (length(tree$children) == 0) {
        ##change!
        if (is.null(mutation))
            mutation <- names(clone)[clone == 1 & tree$states == 0]
        if (length(mutation) > 1) {
            tree$children[[1]] <-
                list(states = clone,
                     children = list(),
                     mutation = mutation[1])            
            tree <- attach2tree(tree, clone, mutation[-1])
        } else {
            tree$children[[1]] <-
                list(states = clone,
                     children = list(),
                     mutation = mutation[1])
        }
        
    } else {
        attached <- FALSE
        for (i in 1:length(tree$children)) {
            if (all(clone[tree$children[[i]]$states == 1] == 1)) {
                tree$children[[i]] <-
                    attach2tree(tree$children[[i]], clone, mutation)
                attached <- TRUE
                break
            }
        }
        
        if (!attached) {
            if (is.null(mutation))
                mutation <- names(clone)[clone == 1 & tree$states == 0]
            
            if (length(mutation) > 1) {
                if (is.null(mutation))
                    mutation <- names(clone)[clone == 1 & tree$states == 0]
                
                tree$children[[length(tree$children) + 1]] <-
                    list(states = clone,
                         children = list(),
                         mutation = mutation[1])
                
                tree <- attach2tree(tree, clone, mutation[-1])
            } else {
                tree$children[[length(tree$children) + 1]] <-
                    list(states = clone,
                         children = list(),
                         mutation = mutation[1])
            }
        }
        
    }
    
    return(tree)
}


addnodesSCITE <- function(tree,
                          current = NULL,
                          scite.tree = scite.tree) {
    if (is.null(current)) {
        current <- 'Root'
    }
    if (current %in% names(tree$states)) {
        tree$states[current] <- 1
    }
    first.child <- which(scite.tree$start %in% current)
    if (length(first.child) == 0) {
        return(list(
            states = tree$states,
            children = list(),
            mutation = current
        ))
    }
    if (length(first.child) == 1) {
        tree$children[[1]] <-
            list(
                states = tree$states,
                children = list(),
                mutation = scite.tree$end[first.child]
            )
    }
    if (length(first.child) > 1) {
        for (i in 1:length(first.child)) {
            tree$children[[i]] <-
                list(
                    states = tree$states,
                    children = list(),
                    mutation = scite.tree$end[first.child[i]]
                )
        }
    }
    if (length(first.child) > 0) {
        for (i in 1:length(tree$children)) {
            tree$children[[i]] <-
                addnodesSCITE(tree$children[[i]],
                              current = scite.tree$end[first.child[i]],
                              scite.tree = scite.tree)
        }
    }
    return(tree)
}


clones2tree <- function(clones) {
    tree <- list(
        states = rep(0, ncol(clones)),
        children = list(),
        mutation = "root"
    )
    names(tree$states) <- colnames(clones)
    while (nrow(clones) > 0) {
        nmuts <- apply(clones, 1, sum)
        mini <- min(nmuts)
        where <- which.min(nmuts)
        
        if (mini == 0) {
            clones <- clones[-where, ]
            next
        }
        
        #find where it attaches
        tree <- attach2tree(tree, unlist(clones[where, ]))
        clones <- clones[-where, ]
    }
    tree
}

topy <- function(tree, init = TRUE) {
    if (is.null(tree$child1) & is.null(tree$child2)) {
        out <- sprintf("Node('%s')", tree$mutation)
    } else if (!is.null(tree$child2) & !is.null(tree$child1)) {
        out <-
            sprintf(
                "Node('%s').%s.%s",
                tree$mutation,
                topy(tree$child1, FALSE),
                topy(tree$child2, FALSE)
            )
    } else if (is.null(tree$child2)) {
        out <-
            sprintf("Node('%s').%s", tree$mutation, topy(tree$child1, FALSE))
    } else {
        stop("Something wrong")
    }
    
    if (!init)
        paste0("addkid(", out , ")")
    else
        out
}

getNodes <- function(tree) {
    out <- tree$mutation
    if (length(tree$children) > 0)
        for (i in 1:length(tree$children))
            out <- c(out, getNodes(tree$children[[i]]))
    out
}

toGraphviz <- function(tree, header = TRUE) {
    if (header)
        out <-
            c("digraph G {",
              "node [color=deeppink4, style=filled, fontcolor=white];")
    else
        out <- c()
    if (length(tree$children) > 0) {
        for (i in 1:length(tree$children)) {
            out <-
                c(out,
                  sprintf("%s -> %s;", tree$mutation, tree$children[[i]]$mutation))
            out <- c(out, toGraphviz(tree$children[[i]], header = FALSE))
        }
    }
    if (header)
        out <- c(out, "}")
    out
}

cutTreeAt <- function(tree, node) {
    if (length(tree$children) > 0) {
        for (i in 1:length(tree$children)) {
            if (tree$children[[i]]$mutation  == node) {
                #if the target node is a child of the current root:
                #Simply return the current root and all its other children as one tree, and a second tree
                #rooted at the target node
                upper <- tree
                upper$children <- upper$children[-i]
                lower <- tree$children[[i]]
                return(list(upper = upper, lower = lower))
            } else {
                #if this is not the case: Search the children
                out <- cutTreeAt(tree$children[[i]], node)
                #in this case, the upper tree is then current root +
                if (!is.null(out)) {
                    newUpper <- tree
                    newUpper$children[[i]] <- out$upper
                    return(list(upper = newUpper, lower = out$lower))
                }
            }
        }
        #nothing found, search in children
    } else
        return(NULL)
}

attachTreeAt <- function(tree, add, node) {
    if (tree$mutation == node) {
        tree$children[[length(tree$children) + 1]] <- add
    } else if (length(tree$children) > 0) {
        for (i in 1:length(tree$children)) {
            tree$children[[i]] <- attachTreeAt(tree$children[[i]], add, node)
        }
    }
    return(tree)
}

permutate <- function(tree) {
    #The basic MCMC move we use is prune and reattach.
    #We sample a node i uniformly from the n available - here try all combinations?
    nodes <- getNodes(tree)
    possibilities <- list()
    for (j in 2:length(nodes)) {
        selected <- nodes[j]
        #and cut the edge leading to this node to remove the subtree from the tree.
        splitTree <- cutTreeAt(tree, selected)
        #Then we sample one of the remaining nodes (including the root) uniformly
        remainingNodes <- getNodes(splitTree$upper)
        
        for (k in 1:length(remainingNodes)) {
            selected <- remainingNodes[k]
            possibilities[[length(possibilities) + 1]] <-
                attachTreeAt(splitTree$upper, splitTree$lower, selected)
        }
    }
    possibilities
}

tree2clones <- function(tree, first = TRUE) {
    out <- tree$states
    if (length(tree$children) > 0) {
        for (i in 1:length(tree$children)) {
            out <- rbind(out, tree2clones(tree$children[[i]], first = FALSE))
        }
    }
    if (first) {
        i <- 2
        while (i < nrow(out)) {
            for (j in 1:(i - 1)) {
                if (identical(out[i, ], out[j, ]))
                    out <- out[-i, ]
            }
            i <- i + 1
        }
        rownames(out) <- as.character(1:nrow(out))
        out <- as.data.frame(out)
    }
    out
}
