% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{getCloneLikelihood}
\alias{getCloneLikelihood}
\alias{getMainClone}
\alias{getConfidence}
\alias{getMut2Clone}
\title{mutationCalls accessors}
\usage{
getCloneLikelihood(mutcall, mainClones = length(mutcall@mut2clone) > 0)

getMainClone(mutcall, mainClones = length(mutcall@mut2clone) > 0)

getConfidence(mutcall, mainClones = length(mutcall@mut2clone) > 0)

getMut2Clone(mutcall)
}
\arguments{
\item{mutcall}{object of class \code{\link{mutationCalls}}.}

\item{mainClones}{Retrieve likelihoods associated with the main
Clones. Defaults to \code{TRUE} if
\code{\link{clusterMetaclones}} has been run.}
}
\value{
Return \code{TRUE} if \code{\link{clusterMetaclones}} has
been run otherwise returns the cell by clone matrix of
likelihood associating each cell to a given clone.
}
\description{
Retrieves the full matrix of likelihoods associating single cells
with clones
}
\section{Functions}{
\itemize{
\item \code{getMainClone()}: Retrieve the most likely clone
associate with each cell.

\item \code{getConfidence()}: Retrieve the likelihood of the most
likely clone for each cell.

\item \code{getMut2Clone()}: Retrieve the assignment of mutations
to clones, once \code{\link{clusterMetaclones}} has been run.

}}
\examples{
load(system.file("extdata/LudwigFig7.Rda",package =
"mitoClone2"))
likelihood_matrix <- getCloneLikelihood(LudwigFig7)
}
