\name{createBifurcate}
\alias{createBifurcate}
\alias{buildApe}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a fully bifurcating object of class 'phylo'.}
\description{createBifurcate and buildApe take a 'paleoPhylo' object and convert it to its fully bifurcating equivalent (i.e. all ancestors terminate if they speciate).  createBifurcate converts the lineage phylogeny into one that only features bifurcating speciations.  buildApe does the conversion to ape.}
\usage{
createBifurcate(pP)
buildApe(pP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}).}
}
\details{
\code{createBifurcate} works by *walking* through the lineage phylogeny until a speciation event occurs.  If the speciation event is a budding-type (where the ancestor lineage persists simultaneously with the descendant), then \code{createBifurcate} splits the ancestor such that the speciation is a bifurcating-type (new, but , unique IDs are generated).  This process is repeated iteratively until extant (or most recently deceased) tips have been dealt with.  The goal is to enable easy conversion to class 'phylo' for use with \code{\link[ape]{read.tree}} formats.

\code{buildApe} converts this into an object of class 'phylo' for use in \code{ape} (e.g. \code{\link[ape]{read.tree}} or \code{\link[ape]{plot.phylo}}), i.e. assigning of node labels and branch lengths, edge matrix and root definitions. The root is usually the first row of the data frame and needs to be unique.

\code{createBifurcate} and \code{buildApe} are kept separate to enable less manipulation before plotting fully bifurcating trees with arguments in \code{\link{drawLinPhylo}}.  
}
\value{
An object of class 'paleoPhylo'; see \code{\link{as.paleoPhylo}}}
\references{\url{www.imperial.ac.uk/research/apurvis/paleoPhylo/}}
\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}
\seealso{\code{\link{as.paleoPhylo}}}
\examples{
data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))
apeTree <- buildApe(createBifurcate(p93))
apeTree

##using ape, do some simple illustratory stuff
##plot a fully bifurcating phylogeny using plot.phylo
plot(reorder(apeTree))
##write apeTree using the Newick (aka New Hampshire) format.
write.tree(apeTree)
}
\keyword{utilities}{createBifurcate}

