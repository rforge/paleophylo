\name{ape2paleoPhylo}
\alias{ape2paleoPhylo}

\title{Convert objects of class 'phylo' to class 'paleoPhylo'.}
\description{
Takes the information contained within objects of class \code{phylo} and converts using \code{\link{as.paleoPhylo}}.
}
\usage{
ape2paleoPhylo(phy,retainNodeLabels,nC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy} {An object of class \code{phylo} WITH KNOWN branch lengths.}
  \item{retainNodeLabels} {OPTIONAL: should internal node labels be retained from ape?  If FALSE, then \code{\link{nchar}} is used to replace node labels with blank names.  Defaults to TRUE.}
  \item{nC} {OPTIONAL: The cut off for \code{\link{nchar}} to overwrite \code{nodelabels}.  Any name with fewer characters than nC is overwritten. Defaults to 5.}
 }
  
\details{
The principal use is to enable plotting with \code{\link{drawLinPhylo}} and its array of axis formats. 

If tip.labels are included, then they are exported as \code{labels} in \code{paleoPhylo}.  Non-tip labels are assumed to be the node label.  Future revisions will provide the option to change this

Objects of class \code{phylo} are defined using edge matrices (see \code{\link[ape]{read.tree}}), whereas \code{\link{paleoPhylo}} objects require a unique root that gave rise to the focal clade.  The unique root is assumed to have geologically infinitessimal duration, such that it is not visible on plots.  Future releases will incorporate analytical tools; at this stage additional functionality will be introduced to adapt this.
}

\value{
An object of class 'paleoPhylo', which is a list containing the vectors defined in \code{\link{as.paleoPhylo}}.
}

\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}

\seealso{\code{\link{as.paleoPhylo}}, \code{\link{drawLinPhylo}}}

\examples{
##load bird order data from ape
data(bird.orders)
bird.orders
brdOrd <- ape2paleoPhylo(bird.orders)
brdOrd

##generate a plot using drawLinPhylo
data(gradstein04)
drawLinPhylo(brdOrd,l2r=TRUE,nmLim=1)

}

