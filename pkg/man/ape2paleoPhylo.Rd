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
  \item{phy}{An object of class phylo.}
  \item{retainNodeLabels}{If FALSE, then \code{\link{nchar}} is used to replace node labels with blank names.  Defaults to TRUE.}
  \item{nC}{number of characters to retain.}
}
  
\details{
The principal use is to enable plotting with \code{\link{drawPhylo}} and its array of axis formats.

If tip.labels are included, then they are exported as \code{labels} in \code{paleoPhylo}.  Non-tip labels are assumed to be the node label.  Future revisions will provide the option to change this.

Objects of class \code{phylo} are defined using edge matrices (see \code{\link[ape]{read.tree}}); \code{\link{paleoPhylo}} objects require a unique root that gave rise to the focal clade.  The unique root is assumed to have geologically infinitessimal duration, such that it is not visible on plots.  Future releases will incorporate analytical tools; at this stage additional functionality will be introduced to adapt this.

}

\value{
An object of class 'paleoPhylo', which is a list containing the vectors defined in \code{\link{as.paleoPhylo}}.
}

\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}
\references{\url{www.bio.ic.ac.uk/research/apurvis/paleoPhylo.htm}}

\seealso{\code{\link{as.paleoPhylo}}, \code{\link{drawPhylo}}}

\examples{
##load bird order data from ape
data(bird.orders)
bird.orders
brdOrd <- ape2paleoPhylo(bird.orders)
brdOrd

##generate a plot using drawPhylo
data(gradstein04)
drawPhylo(brdOrd,l2r=TRUE)

}

