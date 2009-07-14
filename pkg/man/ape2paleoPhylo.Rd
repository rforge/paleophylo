\name{ape2paleoPhylo}
\alias{ape2paleoPhylo}

\title{Convert objects of class 'phylo' to class 'paleoPhylo'.}
\description{
Takes the information contained within objects of class \code{phylo} and converts using \code{\link{as.paleoPhylo}}.
}
\usage{
ape2paleoPhylo(phy,retainNodeLabels=TRUE,nC=0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{phy}{An object of class phylo.}
  \item{retainNodeLabels}{If FALSE, then \code{\link{nchar}} is used to replace node labels with blank names.  Defaults to TRUE.}
  \item{nC}{The boundary of number of characters in tip.labels to retain. Defaults to 0, which means that whole tip labels are retained. Setting \code{nC} to 5 would mean that any label longer than 5 characters would be stripped.}
}
  
\details{
The principal current use is to enable plotting with \code{\link{drawPhylo}} and its array of axis formats.

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
drawPhylo(brdOrd, tmScl=gradstein04, whatTime=c("period","epoch","stage"),
	addTimeLine="c", l2r=TRUE, nmLim=1)

##that's no use though, since all the names are too long for each box
##how do I shorten those names to, say, three letters long?
##generates a new column in gradstein04 containing the shorthand
gradstein04$shStg <- substr(gradstein04$stage,1,3)
str(gradstein04)
drawPhylo(brdOrd, tmScl=gradstein04, whatTime=c("period","epoch","shStg"),
	addTimeLine="c", l2r=TRUE, nmLim=1, cexTime=0.7, cexText=1)

}

