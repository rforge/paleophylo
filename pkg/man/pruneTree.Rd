\name{pruneTree}
\alias{pruneTree}
\alias{route2root}

\title{Prune a tree.}
\description{Prunes a tree to retain only part of it, either the clipping, i.e. tips, or the bush, i.e. the roots. It is intended to work like drop.tip and subclade in \code{ape}.}
\usage{
pruneTree(pP, focLin, keepTips=TRUE, keepFocLin=FALSE, trimFocLin=NULL)
route2root(pP, focLin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}). }
\item{focLin}{The focal lineage to be pruned or to obtain the \code{route2root} from.}
\item{keepTips}{If TRUE (default), then keep the tips from \code{focLin}.}
\item{keepFocLin}{If TRUE (default), then keep the focal lineage, i.e. \code{focLin}, in the pruned tree.}
\item{trimFocLin}{Truncate the focal lineage to start at \code{trimFocLin} if \code{keepTips}==TRUE.}
}
\details{
yes, and also of route2root
}

\value{
An object of class \code{\link{as.paleoPhylo}}.
}

\references{
Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.}
\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}

\seealso{\code{\link{route2root}}, \code{\link{as.paleoPhylo}}, \code{\link{pearson93}}, \code{\link{berggren95}}, \code{\link{drawPhylo}}}

\examples{
data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

par(mfrow=c(3,1))
drawPhylo(p93, l2r=TRUE)

#prune this tree around S. triangularis
p93$nm[which(p93$label=="S.triangularis")]

pTt <- pruneTree(p93, "94") #keep the tips
drawPhylo(pTt, l2r=TRUE)
pTr <- pruneTree(p93, "94", keepTips=FALSE) #keep the rest
drawPhylo(pTr, l2r=TRUE)

}

