\name{pruneTree}
\alias{pruneTree}
\alias{route2root}
\alias{collapseBranches}

\title{Prune a tree.}
\description{Prunes a tree to retain only part of it, either the clipping, i.e. tips, or the bush, i.e. the roots. It is intended to work like drop.tip and subclade in \code{ape}.}

\usage{
pruneTree(pP, focLin=NULL, date=NULL,
  keepTips=TRUE, keepFocLin=TRUE,
  letSpeciate=FALSE, letDie=FALSE, pruneDead=TRUE,
  outPhylo=FALSE, collapseBranches=TRUE)
route2root(pP, focLin)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}). }
\item{focLin}{The focal lineage to be pruned or to obtain the \code{route2root} from.}
\item{date}{asas}
\item{keepTips}{If TRUE (default), then keep the tips from \code{focLin}.}
\item{keepFocLin}{If TRUE (default), then keep the focal lineage, i.e. \code{focLin}, in the pruned tree.}
\item{letSpeciate}{Truncate the focal lineage to start at \code{trimFocLin} if \code{keepTips}==TRUE.}
\item{letDie}{Truncate the focal lineage to start at \code{trimFocLin} if \code{keepTips}==TRUE.}
\item{pruneDead}{Truncate the focal lineage to start at \code{trimFocLin} if \code{keepTips}==TRUE.}
\item{outPhylo}{Convert to an \code{ape} object in the output?}
\item{collapseBranches}{If the pruned tree leaves successive distinct ancestors with persistent but no cladogenesis, then collapse these branches down to a single lineage?}
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

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*example 1: around a focal lineage
par(mfrow=c(2,2))
drawPhylo(p93, l2r=TRUE)
#prune this tree around S. triangularis
p93$nm[which(p93$label=="S.triangularis")]
pTt <- pruneTree(p93, "94") #keep the tips
drawPhylo(pTt, l2r=TRUE)
pTr <- pruneTree(p93, "94", keepTips=FALSE) #keep the rest
drawPhylo(pTr, l2r=TRUE)
#pTrTt <- pruneTree(p93, "94", 56, keepTips=FALSE) #keep the rest
#drawPhylo(pTr, l2r=TRUE)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*example 2: around a focal date
par(mfrow=c(2,4))
eg1 <- pruneTree(p93, date=54.3, keepTips=FALSE)
drawPhylo(eg1, l2r=TRUE, whSpc=0.2, cexLab=1)
eg2 <- pruneTree(p93, date=54.3, keepTips=FALSE, pruneDead=FALSE)
drawPhylo(eg2, l2r=TRUE, whSpc=0.2, cexLab=1)
eg3 <- pruneTree(p93, date=54.3, keepTips=FALSE, letSpeciate=TRUE, pruneDead=TRUE)
drawPhylo(eg3, l2r=TRUE, whSpc=0.2, cexLab=1)
eg4 <- pruneTree(p93, date=54.3, keepTips=FALSE, letSpeciate=FALSE, pruneDead=TRUE)
drawPhylo(eg4, l2r=TRUE, whSpc=0.2, cexLab=1)
eg5 <- pruneTree(p93, date=54.3, keepTips=FALSE, letDie=TRUE, pruneDead=TRUE)
drawPhylo(eg5, l2r=TRUE, whSpc=0.2, cexLab=1)
eg6 <- pruneTree(p93, date=54.3, keepTips=FALSE, letDie=TRUE, pruneDead=FALSE)
drawPhylo(eg6, l2r=TRUE, whSpc=0.2, cexLab=1)
eg7 <- pruneTree(p93, date=54.3, keepTips=FALSE, letSpeciate=TRUE, letDie=TRUE, pruneDead=TRUE)
drawPhylo(eg7, l2r=TRUE, whSpc=0.2, cexLab=1)
eg8 <- pruneTree(p93, date=54.3, keepTips=FALSE, letSpeciate=TRUE, letDie=TRUE, pruneDead=FALSE)
drawPhylo(eg8, l2r=TRUE, whSpc=0.2, cexLab=1)
}

