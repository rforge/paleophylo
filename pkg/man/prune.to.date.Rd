\name{prune.to.date}
\alias{prune.to.date}

\title{
Prunes an object of class "paleoPhylo" to a particular cutoff date
}

\description{
If passed an object of class "paleoPhylo" to a particular cutoff date, then returns the phylogeny of just the species extant at the cutoff. The phylogeny is returned in two formats: paleoPhylo and phylo.
}

\usage{
prune.to.date(pP, date, let.speciate = FALSE, let.die = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}).}
  \item{date}{The cutoff date of choice.}
  \item{let.speciate}{If TRUE, speciations at exactly the cutoff date are taken to have happened.}
  \item{let.die}{IF TRUE, extinctions at exactly the cutoff date are taken to have happened.}
  }
\details{
Enables an event-based protocol to be implemented, through iterative application along a phylogeny.
}

\value{
  \item{paleoPhylo.tree}{The pruned tree as an object of class "paleoPhylo".}
  \item{phylo.tree}{The pruned tree as an object of class "phylo".}
}

\references{Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.}

\author{Andy Purvis \email{a.purvis@imperial.ac.uk}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{as.paleoPhylo}}, \code{\link{pearson93}} }
\examples{

data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))
par(mfrow=c(2,2))
plot(prune.to.date(p93, 54.3)$phylo.tree, cex=0.5,main = "Default")
plot(prune.to.date(p93, 54.3, let.speciate=TRUE)$phylo.tree, cex=0.5,main="Let speciate")
plot(prune.to.date(p93, 54.3, let.die=TRUE)$phylo.tree, cex=0.5,main = "Let die")
plot(prune.to.date(p93, 54.3, let.speciate = TRUE, let.die=TRUE)$phylo.tree, cex=0.5,main = "Live and let die")
par(mfrow=c(1,1))

}
