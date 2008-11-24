\name{paleoPhylo-package}
\alias{paleoPhylo-package}
\alias{paleoPhylo}
\keyword{package}{paleoPhylo}

\docType{package}
\title{Biodiversity Analyses in a Paleontological and Phylogenetic Context}
\description{
paleoPhylo contains a suite of functions to analyze how speciation, extinction & character change contribute to biodiversity as well as tools to display these relationships graphically.  The functions are designed to inform realistic models of clade evolution; examples are drawn from data on planktonic foraminifera or generated randomly.  The package is STILL IN DEVELOPMENT.  It will be added to in due course.
}
\details{
\tabular{ll}{
Package: \tab paleoPhylo\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-09-02\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Thomas H.G. Ezard 

Maintainer: Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}
}
\references{\url{www.bio.ic.ac.uk/research/apurvis/paleoPhylo.htm}}
\keyword{ package }
\seealso{
	\code{\link[ape:ape-package]{ape}}, \code{\link{pearson93}}, \code{\link{buildApe}}, \code{\link{createBifurcate}}, \code{\link{drawLinPhylo}}, \code{\link{ape2paleoPhylo}}
}
\examples{
data(pearson93)
data(berggren95)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))
drawLinPhylo(p93, tmScl=berggren95, whatTime=c("epoch","age"), addTimeLine="t")

apeTree <- buildApe(createBifurcate(p93))
apeTree
}
