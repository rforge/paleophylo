\name{stratUnc}
\alias{stratUnc}

\title{Include stratigraphic uncertainty in paleoPhylo plots.}
\description{
...}
\usage{
stratUnc(uSR=NULL, pP, lwdLin=2, clCd=1:length(pP$nm), ltyLin=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{uSR}{A list (see below).  Defaults to NULL.}
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}). }
  \item{lwdLin}{defaults to 2; base line thickness (is recycled).}
  \item{clCd}{Defaults to black; base color (is recycled).}
  \item{ltyLin}{defaults to 1; base line type (is recycled).}
}
\details{
If a non-null list is passed to \code{stratUnc}, this function checks its suitability.  Specifically, are all types and dates congruent and are all types and styles congruent?
      
If uSR is NULL, then a default list is created for pP.  The default styles are...... These can be edited and added to by any user; an example is given below.
}

\value{
A list containing:}

\references{
To calculate confidence intervals on lineage duration, see: (for example):

Marshall, C. R. 1994. Confidence intervals on stratigraphic ranges: partial relaxation of the assumption of randomly distributed fossil horizons. Paleobiology 20:459-469.

Solow, A. R. 1996. Tests and confidence intervals for a common upper endpoint in fossil taxa. Paleobiology 22:406-410.

Strauss, D., and P. M. Sadler. 1989. Classical confidence intervals and Bayesian probability estimates for ends of local taxon ranges. . Mathematical Geology 21:411-427.

Wagner, P. J. 2000. Likelihood tests of hypothesized durations: determining and accommodating biasing factors. Paleobiology 26:431-449.

Wang, S. C., and C. R. Marshall. 2004. Improved confidence intervals for estimating the position of a mass extinction boundary. Paleobiology 30:5-18.

}
\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}

\seealso{\code{\link{createBifurcate}}, \code{\link{drawPhylo}}}

\examples{
eg1 <- data.frame(nm=paste("id",1:8,sep=""),pn=c("id6","id7","id7","id5","id8",NA,"id8","id6"),st=c(3.1,2,2,1,3,3.2,3,3.1),en=c(0,1.1,0.8,0,0,3.1,2,3),		label=c("E","B","A","C","D","","",""))
pP <- as.paleoPhylo(eg1$nm,eg1$pn,eg1$st,eg1$en,label=eg1$label)
par(mfrow=c(2,1))
##example 1; no uncertainty
drawPhylo(pP,cexLab=1.5, addTimeLine="s", cexText=1.5,hlty=3)
##example 1; hypothetical uncertainty##this is figure 2 in Ezard & Purvis (prov. accpt).
unc <- data.frame(id=c("id1","id1","id3","id4","id4","id2","id2"),st=c(3.1,1,2,1.1,0.1,2,1.7),en=c(2,0,1.1,0.9,0,1.1,1.7),type=c(2,2,3,4,4,2,5))

uSR <- stratUnc(unc,pP)

drawPhylo(pP,uSR,cexLab=1.5, addTimeLine="s", cexText=1.5,hlty=3)

}

