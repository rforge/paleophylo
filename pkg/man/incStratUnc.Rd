\name{incStratUnc}
\alias{incStratUnc}

\title{Include stratigraphic uncertainty in paleoPhylo plots.}
\description{
...}
\usage{
incStratUnc(uSR=NULL, pP, lwdLin=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{uSR}{A list (see below).  Defaults to NULL.}
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}). }
  \item{lwdLin}{defaults to 2; base line thickness.}
}
\details{
If a non-null list is passed to \code{incStratUnc}, this function checks its suitability.  Specifically, are all types and dates congruent and are all types and styles congruent?

If uSR is NULL, then a default list is created for pP.  The default styles are...... These can be edited and added to by any user; an example is given below.
}

\value{
A list containing:}

\references{
Wagner 1995
}
\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}

\seealso{\code{\link{createBifurcate}}, \code{\link{drawLinPhylo}}}

\examples{
eg1 <- data.frame(nm=paste("id",1:8,sep=""),pn=c("id6","id7","id7","id5","id8",NA,"id8","id6"),st=c(3.1,2,2,1,3,3.2,3,3.1),en=c(0,1.1,0.8,0,0,3.1,2,3),		label=c("E","B","A","C","D","","",""))
pP <- as.paleoPhylo(eg1$nm,eg1$pn,eg1$st,eg1$en,label=eg1$label)
par(mfrow=c(2,1))
##example 1; no uncertainty
drawLinPhylo(pP,cexLab=1.5, addTimeLine="s", cexText=1.5,hlty=3)
##example 1; hypothetical uncertainty##this is figure 2 in Ezard & Purvis (prov. accpt).
unc <- data.frame(id=c("id1","id1","id3","id4","id4","id2","id2"),st=c(3.1,1,2,1.1,0.1,2,1.7),en=c(2,0,1.1,0.9,0,1.1,1.7),type=c(2,2,3,4,4,2,5))

uSR <- incStratUnc(unc,pP)

drawLinPhylo(pP,uSR,cexLab=1.5, addTimeLine="s", cexText=1.5,hlty=3)

}

