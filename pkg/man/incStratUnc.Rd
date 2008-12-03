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
eg1 <- data.frame(nm=paste("id",1:8,sep=""), pn=c("id6","id7","id7","id5","id8",NA,"id8","id6"), st=c(3.1,2,2,1,3,3.2,3,3.1), en=c(0,1.1,0.8,0,0,3.1,2,3), label=c(LETTERS[1:5],"","",""))
eg1 <- eg1[rev(order(eg1$st)),]
eg1$label <- c("","","E","","D","A","B","C")
pP <- as.paleoPhylo(eg1$nm,eg1$pn,eg1$st,eg1$en,label=eg1$label)

par(mfrow=c(2,1))
##example 1; no uncertainty
drawLinPhylo(pP,cexLab=1.5, addTimeLine="s", cexText=1.5,hlty=3)
##example 1; hypothetical uncertainty
##this is figure 2 in Ezard & Purvis (prov. accpt).
uSR <- incStratUnc(NULL,pP)
#change A
uSR$dates[[3]] <- c(uSR$dates[[3]][1],2,1,uSR$dates[[3]][2])
uSR$types[[3]] <- c(2,1,2)
#change E
uSR$dates[[6]] <- c(uSR$dates[[6]][1],1.1,uSR$dates[[6]][2])
uSR$types[[6]] <- c(3,1)
#change C
uSR$dates[[8]] <- c(1.05,0.9,0.1,uSR$dates[[8]][2])
uSR$types[[8]] <- c(4,1,4)

drawLinPhylo(pP,uSR=uSR,cexLab=1.5, addTimeLine="s", cexText=1.5,hlty=3)

}

