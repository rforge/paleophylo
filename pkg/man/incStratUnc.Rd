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
par(mfrow=c(3,1))
##example 1; no change
wag <- data.frame(nm=LETTERS[1:8],pn=c("F","G","G","E","H",NA,"H","F"),st=c(3.1,2,2,1,3,3.2,3,3.1),en=c(0,1.1,0.8,0,0,3.1,2,3),label=c(LETTERS[1:5],"","",""))
wag <- wag[rev(order(wag$st)),]
wag$label <- c("","","E","","D","A","B","C")

pP <- as.paleoPhylo(wag$nm,wag$pn,wag$st,wag$en,label=wag$label)

drawLinPhylo(pP,cexLab=1.5,hlty=1)

##example 2; no change
uncStratRange <- incStratUnc(NULL,pP)
#change A

uncStratRange$dates[[3]] <- c(uncStratRange$dates[[3]][1],2,1,uncStratRange$dates[[3]][2])
uncStratRange$types[[3]] <- c(2,1,2)

#change C

uncStratRange$dates[[6]] <- c(uncStratRange$dates[[6]][1],1.1,uncStratRange$dates[[6]][2])
uncStratRange$types[[6]] <- c(3,1)

pP <- as.paleoPhylo(wag$nm,wag$pn,wag$st,wag$en,label=wag$label)
drawLinPhylo(pP,uSR=uncStratRange,cexLab=1.5, addTimeLine="s", cexText=1.5)

##what if that grey is a little light for you?
uncStratRange$styles$possible[2] <- "grey20"

drawLinPhylo(pP,uSR=uncStratRange,cexLab=1.5, addTimeLine="s", cexText=1.5)

}

