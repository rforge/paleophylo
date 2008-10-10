\name{drawLinPhylo}
\alias{drawLinPhylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a lineage phylogeny.}
\description{Generates a cladogram that can distinguishe between budding and bifurcating speciations.  Further details on the distinction can be found in \code{\link{as.paleoPhylo}}.  Additional arguments enable clipping of time axes in variable formats and from different geological scales to the main plot.}
\usage{
drawLinPhylo(pP, 
	addTimeLines = "classic", tmScl = gradstein04, 
	whatTime = c("epoch","stage"), l2r=FALSE, arrowHeads=FALSE,
	nmLim=2, clCd = rep(1, length(pP$nm)), 
	certS = rep(NA, length(pP$nm)), certE = rep(NA, length(pP$nm)), 
	cexText = 0.5, cexTime = 0.5, cexLab=0.5,
	lwdLin = 2, hlty=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}).}
  \item{addTimeLines}{OPTIONAL: add a time axis to the plot.  Three options are accepted: 
  	"none" (plots no time axis); 
  	"classic" (default, plots a set of rectangular boxes with labels on the left-hand side of the cladogram);
  	"tube" (plots grey and white rectangles underneath the cladogram, as London tube maps plot travelcard zones).}
  \item{tmScl}{OTIONAL: some time scale for the time axis.  Defaults to the "Gradstein 2004" Geologic Time Scale (see \code{\link{gradstein04}}) but can take *any* user-read in data frame (see details below).}
  \item{whatTime}{OTIONAL: a vector of column names to be plotted in the time axis on the plot.  Defaults to 'epoch' and 'stage'.}
  \item{l2r}{OTIONAL: should the time axis run from left to right or top to bottom?  Defaults to FALSE, i.e. top to bottom.}
  \item{arrowHeads}{OTIONAL: should arrowheads be drawn on the tips of extant species?  Defaults to FALSE.}
 \item{nmLim}{OPTIONAL: a number of million years that serves as the limit for when the name of the period/epoch/stage/etc... should be printed in the 'classic' time axis.  Defaults to 2.}

  \item{clCd}{OPTIONAL: a vector with color codes for each species lineage.  If the length of clCd is not the same as the number of species lineages, then this will be repeated as outlined in \code{\link{par}}.  Type colors() to see the full range.  Defaults to black.}
  \item{certS}{OPTIONAL: start point of uncertainty in lineage persistence.  Defaults to none.  Use with ....}
  \item{certE}{OPTIONAL: end point of uncertainty in lineage persistence.  Defaults to none.  For an example of these two arguments, see example 3 below.}
  \item{cexText}{OPTIONAL: change the text size of time labels.  See \code{\link{par}}.}
  \item{cexTime}{OTIONAL: change the text size of time numbers.  See \code{\link{par}}.}
  \item{cexLab}{OTIONAL: change the text size of species lineage names.  See \code{\link{par}}.}
  \item{lwdLin}{OTIONAL: change the line thickness.  See \code{\link{par}}.}
  \item{hlty}{OTIONAL: change the horizontal line type.  Takes values for lty in \code{\link{par}}.}
}
\details{
In essence, all this does is generate a blank window and put some lines on top.  As such, its use is not restricted to 'paleoPhylo' objects.  It can be used in conjunction with fully bifurcating trees (\code{\link{ape2paleoPhylo}} can convert an object of class 'phylo' to class 'paleoPhylo' using  \code{\link{as.paleoPhylo}}) to add time axes in multiple forms to plots, and in conjunction with the pdf graphics device (see \code{\link{pdf}}) to output high-quality figures in pdf format (other graphics formats are available; see links from the pdf page).

Currently, any USER READ-IN TIME SCALE MUST HAVE a column called "MA", which contains age estimates for each period/epoch/stage/etc....  See  and copy \code{\link{gradstein04}} or \code{\link{berggren95}} for ways of *ensuring* your data frame will read in correctly.

\code{addTimeLines=="classic"} uses \code{\link{split.screen}}.  The left-hand portion containing the axis has two fixed widths depending on the number of species lineages to be drawn.  In due course, its width will become fully dynamic.
}
\value{
The function opens a graphics window and plots the lineage phylogeny (see \code{\link{as.paleoPhylo}} for details).  It does not return a list, since all the information is contained within input data.}
\references{\url{www.imperial.ac.uk/research/apurvis/paleoPhylo/}}
\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}
\seealso{\code{\link{as.paleoPhylo}}, \code{\link{gradstein04}}, \code{\link{pearson93}}, \code{\link{berggren95}}}
\examples{
data(pearson93)
data(gradstein04)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

##example 1: without a time axis.
drawLinPhylo (p93,addTimeLine="none")

##example 2: with a "classic" time axis.
drawLinPhylo (p93)

##example 3: with a "tube" time axis.
drawLinPhylo (p93,addTimeLine="tube",cexText=3,cexTime=0.8)

##example 4: with a "tube" time axis and different species' names.
drawLinPhylo (p93,addTimeLine="tube",cexText=3,cexTime=0.8)

##example 5: with a "classic" time axis and hypothetical uncertainty.
##suppose (FOR ILLUSTRATION ONLY) that uncertainty exists surrounding lineage "91".
##let this begin 2 years after its origination & ...
##... persist until 4 years before its disapperance/extinction.
certS <- certE <- rep(NA,length(p93$nm))
clCd <- rep("black",length(p93$nm))
certS[which(p93$nm == "91")] <- 2
certE[which(p93$nm == "91")] <- 4
clCd[which(p93$nm == "91")] <- "grey80"
drawLinPhylo (p93,certS=certS,certE=certE,clCd=clCd,cexTime=0.8)

##example 6: with species arranged on the x-axis according ...
## .... to some random trait data from a unifrom distribution.
##Also change the time-axis scale using whatTime.
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End,
	 xx=runif(length(Code),-1,1)))
drawLinPhylo (p93,addTimeLine="tube",
	cexText=3,cexTime=0.8,whatTime=c("period","epoch"))
}
\keyword{utilities}{drawLinPhylo}

