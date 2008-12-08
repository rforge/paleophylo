\name{drawLinPhylo}
\alias{drawLinPhylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a lineage phylogeny.}
\description{Generates a cladogram that can distinguishe between budding and bifurcating speciations.  Further details on the distinction can be found in \code{\link{as.paleoPhylo}}.  Additional arguments enable clipping of time axes in variable formats and from different geological scales to the main plot.}
\usage{
drawLinPhylo(pP, 
	uSR = NULL,
	addTimeLine = "none", tmScl, 
	whatTime = c("epoch","age"), l2r=FALSE,
	nmLim=2, cexText = 0.5, cexTime = 0.5, cexLab=0.5,
	lwdLin = 2, hlty=NULL, barLen=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}). }
  \item{uSR}{A list generated through \code{\link{incStratUnc}} to incorporate stratigraphic uncertainty.  Defaults to NULL. }
  \item{addTimeLine}{OPTIONAL: add a time axis to the plot.  Four options (or only their first letter) are accepted: 
  	"none" (plots no time axis); 
  	"classic" (default, plots a set of rectangular boxes with labels on the left-hand side of the cladogram);
  	"tube" (plots grey and white rectangles underneath the cladogram, as London tube maps plot travelcard zones);
  	"simple" (plots an arrow in the direction of the time axis)}
  \item{tmScl}{OPTIONAL: some time scale for the time axis.  The "Gradstein 2004" (see \code{\link{gradstein04}}) and "Berggren 1995" (see \code{\link{berggren95}}) Geologic Time Scales are included, but this argument can take *any* user-read in data frame (see details below).  Required if addTimeLine is either classic or tube.}
  \item{whatTime}{OPTIONAL: a vector of column names to be plotted in the time axis on the plot.  Defaults to 'epoch' and 'age', as used in \code{\link{berggren95}}.}
  \item{l2r}{OPTIONAL: should the time axis run from left to right or top to bottom?  Defaults to FALSE, i.e. top to bottom.}
 \item{nmLim}{OPTIONAL: a number of million years that serves as the limit for when the name of the period/epoch/stage/etc... should be printed in the 'classic' time axis.  Defaults to 2.}

  \item{cexText}{OPTIONAL: change the text size of time labels.  See \code{\link{par}}.}
  \item{cexTime}{OPTIONAL: change the text size of time numbers.  See \code{\link{par}}.}
  \item{cexLab}{OPTIONAL: change the text size of species lineage names.  See \code{\link{par}}.}
  \item{lwdLin}{OPTIONAL: change the line thickness.  See \code{\link{par}}.}
  \item{hlty}{OPTIONAL: change the horizontal line type if not otherwise specified using \code{\link{incStratUnc}}.  Takes values for lty in \code{\link{par}}.}
 \item{barLen}{The length of the horizontal bars if confidence intervals are being drawn using \code{\link{incStratUnc}}.}
}
\details{
In essence, all this does is generate a blank window and put some lines on top.  As such, its use is not restricted to 'paleoPhylo' objects.  It can be used in conjunction with fully bifurcating trees (\code{\link{ape2paleoPhylo}} can convert an object of class 'phylo' to class 'paleoPhylo' using  \code{\link{as.paleoPhylo}}) to add time axes in multiple forms to plots, and in conjunction with the pdf graphics device (see \code{\link{pdf}}) to output high-quality figures in pdf format (other graphics formats are available; see links from the pdf page).

Currently, any USER READ-IN TIME SCALE MUST HAVE a column called "MA", which contains age estimates for each period/epoch/stage/etc....  See  and copy \code{\link{gradstein04}} or \code{\link{berggren95}} for ways of *ensuring* your data frame will read in correctly.

\code{addTimeLines=="classic"} uses \code{\link{split.screen}}.  The left-hand portion containing the axis has two fixed widths depending on the number of species lineages to be drawn.  In due course, its width will become fully dynamic.
}
\value{
The function opens a graphics window and plots the lineage phylogeny (see \code{\link{as.paleoPhylo}} for details).  It does not return a list, since all the information is contained within input data.}
\references{\url{www.bio.ic.ac.uk/research/apurvis/paleoPhylo.htm}}
\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}
\seealso{\code{\link{as.paleoPhylo}}, \code{\link{gradstein04}}, \code{\link{pearson93}}, \code{\link{berggren95}}}
\examples{
data(pearson93)
data(berggren95)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

##example 1: without a time axis.
drawLinPhylo (p93)

##example 2: with a "classic" time axis.
drawLinPhylo (p93, tmScl=berggren95, addTimeLine="classic")

##example 3: with a "tube" time axis.
drawLinPhylo (p93, tmScl=berggren95, addTimeLine="tube",cexText=3,cexTime=0.8)

##example 4: with a "tube" time axis and different species' names.
drawLinPhylo (p93, tmScl=berggren95, whatTime=c("period","epoch"), addTimeLine="t",cexText=3,cexTime=0.8)

#example 5: with species arranged on the x-axis according ...
## .... to some random trait data from a unifrom distribution.
##Also change the time-axis scale using whatTime.
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End,
	 xx=runif(length(Code),-1,1)))
drawLinPhylo (p93, tmScl=berggren95, addTimeLine="c", cexTime=0.8)
}
\keyword{utilities}{drawLinPhylo}

