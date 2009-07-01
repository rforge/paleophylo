\name{timeAxis}
\alias{timeAxis}

\title{
Plot a time axis to append to, e.g., diversity curve plots.
}

\description{
Works in the same way as the \code{addTimeLine} part of \code{\link{drawPhylo}}, only ....
}
\usage{
timeAxis(tmScl, whatTime, addTimeLine = "classic", l2r = FALSE, nmLim = 2, cexText = 0.5, srtText = 0, cexTime = 0.5, cexLab = 0.5, whSpc = 0.05, lwdLin = 1, dumpLast = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmScl}{The time scale for the time axis.}
  \item{whatTime}{Which columns of the time scale should be plotted}
  \item{addTimeLine}{Two options (or only their first letter) are accepted:   	"classic" (default, plots a set of rectangular boxes with labels on the left-hand side of the cladogram);  	"tube" (plots grey and white rectangles underneath the cladogram, as London tube maps plot travelcard zones).}
  \item{l2r}{Should the time axis run from left to right or top to bottom?  Defaults to FALSE, i.e. top to bottom.} \item{nmLim}{A number of million years that serves as the limit for when the name of the period/epoch/stage/etc... should be printed in the 'classic' time axis.  Defaults to 2.}  \item{cexText}{Change the text size of time labels.  See \code{\link{par}}. Can be a vector for different text sizes for each time classification (see example 2a below).}  \item{srtText}{Change the text direction of time labels.  See \code{\link{par}}. Can be a vector for different rotations for each time classification. It probably makes sense to stick to either 0 (parallel to phylogeny) or 90 (perpendicular to phylogeny).}  \item{cexTime}{Change the text size of time numbers.  See \code{\link{par}}.}  \item{cexLab}{Change the text size of species lineage names.  See \code{\link{par}}.}
 \item{whSpc}{Controls the amount of white space around the phlyogram.  Increasing values of \code{whSpc} will yield more white space, which may be useful when formatting smaller phylograms or if the tip labels are particularly long. It works using the \code{f} in \code{\link{extendrange}}.  Note that \code{mar()} in \code{\link{extendrange}} is not supported.}  \item{lwdLin}{Change the line thickness.  See \code{\link{par}}.}
  \item{dumpLast}{If \code{addTimeLine=="classic" or "tube"} then should the oldest time box be drawn? Removing it might be desired if, eg., an ancestor is shown from an earlier time period but only just which would mean "messy" unclosed time boxes.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\section{Warning }{need to use close.screen()}

\references{\url{www.bio.ic.ac.uk/research/apurvis/paleoPhylo.htm}}\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}\seealso{\code{\link{drawPhylo}}, \code{\link{as.paleoPhylo}}, \code{\link{gradstein04}}, \code{\link{pearson93}}, \code{\link{berggren95}}}

\examples{
##draw a diversity curve for pearson93 phylogeny
data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

data(berggren95)
tmScl <- berggren95
tmScl <- tmScl[tmScl$MA>=min(p93$en),]
tmScl <- tmScl[tmScl$MA<=max(p93$st),]

dvrstyTime <- seq(max(p93$st),min(p93$en), -0.2)
abndncTime <- numeric(length(dvrstyTime))
for (k in 1:length(dvrstyTime))	{abndncTime[k] <- sum(p93$st>dvrstyTime[k] & p93$en<=dvrstyTime[k])}

timeAxis(tmScl, c("epoch", "age"), l2r=TRUE, srtText=c(0,90,0))
lines(-dvrstyTime ,abndncTime/max(abndncTime), lwd=2)

}


