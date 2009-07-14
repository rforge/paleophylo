\name{timeAxis}
\alias{timeAxis}

\title{
Plot a time axis to append to, e.g., diversity curve plots.
}

\description{
Works in the same way as the \code{addTimeLine} part of \code{\link{drawPhylo}}, only ....
}
\usage{
timeAxis(tmScl, whatTime, addTimeLine = "classic", l2r = FALSE, nmLim = 2, cexText = 0.5, srtText = 0, cexTime = 0.5, cexLab = 0.5, whSpc = 0.1, lwdLin = 1, dumpLast = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmScl}{The time scale for the time axis.}
  \item{whatTime}{Which columns of the time scale should be plotted}
  \item{addTimeLine}{Two options (or only their first letter) are accepted: 
  \item{l2r}{Should the time axis run from left to right or top to bottom?  Defaults to FALSE, i.e. top to bottom.}
 \item{whSpc}{Controls the amount of white space around the phlyogram.  Increasing values of \code{whSpc} will yield more white space, which may be useful when formatting smaller phylograms or if the tip labels are particularly long. It works using the \code{f} in \code{\link{extendrange}}.  Note that \code{mar()} in \code{\link{extendrange}} is not supported.}
  \item{dumpLast}{If \code{addTimeLine=="classic" or "tube"} then should the oldest time box be drawn? Removing it might be desired if, eg., an ancestor is shown from an earlier time period but only just which would mean "messy" unclosed time boxes.}
}

\section{Warning }{need to use \code{close.screen()} after the plot if \code{addTimeLine="classic"} as \code{\link{split.screen}} is used to partition up the device.}

\references{\url{www.bio.ic.ac.uk/research/apurvis/paleoPhylo.htm}}

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
for (k in 1:length(dvrstyTime))

timeAxis(tmScl, c("epoch", "age"), l2r=TRUE, srtText=c(0,90,0), 
	cexTime=0.7, cexText=c(1.5,0.7), whSpc=0.12)
lines(-dvrstyTime ,abndncTime/max(abndncTime), lwd=2)
#currently the y-axis scales from 0 to 1, hence the use of proportions here for the diversity curve
abdnLab <- round(seq(0, max(abndncTime), length.out=5),0)
axis(2, at= seq(0,1,.25), labels=abdnLab, line=-3)
mtext("Abundance",2, line=-1, outer=FALSE, cex=1.5)
close.screen(2,all=TRUE)
##the 2 refers to the second panel, which is active and draws the (here) diversity curve

}

