\name{berggren95}
\alias{berggren95}
\alias{brggrnprsn05}

\docType{data}
\title{W.A. Berggren's 1995 geologic timescale}
\description{W.A. Berggren's 1995 geologic timescale for the Cenozoic, giving (sub)tropical planktonic foraminiferal zonation. Additionally, the Berggren and Pearson (2005) revision for the Paleogene in the same format.}
\usage{
data(berggren95)
data(brggrnprsn05)
}
\format{
 A data frame containing the following columns: 
  \describe{
    \item{\code{epoch}}{Epoch name.}
    \item{\code{zone}}{(Sub)tropical planktonic foraminifera zone code.}
    \item{\code{MA}}{Estimated boundary age for the associated interval.}
 }
}

\source{
Berggren, W.A., Kent, D.V., Swisher III, C.C. & Aubry, M.-P., 1995.	A Revised Cenozoic Geochronology and Chronostratigraphy, in Berggren, W.A., Kent, D.V., Swisher III, C.C., Aubry, M.-P. & Hardenbol, J., (eds), Geochronology, Time Scales and Global Stratigraphic Correlation: SEPM (Society for Sedimentary Geology), Special Publication No. 54, p. 129-212.

Berggren, W.A. & Pearson, P.N., 2005. A revised tropical to subtropical Paleogene planktonic foraminiferal zonation. Journal of Foraminiferal Research, 35: 279Ð298
}
 
\section{
}

\seealso{\code{\link{gradstein04}}, \code{\link{drawPhylo}}}

\examples{
data(berggren95)
str(berggren95)

data(brggrnprsn05)
str(brggrnprsn05)
}
\keyword{datasets}
