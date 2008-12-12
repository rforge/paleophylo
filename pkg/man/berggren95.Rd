\name{berggren95}
\alias{berggren95}
\docType{data}
\title{W.A. Berggren's 1995 geologic timescale}
\description{W.A. Berggren's 1995 geologic timescale}
\usage{data(berggren95)}
\format{
 A data frame containing the following columns: 
  \describe{
    \item{\code{eon}}{Eon name.}
    \item{\code{era}}{Era name.}
    \item{\code{period}}{Period name.}
    \item{\code{epoch}}{Epoch name.}
    \item{\code{age}}{Age name.}
    \item{\code{MA}}{Estimated boundary age for the associated interval.}
 }
}
\details{
  ...
  }
\source{
	Berggren, W.A., Kent, D.V., Swisher III, C.C. & Aubry, M.-P., 1995.	A Revised Cenozoic Geochronology and Chronostratigraphy, in Berggren, W.A., Kent, D.V., Swisher III, C.C., Aubry, M.-P. & Hardenbol, J., (eds), Geochronology, Time Scales and Global Stratigraphic Correlation: SEPM (Society for Sedimentary Geology), Special Publication No. 54, p. 129-212.
  
  
  \url{www.chronos.org/}
  }
  
\seealso{\code{\link{gradstein04}}, \code{\link{drawPhylo}}}

\examples{
data(berggren95)
str(berggren95)
}
\keyword{datasets}
