route2root <- function(pP, focLin) 
  {
  if (class(pP) != "paleoPhylo") stop(" objefocLin is not of class 'paleoPhylo'")
  
  x <- which(pP$st==max(pP$st))
  rtDur <- pP$st[x]-pP$en[x]
  rt <- pP$nm[x]
  
  allAnc <- immAnc <- focLin
  dur <- pP$st - pP$en
  tm2rt <- c()
  
  while(immAnc != rt)
    {
    i <- which(pP$nm == immAnc)
    tm2rt <- c(tm2rt, dur[i])
    immAnc <- pP$pn[i]
    allAnc <- c(allAnc,immAnc)
    }
  tm2rt <- c(tm2rt, rtDur)
  if(length(allAnc)!=length(tm2rt)) stop ("Number of lineages does not match number of durations.")
  
  out <- list(path=rev(allAnc), duration=rev(tm2rt), nNode=length(allAnc))
  return(out)
  }