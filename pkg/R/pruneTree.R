pruneTree <- function(pP, focLin, keepTips=TRUE, keepFocLin=FALSE, trimFocLin=NULL)
  {
  if (class(pP) != "paleoPhylo") stop(" object is not of class 'paleoPhylo'")
  
  nn <- length(pP$nm)
  p2r <- vector("list", nn)
  for(n in 1:nn) p2r[[n]] <- route2root(pP, pP$nm[n])$path
  
  kT <- ifelse(keepTips, 1, 0)
  kp   <- function(n) length(intersect(focLin, p2r[[n]]))==kT
  keep <- sapply(1:nn, kp)

  if(keepFocLin) keep[which(pP$nm==focLin)] <- TRUE
  if(!is.null(trimFocLin))
    {
    pP$st[which(pP$nm==focLin)] <- trimFocLin
    trim <- pP$st<=trimFocLin
    keep <- keep & trim
    }
    
  prT  <- with(pP, data.frame(nm, pn, st, en, label, grp))
  prT  <- prT[keep,]
  prT  <- with(prT, as.paleoPhylo(nm, pn, st, en, label=label, grp=grp))

  if(keepTips)
    {#remove the ancestor of the oldest lineage as not in the pruned tree
    whr <- which(prT$st==max(prT$st))
    prT$pn[whr] <- NA 
    }

  return(prT)
  }
  
