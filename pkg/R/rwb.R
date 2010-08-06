rwb <- function(pP, bL=1, st=max(pP$st), en=min(pP$en))
  {
  if(class(pP)!="paleoPhylo") stop("Object is not of class paleoPhylo")
  if (en >= st) stop("End must be later than start.")
  if (bL <= 0) stop("Bin Length must be positive.")
  #check there's an interval within which rates will be calculated
  
  # Lineages persisting to the specified end are not considered to have gone extinct.
  #Events taking place exactly on bin boundaries are taken to occur in the later bin.

  pPcB <- createBifurcate(pP)
  ppdf <- data.frame(nm=pPcB$nm, pn=pPcB$pn, st=pPcB$st, en=pPcB$en)

  if(length(bL)==1) bins <- seq(st, en, -bL)
  if(length(bL)>1) bins <- bL
  nb <- length(bins)
  rwb <- data.frame(binStart=bins, nStart=rep(NA, nb), branchLength=rep(NA, nb), orig=rep(NA, nb), extn=rep(NA, nb), lambda=rep(NA, nb), mu=rep(NA, nb))
 
  for(k in 1:(nb-1))
    {
    inBin <- ppdf[ppdf$st>bins[k+1] & ppdf$en<=bins[k],]
    nn  <- sum(inBin$st>bins[k])
    noDesc <- sapply(1:dim(inBin)[1], function(i) sum(as.character(inBin$pn)==inBin$nm[i], na.rm=TRUE))==0
    noRoot <- !is.na(inBin$pn)
    orig <- sum(inBin$st<=bins[k] & noRoot)/2
    extn <- sum(inBin$en>bins[k+1] & noDesc)

    inBin$st[inBin$st>bins[k]] <- bins[k]
    inBin$en[inBin$en<bins[k+1]] <- bins[k+1]
    brln <- sum(inBin$st-inBin$en)
    rwb[k,2:7] <- c(nn, brln, orig, extn, orig/brln, extn/brln) 
    }
  rwb <- rwb[-nb,]
  rwb <- rwb[rev(order(rwb$binStart)),]
  return(rwb)
  }

