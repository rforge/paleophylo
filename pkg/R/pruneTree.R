pruneTree <- function(pP, focLin=NULL, date=NULL,
  keepTips=TRUE, keepFocLin=TRUE,
  letSpeciate=FALSE, letDie=FALSE, pruneDead=FALSE,
  outPhylo=FALSE, collapseBranches=TRUE)
  {
  #pP <- p93; focLin <- "90" ; date <- 51; keepTips=FALSE; keepFocLin=TRUE
  #pruneDead <- TRUE; letSpeciate <- TRUE; letDie <- TRUE;outPhylo=FALSE
  #test <- pruneTree1(p93, "90", 51, keepTips=FALSE)$paleoPhylo
  #with(test, tapply(nm, pn, length))
  
  if (class(pP) != "paleoPhylo") stop(" object is not of class 'paleoPhylo'")
  nn <- length(pP$nm)
  p2r <- vector("list", nn)
  for(n in 1:nn) p2r[[n]] <- route2root(pP, pP$nm[n])$path
  prT  <- with(pP, data.frame(nm, pn, st, en, label, grp))

  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#*if pruning around a focal lineage
  if(!is.null(focLin))
    {
    kT <- ifelse(keepTips, 1, 0)
    kp   <- function(n) length(intersect(focLin, p2r[[n]]))==kT
    	#could be amended to include all descendants of this lineage, too
    	#that way can prune2end of lineage
    keep <- sapply(1:nn, kp)

    if(keepFocLin) keep[which(pP$nm==focLin)] <- TRUE
    prT  <- prT[keep,]

    if(keepTips)
      {
      #remove the ancestor of the oldest lineage as not in the pruned tree
      whr <- which(prT$st==max(prT$st))
      prT$pn[whr] <- NA 
      }

     if(pruneDead)
        {
        if(letDie) extant <- which(prT$en<date) else extant <- which(prT$en<=date)
        nExt             <- length(extant)
        fromExtantTips   <- vector("list", nExt)  
        for(n in 1:nExt)
          {
          whr <- which(pP$nm==prT$nm[extant[n]])
          if(length(whr)>0)
            {
            fromExtantTips[[n]] <- p2r[[whr]]
            }
          }
          
        unqFromExtant    <- unique(unlist(fromExtantTips))
        nn <- length(unqFromExtant)
        whrs <- numeric(nn)  
        for(n in 1:nn)
          {
          if(length(intersect(prT$nm, unqFromExtant[n]))>0)
            {
            whrs[n] <- which(prT$nm==unqFromExtant[n])
            }
          }  
        whrs <- whrs[whrs!=0]
        prT <- prT[whrs,]
        if(letDie) prT$en[prT$en<date] <- date else prT$en[prT$en<=date] <- date
        }

    }  

  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#*if pruning around a focal date
  if(!is.null(date) & is.null(focLin))
    {
    if(!keepTips)
      {
      if(!is.null(focLin))
        {
        focSt <- prT$st[which(pP$nm==focLin)]
        if(focSt<date)
           stop(paste("the prune2date (", date, ") is before the start date (", focSt, 
              ") of focLin (", as.character(focLin), ")", sep=""))
        }
      if(letSpeciate) prT <- prT[prT$st>=date,] else prT <- prT[prT$st>date,]
          #remove species that started too late
       
      if(pruneDead)
        {
        if(letDie) extant <- which(prT$en<date) else extant <- which(prT$en<=date)
        nExt             <- length(extant)
        fromExtantTips   <- vector("list", nExt)  
        for(n in 1:nExt)
          {
          whr <- which(pP$nm==prT$nm[extant[n]])
          if(length(whr)>0)
            {
            fromExtantTips[[n]] <- p2r[[whr]]
            }
          }
          
        unqFromExtant    <- unique(unlist(fromExtantTips))
        nn <- length(unqFromExtant)
        whrs <- numeric(nn)  
        for(n in 1:nn)
          {
          if(length(intersect(prT$nm, unqFromExtant[n]))>0)
            {
            whrs[n] <- which(prT$nm==unqFromExtant[n])
            }
          }  
        whrs <- whrs[whrs!=0]
        prT <- prT[whrs,]
        if(letDie) prT$en[prT$en<date] <- date else prT$en[prT$en<=date] <- date
        }
     if(!pruneDead) if(letDie) prT$en[prT$en<date] <- date else prT$en[prT$en<=date] <- date
      }
    }
    if(!is.null(focLin) & !is.null(date))
      {
      if(!keepTips)
        {stop("prune2date currently coded only for the 'root-end' of the tree when focal lineage specified.")}
      if(keepTips)
        {
        prT$st[which(prT$nm==focLin)] <- date
        keep <- prT$st<=date
        prT <- prT[keep,]
        }
      } 
    


  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#*output stuff 
  prTpP <- with(prT, as.paleoPhylo(nm, pn, st, en, label=label, grp=grp))
  if(collapseBranches) prTpP <- collapseBranches(prTpP)
  if(outPhylo) prTpa <- reorder(buildApe(createBifurcate(prTpP)))  else prTpa <- NULL
  
  if(!outPhylo) return(prTpP) else return(list(paleoPhylo=prTpP, phylo=prTpa))	
  }