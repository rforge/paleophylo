getXloc <- function (pP) 
  {
  if (class(pP) != "paleoPhylo") stop("object is not of class 'paleoPhylo'")
    {
    dat <- data.frame(nm = pP$nm, pn = pP$pn, st = pP$st, en = pP$en)
    pos <- 0.5
    time <- minTime <- min(-pP$st)
    ids <- pP$nm[1]
    tmFrm <- round(-rev(sort(pP$st)), 5)
    while (time != max(tmFrm))
      {
      time <- tmFrm[(round(tmFrm, 8) > round(time, 8)) == TRUE][1]
      event <- dat[which(-round(dat$st, 5) == round(time, 5)), ]
      parents <- unique(as.character(event$pn))
      inIds <- sapply(1:length(parents), function(i) sum(intersect(ids, parents) == as.list(parents)[[i]]) > 0)
      parents <- c(parents[inIds == TRUE], parents[inIds == FALSE])

      for (i in 1:length(parents))
        {
        singleEvent <- event[as.character(event$pn) == parents[i], ]
        if (length(singleEvent[, 1]) > 1) evnt  <- 1
        if (length(singleEvent[, 1]) == 1) evnt <- 0

        focInd <- as.character(singleEvent$nm)
        pntLoc <- which(ids == parents[i])
        if(time!=minTime)
          {
          fcPrnt <- as.character(unique(singleEvent$pn))
          fcSist <- as.character(unique(singleEvent$nm))
          if(length(pntLoc) == 0) stop(paste("There is a lack of congruence in the tree around",focInd,"\nHave a look at the ancestor",fcPrnt,"but also the sister species",fcSist[1],"and",fcSist[2],"\n"))
          }
           
        parentPos <- pos[pntLoc]
        ids <- c(ids, focInd)
        sortPos <- sort(pos)

        if (length(pos) > 1)
          {
          if (evnt == 0) 
            {
            newPos <- (sortPos[which(as.numeric(sortPos) == as.numeric(parentPos)) - 1] + sortPos[which(as.numeric(sortPos) == as.numeric(parentPos))])/2
            if (length(newPos) == 0 & parentPos == max(pos)) newPos <- extendrange(pos)[2]
            if (length(newPos) == 0 & parentPos == min(pos)) newPos <- extendrange(pos)[1]
            }
          if (evnt == 1)
            {
            newPos1 <- (sortPos[which(as.numeric(sortPos) == as.numeric(parentPos)) - 1] + sortPos[which(as.numeric(sortPos) == as.numeric(parentPos))])/2
            if (parentPos == min(pos)) newPos1 <- extendrange(pos)[1]
            newPos2 <- (sortPos[which(as.numeric(sortPos) == as.numeric(parentPos)) + 1] + sortPos[which(as.numeric(sortPos) == as.numeric(parentPos))])/2
            if (parentPos == max(pos)) newPos2 <- extendrange(pos)[2]
            newPos <- c(newPos1, newPos2)
            if (length(focInd) > 2)
              {
              for (j in 3:length(focInd)) newPos <- c(newPos, (newPos[j - 2] + parentPos)/2)
              }
            }
          pos <- c(pos, newPos)
          }
        if (length(pos) == 1) ifelse(evnt == 0, pos <- seq(0, 1, 1), pos <- c(0.5,  0, 1))

        if (evnt == 1 & length(focInd) == 3 & length(sortPos) == 1) pos <- c(0.5, 0, 1, 0.25)
        if (evnt == 1 & length(focInd) == 4 & length(sortPos) == 1) pos <- c(0.5, 0, 1, 0.25, 0.75)
        if (evnt == 1 & length(focInd) > 4 & length(sortPos) == 1) stop("Cannot yet handle polytomies of order > 4 immediately post root.")
        if (length(pos) != length(ids)) stop(paste("There is not a position for all individuals.\n The problem is something to with", focInd))
      }
    }
  locDat <- data.frame(ids, pos)
  locDat$xx <- as.numeric(as.factor(locDat$pos))/length(locDat$pos)
  if (length(locDat$xx) != length(unique(locDat$xx))) stop("non-unique locations")
  cmbDat <- merge(dat, locDat, by.x = c("nm"), by.y = c("ids"))
  }
  return(cmbDat)
}
