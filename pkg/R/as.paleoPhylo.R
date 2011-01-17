as.paleoPhylo <- function (nm, pn, st, en, xx = NA, label = nm, grp=NA) 
  {
  if(length(nm)!=length(unique(nm))) stop (paste("ID codes are not unique, specifically",nm[duplicated(nm)]))
  if(!is.na(pn[st==max(st)])) 
    {
    warning(paste("The oldest species' ancestor has been changed to 'NA' from", pn[st==max(st)], "."))
    pn[st==max(st)] <- NA
    }

  
  dat<- data.frame(nm,pn,st,en,xx=NA,label,grp)
  dat <- dat[rev(order(dat$st, dat$en)),]
  pP <- list(nm = as.character(dat$nm), pn = as.character(dat$pn), 
  st = dat$st, en = dat$en, xx = dat$xx, label = as.character(dat$label), grp=dat$grp)

  class(pP) <- "paleoPhylo"
  
  xxloc <- getXloc(pP)[, c(1, 6)]
  dummy <- data.frame(nm = pP$nm)
  
  ifelse(is.na(xx),
    pP$xx <- merge(dummy, xxloc, by.x = "nm", by.y = "nm", sort = FALSE)[, 2],
    pP$xx <- (xx - min(xx))/max(xx - min(xx)))
  return(pP)
  }
    