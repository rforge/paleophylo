incStratUnc <- function (uSR = NULL, pP, lwdLin = 2, clCd=1:length(pP$nm), ltyLin=1)	{   	styles <- list(certain = c(lwdLin, "black", 1), lazarus = c(lwdLin, "grey60", 1), extension = c(lwdLin, "black", 2), CI = c(lwdLin * 0.4, "black", 1), point=c(lwdLin * 4, "black", 1))
   	dts <- vector("list", length(pP$nm))
  	for (i in 1:length(pP$st)) {dts[[i]] <- c(pP$st[i], pP$en[i]) }
    usr <- list(id = as.list(pP$nm), dates = dts, types = as.list(rep(1, length(pP$nm))), styles = styles)
  
	if((sum(is.na(pP$grp))==0) && is.null(uSR))
	    {
        for (i in 1:length(pP$nm)) {usr$types[[i]] <- as.numeric(as.factor(as.character(pP$grp)))[i]}
        
        if(length(lwdLin) < length(unique(pP$grp))) lwdLin <- rep(lwdLin, ceiling(length(unique(pP$grp))/length(lwdLin)))
        #if(length(clCd) < length(unique(pP$grp))) clCd <- rep(clCd, ceiling(length(unique(pP$grp))/length(clCd)))
        if(length(ltyLin) < length(unique(pP$grp))) ltyLin <- rep(ltyLin, ceiling(length(unique(pP$grp))/length(ltyLin)))
        
        usr <- list(id=usr$id,dates=usr$dates,types=usr$types, styles=NULL)
        for (i in 1:length(unique(pP$grp))) {usr$styles[[i]] <- c(lwdLin[i], clCd[i], ltyLin[i]) }
        }

    if (!is.null(uSR) & !is.data.frame(uSR)) 
    	{        reqTypes <- unique(unlist(uSR$types))        if (sum(sapply(1:length(reqTypes), function(i) {is.null(uSR$styles[[reqTypes[i]]]) })) != 0)  stop("Not all types have a defined style")        badd <- which(sapply(1:length(uSR$dates), function(i) {length(uSR$dates[[i]]) != (length(uSR$types[[i]]) + 1)        }))        if (length(badd) > 0) {cat(paste("Mismatch at lineage ", badd, sep = ""))}        cat("\n")        if (sum(sapply(1:length(uSR$dates), function(i) { length(uSR$dates[[i]]) != (length(uSR$types[[i]]) +  1)})) != 0)        	stop("The number of dates and types do not match. N(dates)!=(N(types)+1) for all lineages.")        usr <- uSR       	} 
 	if(is.data.frame(uSR)) 		{ 		for (k in 1:length(uSR$id))			{			flag <- FALSE			focID <- which(usr$id==uSR$id[k])			oldDates <- usr$dates[[focID]]			splitLoc <- which(oldDates<=uSR$st[k] & oldDates>=uSR$en[k])[1]			if (is.na(splitLoc)) {splitLoc <- which(oldDates<=uSR$st[k])  ;  flag<-TRUE}			newDates <- rev(sort(c(oldDates[-splitLoc],uSR$st[k],uSR$en[k])))			if (flag==TRUE) {newDates <- c(newDates, oldDates[length(oldDates)])}					oldTypes <- newTypes <- usr$types[[focID]]			if(splitLoc==1) {newTypes <- c(uSR$type[k],oldTypes)}			if(splitLoc==(length(newDates)-1)) {newTypes <- c(oldTypes,uSR$type[k])}			if(splitLoc!=(length(newDates)-1) & splitLoc!=1) 				{				if(uSR$st[k]!=uSR$en[k])  {newTypes <- c(oldTypes[1:(splitLoc-1)],uSR$type[k],oldTypes[(splitLoc+1):length(oldTypes)])}				if(uSR$st[k]==uSR$en[k])  {newTypes <- c(oldTypes[1:(splitLoc-1)],uSR$type[k],oldTypes[(splitLoc-1):(length(oldTypes)-1)])}				}							if((oldDates[1]==uSR$st[k] & oldDates[2]==uSR$en[k] & length(oldDates)==2)==FALSE) {usr$types[[focID]] <- newTypes  ;  usr$dates[[focID]] <- newDates}			if((oldDates[1]==uSR$st[k] & oldDates[2]==uSR$en[k] & length(oldDates)==2)==TRUE)  {usr$types[[focID]] <- uSR$type[k]}					}		}	return(usr)	}
