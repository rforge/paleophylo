 as.paleoPhylo <- function (nm,pn,st,en,xx=NA,label=nm) 
	{
	pP <- list(nm=as.character(nm),pn=as.character(pn),st=st,en=en,xx=xx,label=as.character(label))
	class(pP) <- "paleoPhylo"
	ifelse (is.na(xx), pP$xx<-merge(data.frame(nm=pP$nm),getXloc(pP)[,c(1,6)],by.x="nm",by.y="nm",sort=FALSE)[,2],  pP$xx<-(xx - min(xx))/max(xx - min(xx)))
	return(pP)
	}