`buildApe` <- function(pP)
	{	
	if(class(pP)!="paleoPhylo") stop("object is not of class 'paleoPhylo'")
		{
		get.line<-function (pn,dat,stNode,enNode=stNode)
			{
			tmp <- data.frame(nm=pP$nm[-which(is.na(pP$pn))],pn=pP$pn[-which(is.na(pP$pn))],st=pP$st[-which(is.na(pP$pn))],en=pP$en[-which(is.na(pP$pn))],Speciation=pP$Speciation[-which(is.na(pP$pn))],cd=pP$cd[-which(is.na(pP$pn))])
			tmp<-tmp[as.character(tmp$pn)==pn,]
	
			if(length(tmp[,1])>0) {
					vv<-data.frame(IDcode=tmp$pn, IDname=rep(dat$nm[as.character(dat$nm)==as.character(tmp$pn)[1]],2),descName=tmp$nm,descCode=tmp$cd,duration=tmp$st-tmp$en,
					tip=0+(tmp$Speciation==0),extant= 0 + (tmp$en==min(dat$en)),startNode=rep(stNode,2),endNode=enNode+1:2,startTime=tmp$st,endTime=tmp$en)}
			if(length(tmp[,1])==1) {vv <- vv[1,]}
			vv<-vv[order(rev(vv$startTime)),]
			return(vv)
			}

		pP$Speciation <- as.numeric(sapply(1:length(pP$pn),function(i) sum(match(pP$pn,pP$nm[i],nomatch=0))>0))
		pP$cd <- 1:length(pP$nm)
		stop <- noEvent <- FALSE
		pn <- as.character(pP$nm[pP$st==max(pP$st)])
		clade <- get.line(pn,pP,1000)

		while(stop==FALSE)
			{
			oldLength <- length(clade[,1])
			for(i in 1:length(clade[,1]))
				{
				if (sum(na.omit(as.character(clade$descName)[i])==as.character(clade$IDcode))==0 & clade$tip[i]==0)
					{
					parent <- clade$descName[i]
					nxt <- get.line(parent,pP,clade$endNode[i],max(clade$endNode))
					clade <- rbind(clade,nxt)
					}
				}
			if(oldLength==length(clade[,1])) {stop<-TRUE}
			}
	
		clade$endNode[clade$tip==1]<-clade$descCode[clade$tip==1]
		edge<-matrix(as.numeric(as.factor(c(clade$startNode,clade$endNode))),ncol=2)
		clade$nodeLabel<-edge[,1]

		if (sum(with(clade,tapply(nodeLabel,nodeLabel,length))!=2) !=0) {warning("bifuracting walk through has failed to yield a completely bifurcating tree")}

		tipDetails<-data.frame(cd=as.numeric(as.factor(clade$endNode[clade$tip==1])),nm=as.character(clade$descName)[clade$tip==1])
	    nodeDetails <- data.frame(cd = as.numeric(as.factor(clade$endNode[clade$tip == 0])), nm = as.character(clade$descName)[clade$tip == 0])
	    
		mft <- list(edge = edge, tip.label = as.character(tipDetails[order(tipDetails$cd),2]), 
		               Nnode = length(unique(clade$IDcode)),edge.length=clade$duration, node.label = as.character(nodeDetails[order(nodeDetails $cd),2]))
		class(mft) <- "phylo"	
		return(mft)
		}
	}

