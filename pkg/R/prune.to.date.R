prune.to.date <-
function(pP, date, let.speciate=FALSE, let.die=FALSE){
	#Is passed a paleoPhylo object and a cutoff date, and returns the phylogeny of just the species extant at the cutoff.
	#The phylogenies is returned in two formats: paleoPhylo and phylo.
	#If let.speciate = TRUE, speciations at exactly the cutoff date are taken to have happened.
	#If let.die = TRUE, extinctions at exactly the cutoff date are taken to have happened.

	#check it's a paleoPhylo object
	if (class(pP) !="paleoPhylo")
		stop("Object passed to prune.to.date is not of class paleoPhylo!")
	
	df <- with(pP, data.frame(nm, pn, st, en, xx, label, grp))
	
	#Get rid of lineages that started too late
	if (let.speciate){
		pruned<-subset(df, st >= date)
	} else {
		pruned <- subset(df, st > date)
	}

	if (length(pruned$nm) == 0){
		warning(paste("No lineages at time", date))
		return(NULL)
	}
	
	if (length(pruned$nm) == 1){
		warning(paste("Only one lineage at time", date))
		return(NULL)
	}

	lis <- with(pruned, as.paleoPhylo(nm, pn, st, en, xx, label, grp))	
	cb <-createBifurcate(lis) #Create structures that will let extinct lineages be identified unambiguously
	
	cb$en[cb$en < date] <- 0 #Pad out all terminal edges to the present day
	pruned.ape <-buildApe(cb) #Produce corresponding phylo object
	write.tree(pruned.ape, "tmp.txt") #Write and re-read tree to make it a canonical phylo object
	pruned.ape <- read.tree("tmp.txt")
	
	#Identify lineages extinct before or at cutoff, and prune them (if there are any) from the tree
	if (let.die){
		extinct.by.cutoff <- c(cb$nm[cb$en > date], cb$nm[cb$en == date & !is.element(cb$nm, df$pn)])
	} else {
		extinct.by.cutoff <- cb$nm[cb$en > date]
	}
	
	pruned.extant.ape <- pruned.ape #Initiate phylo object for extant taxa, prior to pruning out those that have gone extinct
	to.drop <- (intersect(pruned.ape$tip.label, extinct.by.cutoff))
	if (length(to.drop)>0){
		pruned.extant.ape <- drop.tip(pruned.extant.ape,to.drop)
	}

	#Now go back to paleoPhylo to trim terminal edges to make tree ultrametric and ending at cutoff
	app <- ape2paleoPhylo(pruned.extant.ape, retainNodeLabels=TRUE, nC=0)
	app$en[app$en < date] <- date #Make tree ultrametric
	pruned.extant.ape2 <-buildApe(app)
	write.tree(pruned.extant.ape2, "tmp.txt") #Write and re-read tree to make it a canonical phylo object
	pruned.extant.ape2 <- read.tree("tmp.txt")
		
	#Set up and fill return structure
	to.return <- list(app, pruned.extant.ape2)
	names(to.return)<- c("paleoPhylo.tree","phylo.tree")
	
	return(to.return)
}

