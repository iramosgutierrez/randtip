phyleticity<- function(tree, genus){

    if(length(genus) != 1){stop("Only one genus accepted.") }

    sp <- tree$tip.label
    taxa.vector <- sp[first.word(sp)==genus]

    if(length(taxa.vector) == 0){
        genus.type <- "Missing"
        return(genus.type)
    }

    if(length(taxa.vector) == 1){
        genus.type <- "Singleton"
        return(genus.type)
    }

    mrca <- phytools::findMRCA(tree = tree, tips = taxa.vector)
    descend <- phytools::getDescendants(tree, mrca)
    desc.tips <- notNA(sp[descend])
    desc.genera <- first.word(desc.tips)

    if(length(unique(desc.genera)) == 1){
        genus.type <- "Monophyletic"
        return(genus.type)
    }
    intruder.tips <- desc.tips[desc.genera != genus]
    if(length(intruder.tips) == 1){
        # there is only one intruder: paraphyletic by singleton)
        genus.type <- "Paraphyletic"
    }else{
        intruder.mrca <- phytools::findMRCA(tree = tree, tips = intruder.tips)
        # if there are more than one intruder...
        desc.intruder.mrca <- phytools::getDescendants(tree, intruder.mrca)
        desc.intruder.tips <- sp[desc.intruder.mrca]
        desc.intruder.tips <- notNA(desc.intruder.tips)
        #intruders grouped: paraphyletic group by monophyletic intruder
        suppressWarnings({
            grouped.intruders <- all(sort(desc.intruder.tips) == sort(intruder.tips))
        })
        if(grouped.intruders){
            genus.type <- "Paraphyletic"
        }else{
            genus.type <- "Polyphyletic"
        }
    }

    return(genus.type)
}

MDCC.phyleticity<-function(input, tree, MDCC.info=list("rank"=NA, "MDCC"=NA),
                           trim=TRUE){

    rank<- MDCC.info$rank
    MDCC <- MDCC.info$MDCC
    input<-input[!is.na(input[,rank]),]

    if(rank=="genus"){
        MDCC.type<-phyleticity(tree, MDCC)
        return(MDCC.type)
    }

    tips<- tree$tip.label[first.word(tree$tip.label) %in% first.word(input$taxon)]
    if(isTRUE(trim) & length(tips)>0){
        tree<- ape::keep.tip(phy = tree, tip = tips)
    }

    species<- input[input[,rank]==MDCC,]
    MDCC.genera<- unique(first.word(species$taxon))
    genera.in.tree<-first.word(tree$tip.label)
    genera.in.tree<-MDCC.genera[MDCC.genera %in% genera.in.tree]
    spp.in.tree<- tree$tip.label[first.word(tree$tip.label) %in% genera.in.tree]

    if(length(spp.in.tree)==0){
        MDCC.type<-"Missing"
        return(MDCC.type)
    }

    if(length(spp.in.tree)==1){
        MDCC.type<-"Singleton"
        return(MDCC.type)
    }

    sp.mrca<- phytools::findMRCA(tree, tips = spp.in.tree)
    descs.num<- phytools::getDescendants(tree,sp.mrca)
    descs.name<-notNA(tree$tip.label[descs.num])

    if(all(first.word(descs.name)%in%MDCC.genera)){
        MDCC.type<-"Monophyletic"
        return(MDCC.type)
    }else{
        intruders<-descs.name[!(first.word(descs.name)%in%MDCC.genera)]
        if(length(intruders)==1){
            MDCC.type<-"Paraphyletic"
            return(MDCC.type)
        }

        intruders.mrca<- ape::getMRCA(phy = tree, tip = intruders)
        intruders.descs.num<- phytools::getDescendants(tree,intruders.mrca)
        intruders.descs.name<-notNA(tree$tip.label[intruders.descs.num])

        if(!any(first.word(intruders.descs.name)%in%MDCC.genera)){
            MDCC.type<-"Paraphyletic"
            return(MDCC.type)
        }else{
            MDCC.type<-"Polyphyletic"
            return(MDCC.type)
        }
    }
}
