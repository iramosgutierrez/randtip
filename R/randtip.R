add.no.exceptions <- function(genus, new.tree, species.table, forbidden.groups,
                              genus.taxa.no.except, taxa.lvl){
    
    taxa.lvls <- c("genus", "family", "order", "class")
    taxa.no.except <- species.table$using.taxa %in% genus.taxa.no.except
    taxa.lvl.use <- unique(species.table[[taxa.lvl]][taxa.no.except])
    genera <- species.table$genus[species.table[[taxa.lvl]] == taxa.lvl.use]
    genus.lab <- stringr::word(new.tree$tip.label,1, sep = "_")
    taxa <- new.tree$tip.label[genus.lab %in% genera]

    if(length(taxa) > 0){
        node <- phytools::findMRCA(tree = new.tree, tips = taxa)
        new.tree<- add.into.node(tree = new.tree, node = node, 
                                 new.tip = genus.taxa.no.except[1], 
                                 exception.list = forbidden.groups)
        if(length(genus.taxa.no.except) > 1){
            new.tips <- genus.taxa.no.except[2:length(genus.taxa.no.except)]
            new.tree <- add.to.singleton(tree = new.tree, 
                                         singleton = genus.taxa.no.except[1],
                                         new.tips = new.tips) 
        }
        
        return(new.tree)
        
    }else{
        if(taxa.lvl == "class"){
            stop(paste0("Genus ", genus, 
                        " has no family, order or class relatives"))
        }
        
        taxa.lvls.next <- taxa.lvls[which(taxa.lvls == taxa.lvl) + 1]        
        return(add.no.exceptions(genus, new.tree, species.table, 
                                 forbidden.groups, genus.taxa.no.except, 
                                 taxa.lvls.next))
    }

}

get.original.names <- function(new.tree, species.table, verbose = FALSE){
    
    for(n in 1:nrow(species.table)){ 
        tip.label <- new.tree$tip.label
        if(species.table$using.taxa[n] == species.table$taxon[n]) next
        if(species.table$using.taxa[n] %in% new.tree$tip.label){
            correct.tip.loc <- new.tree$tip.label == species.table$using.taxa[n]
            new.tree$tip.label[correct.tip.loc] <- species.table$taxon[n]
            if(verbose){
                print(paste0("name correction: ", n, "/", 
                             nrow(species.table), " (", 
                             round(n/nrow(species.table)*100, 2), " %)"))    
            }
        }
    }
    
    return(new.tree)
}

randtip.subsp <- function(new.tree, species.table.dupl, verbose = FALSE){
    
    rep.taxa <- species.table.dupl$taxon 
    rep.taxa.species <- unique(species.table.dupl$using.taxa)

    for(i in 1:length(rep.taxa.species)){
        sp.start<- Sys.time()
        # Select subspecies
        genus <- stringr::word(rep.taxa, 1, sep = "_")
        sp <- stringr::word(rep.taxa, 2, sep= "_")
        ssps <- rep.taxa[paste0(genus, "_", sp) == rep.taxa.species[i]] 
        ssps <- sample(ssps, length(ssps)) 
        genus.tree <- stringr::word(new.tree$tip.label, 1, sep = "_")
        sp.tree <- stringr::word(new.tree$tip.label, 2, sep= "_")
        sp.to.add <- new.tree$tip.label[paste0(genus.tree, "_", sp.tree) == rep.taxa.species[i]]
        if(length(sp.to.add) > 1){
            sp.to.add <- sp.to.add[sp.to.add == rep.taxa.species[i]]
        }
        # Add subspecies as singleton
        new.tree <- add.to.singleton(new.tree, singleton = sp.to.add, 
                                     new.tips = ssps) 
        sp.end<- Sys.time()
#        if(verbose){
#            print(paste0(i, "/", length(rep.taxa.species), " (", 
#                         round(i/length(rep.taxa.species)*100, 2), " %): ",  
#                         rep.taxa.species[i]," (",
#                         round(as.numeric(difftime(sp.end, sp.start, units = "secs")), 2), " sec. out of ",
#                         round(as.numeric(difftime(sp.end, start,     units = "mins")), 2), " mins)"))    
#        }
    }
    return(new.tree)
}

get.forbidden.groups <- function(tree, species.table){
    #forbidden.groups
    genus.types <- species.table$genus.type
    mono_or_para <- (genus.types == "Monophyletic") | 
                        (genus.types =="Paraphyletic")
    forb.genera <- species.table[mono_or_para, "genus"]
    forb.genera <- forb.genera[!duplicated(forb.genera)]
    forbidden.groups<- rep(list(NA),length(forb.genera))
    treelist <- data.frame("taxon" = tree$tip.label, 
                           "genus" = stringr::word(tree$tip.label, 1, 
                                                  sep = "_"))

    if(length(forbidden.groups) > 0){
        for(i in 1:length(forb.genera)){
            forbidden.groups[[i]] <- treelist[treelist$genus == forb.genera[i], 
                                              "taxon"]
        }
    }

    if(any(species.table$genus.type=="Polyphyletic")){
        poly.genera <- species.table[species.table$genus.type == "Polyphyletic",]
        poly.genera <- unique(poly.genera$genus)

        tree.taxa <- tree$tip.label
        genera <- stringr::word(tree.taxa, 1, sep = "_")
        for(genus in poly.genera){            
            genus.tree.taxa <- tree.taxa[genera == genus]
            genus.tree.list<- rep(list(NA), times = length(genus.tree.taxa))

            for(i in 1:length(genus.tree.taxa)){
                sp <- genus.tree.taxa[i]
                if(sp %in% unlist(genus.tree.list)) next
                
                sp.tip<- which(tree$tip.label==sp) 
                par.sib <- get.parent.siblings(tree, sp.tip)

                siblings <- par.sib$siblings
                siblings.genera <- stringr::word(siblings, 1, sep = "_")
                if(sum(siblings.genera != genus) == 1){
                    siblings <- siblings[siblings.genera == genus]
                }else if(sum(siblings.genera != genus) > 1){
                    intruders <- siblings[siblings.genera != genus]
                    intr.mrca <- ape::getMRCA(tree, intruders)
                    intr.desc.tips <- phytools::getDescendants(tree, intr.mrca)
                    intr.desc.names <- tree$tip.label[intr.desc.tips]
                    intr.desc.names <- intr.desc.names[!is.na(intr.desc.names)]

                    intr.genera <- stringr::word(intr.desc.names, 1, sep = "_")
                    if(all(intr.genera != genus)){
                        siblings <- siblings[siblings.genera == genus]
                    }
                }

                siblings.genera <- stringr::word(siblings, 1, sep = "_")
                while(length(unique(siblings.genera)) == 1){ 
                    #tip and parent upstream until they are from different genera
                    sp.tip <- parent
                    parent <- tree$edge[tree$edge[,2]==sp.tip,1]
                    parent.desc <- phytools::getDescendants(tree, parent)
                    siblings <- tree$tip.label[parent.desc]
                    siblings <- siblings[!is.na(siblings)]

                    siblings.genus <- stringr::word(siblings, 1, sep = "_")
                    if(length(siblings[siblings.genus != genus]) == 1){
                        siblings <- siblings[siblings.genus == genus]
                    }else if(length(siblings[siblings.genus != genus]) > 1){
                        intruders <- siblings[siblings.genus != genus]
                        intr.mrca <- ape::getMRCA(tree, intruders)
                        intr.desc.tips <- phytools::getDescendants(tree, intr.mrca)
                        intr.desc.names <- tree$tip.label[intr.desc.tips]
                        intr.desc.names <- intr.desc.names[!is.na(intr.desc.names)]
                        intr.desc.genus <- stringr::word(intr.desc.names, 1, sep = "_")
                        if(all(intr.desc.genus != genus)){
                            
                            siblings <- siblings[siblings.genus == genus] 
                        }
                    }
                }

                descs<- phytools::getDescendants(tree, sp.tip)
                desctips<-tree$tip.label[descs]
                desctips<-desctips[!is.na(desctips)]
                desctips<-desctips[stringr::word(desctips, 1, sep="_") == genus]
                genus.tree.list[[i]]<-desctips
            }

            genus.tree.list<-genus.tree.list[!is.na(genus.tree.list[])]
            genus.tree.list<-genus.tree.list[lengths(genus.tree.list[])>1]
            forbidden.groups<-c(forbidden.groups, genus.tree.list)
        }


    }

    return(forbidden.groups)
}

add.polytomy.types <- function(tree, taxa.table, species.table, type_arg, 
                               insertion, taxa.genera = NULL){
    taxa.lvls <- c("genus", "family", "order", "class")
    taxa.lvl <- stringr::str_extract(type_arg, "\\w+")
    
    if(taxa.lvl == "genus"){
        taxa.lvl.col <- taxa.genera
    }else{
        taxa.lvl.col <- taxa.table[[taxa.lvl]]
    }
    taxa <- unique(taxa.lvl.col)

    for(p in 1:length(taxa)){
        taxon <- taxa[p]
        tree <- unfold.polytomy.types(tree, taxa.table, species.table,
                                      taxa.lvl, insertion, 
                                      taxa.subset = taxon)

    }
    
    return(tree)
}

unfold.polytomy.types <- function(tree, taxa.table, species.table, 
                                  taxa.lvl, insertion, taxa.subset){
    
    taxa.lvls <- c("genus", "family", "order", "class")
       
    taxa.tbl.match <- taxa.table[[taxa.lvl]] == taxa.subset
    #this is redundant, but keeps the structure
    taxa <- taxa.table$taxon[taxa.tbl.match]
    taxa.genera <- species.table$genus[species.table[[taxa.lvl]] == taxa.subset] 
    #print(paste0(p, "/", length(genera), " (",round(p/length(genera)*100, 2), " %). ",
    #           "Adding Gen. ", genus, " (",
    #           length(genus.taxa)," tips). "))
    #species (tips) within class IN ORIGINAL TREE
    tip.genera <- stringr::word(tree$tip.label, 1, sep = "_")
    union.tips <- tree$tip.label[tip.genera %in% taxa.genera]   

    if(length(taxa) == 0) return(tree)

    if(length(union.tips) == 1){
      node <- which(tree$tip.label == union.tips)
      tree <- polytomy.over.node(tree = tree, species = taxa, 
                                 node = node, insertion = insertion)
      return(tree)
    }
    if(length(union.tips)> 1){
      node <- phytools::findMRCA(tree, tips = union.tips)
      tree <- polytomy.over.node(tree = tree, species = taxa, 
                                 node = node, 
                                 insertion = insertion)
      return(tree)
    }
    if(length(union.tips) == 0){
        if(taxa.lvl == "class"){
            ####JOIN TO UPPER TAXONOMIC LEVEL
            message(paste0("ATTENTION: genus ", taxa.genera, 
                           " was not included as no Class coincidences were found"))
            return(tree)
        }
        taxa.lvl.next <- taxa.lvls[which(taxa.lvls == taxa.lvl) + 1]
        taxa.subset <- taxa.table[[taxa.lvl.next]][taxa.tbl.match]
        taxa.subset <- taxa.subset[!duplicated(taxa.subset)]
        
        return(unfold.polytomy.types(tree, taxa.table, species.table, 
                                     taxa.lvl.next, insertion, taxa.subset))
    }
}

get.taxa.to.use <- function(species.table, aggregate.subspecies){
    
    species.table$using.taxa <- NA  
    for(i in 1:nrow(species.table)){          
        taxon <- species.table$taxon[i]
        genus.name <- stringr::word(species.table$taxon[i],1, sep = "_")
        sp.name <- stringr::word(species.table$taxon[i],2, sep = "_")
        if(is.na(species.table$aggregate.subspecies[i])){
            # For NAs there are no exceptions
            if(aggregate.subspecies == TRUE){
                species.table$using.taxa[i] <- paste0(genus.name, "_", 
                                                      sp.name)                                                        
            }else{
                species.table$using.taxa[i] <- taxon
            }
        # Identified with 0 are not clustered. 1 are clustered
        }else if(species.table$aggregate.subspecies[i] == 0){
            species.table$using.taxa[i] <- taxon
        }else if(species.table$aggregate.subspecies[i] == 1){
            species.table$using.taxa[i] <- paste0(genus.name, "_", 
                                                  sp.name)
        }

    }   

    return(species.table$using.taxa)
    
}

#' randtip "MOTHER FUNCTION"
#' @export
randtip <- function(tree, species.table, 
                    type = c("random", "genus.polytomy", "family.polytomy", 
                             "order.polytomy", "class.polytomy"),
                    aggregate.subspecies = TRUE, 
                    insertion = c("random", "middle","long"), 
                    prob = TRUE, verbose = FALSE){
    
    start<- Sys.time()
    
    new.tree <- tree
    species.table.dupl <- NULL
    if(type=="random"){        
        species.table$using.taxa <- get.taxa.to.use(species.table, 
                                                    aggregate.subspecies)
        species.table <- species.table[order(species.table$taxon),]
        is.duplicated <- duplicated(species.table$using.taxa)
        species.table.dupl <- species.table[is.duplicated,]
        species.table <- species.table[!is.duplicated,]

        taxa <- species.table$using.taxa
        taxa <- gsub(" ", "_", taxa)
        taxa <- taxa[!(taxa %in% tree$tip.label)]
        taxa.genera <- stringr::word(taxa, 1, sep = "_")
        taxa.genera <- unique(taxa.genera)
        
        forbidden.groups <- get.forbidden.groups(new.tree, species.table)
        
        for(i in 1: length(taxa.genera)){        #loop 1
            gen.start<- Sys.time()
            genus <- taxa.genera[i]
            genus.match <- species.table$genus==genus
            genus.type <- species.table$genus.type[genus.match]
            genus.type <- unique(genus.type)

            genus.taxa <- taxa[stringr::word(taxa, 1, sep = "_") == genus] 
            genus.taxa <- sample(genus.taxa, length(genus.taxa)) 
            # Tips with no "relatives" information selected (they will be bound as monophyletic)
            grouped.taxa <- species.table$using.taxa[!is.na(species.table$relative.species)]
            grouped.taxa <- grouped.taxa[grouped.taxa %in% genus.taxa]
            genus.taxa <- genus.taxa[!(genus.taxa %in% grouped.taxa)]

            if(verbose){
                print(paste0(i, "/", length(taxa.genera), 
                             " (",round(i/length(taxa.genera)*100, 2), " %). ",
                      "Adding Gen. ", genus, " (", genus.type,", ",
                       length(genus.taxa)," tips). "))
            }
            
            if(length(grouped.taxa) > 0){   
                for(j in 1:length(grouped.taxa)){
                    grouped.taxa.loc <- (species.table$using.taxa == grouped.taxa[j])
                    grouping.taxa <- species.table$relative.species[grouped.taxa.loc]
                    grouping.taxa <- gsub(" ","", grouping.taxa)
                    grouping.taxa <- stringr::str_split(grouping.taxa, 
                                                        split = ",")[[1]]
                    # Select taxa to be grouped
                    grouping.taxa <- grouping.taxa[grouping.taxa %in% tree$tip.label] 
                    # First added as singleton. The rest added as monophyletic
                    if(length(grouping.taxa)==1){
                        new.tree <- add.to.singleton(new.tree, singleton = grouping.taxa, 
                                                    new.tips = grouped.taxa[j])  
                    }else{
                        node <- phytools::findMRCA(new.tree, tips = grouping.taxa)
                        new.tree <- add.into.node(new.tree, new.tip = grouped.taxa[j],
                                                  node = node)} 
                }
                next
            }

            if(genus.type=="Monophyletic"){
                for( j in 1:length(genus.taxa)){
                    new.tree <- add.to.monophyletic(new.tree, 
                                                    new.tip = genus.taxa[j])
                }
            }else if(genus.type=="Paraphyletic"){
                for( j in 1:length(genus.taxa)){
                    new.tree <- add.to.paraphyletic(new.tree, 
                                                    new.tip = genus.taxa[j])
                }
            }else if(genus.type=="Polyphyletic"){
                new.tree<- add.to.polyphyletic(new.tree, species = genus.taxa)
            }else if(genus.type=="Singleton genus"){   
                # All tips added in one step
                genus.match <- (stringr::word(tree$tip.label, 1, sep = "_") == genus)
                new.tree <- add.to.singleton(new.tree, 
                                             singleton = tree$tip.label[genus.match], 
                                             new.tips = genus.taxa)
            }else if(genus.type == "Not included"){
                genus.match <- (species.table$genus == genus)

                # Genera with synonym exceptions
                has.synonyms <- !is.na(species.table$synonim.genus)
                genus.synonyms <- species.table$synonim.genus[genus.match & has.synonyms]
                genus.synonyms.tips <- species.table$using.taxa[genus.match & has.synonyms]
                if(length(genus.synonyms) > 0){
                    new.genus.synonyms.tips <- genus.synonyms.tips
                    for(i in 1:length(genus.synonyms.tips)){
                        new.genus.synonyms.tips[i] <- paste0(genus.synonyms[i], "_", 
                                                         genus.synonyms.tips[i]) 
                        new.genus<- stringr::word(new.genus.synonyms.tips[i], 1, sep = "_")
                        new.genus.match <- species.table$genus == new.genus
                        new.genus.type <- species.table$genus.type[new.genus.match]
                        new.genus.type <- unique(new.genus.type)

                        if(new.genus.type=="Monophyletic"){
                            new.tree <- add.to.monophyletic(new.tree, new.tip = new.genus.synonyms.tips[i])}
                        if(new.genus.type=="Paraphyletic"){
                            new.tree<- add.to.paraphyletic(new.tree, new.tip = new.genus.synonyms.tips[i])}
                        if(new.genus.type=="Polyphyletic"){
                            new.tree <- add.to.polyphyletic(new.tree, species = new.genus.synonyms.tips[i])}
                        if(new.genus.type=="Singleton genus"){
                            genus.tree <- stringr::word(new.tree$tip.label,1,sep= "_")                    
                            new.tree <- add.to.singleton(new.tree,
                                                         singleton = new.tree$tip.label[genus.tree == new.genus],
                                                         new.tips = new.genus.synonyms.tips[i])
                        }

                        new.tree$tip.label <- gsub(pattern = new.genus.synonyms.tips[i],
                                                   replacement = genus.synonyms.tips[i],
                                                   x = new.tree$tip.label)
                    }

                }
                # Genera with siblings exceptions                
                has.siblings <- !is.na(species.table$sibling.genus)
                genus.siblings <- species.table$sibling.genus[genus.match & has.siblings]
                genus.siblings.tips <- species.table$using.taxa[genus.match & has.siblings]
                if(length(genus.siblings) > 0){
                    new.genus <- unique(genus.siblings)
                    if(length(new.genus) >1 ){
                        stop(paste0("Multiple genus siblings for genus ", genus))
                    }
                    genus.siblings.tips <- sample(genus.siblings.tips, length(genus.siblings.tips))
                    genus <- stringr::word(new.tree$tip.label, 1, sep = "_")            
                    adding.node <- phytools::findMRCA(new.tree, 
                                                     tips = new.tree$tip.label[genus == new.genus])
                    new.tree <- add.over.node(new.tree, new.tip = genus.siblings.tips[1], 
                                              node = adding.node)  
                    new.tree <- add.to.singleton(new.tree, singleton = genus.siblings.tips[1],
                                                 new.tips =  genus.siblings.tips[2:length(genus.siblings.tips)])
                }

                # Without exceptions are added randomly within family, order or class
                genus.taxa.no.except<- genus.taxa[!(genus.taxa %in% c(genus.synonyms.tips,
                                                                      genus.siblings.tips))]
                if(length(genus.taxa.no.except) > 0){
                    new.tree <- add.no.exceptions(genus, new.tree, species.table, 
                                                  forbidden.groups,
                                                  genus.taxa.no.except, 
                                                  taxa.lvl = "genus")
                }
            }
            gen.end <- Sys.time()
            if(verbose){
                print(paste0("\U2713", "( done in ",
                           round(as.numeric(difftime(gen.end, gen.start, units = "secs")), 2), " sec. out of ",
                           round(as.numeric(difftime(gen.end, start,     units = "mins")), 2), " mins)"))    
            }
            
        }

        # Phase 2 - subspecies
        if(nrow(species.table.dupl)>0){ 
            new.tree <- randtip.subsp(new.tree, species.table.dupl, verbose)
        }

        new.tree <- get.original.names(new.tree, species.table)        
        
    }else{
        #In polytomy cases, names are not changed
        species.table$using.taxa <- species.table$taxon 

        taxa.table <- species.table[!duplicated(species.table$using.taxa),]
        taxa.table$using.taxa <- gsub(" ", "_", taxa.table$using.taxa)
        taxa.table <- taxa.table[!(taxa.table$using.taxa %in% tree$tip.label),]

        taxa.genera <- taxa.table$genus
        taxa.genera <- sample(taxa.genera, length(taxa.genera))
        
        new.tree <- add.polytomy.types(new.tree, taxa.table, species.table, 
                                       type, insertion, taxa.genera)

    }    

    if(!is.null(species.table.dupl)){
        complete.taxa.list <- c(species.table$taxon, species.table.dupl$taxon)
    }else{
        complete.taxa.list <- species.table$taxon
    }
    
    complete.taxa.list.in.tree <- complete.taxa.list[complete.taxa.list %in% new.tree$tip.label]
    not.included <- complete.taxa.list[!(complete.taxa.list %in% complete.taxa.list.in.tree)]
    if(length(not.included) > 1){
        message(paste0("The following taxa were not included in the tree: ", not.included))
    }

    new.tree <- ape::keep.tip(new.tree, complete.taxa.list.in.tree)
    
    return(new.tree)
}

