get.taxa.to.use <- function(DF1, aggregate.subspecies=TRUE){

  if(isTRUE(aggregate.subspecies)){DF1[is.na(DF1$aggregate.subspecies),"aggregate.subspecies"]<-"1"}else{
                                   DF1[is.na(DF1$aggregate.subspecies),"aggregate.subspecies"]<-"0"}

  DF1$using.taxa <- NA

    for(i in 1:nrow(DF1)){
        taxon <- DF1$taxon[i]
        genus.name <- stringr::word(DF1$taxon[i],1, sep = "_")
        sp.name <-    stringr::word(DF1$taxon[i],2, sep = "_")
        if(DF1$aggregate.subspecies[i]=="1"){
                DF1$using.taxa[i] <- paste0(genus.name, "_", sp.name) }

        if(DF1$aggregate.subspecies[i]=="0"){
                DF1$using.taxa[i] <- taxon}

        # Identified with 0 are not clustered. 1 are clustered
   }

    return(DF1$using.taxa)

}

get.forbidden.groups <- function(tree, DF1){
  #forbidden.groups
  genus.types <- DF1$phyleticity
  mono_or_para <- (genus.types == "Monophyletic") |
    (genus.types =="Paraphyletic")

  forb.genera <- DF1[mono_or_para, "genus"]
  forb.genera <- unique(forb.genera)
  forbidden.groups<- rep(list(NA),length(forb.genera))
  treelist <- data.frame("taxon" = tree$tip.label,
                         "genus" = randtip::firstword(tree$tip.label))

  if(length(forbidden.groups) > 0){
    for(i in 1:length(forb.genera)){
      forbidden.groups[[i]] <- treelist[treelist$genus == forb.genera[i],
                                        "taxon"]
    }
  }

  if(any(DF1$phyleticity=="Polyphyletic")){
    poly.genera <- DF1[DF1$phyleticity == "Polyphyletic",]
    poly.genera <- unique(poly.genera$genus)

    tree.taxa <- tree$tip.label
    genera <- randtip::firstword(tree.taxa)
    for(genus in poly.genera){
      genus.tree.taxa <- tree.taxa[genera == genus]
      genus.tree.mrca<- phytools::findMRCA(tree=tree, tips = genus.tree.taxa)
      genus.tree.mrca.descs<- phytools::getDescendants(tree=tree, node=genus.tree.mrca)

      genus.tree.list<- rep(list(NA), times = length(genus.tree.mrca.descs))

      for(i in 1:length(genus.tree.mrca.descs)){
        node<- genus.tree.mrca.descs[i]
        if(randtip::is.tip(tree, node)){next}
        node.descs<- phytools::getDescendants(tree=tree, node=node)
        tip.descs <- randtip::notNA(tree$tip.label[node.descs])
        tip.descs.gen<- unique(randtip::firstword(tip.descs))
        if(length(tip.descs.gen)==1){
          genus.tree.list[[i]]<-tip.descs
        }


        siblings.genera <- randtip::firstword(siblings)
        while(length(unique(siblings.genera)) == 1){
          if(exists("p")){p<-p+1}else{p=1}#for comprobation, remove afterwards
          print(p)#for comprobation, remove afterwards
          #tip and parent upstream until they are from different genera
          if(exists("parent")){sp.tip <- parent}
          parent <- tree$edge[tree$edge[,2]==sp.tip,1]
          parent.desc <- phytools::getDescendants(tree, parent)
          siblings <- tree$tip.label[parent.desc]
          siblings <- randtip::notNA(siblings)

          siblings.genus <- randtip::firstword(siblings)
          if(length(siblings[siblings.genus != genus]) == 1){
            siblings <- siblings[siblings.genus == genus]
          }else if(length(siblings[siblings.genus != genus]) > 1){
            intruders <- siblings[siblings.genus != genus]
            intr.mrca <- ape::getMRCA(tree, intruders)
            intr.desc.tips <- phytools::getDescendants(tree, intr.mrca)
            intr.desc.names <- tree$tip.label[intr.desc.tips]
            intr.desc.names <- randtip::notNA(intr.desc.names)
            intr.desc.genus <- randtip::firstword(intr.desc.names)
            if(all(intr.desc.genus != genus)){

              siblings <- siblings[siblings.genus == genus]
            }
          }
        }

        descs<- phytools::getDescendants(tree, sp.tip)
        desctips<-tree$tip.label[descs]
        desctips<-desctips[!is.na(desctips)]
        desctips<-desctips[randtip::firstword(desctips) == genus]
        genus.tree.list[[i]]<-desctips

      }


      genus.tree.list<-genus.tree.list[!is.na(genus.tree.list[])]
      forbidden.groups<-c(forbidden.groups, genus.tree.list)
    }


  }

  if(any(DF1$phyleticity=="Singleton genus")){
    sing.genera <- DF1[DF1$phyleticity=="Singleton genus", "genus"]
    sing.genera<- unique(sing.genera)
    sing.groups<- rep(list(NA),length(sing.genera))

    for(i in 1:length(sing.genera)){
      sing.taxa <- treelist[treelist$genus == sing.genera[i], "taxon"]
      if(length(sing.taxa)>1){sing.groups[[i]]<-sing.taxa}else{sing.groups[[i]]<-NA}
    }

    sing.groups<- randtip::notNA(sing.groups)
    forbidden.groups<- c(forbidden.groups, sing.groups)
  }


  return(forbidden.groups)
}

get.original.names <- function(tree, DF1, verbose = FALSE){

    for(n in 1:nrow(DF1)){
        tip.label <- tree$tip.label
        if(DF1$using.taxa[n] == DF1$taxon[n]) next
        if(DF1$using.taxa[n] %in% tree$tip.label){
            correct.tip.loc <- tree$tip.label == DF1$using.taxa[n]
            tree$tip.label[correct.tip.loc] <- DF1$taxon[n]
            if(verbose){
                print(paste0("name correction: ", n, "/",
                             nrow(DF1), " (",
                             round(n/nrow(DF1)*100, 2), " %)"))
            }
        }
    }

    return(tree)
}

randtip.subsp <- function(tree, DF1.dupl, verbose = FALSE){

    rep.taxa <- DF1.dupl$taxon
    rep.taxa.species <- unique(DF1.dupl$using.taxa)

    for(i in 1:length(rep.taxa.species)){
        sp.start<- Sys.time()
        # Select subspecies
        genus <- stringr::word(rep.taxa, 1, sep = "_")
        sp <- stringr::word(rep.taxa, 2, sep= "_")
        ssps <- rep.taxa[paste0(genus, "_", sp) == rep.taxa.species[i]]
        ssps <- sample(ssps, length(ssps))
        genus.tree <- stringr::word(tree$tip.label, 1, sep = "_")
        sp.tree <-    stringr::word(tree$tip.label, 2, sep= "_")
        sp.to.add <- tree$tip.label[rep.taxa.species[i]== paste0(genus.tree, "_", sp.tree)]
        if(length(sp.to.add) > 1){
            sp.to.add <- sp.to.add[sp.to.add == rep.taxa.species[i]]
        }
        if(isFALSE(rep.taxa.species[i]%in%tree$tip.label)){stop("Species ", rep.taxa.species[i], " is not in the tree")}

        # Add subspecies as singleton
        new.tree <- add.to.singleton(tree, singleton = sp.to.add,
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

usingMDCCfinder<- function(DF1, taxon, tree){

  levels<- c("genus","tribe","subfamily","family","order","class")
  MDCC.vect<- vector(mode="character", length = length(taxon))
  MDCC.lev.vect<- vector(mode="character", length = length(taxon))


  for(v in 1:length(taxon)){
  i<- which(DF1$taxon==taxon[v])
  MDCC<-NA
  MDCC.levels<-NA
  for(level in levels){
      if(is.na(MDCC)){
        MDCC<-DF1[i, level]
        phyleticity<-randtip::MDCC.phyleticity(DF1, tree = tree,
                     MDCC.info = list(level=level, MDCC= MDCC))
        if(phyleticity=="Not included"){MDCC<-NA}
        lev<-level
      }else{next}
    }
      MDCC.vect[v]<-MDCC
      MDCC.lev.vect[v]<-lev
    }

  return(list(MDCC=MDCC.vect,MDCC.levels=MDCC.lev.vect) )
}



#' randtip "MOTHER FUNCTION"
#' @export
rand.list <- function(tree, DF1,
                    type = "random",
                    aggregate.subspecies = TRUE,
                    insertion = "random",
                    prob = TRUE, verbose = FALSE, polyphyletic.insertion="freq",
                    trim=TRUE){

    start<- Sys.time()

    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    new.tree <- tree
    DF1.dupl <- NULL
    DF1$taxon <- gsub(" ", "_", DF1$taxon)

    if(type=="random"){
        DF1$using.taxa <- get.taxa.to.use(DF1, aggregate.subspecies)
        DF1 <- DF1[order(DF1$using.taxa),]
        is.duplicated <- duplicated(DF1$using.taxa)
        DF1.dupl <- DF1[ is.duplicated,]
        DF1 <-      DF1[!is.duplicated,]

        DF1$using.MDCC <- usingMDCCfinder(DF1 = DF1, taxon = DF1$taxon, tree = new.tree)[[1]]
        DF1$using.MDCC.lev <- usingMDCCfinder(DF1 = DF1, taxon = DF1$taxon, tree = new.tree)[[2]]

        taxa <- DF1$using.taxa
        taxa <- taxa[!(taxa %in% tree$tip.label)]
        taxa.MDCC <- DF1$using.MDCC[DF1$using.taxa%in%taxa]
        taxa.MDCC <- unique(taxa.MDCC)



        for(i in 1: length(taxa.MDCC)){        #loop 1
            gen.start<- Sys.time()
            forbidden.groups <- get.forbidden.groups(new.tree, DF1)

            MDCC <- taxa.MDCC[i]
            level<- unique(DF1$using.MDCC.lev[DF1$using.MDCC==MDCC])
            genus.match <- DF1$using.MDCC==MDCC
            genus.type <- randtip::MDCC.phyleticity(DF1 = DF1, tree = new.tree,
                          MDCC.info = list(level=level, MDCC=MDCC), trim=F)

            MDCC.taxa <- DF1$using.taxa[DF1$using.MDCC == MDCC]
            MDCC.taxa <- sample(MDCC.taxa, length(MDCC.taxa))


            if(verbose){
              print(paste0(i, "/", length(taxa.genera),
                           " (",round(i/length(taxa.genera)*100, 2), " %). ",
                           "Adding Gen. ", genus, " (", genus.type,", ",
                           length(genus.taxa)," tips). "))
            }




            if(genus.type=="Monophyletic"){
                for( j in 1:length(genus.taxa)){
                  new.tree<-add.to.monophyletic(tree = new.tree, new.tip = genus.taxa[j], prob)
                }
            }else if(genus.type=="Paraphyletic"){#ALGUN ERROR EN PARAPHYLETIC!!!!
                for( j in 1:length(genus.taxa)){
                    new.tree <- add.to.paraphyletic(tree = new.tree,
                                                    new.tip = genus.taxa[j], prob)
                }
            }else if(genus.type=="Polyphyletic"){
                new.tree<- add.to.polyphyletic(tree = new.tree, new.tip = genus.taxa,
                                               polyphyletic.insertion, prob)
            }else if(genus.type=="Singleton genus"){
                # All tips added in one step
                genus.match <- (randtip::firstword(tree$tip.label) == genus)
                new.tree <- add.to.singleton(new.tree,
                                             singleton = tree$tip.label[genus.match],
                                             new.tips = genus.taxa)
            }else if(genus.type == "Not included"){

              otherMDCC.DF<- DF1[!is.na(DF1$other.MDCC)&DF1$genus==genus,]
              MDCC.DF     <- DF1[ is.na(DF1$other.MDCC)&DF1$genus==genus,]

              if(nrow(otherMDCC.DF)>0){
                other.MDCCs<-   unique(otherMDCC.DF$other.MDCC)

                for(other.MDCC in other.MDCCs){
                  other.MDCC.taxa<-  DF1$taxon[ !is.na(DF1$other.MDCC)& DF1$other.MDCC==other.MDCC]
                  other.MDCC.genera<-DF1$genus[ !is.na(DF1$other.MDCC)& DF1$other.MDCC==other.MDCC]
                  other.MDCC.inTree<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)%in%other.MDCC.genera]

                  if(length(other.MDCC.inTree)>1){
                     for(sp in other.MDCC.taxa){
                    other.MDCC.MRCA<- ape::getMRCA(new.tree, tip = other.MDCC.inTree)
                    new.tree<- add.into.node(tree = new.tree, new.tip = sp,
                                             node = other.MDCC.MRCA, prob , exception.list = forbidden.groups)
                  }}

                  if(length(other.MDCC.inTree)==1){
                    new.tree<- add.to.singleton(new.tree, singleton =other.MDCC.inTree, new.tips =other.MDCC.taxa)
                  }

                  if(length(other.MDCC.inTree)==0){
                    MDCC.DF<-rbind(MDCC.DF, otherMDCC.DF[otherMDCC.DF$other.MDCC==other.MDCC,])
                   }}}


              if(nrow(MDCC.DF)>0){

                MDCCs<-   unique(MDCC.DF$MDCC)

                for(MDCC in MDCCs){
                  MDCC.taxa<-  DF1$taxon[ !is.na(DF1$MDCC)& DF1$MDCC==MDCC]
                  MDCC.genera<-DF1$genus[ !is.na(DF1$MDCC)& DF1$MDCC==MDCC]
                  MDCC.inTree<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)%in%MDCC.genera]

                  if(length(MDCC.inTree)>1){

                      sp<- genus.taxa[1]
                      MDCC.MRCA<- ape::getMRCA(new.tree, tip = MDCC.inTree)
                      new.tree<- add.into.node(tree = new.tree, new.tip = sp,
                                               node = MDCC.MRCA, prob , exception.list = forbidden.groups)
                      if(length(genus.taxa)>1){
                        new.tree<-add.to.singleton(new.tree, singleton = genus.taxa[1], new.tips = genus.taxa[-1])
                      }



                    }}

                  if(length(MDCC.inTree)==1){
                    new.tree<- add.to.singleton(new.tree, singleton =MDCC.inTree, new.tips =MDCC.taxa)
                  }

                  if(length(MDCC.inTree)==0){
                    stop("Tips ", paste0(MDCC.taxa, collapse = ", ") ," could not be added as its MDCC was not in represented in the tree")
                  }}
              }

            gen.end <- Sys.time()
            if(verbose){
                print(paste0("\U2713", "( done in ",
                           round(as.numeric(difftime(gen.end, gen.start, units = "secs")), 2), " sec. out of ",
                           round(as.numeric(difftime(gen.end, start,     units = "mins")), 2), " mins)"))
            }

        }

        # Phase 2 - subspecies
        if(nrow(DF1.dupl)>0){
            new.tree <- randtip.subsp(tree = new.tree, DF1.dupl, verbose)
        }

        new.tree <- get.original.names(new.tree, DF1)

    }else{
        #In polytomy cases, names are not changed
        DF1$using.taxa <- DF1$taxon

        taxa.table <- DF1[!duplicated(DF1$using.taxa),]
        taxa.table$using.taxa <- gsub(" ", "_", taxa.table$using.taxa)
        taxa.table <- taxa.table[!(taxa.table$using.taxa %in% tree$tip.label),]

        taxa.genera <- taxa.table$genus
        taxa.genera <- sample(taxa.genera, length(taxa.genera))

        new.tree <- add.polytomy.types(new.tree, taxa.table, species.table,
                                       type, insertion, taxa.genera)

    }

    if(!is.null(DF1.dupl)){
        complete.taxa.list <- c(DF1$taxon, DF1.dupl$taxon)
    }else{
        complete.taxa.list <- DF1$taxon
    }

    complete.taxa.list.in.tree <- complete.taxa.list[complete.taxa.list %in% new.tree$tip.label]
    not.included <- complete.taxa.list[!(complete.taxa.list %in% complete.taxa.list.in.tree)]
    if(length(not.included) > 1){
        message(paste0("The following taxa were not included in the tree: ", not.included))
    }

    if(isTRUE(trim)){new.tree <- ape::keep.tip(new.tree, complete.taxa.list.in.tree)}

    return(new.tree)
}


