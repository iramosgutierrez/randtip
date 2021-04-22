get.taxa.to.use <- function(DF1){
  DF1$using.taxa <- NA

    for(i in 1:nrow(DF1)){
        taxon <- DF1$taxon[i]
        genus.name <- stringr::word(DF1$taxon[i],1, sep = "_")
        sp.name <-    stringr::word(DF1$taxon[i],2, sep = "_")
        if(DF1$agg.ssp[i]=="1"){
                DF1$using.taxa[i] <- paste0(genus.name, "_", sp.name) }

        if(DF1$agg.ssp[i]=="0"){
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

get.permitted.nodes <- function (tree, node){

  nodes<- phytools::getDescendants(tree, node)
  if(length(nodes)==0){stop("Node ", node, " is not present in the given tree")}

  tip.names<- randtip::notNA(tree$tip.label[nodes])
  tip.names.genera<- unique(randtip::firstword(tip.names))
  tip.types.genera<- rep(list(NA), times=length(tip.names.genera))
  names(tip.types.genera)<-tip.names.genera

  tip.nodes.genera<- rep(list(NA), times=length(tip.names.genera))
  names(tip.nodes.genera)<-tip.names.genera

  paraphyletic.intruders <- rep(NA, 1)
  for(g in seq_along(tip.names.genera)){
    gen<- tip.names.genera[g]
    type<- randtip::phyleticity(tree, gen)
    tip.types.genera[g]<-type
    if(type=="Paraphyletic"){
      gen.spp<- tree$tip.label[randtip::firstword(tree$tip.label)==gen]
      mrca<- ape::getMRCA(tree, gen.spp)
      para.descs<- phytools::getDescendants(tree, mrca)
      para.descs<-randtip::notNA(tree$tip.label[para.descs])
      intruders<-para.descs[randtip::firstword(para.descs)!=gen]
      paraphyletic.intruders <-c(paraphyletic.intruders, intruders)
    }
    paraphyletic.intruders<- randtip::notNA(paraphyletic.intruders)
  }

  permitted.nodes<- rep(NA, times= length(nodes))


  for(i in seq_along(nodes)){
    node.i <- nodes[i]


    if(randtip::is.tip(tree, node.i)){
      tip<- tree$tip.label[node.i]
      tip.genus<- randtip::firstword(tip)
      sp.in.tree<-randtip::sp.genus.in.tree(tree, tip.genus)

      if(length(sp.in.tree)!=1){permitted.nodes[i]<-FALSE}
      if(tip%in%paraphyletic.intruders){permitted.nodes[i]<-FALSE}

    }else{

      descs<- phytools::getDescendants(tree, node.i)
      descs<-randtip::notNA (tree$tip.label[descs])
      descs.genera<- unique(randtip::firstword(descs))

      siblings<-randtip::get.parent.siblings(tree, node.i)$siblings
      siblings.genera<- unique(randtip::firstword(siblings))

      if(all(siblings.genera %in% descs.genera)){permitted.nodes[i]<-FALSE}


         }


  }

  paraphyletic.genera<-as.vector(which(tip.types.genera=="Paraphyletic") )
  for(g in paraphyletic.genera){
    gen<-names(tip.types.genera[g])
    gentips<- tree$tip.label[randtip::firstword(tree$tip.label)==gen]
    mrca<- ape::getMRCA(tree, gentips)
    mrcadescs<-phytools::getDescendants(tree, mrca)

    permitted.nodes[which(nodes %in% mrcadescs)]<- FALSE
  }


  permitted.nodes[is.na(permitted.nodes)]<-TRUE

  if(length(nodes[(permitted.nodes)==TRUE])>0){return(nodes[(permitted.nodes)==TRUE])}else{return(node)}
  }

get.original.names <- function(tree, DF1, verbose = FALSE){
  if(verbose){
    cat("Starting name checking...", "\n")
  }

    for(n in 1:nrow(DF1)){
        tip.label <- tree$tip.label
        if(DF1$using.taxa[n] == DF1$taxon[n]) next
        if(DF1$using.taxa[n] %in% tree$tip.label&
           !(DF1$using.taxa[n] %in% DF1$taxon)){
            correct.tip.loc <- tree$tip.label == DF1$using.taxa[n]
            tree$tip.label[correct.tip.loc] <- DF1$taxon[n]

        }
    }
  if(verbose){
    cat("\U2713", "Names checking completed", "\n")
  }
    return(tree)
}

randtip.subsp <- function(tree, DF1.dupl, verbose = FALSE){
    new.tree<-tree
    DF1.dupl<- DF1.dupl[!(DF1.dupl$taxon%in%new.tree$tip.label),]

    rep.taxa <- DF1.dupl$taxon
    rep.taxa.species <- unique(DF1.dupl$using.taxa)

    genus <- stringr::word(rep.taxa, 1, sep = "_")
    sp <- stringr::word(rep.taxa, 2, sep= "_")

    for(i in 1:length(rep.taxa.species)){

        # Select subspecies
        ssps <- rep.taxa[paste0(genus, "_", sp) == rep.taxa.species[i]]
        ssps <- sample(ssps, length(ssps))
        genus.tree <- stringr::word(new.tree$tip.label, 1, sep = "_")
        sp.tree <-    stringr::word(new.tree$tip.label, 2, sep= "_")
        sp.to.add <- new.tree$tip.label[rep.taxa.species[i]== paste0(genus.tree, "_", sp.tree)]
        if(length(sp.to.add) > 1){
            sp.to.add <- sp.to.add[sp.to.add == rep.taxa.species[i]]
        }
        if(!(rep.taxa.species[i]%in%new.tree$tip.label)){stop("Species ", rep.taxa.species[i], " is not in the tree")}

        # Add subspecies as singleton
        new.tree <- add.to.singleton(new.tree, singleton = sp.to.add,
                                     new.tips = ssps)
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

usingMDCCfinder<- function(DF1, taxon, tree, verbose=F){

  levels<- c("genus","tribe","subfamily","family","order","class")
  MDCC.vect<- vector(mode="character", length = length(taxon))
  MDCC.lev.vect<- vector(mode="character", length = length(taxon))

  if(verbose){
    cat(paste0("Searching MDCCs... ", "\n"))
    }
  for(v in 1:length(taxon)){
  i<- which(DF1$taxon==taxon[v])
  MDCC<-as.character(NA)
  MDCC.levels<-as.character(NA)
  for(level in levels){
      if(is.na(MDCC)){
        MDCC<-as.character(DF1[i, level])
        if(!is.na(MDCC)){phyleticity<-randtip::MDCC.phyleticity(DF1, tree = tree,
                     MDCC.info = list(level=level, MDCC= MDCC))
        if(phyleticity=="Not included"){MDCC<-NA}
        }

        lev<-level
      }else{next}
    }
      MDCC.vect[v]<-as.character(MDCC)
      MDCC.lev.vect[v]<-as.character(lev)
      if(is.na(MDCC)){MDCC.lev.vect[v]<-NA}

      if(verbose){
        if(v %in% c(seq(0, length(taxon), 10), length(taxon))){
          cat(paste0(round((v/length(taxon)*100),2), "% completed.\n"))
        }

      }
    }

  return(list(MDCC=MDCC.vect,MDCC.levels=MDCC.lev.vect) )
}



#' randtip "MOTHER FUNCTION"
#' @export
rand.list <- function(tree, DF1,type = "random",agg.ssp = TRUE,
                    prob = TRUE, verbose = FALSE,
                    poly.ins="all",trim=TRUE){

    start<- Sys.time()

    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    DF1<- correct.DF(DF1)
    DF1$taxon <- gsub(" ", "_", DF1$taxon)

    new.tree <- tree

    DF1.rand <- NULL
    DF1.poly <- NULL
    DF1.dupl <- NULL


    DF1_search<- usingMDCCfinder(DF1 = DF1, taxon = DF1$taxon, tree = new.tree, verbose)
    DF1$using.MDCC     <- DF1_search[[1]]
    DF1$using.MDCC.lev <- DF1_search[[2]]
    if(verbose){cat("\n")}

    if(trim){
      trimming.species<- rep(NA, 1)
      for(using.mdcc in unique(DF1$using.MDCC)){
        spp.df<-DF1[DF1$using.MDCC==using.mdcc,]
        using.level<- randtip::notNA(unique(spp.df$using.MDCC.lev))

        mdcc.genera<-randtip::firstword(DF1$taxon[DF1[,using.level]==using.mdcc])
        mdcc.species<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)%in%mdcc.genera]
        trimming.species<- c(trimming.species, mdcc.species)
      }
        trimming.species<- randtip::notNA(trimming.species)
        new.tree <- ape::keep.tip(new.tree, trimming.species)}

    if(isTRUE(agg.ssp)){DF1$agg.ssp[is.na(DF1$agg.ssp)]<-"1"}else{DF1$agg.ssp[is.na(DF1$agg.ssp)]<-"0"}
    if(type=="random"){DF1$rand.type[is.na(DF1$rand.type)]<-"0"}else{DF1$rand.type[is.na(DF1$rand.type)]<-"1"}
    {DF1$poly.ins[is.na(DF1$poly.ins)]<-poly.ins}


    DF1.poly    <- DF1[DF1$rand.type=="1",]
    DF1.nonpoly <- DF1[DF1$rand.type=="0",]


    if(nrow(DF1.nonpoly)>0){
    DF1.nonpoly$using.taxa <- get.taxa.to.use(DF1.nonpoly)
    DF1.nonpoly <- DF1.nonpoly[order(DF1.nonpoly$using.taxa),]
    is.duplicated <- duplicated(DF1.nonpoly$using.taxa)
    DF1.dupl <- DF1.nonpoly[ is.duplicated,]
    DF1.rand <- DF1.nonpoly[!is.duplicated,]
    }


    #Phase 1. Random insertions, non-aggregated
    if(!is.null(DF1.rand)){if(nrow(DF1.rand)>0){
        taxa <- DF1.rand$using.taxa
        taxa <- taxa[!(taxa %in% tree$tip.label)]
        taxa.genera <- randtip::firstword(taxa)
        taxa.genera <- unique(taxa.genera)

        if(verbose){
          cat(paste0("Starting randomization","\n")) }

        for(i in seq_along(taxa.genera)){
            genus <- taxa.genera[i]
            MDCC  <- unique(DF1.rand$using.MDCC    [randtip::firstword(DF1.rand$using.taxa)==genus])
            level <- unique(DF1.rand$using.MDCC.lev[randtip::firstword(DF1.rand$using.taxa)==genus])
if(length(MDCC)>1){stop("Several MDCCs recognised for genus ", genus, ". Please correct your DF1")}

            MDCC.type <- randtip::MDCC.phyleticity(DF1 = DF1.rand, tree = new.tree,
                          MDCC.info = list(level=level, MDCC=MDCC), trim)

            genus.taxa <- taxa[randtip::firstword(taxa)==genus]
            genus.taxa <- sample(genus.taxa, length(genus.taxa))


            if(verbose){
              cat(paste0(i, "/", length(taxa.genera),
                           " (",round(i/length(taxa.genera)*100, 2), "%). ",
                           "Adding ", genus, " to ", MDCC ," (", MDCC.type,", ",
                           length(taxa.genera)," tips).\n")) }




            if(level=="genus"){
            if(MDCC.type=="Monophyletic"){
                for( j in 1:length(genus.taxa)){
                  new.tree<-add.to.monophyletic(tree = new.tree, new.tip = genus.taxa[j], prob=prob)
                  }
            }else if(MDCC.type=="Paraphyletic"){#ALGUN ERROR EN PARAPHYLETIC!!!!
                for( j in 1:length(genus.taxa)){
                    new.tree <- add.to.paraphyletic(tree = new.tree,
                                                    new.tip = genus.taxa[j], prob)
                }
            }else if(MDCC.type=="Polyphyletic"){
              poly.ins<- unique(DF1.rand$poly.ins[DF1.rand$using.taxa %in% genus.taxa])
if(length(poly.ins)>1){stop("Several Polyphyletic insertions recognised for genus ", genus, ". Please correct your DF1")}
                new.tree<- add.to.polyphyletic(tree = new.tree, new.tip = genus.taxa,
                                               poly.ins , prob)
            }else if(MDCC.type=="Singleton MDCC"){
                # All tips added in one step
                singleton <- tree$tip.label[(randtip::firstword(tree$tip.label) == genus)]
                new.tree <- add.to.singleton(tree = new.tree,
                                             singleton = singleton,
                                             new.tips = genus.taxa)}
            }else{

              MDCC.taxa<- DF1$taxon[DF1[,level]==MDCC]
              MDCC.genera<- randtip::notNA(unique(randtip::firstword(MDCC.taxa)))
              MDCC.intree<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)%in%MDCC.genera]

              MDCC.mrca<- ape::getMRCA(new.tree, MDCC.intree)
              permitted.nodes<-get.permitted.nodes(tree=new.tree, node = MDCC.mrca)
              if((length(new.tree$tip.label)+1)%in%permitted.nodes){
                permitted.nodes<-permitted.nodes[-which(permitted.nodes==(length(new.tree$tip.label)+1))]}

              if(length(permitted.nodes)>1){permitted.nodes<-sample(x = permitted.nodes, size = 1, replace = F)}

              new.tree <- randtip::add.over.node(new.tree, new.tip = genus.taxa[1], node = permitted.nodes)
              if(length(genus.taxa)>1){
              new.tree <- randtip::add.to.singleton(tree = new.tree,singleton = genus.taxa[1],
                                        new.tips = genus.taxa[2:length(genus.taxa)] )}
            }



        }
      }}
    #Phase 2 - Random insertions, aggregated subspecies
    if(!is.null(DF1.dupl)){if(nrow(DF1.dupl)>0){
      DF1.masdupl<- DF1.rand[(DF1.rand$using.taxa%in%tree$tip.label) &
                             !(DF1.rand$taxon%in%tree$tip.label) , ]
      DF1.dupl<- rbind(DF1.dupl, DF1.masdupl)
      new.tree <- randtip.subsp(tree = new.tree, DF1.dupl, verbose)
            }}

    if(verbose){cat("\n")}
    if(nrow(DF1.nonpoly)>0){new.tree <- get.original.names(tree = new.tree, DF1 = DF1.nonpoly, verbose )}


    #Polytomies
    if(nrow(DF1.poly)>0){
      if(verbose){
        cat(paste0("\n","Starting polytomies addition \n"))
      }
      DF1.poly<-DF1.poly[!(DF1.poly$taxon %in% new.tree$tip.label),]


      MDCCs<- randtip::notNA(unique(DF1.poly$using.MDCC))
      for(MDCCs.i in MDCCs){

        MDCCs.i.level<- unique(DF1.poly$using.MDCC.lev[DF1.poly[,"using.MDCC"]==MDCCs.i])
        MDCCs.i.level<- randtip::notNA(MDCCs.i.level)
        NAs.i <- which(is.na(DF1.poly[ , MDCCs.i.level]))
        if(length(NAs.i)>0) {DF1.poly<- DF1.poly[-NAs.i,]}
        MDCC.taxa.toAdd <- DF1.poly$taxon[DF1.poly[,"using.MDCC"]==MDCCs.i]
        MDCC.taxa.inDF1 <- DF1$taxon[DF1[, MDCCs.i.level]==MDCCs.i]
        MDCC.genera <- unique(randtip::firstword(MDCC.taxa.inDF1))
        MDCC.taxa.inTree<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)%in%MDCC.genera]

        if(verbose){
          cat(paste0(which(MDCCs==MDCCs.i), "/", length(MDCCs),
                     " (",round(which(MDCCs==MDCCs.i)/length(MDCCs)*100, 2), " %). ",
                     "Adding ", MDCCs.i ," (",
                     length(MDCC.taxa.toAdd)," tips).\n")) }


        if(length(MDCC.taxa.inTree)==1){
          new.tree<- randtip::polytomy.to.singleton(new.tree, MDCC.taxa.inTree,
                                            MDCC.taxa.toAdd, insertion = "long")}

        if(length(MDCC.taxa.inTree)>1) {
          MDCC.mrca<- ape::getMRCA(new.tree, MDCC.taxa.inTree)
          new.tree<- randtip::polytomy.into.node(tree=new.tree,
                            new.tip =MDCC.taxa.toAdd, node =  MDCC.mrca)}
            }



    }




    complete.taxa.list <- DF1$taxon
    complete.taxa.list.in.tree <- complete.taxa.list[complete.taxa.list %in% new.tree$tip.label]
    not.included <- complete.taxa.list[!(complete.taxa.list %in% complete.taxa.list.in.tree)]
    if(length(not.included) > 0){
        message("The following taxa were not included in the tree:\n", paste0(not.included, "\n"))}

    if(isTRUE(trim)){new.tree <- ape::keep.tip(new.tree, complete.taxa.list.in.tree)}


    end <- Sys.time()
    if(verbose){
      cat(paste0("\n","\n","\U2713", " Tip Randomization completed in ",
                 round(as.numeric(difftime(end, start,units = "mins")), 2), " mins\n"))
    }

    return(new.tree)
}






