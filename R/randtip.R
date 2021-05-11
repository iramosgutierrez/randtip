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

get.forbidden.MDCC.nodes <- function(tree,DF1, level, MDCC){
  DF1<- randtip::correct.DF(DF1)
  DF1.mdcc<- DF1[!is.na(DF1[,level]),]
  DF1.mdcc<- DF1.mdcc[DF1.mdcc[,level]==MDCC,]
  forbidden.nodes<- as.numeric(NULL)
  if(level=="genus"){return(NULL)}
  for(lv in 4:(which(colnames(DF1)==level)-1)){
    col.lev<-colnames(DF1.mdcc)[lv]
    MDCCs<-unique(DF1.mdcc[,col.lev])
    MDCCs<-randtip::notNA(MDCCs)
    for(lev in MDCCs){
     phylstat<-randtip::MDCC.phyleticity(DF1, tree = tree,
                                         MDCC.info = list(level=col.lev, MDCC=lev), trim=F )
     if(phylstat%in%c("Monophyletic","Paraphyletic")){
       mdcc.gen<-DF1.mdcc[DF1.mdcc[,col.lev]==lev,"taxon"]
       mdcc.gen<-randtip::firstword(randtip::notNA(mdcc.gen))
       mdcc.sppintree<-sp.genus.in.tree(tree, mdcc.gen)
       mdcc.sppintree.mrca<-ape::getMRCA(tree, tip = mdcc.sppintree)
       mdcc.sppintree.descs<- phytools::getDescendants(tree, mdcc.sppintree.mrca)
       forbidden.nodes<- c(forbidden.nodes, mdcc.sppintree.descs)
    }
      }
    }

  return(unique(forbidden.nodes))
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

usingMDCCfinder<- function(DF1, taxon=NULL, tree, verbose=F){

  if(is.null(taxon)){taxon<- DF1$taxon}
  DF1<-randtip::correct.DF(DF1)
  MDCC.vect<- vector( mode="character", length = length(taxon))
  MDCC.lev.vect<- vector(mode="character", length = length(taxon))
  MDCC.phyletictype.vect<- vector(mode="character", length = length(taxon))

  if(verbose){
    cat(paste0("Searching MDCCs... ", "\n"))
  }

  #manual MDCC search
  taxa<- DF1[!is.na(DF1$taxon1)|!is.na(DF1$taxon2),]
  if(nrow(taxa)>0){
    for(tx in seq_along(taxa$taxon)){
      if(isTRUE(length(strsplit( taxa$taxon1[tx], "_")[[1]])==1 &
         length(strsplit( taxa$taxon2[tx], "_")[[1]])==1 &
         taxa$taxon1[tx]==taxa$taxon2[tx] &
         length(randtip::sp.genus.in.tree(tree, taxa$taxon1[tx]))>0)){

        pos<- which(taxon==taxa$taxon)
        MDCC.vect[pos] <- taxa$taxon1[tx]
        MDCC.lev.vect[pos] <- "Sister genus"
        next
      }

      if(!(taxa$taxon1[tx]%in%tree$tip.label)){taxa$taxon1[tx]<-NA}
      if(!(taxa$taxon2[tx]%in%tree$tip.label)){taxa$taxon2[tx]<-NA}

      if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
         taxa$taxon2[tx] %in% tree$tip.label &
         taxa$taxon1[tx]==taxa$taxon2[tx])){

        pos<- which(taxon==taxa$taxon)
        MDCC.vect[pos] <- taxa$taxon1[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        next
      }

      if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
                taxa$taxon2[tx] %in% tree$tip.label &
                taxa$taxon1[tx]!=taxa$taxon2[tx])){

        pos<- which(taxon==taxa$taxon)
        MDCC.vect[pos] <- paste0("Clade (", taxa$taxon1[tx], "-", taxa$taxon2[tx], ")")
        MDCC.lev.vect[pos] <- "Manual clade"
        next
      }

       if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
         is.na(taxa$taxon2[tx]))){

        pos<- which(taxon==taxa$taxon)
        MDCC.vect[pos] <- taxa$taxon1[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        next
       }

      if(isTRUE(taxa$taxon2[tx] %in% tree$tip.label &
         is.na(taxa$taxon1[tx]))){

        pos<- which(taxon==taxa$taxon)
        MDCC.vect[pos] <- taxa$taxon2[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        next
      }

         }
  }

  #automatic MDCC search
  levels<- c("genus","tribe","subfamily","family","order","class")
  for(v in 1:length(taxon)){

  if(verbose){
      if(v %in% c(seq(0, length(taxon), 10), length(taxon))){
        cat(paste0(round((v/length(taxon)*100),2), "% completed.\n"))
      }

    }

  i<- which(DF1$taxon==taxon[v])
  if((MDCC.vect[v])==""){

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


  }
  if(is.na(MDCC.vect[v])|is.na(MDCC.lev.vect[v])){MDCC.phyletictype.vect[v]<-NA}else{
  MDCC.phyletictype.vect[v]<-randtip::MDCC.phyleticity(DF1 = DF1, tree =  tree,
          MDCC.info = list(level=MDCC.lev.vect[v], MDCC=MDCC.vect[v]), trim = T)
  }
  }


  return(list(MDCC=MDCC.vect,MDCC.levels=MDCC.lev.vect, MDCC.phylstat=MDCC.phyletictype.vect) )
}



#' randtip "MOTHER FUNCTION"
#' @export
rand.list <- function(tree, DF1,type = "random",agg.ssp = TRUE,
                    prob = TRUE, verbose = FALSE,
                    poly.ins="all",trim=TRUE, forceultrametric=F){

    start<- Sys.time()

    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    DF1<- correct.DF(DF1)
    DF1$taxon <- gsub(" ", "_", DF1$taxon)

    new.tree <- tree

    tree$tip.label <- gsub("_x_", "_x-", tree$tip.label)
    tree$tip.label <- gsub("_X_", "_x-", tree$tip.label)

    if(forceultrametric & !ape::is.ultrametric(new.tree)){new.tree<- phytools::force.ultrametric(new.tree)}
    if(isFALSE(forceultrametric) & !ape::is.ultrametric(new.tree)){
      message("Specified tree is not ultrametric. \nTo force the randomization as an ultrametric tree plase set forceultrametric=TRUE")}

    DF1.rand <- NULL
    DF1.poly <- NULL
    DF1.dupl <- NULL


    DF1_search<- randtip::usingMDCCfinder(DF1 = DF1, taxon = DF1$taxon, tree = new.tree, verbose)
    DF1$using.MDCC     <- DF1_search[[1]]
    DF1$using.MDCC.lev <- DF1_search[[2]]
    DF1$using.MDCC.phylstat <- DF1_search[[3]]

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
    DF1.nonpoly <- DF1.nonpoly[order(DF1.nonpoly$taxon),]
    is.duplicated <- duplicated(DF1.nonpoly$using.taxa)
    DF1.dupl <- DF1.nonpoly[ is.duplicated,]
    DF1.rand <- DF1.nonpoly[!is.duplicated,]
    }


    #Phase 1. Random insertions, non-aggregated
    if(!is.null(DF1.rand)){if(nrow(DF1.rand)>0){

        DF1.rand.bind<- DF1.rand[!(DF1.rand$using.taxa %in% new.tree$tip.label),]
        DF1.rand.bind<- DF1.rand.bind[!is.na(DF1.rand.bind$using.MDCC),]

        binding.units<- unique(paste0(randtip::firstword(DF1.rand.bind$using.taxa), "-",DF1.rand.bind$using.MDCC))

        if(verbose){
          cat(paste0("Starting randomization","\n")) }

        for(i in seq_along(binding.units)){
            binding.unit <- binding.units[i]
            genus <- strsplit(binding.unit, "-")[[1]][1]
            MDCC  <- strsplit(binding.unit, "-")[[1]][2]
            level <- unique(DF1.rand.bind$using.MDCC.lev[DF1.rand.bind$using.MDCC==MDCC])
if(length(MDCC)>1){stop("Several MDCC levels recognised for MDCC ", MDCC, ". Please correct your DF1")}

            MDCC.type <- unique(DF1.rand.bind$using.MDCC.phylstat[DF1.rand.bind$using.MDCC==MDCC])
if(length(MDCC.type)>1){stop("Several MDCC phyletic statuses recognised for MDCC ", MDCC, ". Please correct your DF1")}

            unit.taxa <- DF1.rand.bind$using.taxa[randtip::firstword(DF1.rand.bind$using.taxa)==genus&
                                                   DF1.rand.bind$using.MDCC==MDCC]
            unit.taxa <- sample(unit.taxa, length(unit.taxa))


            if(verbose){
              cat(paste0(i, "/", length(binding.units),
                           " (",round(i/length(binding.units)*100, 2), "%). ",
                           "Adding ", genus, " to ", MDCC ," (", MDCC.type, " ", level,", ",
                           length(unit.taxa)," PUTs).\n")) }

            #Manual additions
            if(level=="Sister species"){
              new.tree <- add.to.singleton(tree=new.tree, singleton = MDCC ,new.tips =  unit.taxa)
            }

            if(level=="Sister genus"){
              sister.genus.tips<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)==MDCC]
              sister.genus.mrca<- ape::getMRCA(new.tree, sister.genus.tips)

                new.tree<-add.over.node(tree= new.tree, new.tip = unit.taxa[j], node = sister.genus.mrca)}

            if(level=="Manual clade"){
            table<- DF1.rand.bind[DF1.rand.bind$using.taxa%in%unit.taxa,]
            sp1<- unique(table$taxon1)
            sp2<- unique(table$taxon2)
            clade.mrca<- ape::getMRCA(new.tree, c(sp1, sp2))
            for(j in seq_along(unit.taxa)){
              new.tree <- add.into.node(tree = new.tree, node = clade.mrca,
                                        new.tip = unit.taxa[j], prob = prob )}
              }



            #Automatically searched MDCCs additions
            #Add to genus
            if(level=="genus"){
            if(MDCC.type=="Monophyletic"){
                for( j in 1:length(unit.taxa)){
                  new.tree<-add.to.monophyletic(tree = new.tree, new.tip = unit.taxa[j], prob=prob)
                  }
            }else if(MDCC.type=="Paraphyletic"){#ALGUN ERROR EN PARAPHYLETIC!!!!
                for( j in 1:length(unit.taxa)){
                    new.tree <- add.to.paraphyletic(tree = new.tree,
                                                    new.tip = unit.taxa[j], prob = prob)
                }
            }else if(MDCC.type=="Polyphyletic"){
              poly.ins<- unique(DF1.rand$poly.ins[DF1.rand$using.taxa %in% unit.taxa])
if(length(poly.ins)>1){stop("Several Polyphyletic insertions recognised for genus ", genus, ". Please correct your DF1")}
                new.tree<- add.to.polyphyletic(tree = new.tree, new.tip = unit.taxa,
                                               poly.ins , prob)
            }else if(MDCC.type=="Singleton MDCC"){
                # All tips added in one step
                singleton <- tree$tip.label[(randtip::firstword(tree$tip.label) == genus)]
                new.tree <- add.to.singleton(tree = new.tree,
                                             singleton = singleton,
                                             new.tips = unit.taxa)}
            }
            #Add to other taxonomic MDCC
            if(level%in% c("tribe", "subfamily", "family", "order", "class")){
              DF1.mdcc<-  DF1[!is.na(DF1[,level]),]
              MDCC.taxa<- DF1.mdcc$taxon[DF1.mdcc[,level]==MDCC]
              MDCC.genera<- randtip::notNA(unique(randtip::firstword(MDCC.taxa)))
              MDCC.intree<- sp.genus.in.tree(new.tree, MDCC.genera)

              if(length(MDCC.intree)==1){new.tree<- randtip::add.to.singleton(new.tree, MDCC.intree, unit.taxa)}else

              if(length(MDCC.intree)>1){
              MDCC.mrca<- ape::getMRCA(new.tree, MDCC.intree)
              permitted.nodes<-randtip::get.permitted.nodes(tree=new.tree, node = MDCC.mrca)
              forbidden.nodes<- randtip::get.forbidden.MDCC.nodes(new.tree, DF1, level, MDCC)
              permitted.nodes<- permitted.nodes[!(permitted.nodes%in%forbidden.nodes)]
              if(length(permitted.nodes)==0){
                new.tree<-add.over.node(new.tree, new.tip = unit.taxa, node = MDCC.mrca)
              }else{
              if((length(new.tree$tip.label)+1)%in%permitted.nodes){
                permitted.nodes<-permitted.nodes[-which(permitted.nodes==(length(new.tree$tip.label)+1))]}

              if(length(permitted.nodes)>1){permitted.nodes<-sample(x = permitted.nodes, size = 1, replace = F)}

              new.tree <- randtip::add.over.node(new.tree, new.tip = unit.taxa[1], node = permitted.nodes)
              if(length(unit.taxa)>1){
              new.tree <- randtip::add.to.singleton(tree = new.tree,singleton = unit.taxa[1],
                                        new.tips = unit.taxa[2:length(unit.taxa)] )}}}
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

        DF1.poly<- DF1.poly[!is.na(DF1.poly$using.MDCC),]
        MDCCs.i.level<- unique(DF1.poly$using.MDCC.lev[DF1.poly[,"using.MDCC"]==MDCCs.i])


        MDCC.taxa.toAdd <- DF1.poly$taxon[DF1.poly[,"using.MDCC"]==MDCCs.i]
        MDCC.taxa.inDF1 <- DF1$taxon[DF1[, MDCCs.i.level]==MDCCs.i]
        MDCC.genera <- unique(randtip::firstword(MDCC.taxa.inDF1))
        MDCC.genera <- randtip::notNA(MDCC.genera)
        MDCC.taxa.inTree<- randtip::sp.genus.in.tree(new.tree, MDCC.genera)

        if(verbose){
          cat(paste0(which(MDCCs==MDCCs.i), "/", length(MDCCs),
                     " (",round(which(MDCCs==MDCCs.i)/length(MDCCs)*100, 2), " %). ",
                     "Adding ", MDCCs.i ," (",
                     length(MDCC.taxa.toAdd)," tips).\n")) }


        if(length(MDCC.taxa.inTree)==1){
          new.tree<- randtip::polytomy.to.singleton(tree = new.tree, singleton = MDCC.taxa.inTree,
                                        new.tip =MDCC.taxa.toAdd, insertion = "long")}

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
    tree$tip.label <- gsub("_x-", "_x_", tree$tip.label)
    return(new.tree)
}






