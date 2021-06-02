# Return a vector of all species in the tree that match a given genus
sp.genus.in.tree <- function(tree, genus){
    sp <- tree$tip.label
    taxa.vector <- sp[randtip::firstword(sp)%in%genus]

    return(taxa.vector)
}

get.position <- function(tree, node){
if(randtip::isRoot(tree, node)){position=0
}else{
    df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                       "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )
  edge.length <- df[df$node==node,"length"]

        # Bind at a random point of the branch
        position <- edge.length * stats::runif(1, 0, 1)
}

    return(position)
}

binding.position<- function(tree, node, df=NULL, insertion,  prob){
  position<-list("length"=NA, "where"=NA, "position"=NA)

  if(is.null(df)){df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2], "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )}

  if(ape::is.ultrametric(tree)){position$length<-NULL}else{
    position$length<-abs(rnorm(1, mean=mean(tree$edge.length), sd= sd(tree$edge.length) ))}

    df<- df[df$node==node,]
    position$where<- node

  if(insertion=="polytomy"){
    position$position<- 0
    position$where <- node
  }

  if(insertion=="random"){

    index<- df[df$node==node,"id"]
    pos <- randtip::get.position(tree = tree, node = node)

    position$position<- pos
    position$where <- df[df$id==index,"node"]
  }



  return(position)






}

get.parent.siblings <- function(tree, tip){
    tree.sp <- tree$tip.label
    # Direct ancestor
    parent <- tree$edge[tree$edge[,2] == tip, 1]
    # Ancestor's descendants
    parent.desc <- phytools::getDescendants(tree, parent,curr = NULL)
    siblings <- tree.sp[parent.desc]
    siblings <- siblings[!is.na(tree.sp[parent.desc])]

    return(list(parent = parent,
                siblings = siblings))
}

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

  nodes<- phytools::getDescendants(tree, node,curr = NULL)
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
      para.descs<- phytools::getDescendants(tree, mrca,curr = NULL)
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

      descs<- phytools::getDescendants(tree, node.i,curr = NULL)
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
    mrcadescs<-phytools::getDescendants(tree, mrca,curr = NULL)

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
  for(lv in randtip::randtip_levels()[-1]){
    col.lev<- as.character(lv)
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
        mdcc.sppintree.descs<- phytools::getDescendants(tree, mdcc.sppintree.mrca,curr = NULL)
        forbidden.nodes<- c(forbidden.nodes, mdcc.sppintree.descs)
      }
    }
  }

  return(unique(forbidden.nodes))
}

get.intruder.nodes <- function (tree, DF1, level, MDCC){
  DF1<- randtip::correct.DF(DF1)
  DF1.mdcc<- DF1[!is.na(DF1[,level]),]
  DF1.mdcc<- DF1.mdcc[DF1.mdcc[,level]==MDCC,]

  mdcc.genera<- randtip::firstword(DF1.mdcc$taxon)
  mdcc.intree<- randtip::sp.genus.in.tree(tree, mdcc.genera)
  mdcc.mrca<- ape::getMRCA(tree, mdcc.intree)

  descendants<- phytools::getDescendants(tree, node=mdcc.mrca,curr = NULL)
  descendants<- tree$tip.label[descendants]
  descendants<- randtip::notNA(descendants)

  non.mdcc <- DF1[!is.na(DF1[,level]),]
  non.mdcc <- non.mdcc[non.mdcc[,level]!=MDCC,]
  non.mdcc.genera<- unique(randtip::firstword(non.mdcc$taxon))

  intruders<- descendants[randtip::firstword(descendants)%in%non.mdcc.genera ]

  if(length(intruders)==0){return(NULL)}
  if(length(intruders)==1){return(which(tree$tip.label==intruders))}
  if(length(intruders)> 1){
    intruder.mrca<- ape::getMRCA(tree, intruders)
    intruder.descs<-phytools::getDescendants(tree, node=intruder.mrca,curr = NULL)

    forbidden.nodes<- as.numeric(NULL)
    for(int in intruder.descs){
      des.i.nodes<- phytools::getDescendants(tree, node=int,curr = NULL)
      des.i<- tree$tip.label[des.i.nodes]
      des.i<- randtip::notNA(des.i)
      if(!any(mdcc.intree%in%des.i)){forbidden.nodes<- c(forbidden.nodes,des.i.nodes )}
    }
    return(unique(forbidden.nodes))

  }
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
    if(!(rep.taxa.species[i]%in%new.tree$tip.label)){sp.to.add <- sample(sp.to.add,1)}

    # Add subspecies as singleton
    for(ssp in ssps){
      singleton<-c(sp.to.add,ssps)
      singleton<- singleton[singleton%in%new.tree$tip.label]
      new.tree <- add.to.singleton(new.tree, singleton = singleton,
                                   new.tips = ssp, resp.sing = T, resp.mono = T)}
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
  levels<- randtip::randtip_levels()
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

DF1finder<- function(DF1, taxon, column){
  return(as.character(DF1[DF1$taxon==taxon, column]))
}

firstword<- function(string, sep="_"){
  word<- stringr::word(string, 1, sep=sep)
  return(word)
}

notNA <- function(x){
  vect<- x[!is.na(x)]
  return(vect)
}

get.groups <- function(tree, genus){
  species <- randtip::sp.genus.in.tree(tree, genus)
  sp.mrca<- ape::getMRCA(tree, species)
  mrca.descs <- phytools::getDescendants(tree, sp.mrca,curr = NULL)

  node.descs<- rep(list(NA), times=length(mrca.descs))
  names(node.descs)<- mrca.descs

  node.types<- rep(list(NA), times=length(mrca.descs))
  names(node.types)<- mrca.descs

  for(n in seq_along(mrca.descs)){
     nd<- mrca.descs[n]

    if(randtip::is.tip(tree = tree,node = nd)){
      node.descs[[n]]<- tree$tip.label[nd]
      node.types[[n]]<- "tip"
      next
    }

    nd.descs<- phytools::getDescendants(tree, nd,curr = NULL)
    nd.descs <- randtip::notNA(tree$tip.label[nd.descs])
    nd.genera<- randtip::firstword(nd.descs)

    if(!(genus%in%nd.genera)){
      node.descs[[n]]<- NA
      node.types[[n]]<- NA
      next
    }


    siblings<-randtip::get.parent.siblings(tree, tip = nd)$siblings
    siblings<-siblings[-which(siblings%in%nd.descs)]
    siblings.genus<- randtip::firstword(siblings)
    if(all(nd.genera==genus)){

    # if(genus%in%siblings.genus){
    #   node.descs[[n]]<- NA
    #   node.types[[n]]<- NA}else{
      node.descs[[n]]<- nd.descs
      node.types[[n]]<- "monophyletic"#}
    next}

      intruders<- nd.descs[nd.genera!=genus]

      if(length(intruders)==1){
        node.descs[[n]]<- nd.descs
        node.types[[n]]<- "paraphyletic"
        next}

        intruders.mrca<- ape::getMRCA(tree, intruders)

        intruders.descs<- phytools::getDescendants(tree, intruders.mrca,curr = NULL)
        intruders.descs <- randtip::notNA(tree$tip.label[intruders.descs])
        intruders.genera<- randtip::firstword(intruders.descs)

        if(genus %in% intruders.genera){
          node.descs[[n]]<- nd.descs
          node.types[[n]]<- "polyphyletic"
          next}
        if(length(nd.genera[nd.genera==genus])==1){
            node.descs[[n]]<- nd.descs[nd.genera==genus]
            node.types[[n]]<- "singleton"
            next}

          group<- nd.descs[nd.genera==genus]
          group.mrca<- ape::getMRCA(tree, group)
          group.descs<- phytools::getDescendants(tree, group.mrca,curr = NULL)
          group.descs<- randtip::notNA(tree$tip.label[group.descs])
          group.descs.gen <- randtip::firstword(group.descs)
          if(!all(group.descs.gen==genus)){
            node.descs[[n]]<- nd.descs
            node.types[[n]]<- "paraphyletic"
          }else{
            node.descs[[n]]<- NA
            node.types[[n]]<- NA
          }
        }






  node.types<-node.types[as.vector(which(!is.na(node.types)))]
  node.descs<-node.descs[as.vector(which(!is.na(node.descs)))]

  permitted.type<- as.vector(unlist(node.types))
  permitted.type<- which(!(permitted.type%in% c("polyphyletic")))

  permitted.gen<- grep(paste0(genus,"_"), (node.descs))

  permitted<- intersect(permitted.type,permitted.gen)

  node.types<-node.types[permitted]
  node.descs<-node.descs[permitted]




  return(list(species=node.descs, type=node.types))
}

correct.DF<- function(DF){
  for(i in 1:(ncol(DF))){
    vect<- which(DF[,i]=="")
    if(length(vect)>0){DF[vect,i]<-NA}
    }
  return(DF)
}

isRoot<- function(tree, node){
  tips<- length(tree$tip.label)
  if(node==(tips+1L)){root<-TRUE}else{root<-FALSE}
  return(root)
}

findRoot<- function(tree){
  tips<- length(tree$tip.label)
  return(tips+1L)
}

#' @export
PUT_TIP_col<- function(newtree, oldtree, TIPcol="black", PUTcol="red"){
  col<- vector("character", length(newtree$tip.label))
  col[1:length(col)]<- PUTcol
  col[newtree$tip.label%in%oldtree$tip.label]<-TIPcol
  return(col)

}

randtip_levels<- function(){
  return(as.vector(c("genus","subtribe","tribe","subfamily","family","superfamily","order","class")))
}

sharingtaxa.descs<-function(tree, nodes, MDCC.genera){
  table<- data.frame("node"=as.numeric(nodes), "number"=rep(NA, length(nodes)))
  for(i in seq_along(table$node)){
    node<- table$node[i]
    if(randtip::is.tip(tree,node)){table$number[i]<-1; next}
    descs<- phytools::getDescendants(tree, node,curr = NULL)
    descs<- tree$tip.label[descs]
    descs<- randtip::notNA(descs)
    descs<- descs[randtip::firstword(descs)%in%MDCC.genera]
    table$number[i]<-length(descs)
  }
  table<- table[table$number>0,]
  return(table)
}
