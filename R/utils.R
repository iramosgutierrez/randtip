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
    position$length<-abs(runif(1, 0, max(tree$edge.length)))}

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

get.permitted.nodes <- function (tree, node, respect.para=F){

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

#extract monophyletic groups and tips
  for(i in seq_along(nodes)){
    node.i <- nodes[i]


    if(randtip::is.tip(tree, node.i)){
      tip<- tree$tip.label[node.i]
      tip.genus<- randtip::firstword(tip)
      siblings<-randtip::get.parent.siblings(tree, node.i)$siblings
      siblings.genera<- unique(randtip::firstword(siblings))
      if(all(siblings.genera %in% tip.genus)){permitted.nodes[i]<-FALSE}
      }else{

        descs<- phytools::getDescendants(tree, node.i,curr = NULL)
        descs<-randtip::notNA (tree$tip.label[descs])
        descs.genera<- unique(randtip::firstword(descs))

        if(length(descs.genera)>1){next}
        siblings<-randtip::get.parent.siblings(tree, node.i)$siblings
        siblings.genera<- unique(randtip::firstword(siblings))

        if(all(siblings.genera %in% descs.genera)){permitted.nodes[i]<-FALSE}
      }


  }


#extract paraphyletic groups and tips
  if(respect.para){
  paraphyletic.genera<-as.vector(which(tip.types.genera=="Paraphyletic") )
  paraphyletic.tips <- which(tree$tip.label%in%paraphyletic.intruders)
  for(g in paraphyletic.genera){
    gen<-names(tip.types.genera[g])
    gentips<- tree$tip.label[randtip::firstword(tree$tip.label)==gen]
    mrca<- ape::getMRCA(tree, gentips)
    mrcadescs<-phytools::getDescendants(tree, mrca,curr = NULL)
    mrcadescs<- mrcadescs[!(mrcadescs%in%paraphyletic.tips)]

    permitted.nodes[which(nodes %in% mrcadescs)]<- FALSE
  }}


  permitted.nodes[is.na(permitted.nodes)]<-TRUE

  if(length(nodes[(permitted.nodes)==TRUE])>0){return(nodes[(permitted.nodes)==TRUE])}else{return(node)}
}

get.forbidden.MDCC.nodes <- function(tree,input, rank, MDCC){
  input<- randtip::correct.DF(input)
  input.mdcc<- input[!is.na(input[,rank]),]
  input.mdcc<- input.mdcc[input.mdcc[,rank]==MDCC,]
  forbidden.nodes<- as.numeric(NULL)
  if(rank=="genus"){return(NULL)}
  for(i in 2:which(randtip::randtip_ranks()==rank)){
    lv<-randtip::randtip_ranks()[i]
    col.lev<- as.character(lv)
    MDCCs<-unique(input.mdcc[,col.lev])
    MDCCs<-randtip::notNA(MDCCs)
    for(lev in MDCCs){
      phylstat<-randtip::MDCC.phyleticity(input, tree = tree,
                                          MDCC.info = list(rank=col.lev, MDCC=lev), trim=F )
      if(phylstat%in%c("Monophyletic","Paraphyletic")){
        mdcc.gen<-input.mdcc[input.mdcc[,col.lev]==lev,"taxon"]
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

get.intruder.nodes <- function (tree, input, rank, MDCC){
  input<- randtip::correct.DF(input)
  input.mdcc<- input[!is.na(input[,rank]),]
  input.mdcc<- input.mdcc[input.mdcc[,rank]==MDCC,]

  mdcc.genera<- randtip::firstword(input.mdcc$taxon)
  mdcc.intree<- randtip::sp.genus.in.tree(tree, mdcc.genera)
  mdcc.mrca<- ape::getMRCA(tree, mdcc.intree)

  descendants<- phytools::getDescendants(tree, node=mdcc.mrca,curr = NULL)
  descendants<- tree$tip.label[descendants]
  descendants<- randtip::notNA(descendants)

  non.mdcc <- input[!is.na(input[,rank]),]
  non.mdcc <- non.mdcc[non.mdcc[,rank]!=MDCC,]
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

bind.clump<- function(newtree, tree, input, new.species){

  sp<- paste0(stringr::word(new.species, 1,sep = "_"), "_", stringr::word(new.species, 2,sep = "_"))
  treespp<- paste0(stringr::word(tree$tip.label, 1,sep = "_"), "_", stringr::word(tree$tip.label, 2,sep = "_"))
  newtreespp<- paste0(stringr::word(newtree$tip.label, 1,sep = "_"), "_", stringr::word(newtree$tip.label, 2,sep = "_"))

  if(sp%in%treespp){
    DFspp<- paste0(stringr::word(input$taxon, 1,sep = "_"), "_", stringr::word(input$taxon, 2,sep = "_"))
    clump<- tree$tip.label[treespp==sp]
    clumpDF<-input[input$clump.PUTs==TRUE,]
    clumpDF<-input$taxon[DFspp==sp]
    clump<- unique(c(clump, clumpDF))
    clump<- clump[clump%in%newtree$tip.label]
    if(length(clump)>1){
      mrca<- ape::getMRCA(newtree, clump)
      descs<- newtree$tip.label[phytools::getDescendants(newtree, mrca, curr=F)]
      descs<-randtip::notNA(descs)
      if(any(!(descs%in%clump))){clump<- sample(clump, 1)}
    }
    return(clump)
  }




  using.MDCC<-input[input$taxon==new.species,"using.MDCC"]
  using.MDCC.lev<- input[input$taxon==new.species,"using.MDCC.lev"]

  spp<-input$taxon[input[,using.MDCC.lev]==using.MDCC]
  spp.in.newtree<- spp[spp%in%newtree$tip.label ]
  spp.in.tree<- spp[spp  %in%  tree$tip.label ]

  if(length(spp.in.tree)==0 & length(spp.in.newtree)>0){return(spp.in.newtree)}

}

usingMDCCfinder<- function(input, taxon=NULL, tree, verbose=F){

  if(is.null(taxon)){taxon<- input$taxon}
  input<-randtip::correct.DF(input)
  MDCC.vect<- vector( mode="character", length = length(taxon))
  MDCC.lev.vect<- vector(mode="character", length = length(taxon))
  MDCC.phyletictype.vect<- vector(mode="character", length = length(taxon))

  if(verbose){
    cat(paste0("Searching MDCCs... ", "\n"))
  }

  #manual MDCC search
  taxa<- input[!is.na(input$taxon1)|!is.na(input$taxon2),]
  if(nrow(taxa)>0){
    for(tx in seq_along(taxa$taxon)){
      if(isTRUE(length(strsplit( taxa$taxon1[tx], "_")[[1]])==1 &
                length(strsplit( taxa$taxon2[tx], "_")[[1]])==1 &
                taxa$taxon1[tx]==taxa$taxon2[tx] &
                length(randtip::sp.genus.in.tree(tree, taxa$taxon1[tx]))>0)){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- taxa$taxon1[tx]
        MDCC.lev.vect[pos] <- "Sister genus"
        MDCC.phyletictype.vect[pos]<-"-"
        next
      }

      if(!(taxa$taxon1[tx]%in%tree$tip.label)){taxa$taxon1[tx]<-NA}
      if(!(taxa$taxon2[tx]%in%tree$tip.label)){taxa$taxon2[tx]<-NA}

      if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
                taxa$taxon2[tx] %in% tree$tip.label &
                taxa$taxon1[tx]==taxa$taxon2[tx])){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- taxa$taxon1[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        MDCC.phyletictype.vect[pos]<-"-"
        next
      }

      if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
                taxa$taxon2[tx] %in% tree$tip.label &
                taxa$taxon1[tx]!=taxa$taxon2[tx])){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- paste0("Clade (", taxa$taxon1[tx], "-", taxa$taxon2[tx], ")")
        MDCC.lev.vect[pos] <- "Manual clade"
        MDCC.phyletictype.vect[pos]<-"-"
        next
      }

      if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
                is.na(taxa$taxon2[tx]))){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- taxa$taxon1[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        MDCC.phyletictype.vect[pos]<-"-"
        next
      }

      if(isTRUE(taxa$taxon2[tx] %in% tree$tip.label &
                is.na(taxa$taxon1[tx]))){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- taxa$taxon2[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        MDCC.phyletictype.vect[pos]<-"-"
        next
      }

    }
  }

  #automatic MDCC search
  ranks<- randtip::randtip_ranks()
  taxa<- input[!(!is.na(input$taxon1)|!is.na(input$taxon2)),]

  if(nrow(taxa)>0){
    vect<- which(taxon%in%taxa$taxon)
    for(v in vect){

    if(verbose & (v %in% c(seq(0, length(taxon), 10), length(taxon)))){
        cat(paste0(round((v/length(taxon)*100),2), "% completed.\n"))}




    if(taxon[v]%in%tree$tip.label){
      MDCC.vect[v]<- "Tip"
      MDCC.lev.vect[v]<-"Tip"
      MDCC.phyletictype.vect[v]<-"Tip"
      next
      }
    i<- which(input$taxon==taxon[v])
    if((MDCC.vect[v])==""){

      MDCC<-as.character(NA)
      MDCC.ranks<-as.character(NA)

      for(rank in ranks){
        if(is.na(MDCC)){
          MDCC<-as.character(input[i, rank])
          if(!is.na(MDCC)){phyleticity<-randtip::MDCC.phyleticity(input, tree = tree,
                                         MDCC.info = list(rank=rank, MDCC= MDCC))
          if(phyleticity=="Missing"){MDCC<-NA}
          }

          lev<-rank
        }else{next}
      }
      MDCC.vect[v]<-as.character(MDCC)
      MDCC.lev.vect[v]<-as.character(lev)
      if(is.na(MDCC)){MDCC.lev.vect[v]<-NA}



    if(is.na(MDCC.vect[v])|is.na(MDCC.lev.vect[v])){MDCC.phyletictype.vect[v]<-NA}else{

      if(MDCC%in%taxa$using.MDCC[-which(taxa$taxon==taxon[v])]){
        ps<-taxa$using.MDCC.phylstat[taxa$using.MDCC==MDCC][1]
        MDCC.phyletictype.vect[v]<-MDCC.phyletictype.vect[ps]
      }else{
      MDCC.phyletictype.vect[v]<-randtip::MDCC.phyleticity(input = input, tree =  tree,
               MDCC.info = list(rank=MDCC.lev.vect[v], MDCC=MDCC.vect[v]), trim = T)}
    }
  }}
}

  return(list(MDCC=MDCC.vect,MDCC.ranks=MDCC.lev.vect, MDCC.phylstat=MDCC.phyletictype.vect) )
}

inputfinder<- function(input, taxon, column){
  return(as.character(input[input$taxon==taxon, column]))
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

  for(col in names(DF)){
    DF[,col]<- as.character(DF[,col])
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

randtip_ranks<- function(){
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

combineDF<- function(largeDF, smallDF){
  largeDF<- randtip::correct.DF(largeDF)
  smallDF<- randtip::correct.DF(smallDF)
  largeDF.names<-names(largeDF)
  smallDF.names<-names(smallDF)
  smallDF<- smallDF[,smallDF.names[smallDF.names%in%largeDF.names]]
  notincl.names<- largeDF.names[!(largeDF.names%in%smallDF.names)]
  DF.out<- smallDF
  for(col in notincl.names){
    DF.out[,col]<-NA
  }
  DF.out<-DF.out[,names(largeDF)]
  DF.out<- DF.out[!(DF.out$taxon%in%largeDF$taxon),]
  DF.out<- rbind(largeDF, DF.out)
  return(DF.out)
}

#' this function reveals if a given node at a given tree is an internal node
#' @export
is.node<-function(tree, node){
  if(!(node %in% tree$edge)){
    stop("Node number is not in your tree")
  }
  if(length(phytools::getDescendants(tree = tree, node = node )) > 1){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#' this function reveals if a given node at a given tree is a tip
#' @export
is.tip <-function(tree, node){
  if(!(node %in% tree$edge)){
    stop("Node number is not in your tree")
  }
  if(length(phytools::getDescendants(tree = tree, node = node)) == 1){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
