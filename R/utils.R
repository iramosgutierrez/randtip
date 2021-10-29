#GENERAL FUNCTIONS####

firstword<- function(string, sep="_"){
  word<- stringr::word(string, 1, sep=sep)
  return(word)
}

notNA <- function(x){
  vect<- x[!is.na(x)]
  return(vect)
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

randtip_ranks<- function(){
  return(as.vector(c("genus","subtribe","tribe","subfamily","family","superfamily","order","class")))
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

#' @export
put.tip.col<- function(newtree, oldtree, placed.col="#C23B23", put.col="#3d3d3d"){
  col<- vector("character", length(newtree$tip.label))
  col[1:length(col)]<- put.col
  col[newtree$tip.label%in%oldtree$tip.label]<-placed.col
  return(col)

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


#SPECIFIC FUNCTIONS####
# Return a vector of all species in the tree that match a given genus
sp.genus.in.tree <- function(tree, genus){
    sp <- tree$tip.label
    taxa.vector <- sp[randtip::firstword(sp)%in%genus]

    return(taxa.vector)
}

inputfinder<- function(input, taxon, column){
  return(as.character(input[input$taxon==taxon, column]))
}

usingMDCCfinder<- function(input, taxon=NULL, tree){

  if(is.null(taxon)){taxon<- input$taxon}
  input<-randtip::correct.DF(input)
  MDCC.vect<- vector( mode="character", length = length(taxon))
  MDCC.lev.vect<- vector(mode="character", length = length(taxon))
  #MDCC.phyletictype.vect<- vector(mode="character", length = length(taxon))


    cat(paste0("Searching MDCCs... ", "\n"))


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
        #MDCC.phyletictype.vect[pos]<-"-"
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
        #MDCC.phyletictype.vect[pos]<-"-"
        next
      }

      if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
                taxa$taxon2[tx] %in% tree$tip.label &
                taxa$taxon1[tx]!=taxa$taxon2[tx])){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- paste0("Clade (", taxa$taxon1[tx], "-", taxa$taxon2[tx], ")")
        MDCC.lev.vect[pos] <- "Manual clade"
        #MDCC.phyletictype.vect[pos]<-"-"
        next
      }

      if(isTRUE(taxa$taxon1[tx] %in% tree$tip.label &
                is.na(taxa$taxon2[tx]))){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- taxa$taxon1[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        #MDCC.phyletictype.vect[pos]<-"-"
        next
      }

      if(isTRUE(taxa$taxon2[tx] %in% tree$tip.label &
                is.na(taxa$taxon1[tx]))){

        pos<- which(taxon==taxa$taxon[tx])
        MDCC.vect[pos] <- taxa$taxon2[tx]
        MDCC.lev.vect[pos] <- "Sister species"
        #MDCC.phyletictype.vect[pos]<-"-"
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

      if(v==vect[1]){
       cat(paste0("0%       25%       50%       75%       100%", "\n",
                  "|---------|---------|---------|---------|",   "\n","*"))}




        vec<- seq(from=0, to=40, by=40/length(vect))
        vec<-ceiling(vec)
        vec<- diff(vec)
        cat(strrep("*", times=vec[which(vect==v)]))

       if(v ==vect[length(vect)]){cat("\n")}




    if(taxon[v]%in%tree$tip.label){
      MDCC.vect[v]<- "Tip"
      MDCC.lev.vect[v]<-"Tip"
      #MDCC.phyletictype.vect[v]<-"Tip"
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



    #if(is.na(MDCC.vect[v])|is.na(MDCC.lev.vect[v])){MDCC.phyletictype.vect[v]<-NA}else{
    #
    # if(MDCC%in%taxa$using.MDCC[-which(taxa$taxon==taxon[v])]){
    #   ps<-taxa$using.MDCC.phylstat[taxa$using.MDCC==MDCC][1]
    #   MDCC.phyletictype.vect[v]<-MDCC.phyletictype.vect[ps]
    # }else{
    # MDCC.phyletictype.vect[v]<-randtip::MDCC.phyleticity(input = input, tree =  tree,
    #          MDCC.info = list(rank=MDCC.lev.vect[v], MDCC=MDCC.vect[v]), trim = T)}
    #}
  }}
}

  return(list(MDCC=MDCC.vect,MDCC.ranks=MDCC.lev.vect) )#, MDCC.phylstat=MDCC.phyletictype.vect
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

sharingtaxa.descs<-function(tree, nodes, MDCC.genera){
  table<- data.frame("node"=as.numeric(nodes), "number"=rep(NA, length(nodes)), "tot.number"=rep(NA, length(nodes)))
  for(i in seq_along(table$node)){
    node<- table$node[i]
    if(randtip::is.tip(tree,node)){table$number[i]<-1;table$tot.number[i]<-1; next}
    descs<- phytools::getDescendants(tree, node,curr = NULL)
    table$tot.number[i]<-length(descs)
    descs<- tree$tip.label[descs]
    descs<- randtip::notNA(descs)
    descs<- descs[randtip::firstword(descs)%in%MDCC.genera]
    table$number[i]<-length(descs)
  }
  table<- table[table$number>0,]
  return(table)
}

get.permitted.nodes <- function (tree, input, MDCC, rank, MDCC.type,
                                 polyphyly.scheme, use.paraphyletic, use.singleton){

  input.mdcc<-  input[!is.na(input[,rank]),]
  MDCC.taxa<- input.mdcc$taxon[input.mdcc[,rank]==MDCC]
  MDCC.genera<- randtip::notNA(unique(randtip::firstword(MDCC.taxa)))
  MDCC.intree<- sp.genus.in.tree(tree, MDCC.genera)

  if(MDCC.type=="Monophyletic"){
    MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
    nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
    return(nodes)
  }

  if(MDCC.type=="Paraphyletic"){
    if(!use.paraphyletic){
      MDCC.mrca<- ape::getMRCA(new.tree, MDCC.intree)
      nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
      return(nodes)
    }
    if(use.paraphyletic){
      MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
      nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
      int.nodes<- randtip::get.intruder.nodes(tree, input, rank, MDCC)
      nodes <- nodes[!(nodes%in%int.nodes)]
      return(nodes)
    }

  }

  if(MDCC.type%in% c("Singleton", "Singleton MDCC")){
    if (use.singleton){
      nodes<- which(tree$tip.label==MDCC.intree)
    }
    if(!use.singleton){
      tip<- which(tree$tip.label==MDCC.intree)
      nodes<- phytools::getDescendants(tree, get.parent.siblings(tree, tip)$parent, curr=NULL)
      return(nodes)
    }

  }

  if(MDCC.type=="Polyphyletic"){
    if(polyphyly.scheme=="complete"){
      MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
      nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
      return(nodes)
    }

    if(polyphyly.scheme == "frequentist"){
      MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
      nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
      table<- data.frame("node"=as.numeric(nodes), "descs"=NA, "total.descs"=NA,
                         "sharing.descs"=NA, "eligible"=NA)
      for(i in seq_along(table$node)){
        node<- table$node[i]
        if(randtip::is.tip(tree,node)){
          table$descs[i]   <- NA
          table$total.descs[i]<-1
          if ( tree$tip.label[node]%in% input[input[,rank]==MDCC,"taxon"]){
            table$sharing.descs[i]<-1;table$eligible[i]<- "TRUE"}else{
              table$sharing.descs[i]<-0;table$eligible[i]<- "FALSE"}



        }
        if(randtip::is.node(tree,node)){
          descs<- phytools::getDescendants(tree, node, curr=NULL)
          table$descs[i] <- paste0(descs, collapse = ",")
          descs<-randtip::notNA(tree$tip.label[descs])
          table$total.descs[i]   <- length(descs)
          table$sharing.descs[i] <- length(descs[descs %in% input[input[,rank]==MDCC,"taxon"]])
          if(length(descs[descs %in% input[input[,rank]!=MDCC,"taxon"]])>0){
            table$eligible[i]<-"FALSE"}else{table$eligible[i]<-"TRUE"}
        }
      }
      table<- table[table$eligible=="TRUE",]
      for(i in seq_along(table$node)){
        node<- table$node[i]
        descs<- unique(unlist(strsplit(table$descs, split=",")))
        if(node%in%descs ){table$eligible[i]<- "FALSE"}
      }
      table<- table[table$eligible=="TRUE",]
      if(nrow(table)==1){nodes=table$node}else{
        nodes<- sample(table$node, size = 1,prob=table$total.descs)}
      nodes <- phytools::getDescendants(tree, nodes)




    }

    if(polyphyly.scheme == "largest"){
      MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
      nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
      table<- data.frame("node"=as.numeric(nodes), "descs"=NA, "total.descs"=NA,
                         "sharing.descs"=NA, "eligible"=NA)
      for(i in seq_along(table$node)){
        node<- table$node[i]
        if(randtip::is.tip(tree,node)){
          table$descs[i]   <- NA
          table$total.descs[i]<-1
          if ( tree$tip.label[node]%in% input[input[,rank]==MDCC,"taxon"]){
            table$sharing.descs[i]<-1;table$eligible[i]<- "TRUE"}else{
              table$sharing.descs[i]<-0;table$eligible[i]<- "FALSE"}



        }
        if(randtip::is.node(tree,node)){
          descs<- phytools::getDescendants(tree, node, curr=NULL)
          table$descs[i] <- paste0(descs, collapse = ",")
          descs<-randtip::notNA(tree$tip.label[descs])
          table$total.descs[i]   <- length(descs)
          table$sharing.descs[i] <- length(descs[descs %in% input[input[,rank]==MDCC,"taxon"]])
          if(length(descs[descs %in% input[input[,rank]!=MDCC,"taxon"]])>0){
            table$eligible[i]<-"FALSE"}else{table$eligible[i]<-"TRUE"}
        }
      }
      table<- table[table$eligible=="TRUE",]
      for(i in seq_along(table$node)){
        node<- table$node[i]
        descs<- unique(unlist(strsplit(table$descs, split=",")))
        if(node%in%descs ){table$eligible[i]<- "FALSE"}
      }
      table<- table[table$eligible=="TRUE",]
      table<-table[table$sharing.descs== max(table$sharing.descs),]
      if(nrow(table)>1){table<-table[sample(1:nrow(table), size=1),]}

      nodes<- table$node
      nodes <- phytools::getDescendants(tree, nodes)




    }
  }
  return(nodes)

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


#WORKING SPACE####



perm.nodes<- get.permitted.nodes(tree, input, MDCC, rank, MDCC.type, polyphyly.scheme, use.paraphyletic , use.singleton )

get.forbidden.nodes <- function(tree,input, MDCC, rank, perm.nodes, respect.mono, respect.para){


if(rank=="genus"){return(NULL)}

forbidden.nodes<- vector("numeric")
if(respect.mono){
  for(i in seq_along(perm.nodes)){

    nd<- perm.nodes[i]
    if(nd %in% forbidden.nodes){next}
    if(randtip::is.tip(tree, nd)){next}
    descs.nd<- phytools::getDescendants(tree, nd)
    descs <- randtip::notNA(tree$tip.label[descs.nd])
    genera<- unique(randtip::firstword(descs))
    if(length(genera)==1){forbidden.nodes<- c(forbidden.nodes, descs.nd ); next}

    sub.input <- input[firstword(input$taxon)%in%genera,]
    #if(rank!="genus"){
    for( rk in randtip::randtip_ranks()[2:(which(randtip::randtip_ranks()==rank)-1)]){
      rk.vals<-randtip::notNA(unique(sub.input[,rk]))
      if(length(rk.vals)==1){
        forbidden.nodes<- c(forbidden.nodes, descs.nd ); next}

    #}

}
  }}

if(respect.para){
  for(i in seq_along(perm.nodes)){

    nd<- perm.nodes[i]
    if(nd %in% forbidden.nodes){next}
    if(randtip::is.tip(tree, nd)){next}
    descs.nd<- phytools::getDescendants(tree, nd)
    descs <- randtip::notNA(tree$tip.label[descs.nd])
    genera<- unique(randtip::firstword(descs))


    #account for singleton tips in otherwise monophyletic clusters
    if(length(genera)<3 & any(as.vector(table(firstword(descs)))==1)){
      uniq.gen<-names(table(firstword(descs)))[as.vector(table(firstword(descs)))==1]
      tip<- descs[randtip::firstword(descs)==uniq.gen]
      tip.n<- which(tree$tip.label==tip)
      tip.par<- randtip::get.parent.siblings(tree, tip.n)$parent

      if(tip.par!=nd){forbidden.nodes<- c(forbidden.nodes, descs.nd[which(descs.nd!=tip.n)] ); next}}

    rk.vals.mrca<- vector("numeric", length = length(genera))
    rk.vals.desc<- list( NA)
    for(v in genera){

      ds <- descs[randtip::firstword(descs)%in%v]
      if(length(ds)==1){
        rk.vals.mrca[which(genera==v)]<-which(tree$tip.label==ds)
        rk.vals.desc[[which(genera==v)]]<- NA
        next}

      ds.mrca<- ape::getMRCA(tree, ds)
      rk.vals.mrca[which(genera==v)]<-ds.mrca
      rk.vals.desc[[which(genera==v)]]<- phytools::getDescendants(tree, ds.mrca, curr=NULL)
    }

    nest<- genera[which(rk.vals.mrca==nd)]
    nested<- descs[randtip::firstword(descs)!=nest]
    if(!any(firstword(nested)==nest)){
      para.nodes<- rk.vals.desc[[which(rk.vals.mrca==nd)]]
      intr.nodes <- rk.vals.desc
      intr.nodes[[which(rk.vals.mrca==nd)]]<-NA
      intr.nodes<- randtip::notNA(unlist(intr.nodes))
      para.nodes<- para.nodes[!(para.nodes%in%intr.nodes)]
      forbidden.nodes<- c(forbidden.nodes, para.nodes ); next}

#VOY POR AQUI, HASTA GENERO YA HECHO, FALTA ADECUAR LO SUPERIOR
      sub.input <- input[firstword(input$taxon)%in%genera,]

      for( rk in randtip::randtip_ranks()[2:(which(randtip::randtip_ranks()==rank)-1)]){
        rk.vals<-randtip::notNA(unique(sub.input[,rk]))
         if(length(rk.vals)>1){

           rk.vals.mrca<- vector("numeric", length = length(rk.vals))
           rk.vals.desc<- list( NA)
           for(v in rk.vals){
             gen <- unique(randtip::firstword(sub.input$taxon[sub.input[,rk]==v]))
             ds <- descs[randtip::firstword(descs)%in%gen]
             ds.mrca<- ape::getMRCA(tree, ds)
             rk.vals.mrca[which(rk.vals==v)]<-ds.mrca
             rk.vals.desc[[which(rk.vals==v)]]<- phytools::getDescendants(tree, ds.mrca, curr=NULL)
           }
           if(!any(rk.vals.mrca==nd)){next}


          forbidden.nodes<- c(forbidden.nodes, nd,descs.nd ); next}

      }
  }

  }
}


  return(unique(forbidden.nodes))
}
#nodes<- phytools::getDescendants(tree, node, curr=NULL)
#sub.info<- input[input[,rank]==MDCC,]
#sub.info<-sub.info[randtip::randtip_ranks()[1:which(randtip::randtip_ranks()==rank)-1]]
#
#for (i in 1:length(nodes)){
#  nd <- nodes[i]
#  if(randtip::is.tip(tree, nd)){next}
#  nd.descs<- randtip::notNA(tree$tip.label[phytools::getDescendants(tree, nd, curr=NULL)])
#subranks<- as.vector(NULL, "character")
#for(col in 1:ncol(sub.info)){
#  subranks<- c(subranks, sub.info[,col])
#}
#subranks<- randtip::notNA(subranks)
