# Return a vector of all species in the tree that match a given genus
sp.genus.in.tree <- function(tree, genus){
    sp <- tree$tip.label
    taxa.vector <- sp[randtip::firstword(sp)%in%genus]

    return(taxa.vector)
}

# Convenience function to give bind position for tip
# The returning value is to be used in the position argument
# of phytools bind.tip function.
bind.tip.pos <- function(pos.min, pos.max){
     pos.tip <- stats::runif(1, pos.min, pos.max)
     while(pos.tip == pos.min | pos.tip == pos.max){
        pos.tip <- stats::runif(1, pos.min, pos.max)
     }
     return(pos.tip)
}

# Return index of edge where binding should take place
get.index <- function(tree, how = "sample_simple", node = NULL, df = NULL){

    if(!is.null(node)){if(randtip::isRoot(tree, node)){return(NA)}}

    if(is.null(df)){
        df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                         "length"=tree$edge.length, "id"=1:length(tree$edge[,1]))

    }

    if(!is.null(node)){
        edges <- df[ape::which.edge(tree, node),]
    }else if(how == "sample_simple"){
        edges <- df[sample(x = 1:nrow(df), size = 1),]
    }else if(how == "sample_prob"){
        edges <- df[sample(x = 1:nrow(df), size = 1, prob = abs(df$length)),]
    }else{
        stop("Unrecognized set of arguments to get.index function.")
    }

    to.index <- edges$id

    return(to.index)
}

get.position <- function(tree, node, insertion){
if(randtip::isRoot(tree, node)){position=0
}else{
    edge.length <- tree$edge.length[ape::which.edge(tree, node)]
    if(insertion=="random"){
        # Bind at a random point of the branch
        position <- edge.length * stats::runif(1, 0, 1)
    }else if(insertion=="middle"){
        position <- edge.length*0.5
        # Bind at the middle of the branch
    }else if(insertion=="long"){
        # Bind at the beggining of the branch
        position <- edge.length
    }else{
        stop("Unknown specification of the insertion argument.")
    }}

    return(position)
}


binding.position<- function(tree, node=NULL, df=NULL, insertion,  prob){
  position<-list("length"=NA, "where"=NA, "position"=NA)

  if(is.null(df)){df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2], "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )}

  if(ape::is.ultrametric(tree)){position$length<-NULL}else{
    position$length<-abs(rnorm(1, mean=mean(tree$edge.length), sd= sd(tree$edge.length) ))}

  if(!is.null(node)){
    df<- df[df$node==node,]
    position$where<- node}

  if(insertion=="polytomy"){
    position$position<- 0
    position$where <- node
  }

  if(insertion=="random"){
    if(prob){how <- "sample_prob"}else{how <- "sample_simple"}
    index<- randtip::get.index(tree = tree, df=df , how = how)
    pos <- randtip::get.position(tree = tree, node = df[df$id==index,"node"], insertion = "random")

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
    parent.desc <- phytools::getDescendants(tree, parent)
    siblings <- tree.sp[parent.desc]
    siblings <- siblings[!is.na(tree.sp[parent.desc])]

    return(list(parent = parent,
                siblings = siblings))
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
  mrca.descs <- phytools::getDescendants(tree, sp.mrca)

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

    nd.descs<- phytools::getDescendants(tree, nd)
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

        intruders.descs<- phytools::getDescendants(tree, intruders.mrca)
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
          group.descs<- phytools::getDescendants(tree, group.mrca)
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
  if(node==(tips+1)){root<-TRUE}else{root<-FALSE}
  return(root)
}

findRoot<- function(tree){
  tips<- length(tree$tip.label)
  return(root+1)
}

#' @export
PUT_TIP_col<- function(newtree, oldtree, TIPcol="black", PUTcol="red"){
  col<- vector("character", length(newtree$tip.label))
  col[1:length(col)]<- PUTcol
  col[randtree$tip.label%in%oldtree$tip.label]<-TIPcol
  return(col)

}
