#' Function to add tips over a node. Valid for "sibling" species or genera.
#' @export
add.over.node <- function(tree, new.tip, node){



    if(randtip::isRoot(tree, node)){
      warning("Trying to bind over root node. Addying as polytomy in the root")
      new.tree<-polytomy.into.node(tree=tree, new.tip = new.tip, node=node)
      return(new.tree)

    }
    original.tips<- phytools::getDescendants(tree, node = node)
    original.tips<- randtip::notNA(tree$tip.label[original.tips])

    for(i in seq_along(new.tip)){
      tip<- new.tip[i]
      if(i==1){
    df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2], "length"=tree$edge.length, "id"=1:length(tree$edge.length))


    df <- df[df$node==node,]
    pos<- binding.position(tree = tree, node, insertion = "random")

    new.tree <- phytools::bind.tip(tree, tip, edge.length = pos$length,
                                 where = pos$where, position = pos$position)}else{

       finaltips<- c(original.tips, new.tip[1:i-1])

       original.mrca<- ape::getMRCA(new.tree, original.tips)
       final.mrca <-   ape::getMRCA(new.tree, finaltips)

  permitted.nodes<- c(final.mrca, phytools::getDescendants(new.tree, final.mrca))
  forb.nodes <- phytools::getDescendants(new.tree, original.mrca)
  permitted.nodes<- permitted.nodes[!(permitted.nodes%in%forb.nodes)]

  df <- data.frame("parent"=new.tree$edge[,1], "node"=new.tree$edge[,2], "length"=new.tree$edge.length, "id"=1:length(new.tree$edge.length))

  df <- df[df$node%in%permitted.nodes,]
  pos<- randtip::binding.position(new.tree, df=df, insertion = "random", prob = T)

  new.tree <- phytools::bind.tip(new.tree, tip.label = tip, edge.length = pos$length,
                     where = pos$where, position = pos$position)

                                 }
    }
    return(new.tree)
}

