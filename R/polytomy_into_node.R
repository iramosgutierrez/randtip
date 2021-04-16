#' @export
polytomy.into.node <- function(tree, new.tip, node){

  new.tree <- tree
  node.descs<- phytools::getDescendants(tree=new.tree, node = node)
  tip.descs<-new.tree$tip.label[node.descs]
  tip.descs<-randtip::notNA(tip.descs)

  for(i in 1: length(new.tip)){

    if (i==1) {
    to.index <- randtip::get.index(new.tree, node = node)
    if(is.na(to.index)){
    bind.where <- new.tree$edge[1, 1]}else{
    bind.where <- new.tree$edge[to.index, 2]}
    # Indexing position; node length is 0
    new.tree <- phytools::bind.tip(new.tree,
                                   new.tip[i],
                                   edge.length = NULL,
                                   where = bind.where,
                                   position = 0)
    }else{
     sticksp<-  c(tip.descs,new.tip[1:i-1] )
     node<- phytools::findMRCA(tree=new.tree, tips=sticksp)

     to.index <- randtip::get.index(new.tree, node = node)
     bind.where <- new.tree$edge[to.index, 2]
     # Indexing position; node length is 0
     new.tree <- phytools::bind.tip(new.tree,
                                    new.tip[i],
                                    edge.length = NULL,
                                    where = bind.where,
                                    position = 0)
                                   }

  }


  return(new.tree)
}
