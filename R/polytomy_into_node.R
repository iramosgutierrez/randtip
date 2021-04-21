#' @export
polytomy.into.node <- function(tree, new.tip, node){

  new.tree <- tree
  node.descs<- phytools::getDescendants(tree=new.tree, node = node)
  tip.descs<-new.tree$tip.label[node.descs]
  tip.descs<-randtip::notNA(tip.descs)

  for(i in 1: length(new.tip)){

    if (i==1) {
pos<- binding.position(tree=new.tree, node = node, insertion = "polytomy", prob = T)
    # Indexing position; node length is 0
    new.tree <- phytools::bind.tip(new.tree,
                                   new.tip[i],
                                   edge.length = pos$length,
                                   where = pos$where,
                                   position = 0)
    }else{
     sticksp<-  c(tip.descs,new.tip[1:i-1] )
     node<- phytools::findMRCA(tree=new.tree, tips=sticksp)

     pos<- binding.position(tree=new.tree, node = node, insertion = "polytomy", prob = T)

     # Indexing position; node length is 0
     new.tree <- phytools::bind.tip(new.tree,
                                    new.tip[i],
                                    edge.length = pos$length,
                                    where = pos$where,
                                    position = 0)
                                   }

  }


  return(new.tree)
}
