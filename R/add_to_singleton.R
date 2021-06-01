#' Function to add tip to singleton
#' @export
add.to.singleton <- function(tree, singleton, new.tips, resp.sing=F, resp.mono=F){
  singleton<-gsub(" ", "_", singleton)
  singleton<-singleton[singleton%in%tree$tip.label]

  new.tree <- tree

  if(length(singleton)==1){

    node<-which(new.tree$tip.label==singleton)

   if(isTRUE(resp.sing)){
     pos<- randtip::binding.position(new.tree, node = node, insertion = "random",prob = T)
     new.tree <- phytools::bind.tip(new.tree,
                                       new.tips,
                                       edge.length = pos$length,
                                       where = pos$where,
                                       position = pos$position) }
   if(isFALSE(resp.sing)){
     if(!(randtip::isRoot(new.tree, node))){
     parent<- randtip::get.parent.siblings(new.tree, node)[[1]]
     if(isFALSE(resp.mono)){nodes<- phytools::getDescendants(new.tree, parent)}
     if(isTRUE(resp.mono)) {nodes<- randtip::get.permitted.nodes(new.tree, parent)
     nodes<- nodes[nodes!=parent]
     if(length(nodes)==0){nodes<-node}
     }}else{nodes<-node}
      pos<- randtip::binding.position(new.tree, node = sample(nodes,1), insertion = "random",prob = T)
      new.tree <- phytools::bind.tip(new.tree,
                                     new.tips,
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position) }
  }

  if(length(singleton)> 1){
    node<-which(new.tree$tip.label%in%singleton)
    mrca<- ape::getMRCA(new.tree, singleton)

    if(isTRUE(resp.sing)){
      nodes<- c(mrca, node)
      pos<- randtip::binding.position(new.tree, node = sample(nodes,1), insertion = "random",prob = T)
      new.tree <- phytools::bind.tip(new.tree,
                                     new.tips,
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position) }
    if(isFALSE(resp.sing)){
      if(!(randtip::isRoot(new.tree, mrca))){
      parent<- randtip::get.parent.siblings(new.tree, mrca)[[1]]
      if(isFALSE(resp.mono)){nodes<- phytools::getDescendants(new.tree, parent)}
      if(isTRUE(resp.mono)) {nodes<- randtip::get.permitted.nodes(new.tree, parent)
      nodes<- nodes[nodes!=parent]
      if(length(nodes)==0){nodes<-c(node,mrca)}
      }}else{nodes<-c(node,mrca)}

      pos<- randtip::binding.position(new.tree, node = sample(nodes,1), insertion = "random",prob = T)
      new.tree <- phytools::bind.tip(new.tree,
                                     new.tips,
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position) }

  }





  return(new.tree)
}

