#' @export
polytomy.to.singleton <- function(tree, singleton, new.tip,
                                  insertion = "random" ){ #insertion can be "random", "long" or "middle"

  singleton<-gsub(pattern = " ", "_", singleton)

  new.tree <- tree
  node<-which(new.tree$tip.label==singleton)

  for(i in 1: length(new.tip)){

    if (i==1) {

      pos<- randtip::binding.position(new.tree, node = node, insertion = "random",prob = T)



      new.tree <- phytools::bind.tip(new.tree,
                                     new.tip[i],
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position)
    }else{
      sticksp<-  c(singleton,new.tip[1:i-1] )
      node<- phytools::findMRCA(tree=new.tree, tips=sticksp)
      pos<- randtip::binding.position(new.tree, node = node, insertion = "polytomy",prob = T)



      # Indexing position; node length is 0
      new.tree <- phytools::bind.tip(new.tree,
                                     new.tip[i],
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position)
    }

  }


  return(new.tree)
}
