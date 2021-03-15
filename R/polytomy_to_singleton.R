#' @export
polytomy.to.singleton <- function(tree, singleton, new.tip,
                                  insertion = "random" ){ #insertion can be "random", "long" or "middle"

  singleton<-gsub(pattern = " ", "_", singleton)

  new.tree <- tree
  node<-which(new.tree$tip.label==singleton)

  for(i in 1: length(new.tip)){

    if (i==1) {
      to.index <- get.index(new.tree, node = node)
      bind.where <- new.tree$edge[to.index, 2]
      position<- get.position(tree=new.tree, node = bind.where, insertion = insertion)
      new.tree <- phytools::bind.tip(new.tree,
                                     new.tip[i],
                                     edge.length = NULL,
                                     where = bind.where,
                                     position = position)
    }else{
      sticksp<-  c(singleton,new.tip[1:i-1] )
      node<- phytools::findMRCA(tree=new.tree, tips=sticksp)

      to.index <- get.index(new.tree, node = node)
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
