#' @export
polytomy.to.singleton <- function(tree, singleton, species, 
                                  insertion = c("random", "long", "middle") ){

  new.tree <- tree 
  for(i in 1:length(species)){
    if(i == 1){ 
      sing.node <- which(new.tree$tip.label == singleton)
      
      to.index<- get.index(new.tree, node = sing.node)
      bind.where <- new.tree$edge[to.index, 2]
      position <- get.position(tree, sing.node, insertion)

      new.tree <- phytools::bind.tip(new.tree, 
                                     tip.label =  species[i], 
                                     edge.length = NULL, 
                                     where = bind.where, 
                                     position = position)
    }else{      
      node <- phytools::findMRCA(new.tree, 
                                tips = c(singleton, species[1:(i-1)]))
      new.tree<- polytomy.into.node(tree = new.tree, 
                                    new.tip = species[i],                                     
                                    node = node)
    }
  }

  return(new.tree)
}
