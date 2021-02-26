#' @export
polytomy.into.node <- function(tree, new.tip, node){

  new.tree <- tree 

  to.index <- get.index(new.tree, node = node)
  bind.where <- new.tree$edge[to.index, 2]   
  # Indexing position; node length is 0
  new.tree <- phytools::bind.tip(new.tree, 
                                 new.tip, 
                                 edge.length = NULL, 
                                 where = bind.where, 
                                 position = 0) 
  return(new.tree)
}
