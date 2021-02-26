#' Function to add tips over a node. Valid for "sibling" species or genera.
#' @export
add.over.node <- function(tree, new.tip, node){
   
    to.index <- get.index(tree, node = node)
    bind.where <- tree$edge[to.index,2]   

    edge.length <- tree$edge.length 
    tip.pos <- bind.tip.pos(pos.min = 0, 
                          pos.max = tree$edge.length[to.index])
    new.tree <- phytools::bind.tip(tree, new.tip, edge.length = NULL, 
                                 where = bind.where, position = tip.pos)
    return(new.tree)
}

