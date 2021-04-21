#' Function to add tips over a node. Valid for "sibling" species or genera.
#' @export
add.over.node <- function(tree, new.tip, node){
    df <- data.frame(tree$edge, tree$edge.length, 1:length(tree$edge.length))

    to.index <- get.index(tree, node = node)
    df <- df[to.index,]
    pos<- binding.position(tree = tree, df = df, insertion = "random" , prob=T)

    new.tree <- phytools::bind.tip(tree, new.tip, edge.length = pos$length,
                                 where = pos$where, position = pos$position)
    return(new.tree)
}

