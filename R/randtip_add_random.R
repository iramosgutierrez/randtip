#' modification from phytools "add.random". if prob=TRUE, it remains the same; 
#' else a random node is selected without probability vector
#' 
#' @param tips Character vector with the labels of the tip to add to the 
#'             phylogenetic tree.
#' @param order The addition order the tips. If \dQuote{random}, tips are added
#'              at random. If \dQuote{input}, tips are added sequentially. 
#' @export
add.random2 <- function (tree, n = NULL, tips = NULL, 
                         edge.length = NULL, order = "random", 
                         prob = TRUE){
    if (!inherits(tree, "phylo")){
        stop("tree should be an object of class \"phylo\".")
    }
  
    if(prob == TRUE){
      args <- mget(ls())      
      args$prob <- NULL        
      tree <- do.call(phytools::add.random, args) 
      return(tree)
    }
    
    if(is.null(tips)){
        if(is.null(n)) n <- 1
        tips <- paste0("t", length(tree$tip) + 1:n)
    }else{
        n <- length(tips)  
    } 
    
    if(order != "input"){
        tips <- sample(tips, length(tips))
    }
    
    for(i in 1:n){
        new.tip <- tips[i]
        
        # Position in tree
        to.index <- get.index(tree)
        bind.where <- tree$edge[to.index, 2]   
        
        edge.length <- tree$edge.length
        pos.tip <- bind.tip.pos(pos.min = 0, 
                                pos.max =  edge.length[to.index])
        
        tree <- phytools::bind.tip(tree, new.tip, 
                                   edge.length = NULL, 
                                   where = bind.where, 
                                   position = pos.tip)
    }
    
    return(tree)
}
