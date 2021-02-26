get.interm.edges <- function(tree, tips1, tips2, get.length = TRUE){

    if(get.length){
        df <- data.frame(tree$edge, tree$edge.length)
    }
    length.sum <- NULL
    
    node1 <- phytools::findMRCA(tree, tips = tips1, "node")

    if(any(duplicated(tips2))) {
        child <- which(tree$tip.label == tips2[1])
    } else{
        child <- phytools::findMRCA(tree, tips = tips2, "node")
    }

    parent <- phytools::getParent(tree, child) 
    bind.edges <- matrix(c(parent, child), nrow = 1, ncol = 2)
    colnames(bind.edges) <- c("parent", "child")

    if(tree$edge[ape::which.edge(tree, child), 1] == node1) {
        length.sum <- df[ape::which.edge(tree, child), 3] 
    }
    while(parent != node1) {
        desc.par <- phytools::getDescendants(tree, parent, curr = NULL)
        desc.child <- phytools::getDescendants(tree, child, curr = NULL)
        desc.par <- c(desc.par[!(desc.par %in% desc.child) & desc.par != child], 
                      parent)
        new.edges.loc <- (tree$edge[, 1] %in% desc.par) & (tree$edge[, 2] != child)
        bind.edges <- rbind(bind.edges, tree$edge[new.edges.loc, ])
        child <- parent
        parent <- phytools::getParent(tree, parent)
        bind.edges <- rbind(bind.edges, c(parent, child))
    }
    if(get.length & is.null(length.sum)){
        length.sum <- sum(df[df[,1] %in% bind.edges[,1] & 
                             df[,2] %in% bind.edges[,2], 3])
    }
    
    indexes <- vector(mode = "numeric", length = nrow(bind.edges))
    for(j in 1:nrow(bind.edges)){
        indexes[j] <- ape::which.edge(tree, bind.edges[j,2]) 
    }

    return(list(bind.edges = bind.edges,
                indexes = indexes,
                length.sum = length.sum))

}

#' Function to add species at specified branches. 
#' @param edges matrix with 4 columns character vectors. 
#'        Columns 1 and 2 define the parent node of candidate branches and 
#'        columns 3 and 4 define the corresponding child node 
#'        (mrca of species pairs). Note that a pair of nodes may not 
#'        necessarily define one single branch but a set of them
#' @export 
stick.to.branch <- function(tree, edges, new.tip = new.tip, prob = TRUE){

    df <- data.frame(tree$edge, tree$edge.length)
    
    edges <- cbind(edges, 1)
    if(prob){

        relativ_probs <- tree$edge.length
        for(i in 1:nrow(edges)) {
            tips1.char = edges[i,1:2]
            tips2.char = edges[i,3:4]
            interm.edges <- get.interm.edges(tree, tips1.char, tips2.char,
                                             get.length = TRUE)
            edges[i, 5] <- interm.edges$length.sum
        }
    }else{
        relativ_probs <- rep(1, length(tree$edge.length))
    }

    edge.selec <- sample(1:nrow(edges), 1, prob = edges[,5])
    tips1.char <- edges[edge.selec,1:2]
    tips2.char <- edges[edge.selec,3:4]

    interm.edges <- get.interm.edges(tree, tips1.char, tips2.char)
    
    if(length(interm.edges$indexes) == 1){
        to.index <- interm.edges$indexes
    }else{
        to.index <- sample(interm.edges$indexes, 1, 
                           prob = relativ_probs[interm.edges$indexes])      
    }

    bind.where <- tree$edge[to.index, 2]

    tip.pos <- bind.tip.pos(pos.min = 0, 
                            pos.max = tree$edge.length[to.index])
    
    tree <- phytools::bind.tip(tree, 
                               new.tip, 
                               edge.length = NULL, 
                               where = bind.where, 
                               position = tip.pos)

    return(tree = tree)
}


#' Function to help users visualize the candidate branches for tip insertion with
#' \code{\link{stick.to.branch}} function
#' @export
#' @examples
#' set.seed(1)
#' tr <- ape::rcoal(25)
#' edges <- rbind(c("t22", "t14", "t22", "t22"),
#'                c("t7", "t10", "t11", "t3"))
#' candidate.edges <- get.candidate.edges(tr, edges)
#' # Now plot tree highlighting the candidate edges
#' edge_col <- ifelse(1:nrow(tr$edge) %in% candidate.edges, "red", "black")
#' edge_wid <- ifelse(1:nrow(tr$edge) %in% candidate.edges, 2, 1)     
#' plot(tr, edge.color = edge_col, edge.width = edge_wid)
#' # Using the same tree object tr we can select multiple candidate branches
#' edges <- rbind(c("t22", "t14", "t1", "t11"))
#' candidate.edges <- get.candidate.edges(tr, edges)
#' # Now plot tree highlighting the candidate edges
#' edge_col <- ifelse(1:nrow(tr$edge) %in% candidate.edges, "red", "black")
#' edge_wid <- ifelse(1:nrow(tr$edge) %in% candidate.edges, 2, 1)     
#' plot(tr, edge.color = edge_col, edge.width = edge_wid)
get.candidate.edges <- function(tree, edges){
    candidate.edges <- c()
    for(i in 1:nrow(edges)) {            
        tips1.char = edges[i,1:2]
        tips2.char = edges[i,3:4]
        interm.edges <- get.interm.edges(tree, tips1.char, tips2.char,
                                         get.length = TRUE)
        
        candidate.edges <- c(candidate.edges, interm.edges$indexes)
    }
    
    return(unique(candidate.edges))
}
