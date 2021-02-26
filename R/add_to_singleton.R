#' Function to add tip to singleton 
#' @export
add.to.singleton <- function(tree, singleton, new.tips){
    new.tree <- tree

    singleton.match <- new.tree$tip.label == singleton
    if(!any(singleton.match)){
        stop("Singleton tip not found")
    }    

    to.add <- gsub(" ", "_", new.tips)
    added <- vector("character", length(new.tips) + 1)
    added[1] <- singleton
    for(i in 1:length(to.add)){
        tips.labs <- new.tree$tip.label 
        edge.length <- new.tree$edge.length
        # Add with random branch length
        if(i == 1){    
            new.tip <- to.add[i]
            
            to.index <- get.index(new.tree, 
                                  node = which(tips.labs == added[1]))
            bind.where <- new.tree$edge[to.index,2]
            tip.pos <- bind.tip.pos(pos.min = 0,
                                  pos.max = edge.length[to.index])

            new.tree <- phytools::bind.tip(new.tree, 
                                           new.tip, 
                                           edge.length = NULL, 
                                           where = bind.where, 
                                           position = tip.pos)
            added[i+1] <- new.tip
            
        }else{
            new.tip <- to.add[i]
            species.loc <- match(added, new.tree$tip.label)
            species.loc <- species.loc[!is.na(species.loc)]
            
            df <- data.frame(new.tree$edge,
                             edge.length,
                             1:length(edge.length))
            edges <- df[ape::which.edge(new.tree, species.loc),]
            parent.min <- phytools::getParent(new.tree, min(edges[,1]))
            edges2 <- df[df[,1] == parent.min & df[,2] == min(edges[,1]),]
            edges <- rbind(edges2, edges)
            to.index <- sample(edges[,4], 1, prob = edges[,3])
            bind.where <- new.tree$edge[to.index, 2]
            tip.pos <- bind.tip.pos(pos.min = 0,
                                    pos.max = edge.length[to.index])

            new.tree <- phytools::bind.tip(new.tree, 
                                           new.tip, 
                                           edge.length = NULL, 
                                           where = bind.where, 
                                           position = tip.pos)
            added[i+1] <- new.tip
        }
  }
    
  return(new.tree)
}

