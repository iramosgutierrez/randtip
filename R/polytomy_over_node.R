#' @export
polytomy.over.node<- function(tree, node, species, 
                              insertion = c("random", "long", "middle")){

    new.tree <- tree 
    for(i in 1:length(species)){
        # First species is added to singleton
        if(i==1){   
            to.index <- get.index(new.tree, node = node)
            bind.where <- new.tree$edge[to.index,2]
            
            position <- get.position(new.tree, node, insertion)
    
            new.tree <- phytools::bind.tip(new.tree, 
                                           tip.label =  species[i], 
                                           edge.length = NULL, 
                                           where = bind.where, 
                                           position = position)
        }else{
            # The rest of the species are added to the node
            child.node <- which(new.tree$tip.label == species[i-1])
            node <- phytools::getParent(tree = new.tree, 
                                        node = child.node)
                
            new.tree <- polytomy.into.node(tree = new.tree, 
                                          new.tip = species[i], 
                                          node = node)
        }
    }
    return(new.tree)}

