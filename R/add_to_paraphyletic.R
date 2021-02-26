#' Same process than paraphyletic node, bt detecting automatically the deep and shallow nodes
#' @param group.node mrca node of the paraphyletic clade. If no intern.node is supplied, a random 
#'                    internal node is selected
#' @param intern.node node of the clade to insert the new tip. It must be used in conjuction 
#'                     with group node
#' @export
#' @examples
#' #add.to.paraphyletic(Tree, "Achillea_3_inventatus")
add.to.paraphyletic <- function(tree, new.tip, prob = TRUE, intern.node = NULL,
                                group.node = NULL){
  
  if(is.null(group.node) && !is.null(intern.node)){
      stop("The group.node argument must be supplied along with intern.node argument.")
  }
    
  new.tip <- gsub(" ", "_", new.tip) 
  
  genus <- stringr::word(new.tip, 1, sep="_")
  
  taxa.vector<- sp.genus.in.tree(tree, genus)
  if(length(taxa.vector) == 0){
      stop(paste0("Genus ", genus, " is no included in yout tree."))
  }
    
  sp <- tree$tip.label
  if(!is.null(group.node)){
      mrca <- group.node
      descs <- phytools::getDescendants(tree, group.node)
  }else{
      mrca <- phytools::findMRCA(tree = tree, tips = taxa.vector)
      descs <- phytools::getDescendants(tree, mrca)
  }
  descs.tips <- sp[phytools::getDescendants(tree, mrca)]
  descs.tips <- descs.tips[!is.na(descs.tips)]
  intruder.tips <- descs.tips[!stringr::str_detect(descs.tips, genus)]
    
  if((length(intruder.tips) == 0) & is.null(intern.node)){
      stop(paste0("Genus ", genus, " is not paraphyletic. It is monophyletic."))
  }
  
  if((length(intruder.tips) == 1) & is.null(intern.node)){
      new.tree <- add.into.node(tree = tree, node = mrca, new.tip = new.tip)
  }
  
  if(length(intruder.tips) > 1 | !is.null(intern.node)){
    if(!is.null(intern.node)){
        intruder.descs <- phytools::getDescendants(tree, intern.node)
        intruder.descs <- intruder.descs[!is.na(intruder.descs)]
    }else{
        intruder.mrca <- phytools::findMRCA(tree, tips = intruder.tips)
        intruder.descs <- phytools::getDescendants(tree, intruder.mrca)
        intruder.descs <- intruder.descs[!is.na(intruder.descs)]
        if(any(stringr::str_detect(sp[intruder.descs], genus))){
            stop(paste0("Genus ", genus, 
                        " is not paraphyletic. It is polyphyletic."))
        }
    }

    edge.length <- tree$edge.length
    df <- data.frame(tree$edge, edge.length, 1:length(edge.length)) 
    df <- df[df[,2] %in% descs,]
    df <- df[!(df[,2] %in% intruder.descs),]
    if(prob == TRUE){
        to.index <- get.index(tree, how = "sample_prob", df = df)
    }else{
        to.index <- get.index(tree)
    }
    where <- tree$edge[to.index, 2] 
      
      
    tip.pos <- bind.tip.pos(pos.min = 0, 
                            pos.max = edge.length[to.index])
   
    new.tree <- phytools::bind.tip(tree, new.tip, edge.length = NULL, 
                         where = where, position = tip.pos)
   }
  
    return(new.tree)
}
