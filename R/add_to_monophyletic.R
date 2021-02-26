#' Function to add tips at random monophyletic clade
#' @export
#' @examples
#' #plot(add_to_monophyletic(Tree, "Invent_inventata"))
add.to.monophyletic <- function(tree, new.tip){
    
  new.tip <- gsub(" ", "_", new.tip)    
  
  genus <- stringr::word(new.tip, 1, sep = "_")     
  
  sp <- tree$tip.label
  taxa.vector <- sp[stringr::str_detect(sp, genus)]    
  if(length(taxa.vector) == 0){
      stop(paste0("Genus ", genus, " is not included in yout tree."))
  }
    
  mrca <- phytools::findMRCA(tree = tree, tips = taxa.vector) 
  new.tree <- add.into.node(tree = tree, node = mrca, new.tip = new.tip) 
  return(new.tree)
}
