#' Function to add tips at random monophyletic clade
#' @export
#' @examples
#' #plot(add_to_monophyletic(Tree, "Invent_inventata"))
add.to.monophyletic <- function(tree, new.tip){
  new.tree<- tree
  new.tip <- gsub(" ", "_", new.tip)

  genus <- randtip::firstword(new.tip)

  sp <- new.tree$tip.label
  taxa.vector <- sp[randtip::firstword(sp)==genus]
  if(length(taxa.vector) == 0){
      stop(paste0("Genus ", genus, " is not included in your tree."))
  }

  mrca <- phytools::findMRCA(tree = new.tree, tips = taxa.vector)
  new.tree <- add.into.node(tree = new.tree, node = mrca, new.tip = new.tip, prob )
  return(new.tree)
}
