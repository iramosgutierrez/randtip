#' this function reveals if a given node at a given tree is an internal node
#' @export
is.node<-function(tree, node){
  if(!(node %in% tree$edge)){
      stop("Node number is not in your tree")
  }
  if(length(phytools::getDescendants(tree = tree, node = node )) > 1){
      return(TRUE)
  }else{
      return(FALSE)
  }
}

#' this function reveals if a given node at a given tree is a tip
#' @export
is.tip <-function(tree, node){
  if(!(node %in% tree$edge)){
      stop("Node number is not in your tree")
  }
  if(length(phytools::getDescendants(tree = tree, node = node)) == 1){
      return(TRUE)
  }else{
      return(FALSE)
  }
}
