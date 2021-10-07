#' Function to add species at specified branches.
#' #' @param tree "phylo" object used as backbone tree.
#'    @param edges matrix with 4 columns character vectors.\n
#'        Columns 1 and 2 define the stem node of candidate branches and
#'        columns 3 and 4 define the crown node
#'        (mrca of species pairs). Note that a pair of nodes may not
#'        necessarily define one single branch but a set of them.
#'   @param new.tip Name of the PUT to bind to the specified candidate branches.
#'   @param rand.type "random" or "polytomy"
#'
#' @export
custom.branch <- function(tree, edges, new.tip, rand.type="random", prob=T){

    df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                     "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )
    permittednodes<- as.numeric(NULL)
    root<- randtip::findRoot(tree)
    new.tip  <- gsub(" ", "_", new.tip)
    edges[,1]<- gsub(" ", "_", edges[,1])
    edges[,2]<- gsub(" ", "_", edges[,2])
    edges[,3]<- gsub(" ", "_", edges[,3])
    edges[,4]<- gsub(" ", "_", edges[,4])


    for(i in 1:nrow(edges)){
      if(!all(c(edges[i,1],edges[i,2],edges[i,3],edges[i,4])%in%tree$tip.label)){
        message("Row ", i, " has species not included in the tree and will not be used.")
        next}

    if(edges[i,1]==edges[i,2]){basenode<- which(tree$tip.label==edges[i,1])}else{
      basenode<- ape::getMRCA(tree, c(edges[i,1], edges[i,2]))}

    if(edges[i,3]==edges[i,4]){parnode<- which(tree$tip.label==edges[i,3])}else{
      parnode<- ape::getMRCA(tree, c(edges[i,3], edges[i,4]))}

    if(parnode==basenode){
      message("Row ", i, " specifies a unique node and not a branch, so it will not be used.")
      next}

    perm.nodes.i<- as.numeric(NULL)
    n<- basenode
    p<- parnode
    while(n!=p){
      if(n==root){
        message("Row ", i, " does not reflect a set of branches, so it will not be used.")
        perm.nodes.i<- NULL
        break}
      perm.nodes.i<- c(perm.nodes.i, n)
      n<-df$parent[df$node==n]


    }

    permittednodes<- c(permittednodes, perm.nodes.i)
    }
    if(length(permittednodes)==0){stop("No branches could be selected")}
    permittednodes<- unique(permittednodes)

    df<- df[df$node%in%permittednodes,]

    if(prob){nd<-sample(df$node, 1, prob=df$length)}else{nd<-sample(df$node, 1)}
    if(rand.type=="polytomy"){if(prob){nd<-sample(df$parent, 1, prob=df$length)}else{nd<-sample(df$parent, 1)}}
    bp<-randtip::binding.position(tree, node = nd, insertion = rand.type)
    new.tree<- phytools::bind.tip(tree, new.tip, bp$length, bp$where, bp$position)

    return(new.tree)
}






#' Function to help users visualize the candidate branches for tip insertion with
#' \code{\link{stick.to.branch}} function
#' @export
#' @examples
#' set.seed(1)
plot.custom.branch<- function(tree, edges,
                              candidate.col="#4a8a21", forbidden.col="#3d3d3d",
                              candidate.lwd=1, forbidden.lwd=1,...){

  df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                   "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )
  permittednodes<- as.numeric(NULL)
  root<- randtip::findRoot(tree)
  edges[,1]<- gsub(" ", "_", edges[,1])
  edges[,2]<- gsub(" ", "_", edges[,2])
  edges[,3]<- gsub(" ", "_", edges[,3])
  edges[,4]<- gsub(" ", "_", edges[,4])


  for(i in 1:nrow(edges)){
    if(!all(c(edges[i,1],edges[i,2],edges[i,3],edges[i,4])%in%tree$tip.label)){
      message("Row ", i, " has species not included in the tree and will not be used.")
      next}

    if(edges[i,1]==edges[i,2]){basenode<- which(tree$tip.label==edges[i,1])}else{
      basenode<- ape::getMRCA(tree, c(edges[i,1], edges[i,2]))}

    if(edges[i,3]==edges[i,4]){parnode<- which(tree$tip.label==edges[i,3])}else{
      parnode<- ape::getMRCA(tree, c(edges[i,3], edges[i,4]))}

    if(parnode==basenode){
      message("Row ", i, " specifies a unique node and not a branch, so it will not be used.")
      next}

    perm.nodes.i<- as.numeric(NULL)
    n<- basenode
    p<- parnode
    while(n!=p){
      if(n==root){
        message("Row ", i, " does not reflect a set of branches, so it will not be used.")
        perm.nodes.i<- NULL
        break}
      perm.nodes.i<- c(perm.nodes.i, n)
      n<-df$parent[df$node==n]


    }

    permittednodes<- c(permittednodes, perm.nodes.i)
  }
  if(length(permittednodes)==0){stop("No branches could be selected")}
  permittednodes<- unique(permittednodes)

  col <- vector(mode = "character", length(tree$edge[,1]))
  col[1:length(col)] <- forbidden.col
  col[df$node%in%permittednodes]<- candidate.col

  lwd <- vector(mode = "character", length(tree$edge[,1]))
  lwd[1:length(lwd)] <- forbidden.lwd
  lwd[df$node%in%permittednodes]<- candidate.lwd



  return(ape::plot.phylo(tree, edge.color = col, edge.width = lwd, ...))
}

