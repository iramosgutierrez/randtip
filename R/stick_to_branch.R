#' Function to add species at specified branches.
#' #' @param tree "phylo" object used as backbone tree.
#'    @param edge matrix with 4 columns character vectors.\n
#'        Columns 1 and 2 define the stem node of candidate branches and
#'        columns 3 and 4 define the crown node
#'        (mrca of species pairs). Note that a pair of nodes may not
#'        necessarily define one single branch but a set of them.
#'   @param new.tip Name of the PUT to bind to the specified candidate branches.
#'   @param rand.type "random" or "polytomy"
#'
#' @export
custom.branch <- function(tree, edge, rand.type="random", forceultrametric=F, prob=T){

  if(rand.type == "r"){rand.type <- "random"}
  if(rand.type == "p"){rand.type <- "polytomy"}

  new.tree<- tree

  if(forceultrametric & !ape::is.ultrametric(new.tree)){new.tree<- phytools::force.ultrametric(new.tree)}
  if(!(rand.type %in% c("random", "polytomy"))) {stop("Argument 'rand.type' must be \"random\" or \"polytomy\" ")}
  if(isFALSE(forceultrametric) & !ape::is.ultrametric(new.tree)){
    message("The backbone tree is not ultrametric. \nPlease, set the argument 'forceultrametric' to TRUE if the tree is genuinely ultrametric.")}

  PUTs <- unique(edge[,1])

  for(PUT in PUTs){

    edge.i   <- edge[edge[,1]==PUT,]

    df <- data.frame("parent"=new.tree$edge[,1], "node"=new.tree$edge[,2],
                     "length"= new.tree$edge.length, "id"=1:length(new.tree$edge[,1]) )

    permittednodes<- as.numeric(NULL)
    root<- randtip::findRoot(new.tree)
    PUT  <- gsub(" ", "_", PUT)

    edge.i[,2]<- gsub(" ", "_", edge.i[,2])
    edge.i[,3]<- gsub(" ", "_", edge.i[,3])
    edge.i[,4]<- gsub(" ", "_", edge.i[,4])
    edge.i[,5]<- gsub(" ", "_", edge.i[,5])

    for(i in 1:nrow(edge.i)){
      if(!all(c(edge.i[i,2],edge.i[i,3],edge.i[i,4],edge.i[i,5])%in%new.tree$tip.label)){
        message("Row ", i, " of 'edge' is not defining a phylogenetic branch.")
        next}

    if(edge.i[i,2]==edge.i[i,3]){parnode<- which(new.tree$tip.label==edge.i[i,2])}else{
      parnode<- ape::getMRCA(new.tree, c(edge.i[i,2], edge.i[i,3]))}

    if(edge.i[i,4]==edge.i[i,5]){basenode<- which(new.tree$tip.label==edge.i[i,4])}else{
      basenode<- ape::getMRCA(new.tree, c(edge.i[i,4], edge.i[i,5]))}

    if(edge.i[i,2]==edge.i[i,3] & edge.i[i,2]== edge.i[i,4] &
       edge.i[i,2]==edge.i[i,5]){
      parnode <- randtip::get.parent.siblings(tree, which(tree$tip.label==edge[i,2]))$parent
    }





    if(parnode==basenode){
      message("Row ", i, " of 'edge' is not defining a phylogenetic branch.")
      next}

    perm.nodes.i<- as.numeric(NULL)
    n<- basenode
    p<- parnode
    while(n!=p){
      if(n==root){
        message("Row ", i, " of 'edge' is not defining a phylogenetic branch.")
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

    if(nrow(df)==1){nd <- df$node}
    if(nrow(df)>1 & prob) {nd<-sample(df$node, 1, prob=df$length)}
    if(nrow(df)>1 & !prob){nd<-sample(df$node, 1)}

    if(rand.type=="polytomy"){
      if(nrow(df)==1){nd <- df$parent}
      if(nrow(df)>1 & prob ){nd<-sample(df$parent, 1, prob=df$length)}
      if(nrow(df)>1 & !prob){nd<-sample(df$parent, 1)}
      }
    bp<-randtip::binding.position(new.tree, node = nd, insertion = rand.type)
    new.tree<- phytools::bind.tip(new.tree, PUT, bp$length, bp$where, bp$position)
}
    return(new.tree)
}






#' Function to help users visualize the candidate branches for tip insertion with
#' \code{\link{stick.to.branch}} function
#' @export
#' @examples
#' set.seed(1)
plot.custom.branch<- function(tree, edge, PUT=NULL,
                              candidate.col="#bf2828", forbidden.col="#3d3d3d",
                              candidate.lwd=2, forbidden.lwd=1,...){

  df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                   "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )

  if(length(unique(edge[,1]))==1 & is.null(PUT)){PUT<- unique(edge[,1])}
  if(!is.null(PUT)){
    if(!(PUT %in% edge[,1])){stop("The specified PUT is not contained in the 'edge' data frame.")}
    edge<- edge[edge[,1]==PUT]
  }
  if(length(unique(edge[,1]))>1){stop("Your 'edge' data frame contains more than one PUT. Please specify only one to plot its candidate branches.")}

  permittednodes<- as.numeric(NULL)

  root<- randtip::findRoot(tree)
  edge[,2]<- gsub(" ", "_", edge[,2])
  edge[,3]<- gsub(" ", "_", edge[,3])
  edge[,4]<- gsub(" ", "_", edge[,4])
  edge[,5]<- gsub(" ", "_", edge[,5])


  for(i in 1:nrow(edge)){
    if(!all(c(edge[i,2],edge[i,3],edge[i,4],edge[i,5])%in%tree$tip.label)){
      message("Row ", i, " has species not included in the tree and will not be used.")
      next}

    if(edge[i,4]==edge[i,5]){basenode<- which(tree$tip.label==edge[i,4])}else{
      basenode<- ape::getMRCA(tree, c(edge[i,4], edge[i,5]))}

    if(edge[i,2]==edge[i,3]){parnode<- which(tree$tip.label==edge[i,2])}else{
      parnode<- ape::getMRCA(tree, c(edge[i,2], edge[i,3]))}

    if(edge[i,2]==edge[i,3] & edge[i,2]== edge[i,4] &
       edge[i,2]==edge[i,5]){
      parnode <- randtip::get.parent.siblings(tree, which(tree$tip.label==edge[i,2]))$parent}

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

