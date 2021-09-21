#' Function to obtain tree cut from MDCC
#' @export
get.clade<- function(tree, info, clade){


  rankDF<-info[,c("taxon", randtip::randtip_ranks())]

  rankDF.withclade<-as.data.frame(rankDF[,]==clade)
  rankDF.withclade<- as.vector(colSums(rankDF.withclade, na.rm = T))
  ranks<- which(rankDF.withclade>0)
  if(length(ranks)==0){stop("Specified clade is not reflected in the tree!")}
  if(length(ranks)> 1){stop("Specified clade reflects several ranks. Please correct your 'info' data frame!")}
  rank<-names(rankDF)[ranks]


  spss<- info[which(info[,rank]==clade),]
  genera<- unique(spss$genus)

  cut.list<- tree$tip.label[randtip::firstword(tree$tip.label)%in% genera]
  if(length(cut.list)==0){stop("Specified clade is not reflected in the tree!")}
  if(length(cut.list)==1){stop("Specified clade is related to just 1 tree tip!")}
  cut.node<- ape::getMRCA(tree, tip =cut.list )
  if(length(cut.list)==1){stop("Specified clade is related to just 1 tree tip!")}
  if(cut.node==randtip::findRoot(tree)){
    return(list("Tree"=tree, "info"=info, "rank"=rank, "clade"=clade))}

  subtree<-  phytools::splitTree(tree, split = list("node"=cut.node, "bp"=0))[[2]]
  return(list("Tree"=subtree, "info"=info, "rank"=rank, "clade"=clade))
}

#' Function to plot subtree
#' @export
plot.clade<- function(get.clade.out, coranked.col="#4a8a21",
                      noncoranked.col="#48bce0",unknown.col="#d4c744", ...){
  if(!(is.list(get.clade.out)|all(names(get.clade.out)==c("Tree","info","rank","clade")))){
    stop("Please feed this function with the returned object from get.clade function")
  }

  tipcol<- randtip::clade.col(get.clade.out, coranked.col=coranked.col,
                               noncoranked.col=noncoranked.col, unknown.col=unknown.col)
  return(ape::plot.phylo(get.clade.out$Tree, tip.color = tipcol, ...))
}


clade.col <- function(get.clade.out, coranked.col,
                      noncoranked.col,unknown.col){

  CladeTree<-get.clade.out$Tree
  rank <- get.clade.out$rank
  clade <- get.clade.out$clade
  info   <- get.clade.out$info

  spss<- info[which(info[,rank]==clade),]
  genera<- unique(spss$genus)


  intruders<- info[which(info[,rank]!=clade),]
  intrudergenera<- unique(intruders$genus)

  colours<- vector("character", length(CladeTree$tip.label))
  colours[]<- unknown.col
  colours[randtip::firstword(CladeTree$tip.label)%in%intrudergenera]<- noncoranked.col
  colours[randtip::firstword(CladeTree$tip.label)%in%        genera]<- coranked.col

  return(colours)

}

