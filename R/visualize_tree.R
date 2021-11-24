#' Function to obtain tree cut from MDCC
#' @export
get.clade<- function(info, tree, clade){


  rankDF<-info[,c("taxon", randtip::randtip_ranks())]

  rankDF.withclade<-as.data.frame(rankDF[,]==clade)
  rankDF.withclade<- as.vector(colSums(rankDF.withclade, na.rm = T))
  ranks<- which(rankDF.withclade>0)
  if(length(ranks)==0){stop("Specified clade is not reflected in the 'info' data frame")}
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
plot.clade<- function(get.clade.out, ppcr.col="#4a8a21",
                      nonppcr.col="#48bce0",unknown.col="#adadad", ...){
  if(!(is.list(get.clade.out)|all(names(get.clade.out)==c("Tree","info","rank","clade")))){
    stop("Please feed this function with the returned object from get.clade function")
  }

  tipcol<- randtip::clade.col(get.clade.out, ppcr.col=ppcr.col,
                               nonppcr.col=nonppcr.col, unknown.col=unknown.col)
  return(ape::plot.phylo(get.clade.out$Tree, tip.color = tipcol, ...))
}


clade.col <- function(get.clade.out, ppcr.col="#4a8a21",
                      nonppcr.col="#48bce0",unknown.col="#adadad"){

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
  colours[randtip::firstword(CladeTree$tip.label)%in%intrudergenera]<- nonppcr.col
  colours[randtip::firstword(CladeTree$tip.label)%in%        genera]<- ppcr.col

  return(colours)

}


#' Get PUT or placed color pattern
#'
#' Get a color vector to use in \'tip.color\' argument within plot.phylo function
#'
#' @param newtree An expanded phylogenetic tree.
#' @param oldtree The original backbone tree where the PUTs have been bound.
#' @param placed.col Color to plot phylogenetic tips which were already placed in the original backbone tree. Default value is grey.
#' @param put.col Color to plot bound PUTs in \code{new tree}. Default value is red.
#' @return A vector of length equal to the number of tips in newtree, to be used after \'tip.color\' in plot.phylo function.
#'
#' @author Ignacio Ramos-Gutiérrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @export
put.tip.col<- function(newtree, oldtree, placed.col="#adadad", put.col="#C23B23"){
  col<- vector("character", length(newtree$tip.label))
  col[1:length(col)]<- put.col
  col[newtree$tip.label%in%oldtree$tip.label]<-placed.col
  return(col)

}

