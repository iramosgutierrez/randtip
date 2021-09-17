#' Function to obtain tree cut from MDCC
#' @export
get.clade<- function(tree, DF0, level, clade){


  spss<- DF0[which(DF0[,level]==clade),]
  genera<- unique(spss$genus)

  cut.list<- tree$tip.label[randtip::firstword(tree$tip.label)%in% genera]
  if(length(cut.list)==0){stop("Specified clade is not reflected in the tree!")}
  if(length(cut.list)==1){stop("Specified clade is related to just 1 tree tip!")}
  cut.node<- ape::getMRCA(tree, tip =cut.list )
  subtree<-  phytools::splitTree(tree, split = list("node"=cut.node, "bp"=0))[[2]]
  return(list("Tree"=subtree, "DF0"=DF0, "level"=level, "clade"=clade))
}

#' Function to plot subtree
#' @export
plot.clade<- function(get.clade.out, corankedtaxa.col="#03C03C",
                      intruder.col="#C23B23",stowaway.col="black", ...){
  if(!(is.list(get.clade.out)|all(names(get.clade.out)==c("Tree","DF0","level","clade")))){
    stop("Please feed this function with the returned object from get.clade function")
  }

  edgecol<- randtip::clade.col(get.clade.out, corankedtaxa.col=corankedtaxa.col,
                               intruder.col=intruder.col, stowaway.col=stowaway.col)
  return(ape::plot.phylo(get.clade.out$Tree, edge.color = edgecol, ...))
}


clade.col <- function(get.clade.out, corankedtaxa.col,
                      intruder.col,stowaway.col){

  CladeTree<-get.clade.out$Tree
  level <- get.clade.out$level
  clade <- get.clade.out$clade
  DF0   <- get.clade.out$DF0

  spss<- DF0[which(DF0[,level]==clade),]
  genera<- unique(spss$genus)


  intruders<- DF0[which(DF0[,level]!=clade),]
  intrudergenera<- unique(intruders$genus)

  colours<- vector("character", length(CladeTree$tip.label))
  colours[]<- stowaway.col
  colours[randtip::firstword(CladeTree$tip.label)%in%intrudergenera]<- intruder.col
  colours[randtip::firstword(CladeTree$tip.label)%in%        genera]<- corankedtaxa.col

  return(colours)

}

