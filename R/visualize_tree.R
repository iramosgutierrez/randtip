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
plot.clade<- function(get.clade.out, corankedtaxa.col="#03C03C",
                      intruder.col="#C23B23",stowaway.col="black", ...){
  if(!(is.list(get.clade.out)|all(names(get.clade.out)==c("Tree","info","rank","clade")))){
    stop("Please feed this function with the returned object from get.clade function")
  }

  tipcol<- randtip::clade.col(get.clade.out, corankedtaxa.col=corankedtaxa.col,
                               intruder.col=intruder.col, stowaway.col=stowaway.col)
  return(ape::plot.phylo(get.clade.out$Tree, tip.color = tipcol, ...))
}


clade.col <- function(get.clade.out, corankedtaxa.col,
                      intruder.col,stowaway.col){

  CladeTree<-get.clade.out$Tree
  rank <- get.clade.out$rank
  clade <- get.clade.out$clade
  info   <- get.clade.out$info

  spss<- info[which(info[,rank]==clade),]
  genera<- unique(spss$genus)


  intruders<- info[which(info[,rank]!=clade),]
  intrudergenera<- unique(intruders$genus)

  colours<- vector("character", length(CladeTree$tip.label))
  colours[]<- stowaway.col
  colours[randtip::firstword(CladeTree$tip.label)%in%intrudergenera]<- intruder.col
  colours[randtip::firstword(CladeTree$tip.label)%in%        genera]<- corankedtaxa.col

  return(colours)

}

