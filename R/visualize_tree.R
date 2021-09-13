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

clade.col <- function(get.clade.output, sharingtaxa.col="green",
                      intruder.col="red",stowaway.col="black"){

  CladeTree<-get.clade.output$Tree
  level <- get.clade.output$level
  clade <- get.clade.output$clade
  DF0   <- get.clade.output$DF0

  spss<- DF0[which(DF0[,level]==clade),]
  genera<- unique(spss$genus)


  intruders<- DF0[which(DF0[,level]!=clade),]
  intrudergenera<- unique(intruders$genus)

  colours<- vector("character", length(CladeTree$tip.label))
  colours[]<- stowaway.col
  colours[randtip::firstword(CladeTree$tip.label)%in%intrudergenera]<- intruder.col
  colours[randtip::firstword(CladeTree$tip.label)%in%        genera]<- sharingtaxa.col

  return(colours)

}

