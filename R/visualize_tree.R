#' Function to obtain tree cut from MDCC
#' @export
getClade<- function(tree, DF0, cladeinfo=list("level"=NA, "clade"=NA)){
if(is.null(names(cladeinfo))){names(cladeinfo)<- c("level", "clade")}

  level<-cladeinfo$level
  clade<- cladeinfo$clade

  spss<- DF0[which(DF0[,level]==clade),]
  genera<- unique(spss$genus)

  cut.list<- tree$tip.label[randtip::firstword(tree$tip.label)%in% genera]
  if(length(cut.list)==0){stop("Specified clade is not reflected in the tree!")}
  if(length(cut.list)==1){stop("Specified clade is related to just 1 tree tip!")}
  cut.node<- ape::getMRCA(tree, tip =cut.list )
  subtree<-  phytools::splitTree(tree, split = list("node"=cut.node, "bp"=0))[[2]]
  return(list("Tree"=subtree, "DF0"=DF0, "level"=level, "clade"=clade))
}

clade_col <- function(getClade.output, sharingtaxa.col="green",
                      intruder.col="red",stowaway.col="black"){

  CladeTree<-getClade.output$Tree
  level <- getClade.output$level
  clade <- getClade.output$clade
  DF0   <- getClade.output$DF0

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

#a<-getClade(tree, DF0, list("genus","Dendrolagus"))
#plot(a[[1]], tip.color = clade_col(a))
