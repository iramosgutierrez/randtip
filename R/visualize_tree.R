#' Function to plot tree cut from MDCC
#' @export
plotTree.MDCC<- function(tree, DF1, MDCC, DF2=NULL){

  spss<- DF1[c(which(MDCC==DF1$MDCC), which(MDCC==DF1$other.MDCC), which(MDCC==DF1$genus)),]
  if(!is.null(DF2)){spss2<- DF2[which(MDCC==DF2$MDCC),]}
  genera<- unique(spss$genus)
  if(!is.null(DF2)){genera<- unique(c(spss$genus, spss2$genus)) }
  cut.list<- tree$tip.label[stringr::word(tree$tip.label, 1, sep="_")%in% genera]
  if(length(cut.list)==0){stop("Specified MDCC is not reflected in the tree!")}
  if(length(cut.list)==1){stop("Specified MDCC is related to just 1 tree tip!")}
  cut.node<- ape::getMRCA(tree, tip =cut.list )
  subtree<-  phytools::splitTree(tree, split = list("node"=cut.node, "bp"=0))[[2]]
  return(plot(subtree))
}
#plotTree.MDCC(tree25, DF1 =DF1, MDCC="Pinaceae")


#' Function to obtain tree cut from MDCC
#' @export
cutTree.MDCC<- function(tree, DF1, MDCC, DF2=NULL){

  spss<- DF1[c(which(MDCC==DF1$MDCC), which(MDCC==DF1$other.MDCC), which(MDCC==DF1$genus)),]
  if(!is.null(DF2)){spss2<- DF2[which(MDCC==DF2$MDCC),]}
  genera<- unique(spss$genus)
  if(!is.null(DF2)){genera<- unique(c(spss$genus, spss2$genus)) }
  cut.list<- tree$tip.label[stringr::word(tree$tip.label, 1, sep="_")%in% genera]
  if(length(cut.list)==0){stop("Specified MDCC is not reflected in the tree!")}
  if(length(cut.list)==1){stop("Specified MDCC is related to just 1 tree tip!")}
  cut.node<- ape::getMRCA(tree, tip =cut.list )
  subtree<-  phytools::splitTree(tree, split = list("node"=cut.node, "bp"=0))[[2]]
  return((subtree))
}
#cutTree.MDCC(tree, DF1, MDCC="Clinopodium")
