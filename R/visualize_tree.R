#' Function to obtain tree cut from MDCC
#' @export
split.MDCC<- function(tree, DF1, MDCC.info=list("level"=NA, "MDCC"=NA), DF2=NULL){

  level<-MDCC.info$level
  MDCC<- MDCC.info$MDCC

  spss<- DF1[which(DF1[,level]==MDCC),]
  if(!is.null(DF2)){spss2<- DF2[which(DF2[,level]==MDCC),]}
  genera<- unique(spss$genus)
  if(!is.null(DF2)){genera<- unique(c(spss$genus, spss2$genus)) }
  cut.list<- tree$tip.label[randtip::firstword(tree$tip.label)%in% genera]
  if(length(cut.list)==0){stop("Specified MDCC is not reflected in the tree!")}
  if(length(cut.list)==1){stop("Specified MDCC is related to just 1 tree tip!")}
  cut.node<- ape::getMRCA(tree, tip =cut.list )
  subtree<-  phytools::splitTree(tree, split = list("node"=cut.node, "bp"=0))[[2]]
  return((subtree))
}

