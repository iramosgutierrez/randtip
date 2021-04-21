#' Insert in paraphyletic. i no nodes are given, they will be automatically searched.
#' @param group.node mrca node of the paraphyletic clade. It must be used in conjuction
#'                    with intern node
#' @param intern.node mrca node of the clade where not to insert the new tip. It must be used in conjuction
#'                     with group node
#' @export
#' @examples
#' #add.to.paraphyletic(Tree, "Achillea_3_inventatus")

add.to.paraphyletic <- function(tree, new.tip, prob = TRUE, intern.node = NULL,
                                group.node = NULL){

  if((is.null(group.node) && !is.null(intern.node)) | (!is.null(group.node) && is.null(intern.node))){
    stop("The group.node argument must be supplied along with intern.node argument.")
  }


  if(is.null(group.node) && is.null(intern.node)){

  new.tip <- gsub(" ", "_", new.tip)
  genus <- randtip::firstword(new.tip)

  taxa.vector<- sp.genus.in.tree(tree, genus)
  if(length(taxa.vector) == 0){
    stop(paste0("Genus ", genus, " is no included in your tree."))
    }

  sp <- tree$tip.label

    mrca <- phytools::findMRCA(tree = tree, tips = taxa.vector)
    descs <- phytools::getDescendants(tree, mrca)

  descs.tips <- sp[phytools::getDescendants(tree, mrca)]
  descs.tips <- randtip::notNA(descs.tips)
  intruder.tips <- descs.tips[randtip::firstword(descs.tips)!=genus]

  if(length(intruder.tips) == 0){
    stop(paste0("Genus ", genus, " is not paraphyletic. It is monophyletic."))
  }

  if(length(intruder.tips) == 1 ){
    new.tree <- add.into.node(tree = tree, node = mrca, new.tip = new.tip)
  }

  if(length(intruder.tips) > 1 ){

      intruder.mrca <- phytools::findMRCA(tree, tips = intruder.tips)
      intruder.descs <- phytools::getDescendants(tree, intruder.mrca)
      intruder.desc.names<- tree$tip.label[intruder.descs]
      intruder.desc.names <- randtip::notNA(intruder.desc.names)
      if(any(randtip::firstword(intruder.desc.names)==genus)){
        stop(paste0("Genus ", genus,
                    " is not paraphyletic. It is polyphyletic."))
      }

  }}else{
      descs<-phytools::getDescendants(tree = tree, node = group.node)
      intruder.descs<-phytools::getDescendants(tree = tree, node = intern.node)
    }


    df <- data.frame("parent.node"=tree$edge[,1],"child.node"= tree$edge[,2],"length"=tree$edge.length, "id" =1:length(tree$edge.length))
    df <- df[df[,2] %in% descs,]
    if(exists("intruder.descs")){df <- df[!(df[,2] %in% intruder.descs),]}
    pos<- randtip::binding.position(tree = tree, df = df, insertion = "random", prob=prob )

    new.tree <- phytools::bind.tip(tree, new.tip, edge.length = pos$length,
                                   where = pos$where, position = pos$position)


  return(new.tree)
}
