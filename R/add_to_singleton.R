#' Function to add tip to singleton
#' @export
add.to.singleton <- function(tree, singleton, new.tips){
  singleton<-gsub(pattern = " ", "_", singleton)

  new.tree <- tree
  node<-which(new.tree$tip.label==singleton)

  for(i in 1: length(new.tips)){

    if (i==1) {

      pos<- randtip::binding.position(new.tree, node = node, insertion = "random",prob = T)



      new.tree <- phytools::bind.tip(new.tree,
                                     new.tips[i],
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position)
    }else{
      sticksp<-  c(singleton,new.tips[1:i-1] )
      node<- phytools::findMRCA(tree=new.tree, tips=sticksp)
      permitted<- c(node, phytools::getDescendants(new.tree, node))
      df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2], "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )
      df<- df[df$node%in%permitted,]
      pos<- randtip::binding.position(new.tree, df = df, insertion = "random",prob = T)

      new.tree <- phytools::bind.tip(new.tree,
                                     new.tips[i],
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position)
    }

  }


  return(new.tree)
}

