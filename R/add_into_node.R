#' Function to add tips at random from a given node
#' @export
#' @examples
#' exception.list<-list(c("Abies_alba","Abies_pinsapo"),
#'                      c("Achillea_pyrenaica","Achillea_santolinoides",
#'                        "Achillea_chamaemelifolia"))


add.into.node <- function(tree, node, new.tip, exception.list = NULL, prob){

  new.tip <- gsub(" ", "_", new.tip)

  if(is.null(exception.list)){
    descs<-phytools::getDescendants(tree = tree, node = node, curr = T)
    descs.DF<- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2], "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )
    descs.DF<-descs.DF[descs.DF$node%in%descs,]

    if(isTRUE(prob)){ rand.desc<- sample(descs.DF$node, size = 1, prob = descs.DF$length)}else{rand.desc<- sample(descs.DF$node, size = 1, prob = NULL)}
    new.tree<-randtip::add.over.node(tree = tree, node=rand.desc, new.tip = new.tip)

    return(new.tree)
       }else{



  forbidden.taxa<- unlist(exception.list)
  descs <- phytools::getDescendants(tree = tree, node = node)
  descs.tips <- tree$tip.label[descs]
  descs.tips <- descs.tips[!is.na(descs.tips)]

  if(all(descs.tips %in% forbidden.taxa)){
    warning(paste0("All descendant tips coincide with the exception list! ",
                   "Adding over the node of the most recent common ancestor"))
    rand.sp <- sample(descs.tips, 1)
    rand.grp <- exception.list[[grep(rand.sp, exception.list )]]
    new.tree <- add.over.node(tree, new.tip = new.tip,
                              node = ape::getMRCA(tree, tip = rand.grp))
    return(new.tree)
  }else{


  df <- data.frame(tree$edge, 1:nrow(tree$edge), tree$edge.length)
  colnames(df) <- c("parent.node", "child.node", "edge.id", "edge.length")

  # List to store identity of forbidden branches
  forbid.edge <- vector(mode = "list",
                        length = length(exception.list))

  for(i in 1:length(exception.list)){
    except.grp <- exception.list[[i]]
    id.tip <- NA

    for(j in 1:length(except.grp)){
      id.tip <- c(id.tip,
                  which(tree$tip.label == except.grp[j]))

    }
    if(!all(is.na(id.tip))){
      id.tip <- id.tip[!(is.na(id.tip))]
    }else{
      next
    }

    pare.tip <- NA
    # Identify group terminal branches and identity of internal nodes
    for(j in 1:length(id.tip)) {
      pare.tip <- c(pare.tip,
                    phytools::getParent(tree, id.tip[j]))

    }
    pare.tip <- unique(pare.tip[!is.na(pare.tip)])
    id.tip.mrca <- ape::getMRCA(tree,id.tip)
    is.id.tip <- (df$child.node %in% id.tip)
    forbid.nodes <- (df$child.node <= max(pare.tip)) & (df$child.node > id.tip.mrca)
    forbid.edge[[i]] <- df[(forbid.nodes | is.id.tip), "edge.id"]

  }

  # Subset df to include only the edges allowed for tip insertion
  df.to.stick <- df[!(df$edge.id %in% unlist(forbid.edge)), ]
  df.to.stick <- df.to.stick[df.to.stick$child.node %in% descs,]

  # Sampling probability is given by edge lenght
  to.index <- sample(df.to.stick$edge.id, 1, prob = df.to.stick$edge.length)
  bind.where <- df[to.index, "child.node"]

  tip.pos <- bind.tip.pos(pos.min = 0,
                          pos.max = tree$edge.length[to.index])
  new.tree <- phytools::bind.tip(tree, new.tip, edge.length = NULL,
                                 where = bind.where, position = tip.pos)

  return(new.tree)
}}}
