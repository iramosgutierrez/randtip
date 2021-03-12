# Return a vector of all species in the tree that match a given genus
sp.genus.in.tree <- function(tree, genus){
    sp <- tree$tip.label
    taxa.vector <- sp[randtip::firstword(sp)==genus]

    return(taxa.vector)
}

# Convenience function to give bind position for tip
# The returning value is to be used in the position argument
# of phytools bind.tip function.
bind.tip.pos <- function(pos.min, pos.max){
     pos.tip <- stats::runif(1, pos.min, pos.max)
     while(pos.tip == pos.min | pos.tip == pos.max){
        pos.tip <- stats::runif(1, pos.min, pos.max)
     }
     return(pos.tip)
}

# Return index of edge where binding should take place
get.index <- function(tree, how = "sample_simple", node = NULL, df = NULL){
    if(is.null(df)){
        df <- data.frame(tree$edge, tree$edge.length, 1:length(tree$edge.length))
    }

    if(!is.null(node)){
        edges <- df[ape::which.edge(tree, node),]
    }else if(how == "sample_simple"){
        edges <- df[sample(x = 1:nrow(df), size = 1),]
    }else if(how == "sample_prob"){
        edges <- df[sample(x = 1:nrow(df), size = 1, prob = df[,3]),]
    }else{
        stop("Unrecognized set of arguments to get.index function.")
    }

    to.index <- edges[,4]

    return(to.index)
}

get.position <- function(tree, node, insertion){

    edge.length <- tree$edge.length[ape::which.edge(tree, node)]
    if(insertion=="random"){
        # Bind at a random point of the branch
        position <- edge.length * stats::runif(1, 0, 1)
    }else if(insertion=="middle"){
        position <- edge.length*0.5
        # Bind at the middle of the branch
    }else if(insertion=="long"){
        # Bind at the beggining of the branch
        position <- edge.length
    }else{
        stop("Unknown specification of the insertion argument.")
    }

    return(position)
}


get.parent.siblings <- function(tree, tip){
    tree.sp <- tree$tip.label
    # Direct ancestor
    parent <- tree$edge[tree$edge[,2] == tip, 1]
    # Ancestor's descendants
    parent.desc <- phytools::getDescendants(tree, parent)
    siblings <- tree.sp[parent.desc]
    siblings <- siblings[!is.na(tree.sp[parent.desc])]

    return(list(parent = parent,
                siblings = siblings))
}

par.sib.diff.genera <- function(tree, siblings.genera, tip){
    # At least one sibling is from the same genus
    while(length(unique(siblings.genera)) == 1){
        # Tip and parent upstream until they are from different genera
        par.sib <- get.parent.siblings(tree, tip)
        siblings.genera <- stringr::word(par.sib$siblings, 1, sep = "_")
        tip <- par.sib$parent
    }

    return(par.sib)
}

sib.intruders <- function(tree, genus, siblings){
    sibling.genera <- stringr::word(siblings, 1, sep = "_")

    if(length(siblings[sibling.genera != genus]) == 1){
        siblings <- siblings[sibling.genera == genus]
    }else if(length(siblings[sibling.genera != genus]) > 1){
        intruders <- siblings[sibling.genera != genus]
        intr.mrca <- ape::getMRCA(tree, intruders)
        intr.desc.tips <- phytools::getDescendants(tree, intr.mrca)
        intr.desc.names <- tree$tip.label[intr.desc.tips]
        intr.desc.names<- intr.desc.names[!is.na(intr.desc.names)]
        intr.genera <- stringr::word(intr.desc.names, 1, sep = "_")
        if(all(intr.genera != genus)){
            siblings<-siblings[sibling.genera==genus]
        }
    }

    return(siblings)
}


get.grouped <- function(tree, siblings.genera, genus, tip){

    par.sib <- par.sib.diff.genera(tree, siblings.genera, tip)
    siblings <- par.sib$siblings
    siblings.genera <- stringr::word(siblings, 1, sep = "_")
    gen.mrca <- phytools::findMRCA(tree = tree,
                                  tips = siblings[siblings.genera == genus])
    gen.mrca.desc <- phytools::getDescendants(tree, gen.mrca)
    # Non-node descendants
    grouped <- tree$tip.label[gen.mrca.desc][!is.na(gen.mrca.desc)]
    # MRCA descendants
    grouped.gen <- stringr::word(grouped, 1, sep = "_")

    return(list(grouped = grouped,
                grouped.gen = grouped.gen,
                gen.mrca = gen.mrca))
}

# Convenience function to add species to node according to group and given
# a vector (group.gen) given by get.grouped function
# Convenience function to add species to node according to phyletic group
add.phyletic.group <- function(tree, sp, genus, grouped, grouped.gen, mrca){

    if(length(unique(grouped.gen))==1){
        # Monophyletic subgroup
        new.tree <- add.into.node(tree, mrca, sp)
    }else if(length(unique(grouped.gen)) > 1){
        # Paraphyletic group
        intruders <- grouped[grouped.gen != genus]
        if(length(intruders) == 1){
            # Singleton intruders; added as monophyletic
            new.tree <- add.into.node(tree, mrca, sp)
        }else{
            # Monophyletic intruders; added as paraphyletic
            intruders.mrca <- phytools::findMRCA(tree = tree, tips = intruders)
            new.tree <- add.to.paraphyletic(tree,
                                            new.tip = sp,
                                            group.node = mrca,
                                            intern.node = intruders.mrca)
        }
    }

    return(new.tree)
}


firstword<- function(string, sep="_"){
  word<- stringr::word(string, 1, sep=sep)
  return(word)
}

notNA <- function(x){
  vect<- x[!is.na(x)]
  return(vect)
}
