#' A function to stick species at random within a polyphyletic clade
#' @export
add.to.polyphyletic <- function(tree, new.tip, polyphyly.scheme = "freq", prob=T, respect.mono=F, respect.para=F){

    new.tip <- gsub(" ", "_", new.tip)
    genus <- unique(randtip::firstword(new.tip))
    if(length(genus) > 1){
        stop("Species from more than 1 genera are included.")
    }
    taxa.vector <- sp.genus.in.tree(tree, genus = genus)
    if(length(taxa.vector)==0){
      stop(paste0("Genus ", genus, " is not included in your tree."))
    }

    new.tree <- tree
    new.tree.sp <- new.tree$tip.label

    if(polyphyly.scheme == "complete"){
        for(sp in new.tip){
            mrca <- phytools::findMRCA(new.tree, taxa.vector)
            if(isFALSE(respect.mono)) {nodes<- phytools::getDescendants(new.tree, mrca)}
            if(isTRUE(respect.mono)& isFALSE(respect.para)){nodes<-get.permitted.nodes(new.tree, mrca, respect.para = F)}
            if(isTRUE(respect.mono)& isTRUE(respect.para)) {nodes<-get.permitted.nodes(new.tree, mrca, respect.para = T)}
            bpos<-binding.position(new.tree, node = sample(nodes,1), insertion="random", prob=prob)

            new.tree <- phytools::bind.tip(tree = new.tree, tip.label = sp, edge.length = bpos$length,
                                 where = bpos$where, position = bpos$position)
        }

        return(new.tree)
    }

    if(polyphyly.scheme == "freq"){



      groups <- get.groups(tree = new.tree, genus=genus)

      for(sp in new.tip){
        randsp<- sample(x = taxa.vector, size = 1)
        slot<- grep(randsp, x = groups$species)
        if(length(slot)>1){slot<-sample(x=slot, size = 1)}

        group <- groups$species[[slot]]
        group.type <- groups$type[[slot]]

        if(group.type == "singleton" | group.type == "tip"){
           new.tree <- randtip::add.to.singleton(tree = new.tree, singleton = group, new.tips = sp)}
        if(group.type == "monophyletic"){
          mrca <- phytools::findMRCA(new.tree, tips = group)
          new.tree <- randtip::add.into.node(tree = new.tree, node = mrca,
                                             new.tip = sp, prob=prob)}
        if(group.type == "paraphyletic"){
          mrca <- phytools::findMRCA(new.tree, tips = group)
          mrca.desc <- phytools::getDescendants(new.tree, mrca)
          mrca.desc <- mrca.desc[!is.na(mrca.desc)]
          mrca.desc <- randtip::notNA(new.tree$tip.label[mrca.desc])
          mrca.desc.gen <- randtip::firstword(mrca.desc)
          intruders <- mrca.desc[mrca.desc.gen != genus]
          if(length(intruders)==1){
            new.tree<-add.into.node(new.tree, node = mrca, new.tip = sp, prob = prob)
          }else{
          intruders.mrca <- phytools::findMRCA(new.tree, tips = intruders)
          new.tree <- add.to.paraphyletic(tree = new.tree,
                                          new.tip = sp,
                                          group.node = mrca,
                                          intern.node = intruders.mrca,
                                          prob=prob)}
          }
        if(group.type == "polyphyletic"){
          stop("A polyphyletic subgroup can't be selected. Code Error")
          }
      }
      }

    if(polyphyly.scheme == "large"){


        for(sp in new.tip){

          groups <- get.groups(tree = new.tree, genus=genus)

          groups.sz <- vector(mode="numeric", length = length(groups$species))
          for(g in seq_along(groups$species)){
            vect<- groups$species[[g]]
            groups.sz[g]<-length(vect[randtip::firstword(vect)==genus])}

          max.sz <- max(as.numeric(groups.sz))
          max.id <- which(groups.sz == max.sz)
          if(length(max.id)>1){max.id<-sample(max.id, size = 1)}

          group <- groups$species[[max.id]]
          group.type <- groups$type[[max.id]]


          if(group.type == "singleton" | group.type == "tip"){
            new.tree <- randtip::add.to.singleton(tree = new.tree, singleton = group, new.tips = sp)}
          if(group.type == "monophyletic"){
            mrca <- phytools::findMRCA(new.tree, tips = group)
            new.tree <- randtip::add.into.node(tree = new.tree, node = mrca,
                                               new.tip = sp, prob=prob)}
          if(group.type == "paraphyletic"){
            mrca <- phytools::findMRCA(new.tree, tips = group)
            mrca.desc <- phytools::getDescendants(new.tree, mrca)
            mrca.desc <- mrca.desc[!is.na(mrca.desc)]
            mrca.desc <- randtip::notNA(new.tree$tip.label[mrca.desc])
            mrca.desc.gen <- randtip::firstword(mrca.desc)
            intruders <- mrca.desc[mrca.desc.gen != genus]
            if(length(intruders)==1){
              new.tree<-add.into.node(new.tree, node = mrca, new.tip = sp, prob = prob)
            }else{
              intruders.mrca <- phytools::findMRCA(new.tree, tips = intruders)
              new.tree <- add.to.paraphyletic(tree = new.tree,
                                              new.tip = sp,
                                              group.node = mrca,
                                              intern.node = intruders.mrca,
                                              prob=prob)}
          }
          if(group.type == "polyphyletic"){
            stop("A polyphyletic subgroup can´t be selected. Code Error")
          }
            }

    }

    return(new.tree)

}

