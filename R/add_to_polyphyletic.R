#' A function to stick species at random within a polyphyletic clade
#' @export
add.to.polyphyletic <- function(tree, new.tip, poly.ins = "freq", prob=T){

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

    if(poly.ins == "all"){
        for(sp in new.tip){
            mrca <- phytools::findMRCA(new.tree, taxa.vector)
            new.tree <- add.into.node(new.tree, node = mrca, new.tip = sp, prob = T)
        }

        return(new.tree)
    }

    if(poly.ins == "freq"){


      groups <- rep(list(NA), times = length(taxa.vector))
      names(groups) <- taxa.vector
      group.types <- rep(list(NA), times = length(taxa.vector))
      names(group.types) <- taxa.vector

      for(t in 1:length(taxa.vector)){
        taxon <- taxa.vector[t]
        taxon.tip <- which(new.tree$tip.label == taxon)
        par.sib <- get.parent.siblings(new.tree, taxon.tip)

        siblings.genus <- randtip::firstword(par.sib$siblings)

        if(sum(siblings.genus == genus) == 1){
          groups[[t]]<- taxon
          group.types[[t]]<- "singleton"
        }
        if(sum(siblings.genus == genus) > 1){
          grouped.list <- get.grouped(new.tree, siblings.genus,
                                      genus, tip = par.sib$parent)

          grouped <-randtip::notNA(grouped.list$grouped)
          grouped.gen <- randtip::notNA(grouped.list$grouped.gen)
          grouped.gen.uniq <- randtip::notNA(unique(grouped.gen))

          if(length(grouped.gen.uniq) == 1){
            groups[[t]] <- grouped
            group.types[[t]] <- "monophyletic"
          }
          if(length(grouped.gen.uniq) > 1){
            groups[[t]] <- grouped[grouped.gen == genus]

            group.mrca <- ape::getMRCA(new.tree, grouped )
            group.descs<- phytools::getDescendants(new.tree, group.mrca)
            group.descs<- randtip::notNA(new.tree$tip.label[group.descs])
            intruders <- group.descs[randtip::firstword(group.descs)!=genus]
            if(length(intruders)==1){group.types[[t]]<- "paraphyletic"
            next}
            intruders.mrca<- ape::getMRCA(new.tree, intruders )
            intruders.descs<- phytools::getDescendants(new.tree, intruders.mrca)
            intruders.descs<- randtip::notNA(new.tree$tip.label[intruders.descs])
            intruders.genera <- randtip::firstword(intruders.descs)
            if(genus %in% intruders.genera){group.types[[t]]<- "polyphyletic"}else{
              group.types[[t]]<- "paraphyletic"
            }

          }
        }
      }

      non.polyphyletic<- which(as.vector(unlist(group.types))!="polyphyletic")
      groups<-groups[non.polyphyletic]
      group.types<-group.types[non.polyphyletic]
      taxa.vector.select<- taxa.vector[non.polyphyletic]

      for(sp in new.tip){
        randsp<- sample(x = taxa.vector.select, size = 1)
        slot<- grep(randsp, x = groups)
        if(length(slot)>1){slot<-sample(x=slot, size = 1)}

        group <- groups[[slot]]
        group.type <- group.types[[slot]]

        if(group.type == "singleton"){
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

    if(poly.ins == "large"){


        for(sp in new.tip){

          groups <- rep(list(NA), times = length(taxa.vector))
          names(groups) <- taxa.vector
          group.types <- rep(list(NA), times = length(taxa.vector))
          names(group.types) <- taxa.vector

          for(t in 1:length(taxa.vector)){
            taxon <- taxa.vector[t]
            taxon.tip <- which(new.tree$tip.label == taxon)
            par.sib <- get.parent.siblings(new.tree, taxon.tip)

            siblings.genus <- randtip::firstword(par.sib$siblings)

            if(sum(siblings.genus == genus) == 1){
              groups[[t]]<- taxon
              group.types[[t]]<- "singleton"
            }
            if(sum(siblings.genus == genus) > 1){
              grouped.list <- get.grouped(new.tree, siblings.genus,
                                          genus, tip = par.sib$parent)

              grouped <-randtip::notNA(grouped.list$grouped)
              grouped.gen <- randtip::notNA(grouped.list$grouped.gen)
              grouped.gen.uniq <- randtip::notNA(unique(grouped.gen))

              if(length(grouped.gen.uniq) == 1){
                groups[[t]] <- grouped
                group.types[[t]] <- "monophyletic"
              }
              if(length(grouped.gen.uniq) > 1){
                groups[[t]] <- grouped[grouped.gen == genus]

                group.mrca <- ape::getMRCA(new.tree, grouped )
                group.descs<- phytools::getDescendants(new.tree, group.mrca)
                group.descs<- randtip::notNA(new.tree$tip.label[group.descs])
                intruders <- group.descs[randtip::firstword(group.descs)!=genus]
                if(length(intruders)==1){group.types[[t]]<- "paraphyletic"
                next}
                intruders.mrca<- ape::getMRCA(new.tree, intruders )
                intruders.descs<- phytools::getDescendants(new.tree, intruders.mrca)
                intruders.descs<- randtip::notNA(new.tree$tip.label[intruders.descs])
                intruders.genera <- randtip::firstword(intruders.descs)
                if(genus %in% intruders.genera){group.types[[t]]<- "polyphyletic"}else{
                  group.types[[t]]<- "paraphyletic"
                }

              }
            }
          }

          non.polyphyletic<- which(as.vector(unlist(group.types))!="polyphyletic")
          groups<-groups[non.polyphyletic]
          group.types<-group.types[non.polyphyletic]
          taxa.vector.select<- taxa.vector[non.polyphyletic]

          groups.sz <- lengths(groups)
          max.sz <- max(as.numeric(groups.sz))
          group.types <- group.types[groups.sz == max.sz]
          groups <- groups[groups.sz == max.sz]

          slot <- sample(1:length(groups), size = 1)
          group <- groups[[slot]]
          group.type <- group.types[[slot]]



          if(group.type == "singleton"){
            new.tree <- randtip::add.to.singleton(tree = new.tree, singleton = group, new.tips = sp)}

            if(group.type == "monophyletic"){
                mrca <- phytools::findMRCA(new.tree, tips = group)
                new.tree <- randtip::add.into.node(tree = new.tree, node = mrca,
                                          new.tip = sp, prob=prob)}
            if(group.type == "paraphyletic"){
                mrca <- phytools::findMRCA(new.tree, tips = group)
                mrca.desc <- phytools::getDescendants(new.tree, mrca)
                mrca.desc <- mrca.desc[!is.na(mrca.desc)]
                mrca.desc.gen <- stringr::word(mrca.desc, 1, sep = "_")
                intruders <- mrca.desc[mrca.desc.gen != genus]
                intruders.mrca <- phytools::findMRCA(new.tree, tips = intruders)
                new.tree <- add.to.paraphyletic(tree = new.tree,
                                                new.tip = sp,
                                                group.node = mrca,
                                                intern.node = intruders.mrca, prob=prob)}
            }

    }

    return(new.tree)

}


