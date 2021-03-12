#' A function to stick species at random within a polyphyletic clade
#' @export
add.to.polyphyletic <- function(tree, new.tip, polyphyletic.insertion = "freq", prob=T){

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

    if(polyphyletic.insertion == "all"){
        for(sp in new.tip){
            mrca <- phytools::findMRCA(new.tree, taxa.vector)
            new.tree <- add.into.node(new.tree, node = mrca, new.tip = sp, prob = T)
        }

        return(new.tree)
    }

    if(polyphyletic.insertion == "freq"){
        for(sticked.species in new.tip){

            # Choose a species at random
            # Allows a probability based on frequencies inside each group
            sticking.species<- sample(taxa.vector, 1)
            sticking.tip <- which(new.tree.sp == sticking.species)

            par.sib <- get.parent.siblings(new.tree, sticking.tip)

            siblings.genus <- stringr::word(par.sib$siblings, 1, sep = "_")

            if(sum(siblings.genus == genus) == 1){
                # Siblings are from a different genus
                # It is a singleton inside another clade
                new.tree <- add.to.singleton(new.tree, sticking.species, new.tips = sticked.species)
            }else if(sum(siblings.genus == genus) > 1){
                grouped.list <- get.grouped(new.tree, siblings.genus,
                                            genus, tip = par.sib$parent)

                gen.mrca <- grouped.list$gen.mrca

                grouped <- grouped.list$grouped
                grouped.gen <- grouped.list$grouped.gen

                # Add to tree according to phyletic group
                new.tree <- add.phyletic.group(new.tree, sticked.species,
                                               genus, grouped,
                                               grouped.gen, gen.mrca)
            }
        }
    }

    if(polyphyletic.insertion == "large"){
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

                grouped <- grouped.list$grouped
                grouped.gen <- grouped.list$grouped.gen
                grouped.gen.uniq <- unique(grouped.gen)

                if(length(grouped.gen.uniq) == 1){
                    groups[[t]] <- grouped
                    group.types[[t]] <- "monophyletic"
                }
                if(length(grouped.gen.uniq) > 1){
                    groups[[t]] <- grouped[grouped.gen.uniq == genus]
                    group.types[[t]]<- "paraphyletic"
                }
            }
        }

        groups.sz <- lengths(groups)
        max.sz <- max(as.numeric(groups.sz))
        group.types <- group.types[groups.sz == max.sz]
        groups <- groups[groups.sz == max.sz]

        slot <- sample(1:length(groups), size = 1)
        group <- groups[[slot]]
        group.type <- group.types[[slot]]

        for(sp in new.tip){
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
