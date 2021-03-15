#' function phyleticity retrieves genus type given a certain tree
#'
#' @export
#' @examples
#' #phyleticity(Tree, "Invent")
phyleticity<- function(tree, genus){

    if(length(genus) != 1){ stop("Only one genus accepted.") }

    sp <- tree$tip.label
    taxa.vector <- sp[randtip::firstword(sp)==genus]

    if(length(taxa.vector) == 0){
        genus.type <- "Not included"
        return(genus.type)
    }

    if(length(taxa.vector) == 1){
        genus.type <- "Singleton genus"
        return(genus.type)
    }

    mrca <- phytools::findMRCA(tree = tree, tips = taxa.vector)
    descend <- phytools::getDescendants(tree, mrca)
    desc.tips <- sp[descend]
    desc.tips <- randtip::notNA(desc.tips)
    desc.genera <- randtip::firstword(desc.tips)
    if(length(unique(desc.genera)) == 1){
        genus.type <- "Monophyletic"
        return(genus.type)
    }

    intruder.tips <- desc.tips[desc.genera != genus]
    intruder.mrca <- phytools::findMRCA(tree = tree, tips = intruder.tips)
    if(length(intruder.tips) == 1){
        #(and there is only one intruder: paraphyletic by singleton)
        genus.type <- "Paraphyletic"
    }else{
        # if there are more than one intruder...
        desc.intruder.mrca <- phytools::getDescendants(tree, intruder.mrca)
        desc.intruder.tips <- sp[desc.intruder.mrca]
        desc.intruder.tips <- notNA(desc.intruder.tips)
        #intruders grouped: paraphyletic group by monophyletic intruder
        suppressWarnings({
            grouped.intruders <- all(sort(desc.intruder.tips) ==
                                     sort(intruder.tips))
        })
        if(grouped.intruders){
            genus.type <- "Paraphyletic"
        }else{
            genus.type <- "Polyphyletic"
        }
    }

    return(genus.type)
}
