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
        genus.type <- "Missing"
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



MDCC.phyleticity<-function(DF1, tree, MDCC.info=list("level"=NA, "MDCC"=NA), trim=T){


  level<- MDCC.info$level
  MDCC <- MDCC.info$MDCC
  DF1<-DF1[!is.na(DF1[,level]),]

   tips<- tree$tip.label[randtip::firstword(tree$tip.label)%in%randtip::firstword(DF1$taxon)]
   if(isTRUE(trim)){tree<- ape::keep.tip(phy =tree, tip = tips)}

   species<- DF1[which(DF1[,level]==MDCC),]
   MDCC.genera<- unique(randtip::firstword(species$taxon))
   genera.in.tree<-randtip::firstword(tree$tip.label)
   genera.in.tree<-MDCC.genera[MDCC.genera%in%genera.in.tree]
   spp.in.tree<- tree$tip.label[randtip::firstword(tree$tip.label)%in%genera.in.tree]

   if(length(spp.in.tree)==0){MDCC.type<-"Missing"
   return(MDCC.type) }

   if(length(spp.in.tree)==1){MDCC.type<-"Singleton MDCC"
   return(MDCC.type) }

  sp.mrca<- phytools::findMRCA(tree, tips = spp.in.tree)
  descs.num<- phytools::getDescendants(tree,sp.mrca)
  descs.name<-randtip::notNA(tree$tip.label[descs.num])




 if(all(randtip::firstword(descs.name)%in%MDCC.genera)){MDCC.type<-"Monophyletic"
 return(MDCC.type) }else{

   intruders<-descs.name[!(randtip::firstword(descs.name)%in%MDCC.genera)]

   if(length(intruders)==1){MDCC.type<-"Paraphyletic"
   return(MDCC.type) }

   intruders.mrca<- ape::getMRCA(phy = tree, tip = intruders)
   intruders.descs.num<- phytools::getDescendants(tree,intruders.mrca)
   intruders.descs.name<-randtip::notNA(tree$tip.label[intruders.descs.num])

   if(!any(randtip::firstword(intruders.descs.name)%in%MDCC.genera)){MDCC.type<-"Paraphyletic"
   return(MDCC.type) }else{
     MDCC.type<-"Polyphyletic"
     return(MDCC.type)}


 }
}
