#' Function to extract clades from the backbone phylogeny.
#'
#' This function serves to extract the clade in the backbone phylogeny that includes all the 
#' species in the specified group (as defined in the info data frame).
#'
#' @param info An 'info' (or 'input') data frame.
#' @param tree The backbone tree to be subsetted.
#' @param clade Taxonomic group of species that defines the clade to extract ## CREO QUE ES MÁS CORRECTO QUE EL ARGUMENTO SEA "GROUP", NO "CLADE", YA QUE LO QUE SE DEFINE ES UN GRUPO, QUE A SU VEZ DEFINE UN CLADO.
#'
#' @return A list with four objects that will be used for automatic plotting
#'         with the \code{\link{plot_clade}} function.
#'  'Tree' is the subtree; 'info' is the handed info
#'         data frame; 'rank' is the taxonomic rank of the specified group,
#'         and 'clade' is the name of the clade. ## SAME AS ABOVE, SHOULD BE 'group'

#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#' catspecies <- c("Lynx_lynx", "Panthera_uncia",
#' "Panthera_onca", "Felis_catus", "Puma_concolor",
#' "Lynx_canadensis", "Panthera_tigris", "Panthera_leo",
#' "Felis_silvestris")
#'
#' cats.info <- build_info(species=catspecies, tree= cats,
#'      find.ranks=TRUE, db="ncbi", mode="backbone")
#'
#' felinae.clade <- get_clade(info=cats.info,
#' tree=cats, group="Felinae")
#'
#' @export
get_clade<- function(info, tree, group){

    rankDF<-info[,c("taxon", randtip_ranks())]
    rankDF.withgroup<-as.data.frame(rankDF[,]==group)
    rankDF.withgroup<- as.vector(colSums(rankDF.withgroup, na.rm = T))
    ranks<- which(rankDF.withgroup>0)
    if(length(ranks)==0){
        stop("Specified taxonomic group not found in 'info'")}
    if(length(ranks)> 1){
        stop("Specified group is duplicated through the taxonomic hierarchy. ",
            "Please, correct your 'info' data frame!")
    }
    rank<-names(rankDF)[ranks]

    spss<- info[which(info[,rank]==group),]
    genera<- unique(spss$genus)

    cut.list<- tree$tip.label[first_word(tree$tip.label)%in% genera]

    if(length(cut.list)==0){stop("Specified taxonomic group not found in the tree!")}
    if(length(cut.list)==1){stop("Specified taxonomic group is represented by one single tip in the phylogeny!")}
    cut.node<- ape::getMRCA(tree, tip = cut.list )
    if(cut.node==findRoot(tree)){
        return(list("Tree"=tree, "info"=info, "rank"=rank, "group"=group))
    }

    subtree<-  ape::extract.clade(tree, node = cut.node, root.edge = 0 )

    return(list("Tree"=subtree, "info"=info, "rank"=rank, "group"=group))
}



#' Function to plot clades extracted from the backbone phylogeny. 
#'

#' @param get.clade.out Output of the \code{\link{get_clade}} function (a phylogenetic clade)
#' @param ppcr.col Color for the tips representing species placed in the specified taxonomic group (default is green).
#' @param nonppcr.col Color for the tips representing other species but those placed in the specified taxonomic groupnt clade (default is blue).
#' @param unknown.col Color for the tips representing species without taxonomic information (default is grey).
#' @param ... Further arguments to pass through \code{\link{plot.phylo}} function.
#'
#' @return Less inclusive clade that includes all the species in the specified taxonomic group. 

#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#' catspecies <- c("Lynx_lynx", "Panthera_uncia",
#' "Panthera_onca", "Felis_catus", "Puma_concolor",
#' "Lynx_canadensis", "Panthera_tigris", "Panthera_leo",
#' "Felis_silvestris")
#'
#' cats.info <- build_info(species=catspecies, tree= cats,
#'      find.ranks=TRUE, db="ncbi", mode="backbone")
#'
#' felinae.clade <- get_clade(info=cats.info,
#' tree=cats, group="Felinae")
#'
#' plot_clade(felinae.clade, ppcr.col="green",
#' nonppcr.col="red",unknown.col="grey" )
#'
#' @export
plot_clade<- function(get.clade.out, ppcr.col="#4a8a21",
                      nonppcr.col="#48bce0",unknown.col="#adadad", ...){


    get.clade.names <- c("Tree", "info", "rank", "clade")
    if(!(is.list(get.clade.out)|all(names(get.clade.out)==get.clade.names))){
        stop("Please, feed this function with the output from  get_clade function")
    }

    tipcol <- clade_col(get.clade.out, ppcr.col=ppcr.col,
                      nonppcr.col=nonppcr.col, unknown.col=unknown.col)
    legcols <- c(ppcr.col, nonppcr.col, unknown.col)
    legnames<- c("PPCR taxa", "Non-PPCR taxa", "Unknown")


    ape::plot.phylo(get.clade.out$Tree, tip.color = tipcol, ...)
    graphics::legend("topright", legend = legnames, border = legcols, fill=legcols,
           cex = 0.7, bty = "n",  title = paste0(get.clade.out$group, " ",
                                                 get.clade.out$rank))


}

#' Get PUT or placed color pattern
#'
#' Set a color pattern to distinguish between phylogenetically placed taxa and PUTs when plotting expanded trees
#' with the \code{plot.phylo} function of ‘ape’ R package
#'
#' @param newtree An expanded phylogenetic tree.
#' @param oldtree The original backbone tree.
#' @param placed.col Color for phylogenetic tips that were already placed in the backbone tree (default is grey).
#' @param put.col Color for phylogenetic tips representing PUTs in \code{new tree} (default is red).
#' @return A vector of length equal to the number of tips in newtree that defines the argument \code{tip.color} in plot.phylo function.
#'         
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#'
#' catspecies <- c("Lynx_lynx", "Panthera_uncia",
#' "Panthera_onca", "Felis_catus", "Puma_concolor",
#' "Lynx_canadensis", "Panthera_tigris", "Panthera_leo",
#' "Felis_silvestris")
#'
#' cats.info <- build_info(species=catspecies, tree= cats,
#'      find.ranks=TRUE, db="ncbi", mode="backbone")
#'
#' cats.input <- info2input(info=cats.info, tree=cats)
#'
#' expanded.cats <- rand_tip(input=cats.input,
#'  tree=cats, rand.type = "polytomy",
#'  forceultrametric = TRUE)
#'
#' #Set color pattern for original tips and PUTs, respectively
#' cats.tip.cols <- put_tip_col(newtree = expanded.cats,
#'  oldtree = cats, placed.col="black", put.col="red")
#'
#' plot.phylo(expanded.cats, tip.color = cats.tip.cols)
#'
#' @export
put_tip_col<- function(newtree, oldtree, placed.col="#adadad", put.col="#C23B23"){
    col<- vector("character", length(newtree$tip.label))
    col[1:length(col)]<- put.col
    col[newtree$tip.label%in%oldtree$tip.label]<-placed.col
    return(col)
}

clade_col <- function(get.clade.out, ppcr.col="#4a8a21",
                      nonppcr.col="#48bce0",unknown.col="#adadad"){


    get.clade.names <- c("Tree", "info", "rank", "clade")
    if(!(is.list(get.clade.out)|all(names(get.clade.out)==get.clade.names))){
        stop("Please, feed this function with the output from get_clade function")
    }

    group.tree<-get.clade.out$Tree
    rank <- get.clade.out$rank
    group <- get.clade.out$group
    info   <- get.clade.out$info

    spss<- info[which(info[,rank]==group),]
    genera<- unique(spss$genus)

    intruders<- info[which(info[,rank]!=group),]
    intrudergenera<- unique(intruders$genus)

    colours<- vector("character", length(group.tree$tip.label))
    colours[]<- unknown.col
    colours[first_word(group.tree$tip.label)%in%intrudergenera]<- nonppcr.col
    colours[first_word(group.tree$tip.label)%in%genera]<- ppcr.col

    return(colours)

}
