#' Function to extract a subtree representing a specified group.
#'
#' Obtain a phylogenetic tree representing a specified group split from an
#' original backbone tree.
#' Note that the tree splitting will be performed at the MRCA node of all the taxa included in the
#' group in the 'info' object.
#'
#' @param info An 'info' (or 'input') object.
#' @param tree The original backbone tree to be split.
#' @param group group to be extracted
#'
#' @return A list of four objects which will be used for automatic plotting
#'         using \code{\link{plot_clade}} function.
#'  'Tree' will contain the splitted tree; 'info' will contain the handed info
#'         file; 'rank' will contain the taxonomic rank of the specified group,
#'         and 'group' will contain the group name.
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
        stop("Specified group is not reflected in the 'info' data frame")}
    if(length(ranks)> 1){
        stop("Specified group reflect several ranks. ",
            "Please correct your 'info' data frame!")
    }
    rank<-names(rankDF)[ranks]

    spss<- info[which(info[,rank]==group),]
    genera<- unique(spss$genus)

    cut.list<- tree$tip.label[first_word(tree$tip.label)%in% genera]
    if(length(cut.list)==0){stop("Specified group is not reflected in the tree!")}
    if(length(cut.list)==1){stop("Specified group is represented by a single tip in the phylogeny!")}
    cut.node<- ape::getMRCA(tree, tip = cut.list )
    if(cut.node==findRoot(tree)){
        return(list("Tree"=tree, "info"=info, "rank"=rank, "group"=group))
    }

    subtree<-  ape::extract.clade(tree, node = cut.node, root.edge = 0 )

    return(list("Tree"=subtree, "info"=info, "rank"=rank, "group"=group))
}


#' Function to plot a subtree representing a specified group.
#'
#' Plot a phylogenetic tree splitted from the backbone tree
#' using \code{get_clade} function.
#'
#' @param get.clade.out Output from \code{\link{get_clade}} function.
#' @param ppcr.col Color to represent tips included in the specified group.
#'                 Default value is green.
#' @param nonppcr.col Color to represent tips included in a different group
#'                    from the specified one (at the same taxonomic rank).
#'                    Default value is blue.
#' @param unknown.col Color to represent tips without taxonomic information
#'                    at the specifed group's taxonomic rank.
#'                    Default value is grey
#' @param ... Arguments to pass through \code{\link{plot.phylo}} function.
#'
#' @return A plot representing the group specified in \code{\link{get_clade}}
#'         function using the selected colors.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#' #First the group information must be obtained
#'
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
#' #Then it can be plotted
#' plot_clade(felinae.clade, ppcr.col="green",
#' nonppcr.col="red",unknown.col="grey" )
#'
#' @export
plot_clade<- function(get.clade.out, ppcr.col="#4a8a21",
                      nonppcr.col="#48bce0",unknown.col="#adadad", ...){

    get.group.names <- c("Tree", "info", "rank", "group")
    if(!(is.list(get.clade.out)|all(names(get.clade.out)==get.group.names))){
        stop("Please feed this function with the returned object from ",
            "get_clade function")
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
#' Set color pattern to distinguish between phylogenetically placed taxa and PUTs when plotting trees
#' with the \code{plot.phylo} function of ‘ape’ R package using the \code{tip.color} argument
#'
#' @param newtree An expanded phylogenetic tree.
#' @param oldtree The original backbone tree where the PUTs have been bound.
#' @param placed.col Color to plot phylogenetic tips which were already placed in the original backbone tree.
#'                   Default value is grey.
#' @param put.col Color to plot bound PUTs in \code{new tree}. Default value is red.
#' @return A vector of length equal to the number of tips in newtree, to be used after \code{tip.color}
#'         in plot.phylo function.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#' #Perform a tree expansion
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
#' #Set the colours for original tips and bound PUTs
#' cats.tip.cols <- put_tip_col(newtree = expanded.cats,
#'  oldtree = cats, placed.col="black", put.col="red")
#'
#' #Plot the resulting tree visualizing original tips and PUTs
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

    get.group.names <- c("Tree", "info", "rank", "group")
    if(!(is.list(get.clade.out)|all(names(get.clade.out)==get.group.names))){
        stop("Please feed this function with the returned object from ",
            "get_clade function")
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
