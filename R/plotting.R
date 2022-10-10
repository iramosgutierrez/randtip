#' Function to extract a subtree representing a specified clade.
#'
#' Obtain a phylogenetic tree representing a specified clade split from an
#' original backbone tree.
#' Note that the tree splitting will be performed at the MRCA node of all the taxa included in the
#' clade in the 'info' object.
#'
#' @param info An 'info' (or 'input') object.
#' @param tree The original backbone tree to be split.
#' @param clade Clade to be extracted
#'
#' @return A list of four objects which will be used for automatic plotting
#'         using \code{\link{plot.clade}} function.
#'  'Tree' will contain the splitted tree; 'info' will contain the handed info
#'         file; 'rank' will contain the taxonomic rank of the specified clade,
#'         and 'clade' will contain the clade name.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#' felinae.clade <- get.clade(info=cats.info,
#' tree=cats, clade="Felinae")
#'
#' @export
get.clade<- function(info, tree, clade){

    rankDF<-info[,c("taxon", randtip.ranks())]
    rankDF.withclade<-as.data.frame(rankDF[,]==clade)
    rankDF.withclade<- as.vector(colSums(rankDF.withclade, na.rm = T))
    ranks<- which(rankDF.withclade>0)
    if(length(ranks)==0){
        stop("Specified clade is not reflected in the 'info' data frame")}
    if(length(ranks)> 1){
        stop("Specified clade reflect several ranks. ",
            "Please correct your 'info' data frame!")
    }
    rank<-names(rankDF)[ranks]

    spss<- info[which(info[,rank]==clade),]
    genera<- unique(spss$genus)

    cut.list<- tree$tip.label[first.word(tree$tip.label)%in% genera]
    if(length(cut.list)==0){stop("Specified clade is not reflected in the tree!")}
    if(length(cut.list)==1){stop("Specified clade is represented by a single tip in the phylogeny!")}
    cut.node<- ape::getMRCA(tree, tip = cut.list )
    if(cut.node==findRoot(tree)){
        return(list("Tree"=tree, "info"=info, "rank"=rank, "clade"=clade))
    }

    subtree<-  ape::extract.clade(tree, node = cut.node, root.edge = 0 )

    return(list("Tree"=subtree, "info"=info, "rank"=rank, "clade"=clade))
}


#' Function to plot a subtree representing a specified clade.
#'
#' Plot a phylogenetic tree splitted from the backbone tree
#' using \code{get.clade} function.
#'
#' @param get.clade.out Output from \code{\link{get.clade}} function.
#' @param ppcr.col Color to represent tips included in the specified clade.
#'                 Default value is green.
#' @param nonppcr.col Color to represent tips included in a different clade
#'                    from the specified one (at the same taxonomic rank).
#'                    Default value is blue.
#' @param unknown.col Color to represent tips without taxonomic information
#'                    at the specifed clade's taxonomic rank.
#'                    Default value is grey
#' @param ... Arguments to pass through \code{\link{plot.phylo}} function.
#'
#' @return A plot representing the clade specified in \code{\link{get.clade}}
#'         function using the selected colors.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#' #First the clade information must be obtained
#' felinae.clade <- get.clade(info=cats.info,
#' tree=cats, clade="Felinae")
#'
#' #Then it can be plotted
#' plot.clade(felinae.clade, ppcr.col="green",
#' nonppcr.col="red",unknown.col="grey" )
#'
#' @export plot.clade
#' @export
plot.clade<- function(get.clade.out, ppcr.col="#4a8a21",
                      nonppcr.col="#48bce0",unknown.col="#adadad", ...){

    get.clade.names <- c("Tree", "info", "rank", "clade")
    if(!(is.list(get.clade.out)|all(names(get.clade.out)==get.clade.names))){
        stop("Please feed this function with the returned object from ",
            "get.clade function")
    }

    tipcol <- clade.col(get.clade.out, ppcr.col=ppcr.col,
                      nonppcr.col=nonppcr.col, unknown.col=unknown.col)

    return(ape::plot.phylo(get.clade.out$Tree, tip.color = tipcol, ...))
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
#' expanded.cats <- rand.tip(input=cats.input,
#'  tree=cats, rand.type = "polytomy",
#'  forceultrametric = T)
#'
#' #Set the colours for original tips and bound PUTs
#' cats.tip.cols <- put.tip.col(newtree = expanded.cats,
#'  oldtree = cats, placed.col="black", put.col="red")
#'
#' #Plot the resulting tree visualizing original tips and PUTs
#' plot.phylo(expanded.cats, tip.color = cats.tip.cols)
#'
#' @export
put.tip.col<- function(newtree, oldtree, placed.col="#adadad", put.col="#C23B23"){
    col<- vector("character", length(newtree$tip.label))
    col[1:length(col)]<- put.col
    col[newtree$tip.label%in%oldtree$tip.label]<-placed.col
    return(col)
}

clade.col <- function(get.clade.out, ppcr.col="#4a8a21",
                      nonppcr.col="#48bce0",unknown.col="#adadad"){

    get.clade.names <- c("Tree", "info", "rank", "clade")
    if(!(is.list(get.clade.out)|all(names(get.clade.out)==get.clade.names))){
        stop("Please feed this function with the returned object from ",
            "get.clade function")
    }

    clade.tree<-get.clade.out$Tree
    rank <- get.clade.out$rank
    clade <- get.clade.out$clade
    info   <- get.clade.out$info

    spss<- info[which(info[,rank]==clade),]
    genera<- unique(spss$genus)

    intruders<- info[which(info[,rank]!=clade),]
    intrudergenera<- unique(intruders$genus)

    colours<- vector("character", length(clade.tree$tip.label))
    colours[]<- unknown.col
    colours[first.word(clade.tree$tip.label)%in%intrudergenera]<- nonppcr.col
    colours[first.word(clade.tree$tip.label)%in%genera]<- ppcr.col

    return(colours)

}
