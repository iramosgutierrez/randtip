#'
#' Bind PUTs at completely customized tree branches
#'
#' @param tree "phylo" object used as backbone tree.
#' @param edges matrix with 5 character vector columns.
#'     Column 1 secifies the PUT to which the candidate branches in the
#'     row refer to.
#'     Columns 2 and 3 must represent two species whose MRCA represents the
#'     stem node of candidate branches and
#'     columns 4 and 5 two species whose MRCA define the crown node.
#'     Note that a pair of nodes may not necessarily define one single branch
#'     but a set of them.
#'     Two identical species will define a tree tip rather than an internal node.
#'     If the pairs of species in columns 2&3 and 4&5 are the same, the set of
#'     branches to select will be defined as the complete clade below ther MRCA.
#' @param rand.type "random" or "polytomy". Default value is "random".
#' @param forceultrametric Whether or not the backbone tree will be forced to be ultrametric,
#'                         only in case it is not. Default value is FALSE.
#' @param prob Whether or not branch selection probability must be proportional
#'             to branch length or equiprobable. Default value is TRUE.
#'
#' @return An expanded phylogeny.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examplesIf interactive()
#' #Create a 'edges' dataframe
#' cats.edges <- data.frame(
#'  "PUT"= "Felis_catus",
#'  "parent1"= "Felis_silvestris",
#'  "parent2"= "Felis_silvestris",
#'  "child1"= "Felis_silvestris",
#'  "child2"= "Felis_silvestris")
#'
#' #Bind the PUT to one of the selected branches
#' cats.expanded <- custom_branch(tree=cats,
#'  edges=cats.edges, forceultrametric=TRUE)
#' @export
custom_branch <- function(tree, edges, rand.type="random",
                          forceultrametric=FALSE, prob=TRUE){

    if(rand.type == "r"){rand.type <- "random"}
    if(rand.type == "p"){rand.type <- "polytomy"}

    new.tree<- tree

    if(forceultrametric & !ape::is.ultrametric(new.tree)){
        new.tree<- phytools::force.ultrametric(new.tree)
    }
    if(!(rand.type %in% c("random", "polytomy"))) {
        stop("Argument 'rand.type' must be \"random\" or \"polytomy\" ")
    }
    if(isFALSE(forceultrametric) & !ape::is.ultrametric(new.tree)){
        message("The backbone tree is not ultrametric. ",
                "\nPlease, set the argument 'forceultrametric' to TRUE ",
                "if the tree is genuinely ultrametric.")
    }

    PUTs <- unique(edges[,1])

    for(PUT in PUTs){
        edges.i <- edges[edges[,1]==PUT,]

        df <- data.frame("parent"=new.tree$edge[,1],
                        "node"=new.tree$edge[,2],
                        "length"= new.tree$edge.length,
                        "id"=1:length(new.tree$edge[,1]) )

        root<- findRoot(new.tree)
        PUT  <- gsub(" ", "_", PUT)

        edges.i[,2]<- gsub(" ", "_", edges.i[,2])
        edges.i[,3]<- gsub(" ", "_", edges.i[,3])
        edges.i[,4]<- gsub(" ", "_", edges.i[,4])
        edges.i[,5]<- gsub(" ", "_", edges.i[,5])

        permittednodes <- get_permitted_nodes_custom(new.tree, df, edges.i, root)

        df<- df[df$node%in%permittednodes,]

        if(nrow(df)==1){nd <- df$node}
        if(nrow(df)>1 & prob) {nd<-sample(df$node, 1, prob=df$length)}
        if(nrow(df)>1 & !prob){nd<-sample(df$node, 1)}

        if(rand.type=="polytomy"){
            if(nrow(df)==1){nd <- df$parent}
            if(nrow(df)>1 & prob ){nd<-sample(df$parent, 1, prob=df$length)}
            if(nrow(df)>1 & !prob){nd<-sample(df$parent, 1)}
        }

        bp<-binding_position(new.tree, node = nd, insertion = rand.type,
                             ultrametric = ape::is.ultrametric(new.tree))
        new.tree<- phytools::bind.tip(new.tree, PUT, bp$length,
                                      bp$where, bp$position)
    }

    return(new.tree)
}


#' plot.custom.branch
#'
#' Function to help users visualize the candidate branches for tip insertion with \code{\link{custom_branch}} function
#'
#' @param tree "phylo" object used as backbone tree.
#' @param edges matrix with 5 character vector columns.
#'        Column 1 secifies the PUT to which the candidate branches in
#'        the row refer to.
#'        Columns 2 and 3 must represent two species whose MRCA represents the
#'        stem node of candidate branches and
#'        columns 4 and 5 two species whose MRCA define the crown node.
#'        Note that a pair of nodes may not necessarily define one single branch
#'        but a set of them.
#'        Two identical species will define a tree tip rather than an internal
#'        node.
#'        If the pairs of species in columns 2&3 and 4&5 are the same, the
#'        set of branches to select will be defined as the complete clade
#'        below their MRCA.
#' @param PUT If the \code{edges} data frame refers to more than one PUT,
#'            which one's set of branches must be plotted.
#' @param candidate.col Color to represent branches defined by the \code{edges}
#'                      data frame as candidates. Default value is red.
#' @param forbidden.col Color to represent branches not defined by the
#'                      \code{edges} data frame as candidates. Default value is
#'                      black.
#' @param candidate.lwd Line width to represent branches defined by the
#'                      \code{edges} data frame as candidates. Default value 2.
#' @param forbidden.lwd Line width to represent branches not defined by the
#'                      \code{edges} data frame as candidates. Default value 1.
#' @param ... Arguments to pass through \code{\link[ape]{plot.phylo}} function.
#'
#' @examplesIf interactive()
#' #Create a 'edges' dataframe
#' cats.edges <- data.frame(
#'  "PUT"= "Felis_catus",
#'  "parent1"= "Felis_silvestris",
#'  "parent2"= "Felis_silvestris",
#'  "child1"= "Felis_silvestris",
#'  "child2"= "Felis_silvestris")
#'
#' #Plot the tree highlighting candidate branches
#' plot_custom_branch(tree=cats, edges=cats.edges)
#'
#' @export
plot_custom_branch<- function(tree, edges, PUT=NULL,
                              candidate.col="#bf2828", forbidden.col="#3d3d3d",
                              candidate.lwd=2, forbidden.lwd=1,...){

    df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                     "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )

    if(length(unique(edges[,1]))==1 & is.null(PUT)){PUT<- unique(edges[,1])}
    if(!is.null(PUT)){
        if(!(PUT %in% edges[,1])){
            stop("The specified PUT is not contained in the 'edges' data frame.")
        }
        edges<- edges[edges[,1]==PUT]
    }
    if(length(unique(edges[,1]))>1){
        stop("Your 'edges' data frame contains more than one PUT. ",
            "Please specify only one PUT to plot its candidate branches.")
    }

    root<- findRoot(tree)
    edges[,2]<- gsub(" ", "_", edges[,2])
    edges[,3]<- gsub(" ", "_", edges[,3])
    edges[,4]<- gsub(" ", "_", edges[,4])
    edges[,5]<- gsub(" ", "_", edges[,5])

    permittednodes <- get_permitted_nodes_custom(tree, df, edges, root)

    col <- vector(mode = "character", length(tree$edge[,1]))
    col[1:length(col)] <- forbidden.col
    col[df$node%in%permittednodes]<- candidate.col

    lwd <- vector(mode = "character", length(tree$edge[,1]))
    lwd[1:length(lwd)] <- forbidden.lwd
    lwd[df$node%in%permittednodes]<- candidate.lwd

    return(ape::plot.phylo(tree, edge.color = col, edge.width = lwd, ...))
}

get_permitted_nodes_custom <- function(tree, df, edges, root){

    permittednodes<- as.numeric(NULL)

    for(i in 1:nrow(edges)){
        if(!all(c(edges[i,2],edges[i,3],edges[i,4],edges[i,5])%in%tree$tip.label)){
            message("Row ", i, " has species not included in the tree ",
                    "and will not be used.")
            next
        }

        if(edges[i,4]==edges[i,5]){
            basenode<- which(tree$tip.label==edges[i,4])
        }else{
            basenode<- ape::getMRCA(tree, c(edges[i,4], edges[i,5]))
        }
        if(edges[i,2]==edges[i,3]){
            parnode<- which(tree$tip.label==edges[i,2])
        }else{
            parnode<- ape::getMRCA(tree, c(edges[i,2], edges[i,3]))
        }
        if(edges[i,2]==edges[i,3] & edges[i,2]== edges[i,4] &
            edges[i,2]==edges[i,5]){

            parnode <- get_parent_siblings(tree,
                                          which(tree$tip.label==edges[i,2]))$parent
        }
        if(edges[i,2]==edges[i,4]|edges[i,2]==edges[i,5] &
            edges[i,3]==edges[i,4]|edges[i,3]==edges[i,5]){
            equal<-TRUE
        }else{
            equal<-FALSE
        }

        if(parnode==basenode & isFALSE(equal) ){
            message("Row ", i, " of 'edges' is not defining a phylogenetic branch.")
            next
        }
        if(parnode==basenode & isTRUE(equal) ){
            mrca<- ape::getMRCA(tree, c(edges[i,2], edges[i,3]))
            perm.nodes.i<-phytools::getDescendants(tree, mrca, curr=NULL)
            permittednodes<- c(permittednodes, perm.nodes.i)
            next
        }

        perm.nodes.i<- as.numeric(NULL)
        n<- basenode
        p<- parnode
        while(n!=p){
            if(n==root){
                message("Row ", i, " does not reflect a set of branches, ",
                        "so it will not be used.")
                perm.nodes.i<- NULL
                break
            }
            perm.nodes.i<- c(perm.nodes.i, n)
            n<-df$parent[df$node==n]
        }
        permittednodes<- c(permittednodes, perm.nodes.i)
    }

    if(length(permittednodes)==0){stop("No branches could be selected")}
    permittednodes<- unique(permittednodes)

    return(permittednodes)

}
