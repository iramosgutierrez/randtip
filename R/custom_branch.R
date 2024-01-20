#'
#' Bind PUTs to customized sets of phylogenetic edges
#'
#' @param tree backbone tree of class phylo.
#' @param edges matrix with five columns (character vectors). Each row defines 
#'              a set of consecutive edges between two nodes, and multiple rows
#'              can be used for the same PUT.
#'     Column 1 includes the PUTs.
#'     Columns 2 and 3 include a pair of species whose MRCA defines the older node, 
#'     whereas the MRCA of the species in columns 4 and 5 defines the younger node.
#'     Using the same species in columns 4 and 5 will define a terminal node as the
#'     younger node. If columns 2-3 are filled with a species and columns 4-5 are
#'     filled with another species, then all edges below the MRCA of the two species
#'     will be set as candidates. If the four columns are filled with the same
#'     species, the PUT will be bound as sister to this species.
#' @param rand.type "random" or "polytomy". Default is "random".
#' @param forceultrametric Whether or not to force the backbone tree to be ultrametric
#'                         (in case it is detected as non-ultrametric). Default is FALSE.
#' @param prob Whether or not the probability for selecting edges must be proportional
#'             to their length. If FALSE (default is TRUE), all edges have equal probability.
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
#' #Bind the PUT to any of the candidate edges as defined in cats.edges 
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
#' Auxiliar function to visualize candidate edges for PUT binding with \code{\link{custom_branch}}
#'
#' @param tree backbone tree of class phylo.
#' @param edges matrix with five columns (character vectors). Each row defines 
#'              a set of consecutive edges between two nodes, and multiple rows
#'              can be used for the same PUT.
#'     Column 1 includes the PUTs.
#'     Columns 2 and 3 include a pair of species whose MRCA defines the older node, 
#'     whereas the MRCA of the species in columns 4 and 5 defines the younger node.
#'     Using the same species in columns 4 and 5 will define a terminal node as the
#'     younger node. If columns 2-3 are filled with a species and columns 4-5 are
#'     filled with another species, then all edges below the MRCA of the two species
#'     will be set as candidates. If the four columns are filled with the same
#'     species, the PUT will be bound as sister to this species.
#' @param PUT If the \code{edges} data frame includes binding information for multiple
#'            PUTs, specifies the PUT whose candidate edges are to be depicted.
#' @param candidate.col Color for candidate edges (default is "red").
#' @param forbidden.col Color for non-candidate edges (default is "black").
#' @param candidate.lwd Line width for candidate edges (default is 2).
#' @param forbidden.lwd ine width for non-candidate edges (default is 1).
#' @param ... further arguments to be passed to \code{\link[ape]{plot.phylo}}.
#'
#' @examplesIf interactive()
#' cats.edges <- data.frame(
#'  "PUT"= "Felis_catus",
#'  "parent1"= "Felis_silvestris",
#'  "parent2"= "Felis_silvestris",
#'  "child1"= "Felis_silvestris",
#'  "child2"= "Felis_silvestris")
#' 
#' #Depict candidate edges on the phylogeny
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
            stop("The specified PUT is included in the 'edges' data frame.")
        }
        edges<- edges[edges[,1]==PUT]
    }
    if(length(unique(edges[,1]))>1){
        stop("Your 'edges' data frame contains more than one PUT. ",
            "Please specify only one PUT to plot candidate edges.")
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
                    "and will not be used.")              ###  ¿Y SI NO SE PUEDE DEFINIR UNO DE LOS NODOS?
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
            message("Row ", i, " of 'edges' is not defining a set of phylogenetic edges.")  ### CORRECT? Un único edge definido por una row es un caso particular, no la norma
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
                message("Row ", i, " does not define a set of edges, ",
                        "and it will not be used.")     ### ¿Y SI ES EL ÚNICO RAW QUÉ PASA?
                perm.nodes.i<- NULL
                break
            }
            perm.nodes.i<- c(perm.nodes.i, n)
            n<-df$parent[df$node==n]
        }
        permittednodes<- c(permittednodes, perm.nodes.i)
    }

    if(length(permittednodes)==0){stop("No candidate edges could be defined")}   # CHECK THAT THIS IS CORRECT
    permittednodes<- unique(permittednodes)

    return(permittednodes)

}
