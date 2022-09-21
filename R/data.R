
#' cats phylogenetic tree
#'
#' Phylogenetic tree used for randtip package examples
#'
#' @format A phylogenetic tree with 14 tips and 13 internal nodes
#' \describe{
#' Phylogenetic tree regarding 14 species of felids (members of the Felidae family).
#' This tree has been extracted from a complete mammal tree by Upham et al. (2019) for example purposes
#' }
"cats"


#' mythology
#'
#' Data used in randtip tutorial (see Supplementary Material in https://doi.org/10.1101/2022.01.04.474924).
#'
#' @format A list with 5 slots
#' 
#' \describe{
#' \item{back.tree}{Imaginary backbone tree where PUTs will be bound.}
#' \item{sp.list}{List of species to include in the final expanded tree}
#' \item{info.list}{'info' file resulting after running build.info using mode="list"}
#' \item{info.backbone}{'info' file resulting after running build.info using mode="backbone"}
#' \item{edges}{'edges' file used for manual selection of branches where to bind the PUT "Monoceros_x_alaricornus"}
#' }
"mythology"