
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
#' \item{info.list}{'info' file resulting after running build_info using mode="list"}
#' \item{info.backbone}{'info' file resulting after running build_info using mode="backbone"}
#' \item{edges}{'edges' file used for manual selection of branches where to bind the PUT "Monoceros_x_alaricornus"}
#' }
"mythology"


#' birds
#'
#' Phylogenetic tree obtained from Kimball et al. (2019).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 655 species of birds.
#' Data obtained from https://doi.org/10.3390/d11070109
#' }
"birds"


#' plants
#'
#' Phylogenetic tree obtained from Jin & Quian (2019), based on Smith & Brown (2018) and Zanne et al. (2014).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 74531 species of vascular plants.
#' Data obtained from https://github.com/jinyizju/V.PhyloMaker/tree/master/data
#' }
"plants"


#' fish
#'
#' Phylogenetic tree obtained from The Fish Tree of Life.
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 11638 species of birds.
#' Data obtained from https://fishtreeoflife.org/downloads/
#' }
"fish"


#' mammals
#'
#' Phylogenetic tree obtained from Upham et al. (2019).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 4099 species of mammals.
#' Data obtained from https://doi.org/10.1371/journal.pbio.3000494
#' }
"mammals"


#' squamates
#'
#' Phylogenetic tree obtained from Pyron et al. (2013).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 4161 species of squamates.
#' Data obtained from https://doi.org/10.1186/1471-2148-13-93
#' }
"squamates"


#' birds.info
#'
#' randtip 'info' file created after birds phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 655 species of birds included on
#' 'birds' object. This object was created after 'ncbi' database taxonomy.
#'
#' }
"birds.info"


#' plants.info
#'
#' randtip 'info' file created after plants phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 74531 species of vascular plants
#'  included on 'plants' object. This object was created after 'ncbi' database taxonomy.
#'
#' }
"plants.info"


#' fish.info
#'
#' randtip 'info' file created after fish phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 11638 species of vascular fish
#'  included on 'fish' object. This object was created after 'ncbi' database taxonomy.
#'
#' }
"fish.info"

#' mammals.info
#'
#' randtip 'info' file created after mammals phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 4099 species of vascular mammals
#'  included on 'mammals' object. This object was created after 'ncbi' database taxonomy.
#'
#' }
"mammals.info"


#' squamates.info
#'
#' randtip 'info' file created after squamates phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 4161 species of vascular squamates
#'  included on 'squamates' object. This object was created after 'ncbi' database taxonomy.
#'
#' }
"squamates.info"
