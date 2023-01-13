
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
#' Phylogenetic tree based in Jetz et al. (2012).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 9993 species of birds.
#' Data obtained from https://github.com/megatrees/bird_20221117
#' }
"birds"


#' plants
#'
#' Phylogenetic tree obtained from Jin & Quian (2022), based on Smith & Brown (2018) and Zanne et al. (2014).
#' Nomenclature stands for World Checklist of Vascular Plants (WCVP).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 73013 species of vascular plants.
#' Data obtained from https://github.com/megatrees/plant_20221214
#' }
"plants"


#' fish
#'
#' Phylogenetic tree obtained from Rabosky et al. (2018).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 31521 species of fish
#' Data obtained from https://github.com/megatrees/fish_20221117
#' }
"fish"


#' mammals
#'
#' Phylogenetic tree obtained from Upham et al. (2019).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 5911 species of mammals.
#' Data obtained from https://github.com/megatrees/mammal_20221117
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


#' amphibians
#'
#' Phylogenetic tree obtained from Jetz & Pyron, (2018).
#'
#' @format a Phylo object
#' \describe{
#' Phylogenetic tree comprising 7238 species of amphibians
#' Data obtained from https://github.com/megatrees/amphibian_20221117
#' }
"amphibians"


#' birds.info
#'
#' randtip 'info' file created after birds phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 9993 species of birds included in
#' 'birds' object. This object was mainly created after 'ncbi' database taxonomy.
#'
#' }
"birds.info"


#' plants.info
#'
#' randtip 'info' file created after plants phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 73013 species of vascular plants
#'  included in 'plants' object. This object was mainly created after 'ncbi' database taxonomy.
#'
#' }
"plants.info"


#' fish.info
#'
#' randtip 'info' file created after fish phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 31521 species of vascular fish
#'  included in 'fish' object. This object was created after The Fish Tree of Life taxonomy.
#'
#' }
"fish.info"

#' mammals.info
#'
#' randtip 'info' file created after mammals phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 5911 species of vascular mammals
#'  included in 'mammals' object. This object was mainly created after 'ncbi' database taxonomy.
#'
#' }
"mammals.info"


#' squamates.info
#'
#' randtip 'info' file created after squamates phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 4161 species of squamates
#'  included in 'squamates' object. This object was created after 'ncbi' database taxonomy.
#'
#' }
"squamates.info"


#' amphibians.info
#'
#' randtip 'info' file created after amphibians phylogeny
#'
#' @format a 'info' data frame
#' \describe{
#' Data frame with taxonomic information on the 7238 species of amphibians
#'  included in 'amphibians' object. This object was created after 'ncbi' database taxonomy.
#'
#' }
"amphibians.info"
