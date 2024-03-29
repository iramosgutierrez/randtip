% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_input.R
\name{check_info}
\alias{check_info}
\title{Function to check randtip's files.}
\usage{
my.check <- check_info(info = my.info, tree = tree, sim = 0.85,
                   search.typos = TRUE, find.phyleticity = TRUE,
                   verbose = TRUE, parallelize = TRUE, ncores = 2 )
}
\arguments{
\item{info}{An 'info' object.}

\item{tree}{The original backbone tree.}

\item{sim}{Name similarity threshold to detect possible misspellings
on tip labels. Default value is 0.85. Similarity is obtained
with \code{stringsim} function from \code{stringdist} package.
See \link[stringdist]{stringsim} for details.}

\item{find.phyleticity}{Logical. Should or not the phyletic nature o the
taxonomic ranks be evaluated.}

\item{search.typos}{Logical. Should or not the function search for possible
misspelling on tip labels. This match lookup will be performed for
all PUTs using all tree tips.}

\item{verbose}{Logical. Should or not progress be printed.}

\item{parallelize}{Logical. If TRUE it allows the function to look for
phyletic status using multiple processing
cores.}

\item{ncores}{Number of cores to use in parallelization. If no number
is provided it defaults to all but one of system logical
cores.}
}
\value{
A data frame containing possible typographic errors,
taxonomic ranks extracted from 'info' and the phyletic
nature of each of them.
}
\description{
Function to account for the PUT status of the species in 'info',
spelling errors, putative MDCCs and the phyletic nature of
groups of PPCR species.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
catspecies <- c("Lynx_lynx", "Panthera_uncia",
"Panthera_onca", "Felis_catus", "Puma_concolor",
"Lynx_canadensis", "Panthera_tigris", "Panthera_leo",
"Felis_silvestris")

cats.info <- build_info(species=catspecies, tree= cats,
     find.ranks=TRUE, db="ncbi", mode="backbone")

cats.checked <- check_info(info=cats.info, tree=cats, sim=0.75, parallelize = FALSE)
\dontshow{\}) # examplesIf}
}
\author{
Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
}
