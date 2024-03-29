
#' Create a 'info' data frame.
#'
#' This function creates an 'info' object for a given list of species.
#'
#' @usage my.info <- build_info(species = species.list, tree = tree, db="ncbi",
#'                    mode="list", find.ranks = TRUE, interactive =FALSE,
#'                    genus = FALSE, prior.info = NULL, verbose = TRUE)
#'
#' @param species Character vector or a single-column data frame including
#'                the species of interest. Word breakers must be blanks (" ")
#'                or underscores ("_").
#' @param tree A 'phylo' object with the backbone tree (set to NULL if
#'             \code{mode} is "list").
#' @param find.ranks Logical. If TRUE, taxonomic information will be retrieved
#'                   from the specified taxonomic repository to identify
#'                   supra-generic MDCCs for the PUTs.
#' @param db Taxonomic repository to query if \code{find.ranks} is
#'           set to TRUE. One of 'ncbi' (default), 'itis', 'gbif'
#'           or 'bold'.
#' @param mode If mode is "list", the info data frame will be filled with
#'             the species provided in the \code{species} argument. If mode
#'             is "backbone", the info data frame will also include all the
#'             tips in the backbone tree.
#' @param interactive Logical. Whether or not ambiguous species names will
#'                    be resolved manually as they appear when retrieving
#'                    taxonomic information. If FALSE, NAs will be returned
#'                    in the corresponding row of info.
#' @param genus Logical. Whether or not the backbone tree is resolved to the
#'              genus level. If TRUE, all the taxa in the tips of the phylogeny
#'              and the species vector must represent genera.
#' @param prior.info A previously created info data frame that is recycled to
#'                   build the new info.
#' @param verbose Logical. Whether or not to print the progress.
#'
#' @return An 'info' data frame
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examplesIf interactive()
#' # Create a list of species to include in the phylogeny
#'  catspecies <- c("Lynx_lynx",
#' "Panthera_uncia",
#' "Panthera_onca",
#' "Felis_catus",
#' "Puma_concolor",
#' "Lynx_canadensis",
#' "Panthera_tigris",
#' "Panthera_leo",
#' "Felis_silvestris")
#'
#' #Create the 'info' data frame
#' cats.info <- build_info(species=catspecies, tree= cats,
#'      find.ranks=TRUE, db="ncbi", mode="backbone")
#' @export
build_info<- function(species, tree=NULL, find.ranks=TRUE, db="ncbi",mode="backbone",
                      interactive=FALSE, genus=FALSE, prior.info=NULL, verbose = T){


    if(is.data.frame(species)){
        if(ncol(species)!=1){
            stop("Species must be provided as a character vector ",
                 "or single-column data frame.")
        }else{species <- species[,1]}
    }
    if(!(is.vector(species))){
        stop("Species must be provided as a character vector ",
             "or single-column data frame.")
    }
    species <- remove_spaces(species)
    duplicated_sp <-  species[duplicated(species)]
    species <- gsub("_x_|_X_", "_x-", species)

    if(length(duplicated_sp)>=1){
        stop("Taxa ", paste0(duplicated_sp, collapse=", ")," are duplicated.")
    }
    if(is.null(tree) & mode=="backbone"){
        stop("Parameter 'mode' is set to 'backbone', ",
            "but the backbone tree is missing. ",
            "Please, provide a backbone tree.")
    }
    if(!is.null(tree) & !inherits(tree, "phylo")){
        stop("Backbone tree must be an object of class 'phylo'.")
    }

    if(!(mode %in% c("list", "backbone"))){
        stop("Parameter 'mode' must be 'list' or 'backbone' ")
    }
    if(!(db %in% c("ncbi", "itis", "gbif", "bold"))){
        stop(paste0(db, " is not a supported taxonomic repository."))
    }

    tree$tip.label <- gsub("_x_|_X_", "_x-", tree$tip.label)
    spp.in.tree<- tree$tip.label

    only.genus<- !grepl("_", species)
    only.tree.genus<- !grepl("_", spp.in.tree)
    if(isTRUE(genus)){
      if(!all(only.genus)){
        stop("Taxa must represent genera only for \"genus\" mode")
      }
      if(!all(only.tree.genus)){
        stop("Phylogenetic tips must represent genera only for \"genus\" mode")
      }

    }else{
    # Put suffix _sp in taxa with only the genus
    for(t in which(only.genus)){
      taxon.suffix.i <- "_sp."
      i <- 2
      while(paste0(species[t], taxon.suffix.i) %in% tree$tip.label){
        taxon.suffix.i <- paste0("_sp", i, ".")
        i <- i + 1
      }
      species[t]<-paste0(species[t], taxon.suffix.i)
    }
    }

     spp.original<- species

    if(any(first_word(spp.in.tree)=="X")|any(first_word(spp.in.tree)=="x")){
        tree$tip.label[first_word(tree$tip.label)=="x"] <-
            gsub("x_", "X-", tree$tip.label[first_word(tree$tip.label)=="x"])
        tree$tip.label[first_word(tree$tip.label)=="X"] <-
            gsub("X_", "X-", tree$tip.label[first_word(tree$tip.label)=="X"])
    }

    if(mode=="backbone"){
         species <- c(species, spp.in.tree[!(spp.in.tree%in%species)])
    }

    names_df <- c("taxon", randtip_ranks(),
                  "taxon1", "taxon2", "rand.type", "polyphyly.scheme",
                  "use.paraphyletic", "use.singleton","use.stem",
                  "respect.mono", "respect.para","clump.puts",
                  "prob","keep.tip")
    info<- as.data.frame(matrix(nrow = length(species),
                                ncol = length(names_df),
                                dimnames = list(NULL, names_df)))
    info$taxon <- species


    if(!is.null(prior.info)){
      if(!all(names(info)==names(prior.info))){
        stop("Invalid column names for prior.info. They should match column names of an info data frame.")

        }else{

          #species in prior must be changed to be same
          prior.info$taxon <- remove_spaces(prior.info$taxon)

          spp.in.prior <- info$taxon[info$taxon %in% prior.info$taxon]

          spp.gen.in.prior <- info$taxon[!(info$taxon%in%prior.info$taxon) &
                                        first_word(info$taxon) %in% prior.info$genus]
          genera.in.prior <- unique(first_word(spp.gen.in.prior))


          if(length(spp.in.prior)>0){
            prior.info.cut <- prior.info[prior.info$taxon%in%spp.in.prior,]
            v.spp<- order(spp.in.prior)

            prior.info.cut <-prior.info.cut[v.spp,]

            info[info$taxon%in% prior.info.cut$taxon,]<- prior.info.cut



            ord1 <- match(spp.original, info$taxon)
            ord2 <- match(info$taxon[!(info$taxon%in%spp.original)], info$taxon)
            info <- info[c(ord1,ord2),]

            }
            #fill in taxonomic information with that included in prior.info  #  REMOVE THIS

          if(length(genera.in.prior)>0){
            for(gen in genera.in.prior){
            gen.sp <- spp.gen.in.prior[first_word(spp.gen.in.prior)==gen]
            prior.info.cut <- prior.info[first_word(prior.info$taxon)==gen,1:9]

            if(nrow(unique(prior.info.cut[,2:9]))==1){
              info[info$taxon%in%gen.sp, 2:9] <- unique(prior.info.cut[,2:9])
            }


          }}



        }

    }#include information from prior.info

    info$genus<- first_word(info$taxon)
    genera <- unique(info$genus)
    if(!is.null(prior.info)){
      unmatched <- info[which(rowSums(is.na(info[,3:9]))==7),]
      genera <- unique(first_word(unmatched$taxon))
    }#search only for unmatched genera  #  REMOVE THIS

    cols.select <- c("taxon1", "taxon2","rand.type", "polyphyly.scheme",
                    "use.paraphyletic", "use.singleton",
                    "use.stem","respect.mono","respect.para",
                    "clump.puts", "prob" )


    if(find.ranks & length(genera)>0){
        info <- search_taxize(info, genera, interactive, db, verbose=verbose)
    }

    info[!(species %in% spp.original), cols.select] <- "-"
    info$keep.tip <- ifelse(species %in% spp.original, 1, 0)
    for(rank in randtip_ranks()){
      info[,rank] <- gsub(" ", "_", info[,rank])
    }

    nonfoundtaxa<-info[is.na(info$subtribe)&is.na(info$tribe)&is.na(info$subfamily)&
                       is.na(info$family)&is.na(info$superfamily)&is.na(info$order)&
                       is.na(info$class),]
    if(nrow(nonfoundtaxa)>0 & isTRUE(find.ranks)){
        nonfoundgenera <- unique(first_word(nonfoundtaxa$taxon))
        message(paste0("The following genera were detected as ambiguous or missing. Please, consider checking them manually.\n"),
                paste0(nonfoundgenera, "\n"))
    }


    return(info)
}



#' Function to check randtip's files.
#'
#' Function to account for the PUT status of the species in 'info',
#' spelling errors, putative MDCCs and the phyletic nature of
#' groups of PPCR species.
#'
#' @usage my.check <- check_info(info = my.info, tree = tree, sim = 0.85,
#'                    search.typos = TRUE, find.phyleticity = TRUE,
#'                    verbose = TRUE, parallelize = TRUE, ncores = 2 )
#'
#' @param info An 'info' object.
#' @param tree The original backbone tree.
#' @param search.typos Logical. Should or not the function search for possible
#'            misspelling on tip labels. This match lookup will be performed for
#'            all PUTs using all tree tips.
#' @param sim Name similarity threshold to detect possible misspellings
#'            on tip labels. Default value is 0.85. Similarity is obtained
#'            with \code{stringsim} function from \code{stringdist} package.
#'            See \link[stringdist]{stringsim} for details.
#' @param find.phyleticity Logical. Should or not the phyletic nature o the
#'            taxonomic ranks be evaluated.
#' @param verbose Logical. Should or not progress be printed.
#' @param parallelize Logical. If TRUE it allows the function to look for
#'                             phyletic status using multiple processing
#'                             cores.
#' @param ncores Number of cores to use in parallelization. If no number
#'                is provided it defaults to all but one of system logical
#'                cores.
#'
#' @return A data frame containing possible typographic errors,
#'         taxonomic ranks extracted from 'info' and the phyletic
#'         nature of each of them.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examplesIf interactive()
#' catspecies <- c("Lynx_lynx", "Panthera_uncia",
#' "Panthera_onca", "Felis_catus", "Puma_concolor",
#' "Lynx_canadensis", "Panthera_tigris", "Panthera_leo",
#' "Felis_silvestris")
#'
#' cats.info <- build_info(species=catspecies, tree= cats,
#'      find.ranks=TRUE, db="ncbi", mode="backbone")
#'
#' cats.checked <- check_info(info=cats.info, tree=cats, sim=0.75, parallelize = FALSE)
#'
#' @export
check_info<- function(info, tree, sim=0.85, find.phyleticity=TRUE,search.typos =TRUE,
                      verbose=TRUE, parallelize = TRUE, ncores = NULL){

    #if(file.exists(info)){

    #  if(grep(getwd(), info)==1){filedir <-  info}else{
    #    filedir <- paste0(getwd(), "/", info)
    #  }

    #  cat(paste0("Reading info file from\n", filedir))
    #  info <- read.table(info)
    #  }

    if(parallelize){
        if(is.null(ncores)){
            cat("\nncores argument was not provided.",
                "Using all system cores minus one.\n\n")
            ncores <- parallel::detectCores(logical = TRUE) - 1
        }else if(ncores > (parallel::detectCores(logical = TRUE) - 1)){
            cat("\nNumber of cores not availble.",
                "Using all system cores minus one.\n\n")
            ncores <- parallel::detectCores(logical = TRUE) - 1
        }
    }



    if(is.null(info)){stop("Data frame 'info' is missing.")}
    if(is.null(tree)){stop("Backbone tree is missing.")}

    info <- correct_DF(info)
    info$keep.tip[is.na(info$keep.tip)] <- "1"
    if(all(info$keep.tip != "1")){
        stop("No species in info with keep.tip equal to '1'.",            # THIS IS WEIRD. THE USER SHOULD NOT EDIT THAT COLUMN. I WOULD SIMPLY SAY SOMETHING LIKE "Cannot drop all species from the tree. At least some of them must have keep.tip = 1 status"
             " Please set keep.tip to '1' for every species to ", #
            "keep in the final tree.")
    }

    info$taxon<- remove_spaces(info$taxon)
    info.taxa <- info$taxon

    tree$tip.label<- remove_spaces(tree$tip.label)
    tree.taxa<- tree$tip.label

    DF <- info[info$keep.tip=="1", c("taxon",randtip_ranks())]
    DF$PUT.status <- NA
    DF$Typo<- FALSE
    DF$Typo.names<- NA

    DF$PUT.status <- ifelse(DF$taxon %in% tree.taxa, "Tip", "PUT")

    # Look for name similarities
    if(search.typos){
        PUTs <- which(DF$PUT.status == "PUT")

        cat(paste0("Searching for possible misspellings\n"))
        putlength <- length(which(DF$PUT.status == "PUT"))

        if(parallelize){
            # Create cluster needed for both search_typos and check_phyletic function
            cl <- parallel::makeCluster(ncores)
            parallel::clusterExport(cl, c("check_phyletic", "notNA",
                                          "MDCC_phyleticity", "phyleticity",
                                          "first_word", "search_typos"),
                                    envir = environment())
            parallel::clusterEvalQ(cl, library(stringdist))
            DF_out <- parallel::parLapply(cl, PUTs,
                                          function(PUT_i, PUTs, putlength, DF, tree.taxa,
                                                   sim, verbose){
                                              search_typos(PUT_i, PUTs, putlength, DF, tree.taxa,
                                                           sim, verbose)
                                          }, PUTs, putlength, DF, tree.taxa, sim, verbose)
            for(i in seq_along(PUTs)){
                DF[PUTs[i], ] <- DF_out[[i]]
            }

        }else{
            for(PUT_i in PUTs){
                DF[PUT_i, ] <- search_typos(PUT_i, PUTs, putlength, DF, tree.taxa, sim, verbose)
            }
        }

    }

    # Taxonomy lookup:
    ranks<-randtip_ranks()
    if(parallelize){
        if(!(search.typos)){
            cl <- parallel::makeCluster(ncores)
            parallel::clusterExport(cl, c("check_phyletic", "notNA",
                                          "MDCC_phyleticity", "phyleticity",
                                          "first_word"),
                                    envir = environment())
        }
        if(verbose & find.phyleticity){cat("Checking phyletic status in parallel.\n")}
        DF_out <- parallel::parLapply(cl, 1:length(ranks),
                                     function(rank_i, ranks, DF, info, tree,
                                              find.phyleticity, verbose){
                                         check_phyletic(ranks, rank_i,
                                                        DF, info, tree,
                                                        find.phyleticity,
                                                        verbose)
                                     }, ranks, DF, info, tree, find.phyleticity,
                                     verbose)
        parallel::stopCluster(cl)

        for(rank_i in seq_along(ranks)){
            DF_col <- paste0(ranks[rank_i], "_phyletic.status")
            DF[[DF_col]] <- DF_out[[rank_i]]
        }

    }else{
        for(rank_i in seq_along(ranks)){
            DF_col <- paste0(ranks[rank_i], "_phyletic.status")
            DF[[DF_col]] <- check_phyletic(ranks, rank_i,
                                           DF, info, tree,
                                           find.phyleticity,
                                           verbose)
        }
    }

    if(length(DF$Typo[DF$Typo==TRUE])>0){
        message("There might be misspelling errors in ",
                "the species list or the phylogenetic tips. ",
                "Please, check the TYPO column in the outputted data frame.\n")
    }

    DF<-DF[,c("taxon", "PUT.status", "Typo", "Typo.names","genus",
              "genus_phyletic.status", "subtribe" ,
              "subtribe_phyletic.status","tribe" , "tribe_phyletic.status",
              "subfamily","subfamily_phyletic.status","family",
              "family_phyletic.status", "superfamily",
              "superfamily_phyletic.status", "order","order_phyletic.status",
              "class","class_phyletic.status")]

    # Check tips
    tips<- tree$tip.label
    if(length(tips[duplicated(tips)])>=1){
        message("Tips ", tips[duplicated(tips)],
                " are duplicated in the phylogeny tips.",
                " Please remove one of them.\n")
    }

    subsp.tips<- tips[sapply(strsplit(tips, "_"), length)>2]
    for(ssp in subsp.tips){
        nomials<-strsplit(ssp, split="_")[[1]]
        if(paste0(nomials[1], "_", nomials[2])%in%tips &
              any(nomials[3:length(nomials)]==nomials[2])){

            message("Tips ", ssp, " and " , paste0(nomials[1], "_", nomials[2]),
                    " may represent the same taxon. Please consider ",
                    "removing one of them.\n" )
        }
    }

    # Tree ultrametricity evaluation
    if(!ape::is.ultrametric(tree)){
        message("The backbone tree is not ultrametric.")
    }

    return(DF)
}



#' Convert 'info' to 'input'.
#'
#' Convert an 'info' data frame into an 'input' object.
#'
#' @param info A definitive info data frame.
#' @param tree Backbone tree.
#' @param parallelize Logical. If TRUE it allows the function to look for
#'                             phyletic status using multiple processing
#'                             cores.
#' @param ncores Number of cores to use in parallelization. If no number
#'                is provided it defaults to all but one of system logical
#'                cores.
#' @param verbose Logical. Whether or not to print the progress.
#'
#' @return An input data frame to feed the \code{rand_tip} function.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examplesIf interactive()
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
#' @export
info2input<- function(info, tree, parallelize = T, ncores = NULL, verbose=T){

    if(parallelize){
        if(is.null(ncores)){
            cat("\nncores argument was not provided.",
                "Using all but one of system cores.\n\n")
            ncores <- parallel::detectCores(logical = TRUE) - 1
        }else if(ncores > (parallel::detectCores(logical = TRUE) - 1)){
            cat("\nNumber of cores not availble.",
                "Using all system cores but one.\n\n")
            ncores <- parallel::detectCores(logical = TRUE) - 1
        }
    }

    input.to.mdcc <- input_to_MDCCfinder(info, tree)
    input <- input.to.mdcc$input
    tree <- input.to.mdcc$tree
    taxon.in.tree <- input.to.mdcc$taxon.in.tree

    if(parallelize){
        if(verbose){cat("Searching MDCCs in parallel\n")}

        taxon = input$taxon[!(taxon.in.tree)]
        silent=T

        cl <- parallel::makeCluster(ncores)
        parallel::clusterExport(cl, c("usingMDCCfinder", "correct_DF",
                                    "randtip_ranks", "first_word",
                                    "sp_genus_in_tree"),
                                envir = environment())

        input_out <- parallel::parLapply(cl, taxon,
                                         function(taxon_i, input,tree, silent){
                                             usingMDCCfinder(input = input,
                                                             taxon = taxon_i,
                                                             tree = tree,
                                                             silent = T)
                                       }, input, tree, silent)
        parallel::stopCluster(cl)

        for(i in seq_along(taxon)){
            pos <- which(input$taxon==taxon[i])
            input[pos,"MDCC" ] <- input_out[[i]][[1]]
            input[pos,"MDCC.rank" ] <- input_out[[i]][[2]]
        }

    }else{
        input_search<- usingMDCCfinder(input = input,
                                       taxon = input$taxon[!(taxon.in.tree)],
                                       tree = tree,
                                       silent = !verbose)

        input$MDCC[!(taxon.in.tree)] <- input_search[[1]]
        input$MDCC.rank[!(taxon.in.tree)] <- input_search[[2]]
    }

    # Taxa with no MDCC
    not.included <- input[is.na(input$MDCC),]
    if(length(not.included$taxon) > 0){
        message("The following taxa were not assigned MDCC and will not be ",
                "bound to the tree:\n",
                paste0(not.included$taxon, "\n"))
    }

    return(input)
}


search_taxize <- function(info, genera, interactive, db, verbose=T){
    searching.categories<- randtip_ranks()[-1]

    for(i in 1:length(genera)){
        tryCatch({
            if(interactive){
                search <- suppressMessages(taxize::classification(as.character(genera[i]),
                                                                  db = db))[[1]]
            }else if (db!="itis"){
                out<-utils::capture.output(suppressMessages(
                  search <- taxize::classification(as.character(genera[i]),
                                                   db = db, rows=Inf)[[1]]))
            }else{
              out<-utils::capture.output(suppressMessages(
                search <- taxize::classification(as.character(genera[i]),
                                                 db = db)[[1]]))
            }
            for(cat in searching.categories){
                if(length(search[which(search$rank==cat), "name"])==0){
                    info[info$genus==genera[i], cat]<-NA
                }else{
                    cats<-search[which(search$rank==cat), "name"]
                    if(length(cats)>1){cats<-cats[1]}
                    info[info$genus==genera[i], cat]<- cats}
            }
          }, error=function(e){
              # Assign NA to fetching errors
              info[info$genus==genera[i], searching.categories] <- NA
        })

        # Avoid ip blocks. Taxize allows only 3 searches per second.
        Sys.sleep(0.33)

        if(!interactive & verbose){
            if(i==1 & verbose){
                cat(paste0("Retrieving taxonomic information from ", db, " database.\n",
                          "0%       25%       50%       75%       100%", "\n",
                          "|---------|---------|---------|---------|", "\n"))
            }

            v<- seq(from=0, to=40, by=40/length(genera))
            v<- diff(ceiling(v))
            cat(strrep("*", times=v[i]))

            if(i ==length(genera)){cat( "*\n")}
        }
    }

    return(info)

}

# Creates an input for MDCC finder and facilitates unit testing of the
# usingMDCCfinder function (utils.R source file)
input_to_MDCCfinder <- function(info, tree){

    input<-info
    input[is.na(input$keep.tip), "keep.tip"]<-"1"
    input<- correct_DF(input)
    input$taxon <- gsub(" ", "_", input$taxon)

    if(any(duplicated(input$taxon))){
        duptax <- unique(input$taxon[duplicated(input$taxon)])

        lines <- which(input$taxon%in%duptax & input$keep.tip=="0")
        input <- input[-lines,]
        rm(lines)
    }

    if(any(duplicated(input$taxon))){
        lines <- which(duplicated(input$taxon))
        input <- input[-lines,]
    }

    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    tree$tip.label <- gsub("_x_|_X_", "_x-", tree$tip.label)

    col.args <- c("rand.type", "polyphyly.scheme")
    col.args.logical <- c("use.paraphyletic", "use.singleton", "use.stem",
                          "respect.mono", "respect.para", "clump.puts", "prob")

    rand.type <- c("random", "polytomy", "-")
    polyphyly.scheme <- c("complete", "largest", "frequentist", "-")
    logical.args <- c("TRUE", "FALSE", "-")
    for(arg in c(col.args, col.args.logical)){
        if(arg %in% col.args){
            allowed.args <- get(arg)
        }else{
            allowed.args <- logical.args
        }

        if(all(notNA(info[[arg]]) %in% allowed.args)){next}

        errortaxa<- info[!is.na(info[[arg]]),]
        errortaxa<- errortaxa$taxon[!(errortaxa[[arg]] %in% allowed.args)]

        if(arg %in% col.args.logical){
            stop(paste0("\"", arg, "\" argument must be logical. ",
                        "Please check your info at the following taxa:\n",
                        paste0(errortaxa, collapse = "\n")))
        }else{
            allowed.args.str <- paste(allowed.args, collapse = "' or '")
            stop(paste0("\"", arg, "\" argument must be \'",
                        allowed.args.str, "'. ",
                        "Please check your info at the following taxa:\n",
                        paste0(errortaxa, collapse = "\n")))
        }
    }

    input$MDCC <- as.character(NA)
    input$MDCC.rank <- as.character(NA)

    taxon.in.tree <- (input$taxon %in% tree$tip.label)
    input$MDCC[taxon.in.tree] <- "Tip"
    input$MDCC.rank[taxon.in.tree] <- "Tip"

    return(list(input = input, tree = tree, taxon.in.tree = taxon.in.tree))
}


# Check phyletic status
check_phyletic <- function(ranks, rank_i, DF, info, tree, find.phyleticity,
                           verbose){

    rank <- ranks[rank_i]

    phyletic_status <- rep(NA, times = nrow(DF))

    groups<- notNA(unique(DF[,rank]))

    if(verbose & length(groups) & isTRUE(find.phyleticity)>0){
        cat(paste0("Checking phyletic status at ", rank, " level...\n"))

        cat(paste0("0%       25%       50%       75%       100%", "\n",
                   "|---------|---------|---------|---------|", "\n"))
    }
    for(group in groups){
        if(isTRUE(find.phyleticity)){
            phyle.type<- MDCC_phyleticity(info, tree,
                                        MDCC.info = list("rank"= rank,
                                                         "MDCC"= group))
        }else{
            phyle.type <- "unknown"
        }

        phyletic_status[which(DF[,rank]==group)] <-phyle.type

        if(verbose & find.phyleticity){
          v<- seq(from=0, to=40, by=40/length(groups))
          v<- diff(ceiling(v))
          cat(strrep("*", times=v[which(groups==group)]))

          if(group == groups[length(groups)]){cat("*\n")}
        }
    }

    return(phyletic_status)
}

search_typos <- function(PUT_i, PUTs, putlength, DF, tree.taxa, sim, verbose){

    tax<-DF$taxon[PUT_i]
    sim.search<-tree.taxa[stringdist::stringsim(tree.taxa,tax)>sim]
    if(length(sim.search)>0){
        DF$Typo[PUT_i]<- TRUE
        sim.search<-paste0(sim.search, collapse = " / ")
        DF$Typo.names[PUT_i]<- sim.search
    }
    if(verbose){
        if (PUT_i == PUTs[1]) {
            cat(paste0("0%       25%       50%       75%       100%",
                        "\n", "|---------|---------|---------|---------|",
                        "\n"))
        }
        v <- seq(from = 0, to = 40, by = 40/putlength)
        v <- diff(ceiling(v))
        pos <- which(PUTs==PUT_i)
        cat(paste0(strrep("*", times = v[pos])))
        if (PUT_i == PUTs[putlength]) {
            cat("*\n")
        }
    }

    return(DF[PUT_i, ])

}

