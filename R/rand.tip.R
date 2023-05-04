
#' Tree expansion function.
#'
#' Expand a phylogeny binding PUTs to a backbone tree.
#'
#' @param input An 'input' data frame obtained with \code{\link{info2input}} function.
#' @param tree A backbone tree.
#' @param rand.type For all PUTs not specified individually in 'input', which randomization type ("random" or
#'                  "polytomy") must be carried out. Default value is "random".
#' @param polyphyly.scheme For all PUTs not specified individually in 'input', which polyphyly
#'                         scheme ("largest", "complete" or "frequentist") must be used. Default value is "largest".
#' @param use.paraphyletic For all PUTs not specified individually in 'input', whether or not should paraphyletic
#'                         clades be taken into account or not. Default value is TRUE.
#' @param use.singleton For all PUTs not specified individually in 'input', should or not singleton MDCCs be
#'                      considered for binding as a sister species, or contrarily binding should be performed
#'                      anywhere below the parent node. Default value is TRUE.
#' @param use.stem For all PUTs not specified individually in 'input', whether or not should the stem branch be
#'                 considered as candidate for binding.  Default value is FALSE.
#' @param respect.mono For all PUTs not specified individually in 'input', whether or not monophyletic groups
#'                     should be respected when binding a PUT. Default value is TRUE.
#' @param respect.para For all PUTs not specified individually in 'input', whether or not paraphyletic groups
#'                     should be respected when binding a PUT. Default value is TRUE.
#' @param clump.puts For all PUTs not specified individually in 'input', whether or not co-ranked PUTs should be
#'                   clumped together in the phylogeny in case their taxonomic group is missing in the tree.
#'                   Will also clump conspecific PUTs. Default value is TRUE.
#' @param prob For all PUTs not specified individually in 'input', whether or not branch selection probability
#'             must be proportional to branch length or equiprobable. Default value is TRUE.
#' @param prune Whether or not the newly expanded tree will include only the species in the user's list.
#'              Default value is TRUE.
#' @param forceultrametric Whether or not the backbone tree will be forced to be ultrametric, only in case it is
#'                         not. Default value is FALSE.
#' @param verbose Whether or not to print information about the flow of the function. Default value is TRUE.
#'
#' @details The optimum binding procedure might be different in each case,
#'          given the phyletic nature of the MDCC to which the PUT should be bound to,
#'          and different reasons might drive the user to choose one or another.
#' \itemize{
#' \item{polyphyly.scheme} {The resulting outputs may be very different depending on the
#'                  selected parameter. If 'complete' procedure is selected, the binding
#'                  will be performed as a monophyletic group using the MRCA of
#'                  all the known PPCR species. The 'largest' procedure, otherwise,
#'                  will use as binding group the largest monophyletic group formed
#'                  by PPCR species. Ultimately, the 'frequentist' option will
#'                  split the group in monophyletic and singleton chunks and
#'                  for every PUT to be bound to the MDCC, one of the chunks will
#'                  be selected weighting the probability by the number of original tips.}
#'
#'
#' \item{use.paraphyletic}{ This parameter should be used as TRUE in case we are certain of
#'                  the paraphyletic nature of one group, but disregarded if the
#'                  "intruder" group within the otherwise monophyletic group is possibly
#'                  due to an incorrect placement in the creation of the phylogeny.}
#'
#' \item{forceultrametric} {It is important to note that forcing phylogenies to be ultrametric
#'                  in this way should not be taken as a formal statistical approach for inferring an
#'                  ultrametric tree but a method to be deployed whenever a genuinely ultrametric
#'                  phylogeny read from file fails due to issues related to numerical precision
#'                  (Revell, 2012). Thus, we strongly recommend the user to visually explore
#'                  phylogenetic trees that fail the ultrametricity test of check.info before
#'                   assuming the failure is due to numerical precision of computer machinery.}
#'}
#' @return An expanded phylogeny.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples
#'
#'  catspecies <- c("Lynx_lynx", "Panthera_uncia",
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
#'  tree=cats, rand.type = "random",
#'   forceultrametric = F)
#'
#' @export
rand_tip <- function(input, tree,rand.type = "random",
                    polyphyly.scheme="largest", use.paraphyletic=TRUE,use.singleton=TRUE, use.stem=FALSE,
                    respect.mono=TRUE, respect.para=TRUE, clump.puts = TRUE, prob=TRUE,
                    prune=TRUE, forceultrametric=FALSE, verbose = TRUE){

  #if(file.exists(input)){
  #  if(grep(getwd(), input)==1){filedir <-  input}else{
  #    filedir <- paste0(getwd(), "/", input)
  #  }

  #  cat(paste0("Reading input file from\n", filedir))
  #  input <- read.table(input)
  #}

    if(rand.type == "r"){rand.type <- "random"}
    if(rand.type == "p"){rand.type <- "polytomy"}

    if(polyphyly.scheme == "c"){polyphyly.scheme <- "complete"}
    if(polyphyly.scheme == "f"){polyphyly.scheme <- "frequentist"}
    if(polyphyly.scheme == "l"){polyphyly.scheme <- "largest"}

    if (!inherits(tree, "phylo")){
        stop("Backbone tree must be an object of class \"phylo\"")
    }
    if(!(rand.type %in% c("random", "polytomy"))){
        stop("Argument 'rand.type' must be \"random\" or \"polytomy\" ")
    }
    if(!(polyphyly.scheme %in% c("frequentist", "complete", "largest"))){
        stop("Argument 'polyphyly.scheme' must be \"frequentist\", \"complete\" or \"largest\" ")
    }

    if(length(tree$tip.label[duplicated(tree$tip.label)])>=1){
        stop("Tips ",
            paste0(tree$tip.label[duplicated(tree$tip.label)], collapse=", "),
            " are duplicated in the backbone tree. Please remove one of them.")
    }

    start<- Sys.time()

    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    input<- correct_DF(input)
    originalinput<- input
    input<- input[!is.na(input$MDCC),]
    input$taxon <- gsub(" ", "_", input$taxon)

    tree$tip.label <- gsub("_x_", "_x-", tree$tip.label)
    tree$tip.label <- gsub("_X_", "_x-", tree$tip.label)

    new.tree <- tree
    if(is.null(tree$edge.length)){new.tree$edge.length<-rep(1, nrow(new.tree$edge))}

    if(forceultrametric & !ape::is.ultrametric(new.tree)){new.tree<- phytools::force.ultrametric(new.tree, method = "extend")}
    if(isFALSE(forceultrametric) & !ape::is.ultrametric(new.tree)){
        message("The backbone tree is not ultrametric.",
                "\nPlease, set the argument 'forceultrametric' to TRUE if the tree is genuinely ultrametric.")
    }

    if(prune){
        if(length(input$taxon[input$MDCC=="Tip"])>0){
            trimming.species<- input$taxon[input$MDCC=="Tip"]
        }else{
            trimming.species<- as.vector(NULL)
        }

        for(using.mdcc in as.character(unique(input$MDCC))){

            if(using.mdcc=="Tip") next

            spp.df<-input[input$MDCC==using.mdcc,]
            using.rank<- as.character(notNA(unique(spp.df$MDCC.rank)))

            if(!(using.rank%in%names(input))){
                mdcc.genera<-first_word(spp.df[,c("taxon1","taxon2")])
            }else{
                mdcc.genera<-first_word(input$taxon[input[,using.rank]==using.mdcc])
            }
            mdcc.species<- new.tree$tip.label[first_word(new.tree$tip.label)%in%mdcc.genera]
            trimming.species<- c(trimming.species, mdcc.species)
            trimming.species<-trimming.species[trimming.species%in%new.tree$tip.label]
        }

        trimming.species<- notNA(trimming.species)
        new.tree <- ape::keep.tip(new.tree, trimming.species)
    }

    input[is.na(input$rand.type), "rand.type"]<-rand.type
    input[is.na(input$polyphyly.scheme), "polyphyly.scheme"]<-polyphyly.scheme
    input[is.na(input$use.paraphyletic) , "use.paraphyletic"] <- use.paraphyletic
    input[is.na(input$use.stem) , "use.stem"] <- use.stem
    input[is.na(input$clump.puts) , "clump.puts"] <- clump.puts
    input[is.na(input$respect.mono), "respect.mono"]<- respect.mono
    input[is.na(input$respect.para), "respect.para"]<- respect.para
    input[is.na(input$use.singleton), "use.singleton"]<- use.singleton
    input[is.na(input$prob) , "prob"] <- prob

    input$use.paraphyletic  <- as.logical(input$use.paraphyletic)
    input$use.singleton <- as.logical(input$use.singleton)
    input$use.paraphyletic <- as.logical(input$use.paraphyletic)
    input$use.stem <- as.logical(input$use.stem)
    input$clump.puts <- as.logical(input$clump.puts)
    input$respect.mono <- as.logical(input$respect.mono)
    input$respect.para <- as.logical(input$respect.para)
    input$prob <- as.logical(input$prob)

    input.bind<- input[!(input$taxon %in% new.tree$tip.label),]
    input.bind<- input.bind[!is.na(input.bind$MDCC),]


    manual.mdcc.taxa<-input.bind$taxon[!is.na(input.bind$taxon1)|!is.na(input.bind$taxon2)]
    rand.PUTs<- input.bind$taxon
    rand.PUTs<-sample(rand.PUTs, length(rand.PUTs), replace = F)
    rand.PUTs<- c(rand.PUTs[rand.PUTs%in%manual.mdcc.taxa], rand.PUTs[!(rand.PUTs%in%manual.mdcc.taxa)])

    for(i in seq_along(rand.PUTs)){
        PUT <- rand.PUTs[i]

        MDCC  <- inputfinder(input.bind,PUT, "MDCC")
        rank <- inputfinder(input.bind,PUT, "MDCC.rank")

        if(rank%in%randtip_ranks()){
            MDCC.type <- MDCC_phyleticity(input, new.tree,
                                          MDCC.info = list(rank=rank,MDCC=MDCC), trim=F)
        }else{
            MDCC.type <- rank
        }

        rand.type <- inputfinder(input.bind, PUT, "rand.type")

        use.singleton <- as.logical(inputfinder(input.bind, PUT, "use.singleton"))
        polyphyly.scheme<- as.character(inputfinder(input.bind,PUT, "polyphyly.scheme"))

        respect.mono <- as.logical(inputfinder(input.bind, PUT, "respect.mono"))
        respect.para <- as.logical(inputfinder(input.bind, PUT, "respect.para"))

        clump.PUT<-as.logical(inputfinder(input.bind, PUT, "clump.puts"))

        prob<-as.logical(inputfinder(input.bind, PUT, "prob"))


        if(isTRUE(clump.PUT)){
            clump<- bind_clump(new.tree, tree, input, PUT)
            if(!is.null(unlist(clump))){
                MDCC<-clump$MDCC
                rank<-clump$rank
                MDCC.type<-clump$MDCC.type
            }
        }

        if(verbose){
            cat(paste0(i, "/", length(rand.PUTs),
                       " (",round(i/length(rand.PUTs)*100, 2), "%). ",
                       "Binding ", PUT, " to ", MDCC ,"\n"))
        }


        perm.nodes<- NULL
        forbidden.nodes<-NULL
        if(rank=="Sister species"){
            perm.nodes <- which(new.tree$tip.label==MDCC)
        }

        if(rank=="Sister genus"){
            sister.genus.tips<- new.tree$tip.label[first_word(new.tree$tip.label)==MDCC]
            sister.genus.mrca<- ape::getMRCA(new.tree, sister.genus.tips)
            perm.nodes <- sister.genus.mrca
        }

        if(rank=="Manual setting"){
            sp1<- inputfinder(input.bind, PUT, "taxon1")
            sp2<- inputfinder(input.bind, PUT, "taxon2")
            clade.mrca<- ape::getMRCA(new.tree, c(sp1, sp2))

            perm.nodes <- phytools::getDescendants(new.tree, clade.mrca, curr=NULL)
            if(use.stem){perm.nodes <- c(clade.mrca, perm.nodes)}

        }

        if(rank=="species"){
            new.tree<-add_to_singleton(new.tree, clump$taxa, PUT, use.singleton=T)
            next
        }

        if(rand.type=="random"){
            if(is.null(perm.nodes)){
                perm.nodes<- get_permitted_nodes(new.tree, input,
                                                 MDCC, rank, MDCC.type,
                                                 polyphyly.scheme, use.paraphyletic,
                                                 use.singleton, use.stem=TRUE)

                forbidden.nodes<- get_forbidden_nodes(new.tree,input, MDCC, rank,
                                                      perm.nodes, respect.mono,
                                                      respect.para)

                if(all(perm.nodes[-1]%in%forbidden.nodes)){
                    perm.nodes <- perm.nodes[1]
                }else{
                    if(!use.stem){
                        perm.nodes<- perm.nodes[-1]
                        perm.nodes<-perm.nodes[!(perm.nodes%in%forbidden.nodes)]
                    }
                    if(use.stem){
                        perm.nodes<-perm.nodes[!(perm.nodes%in%forbidden.nodes)]
                    }
                }
                if(is.null(perm.nodes)){
                    perm.nodes<- get_permitted_nodes(new.tree, input, MDCC, rank, MDCC.type,
                                                     polyphyly.scheme, use.paraphyletic,
                                                     use.singleton, use.stem)
                }
            }

            adding.DF<- data.frame("parent"=new.tree$edge[,1], "node"=new.tree$edge[,2],
                                   "length"= new.tree$edge.length )
            adding.DF<- adding.DF[adding.DF$node %in% perm.nodes,]
            adding.DF$id<- 1:nrow(adding.DF)

            if(nrow(adding.DF)==1){node <- adding.DF$node}
            if(nrow(adding.DF) >1 & prob) {node<-sample(adding.DF$node, 1, prob = adding.DF$length)}
            if(nrow(adding.DF) >1 & !prob){node<-sample(adding.DF$node, 1)}

            bind.pos<- binding_position(new.tree, node,  insertion = "random",  prob)

            new.tree <- phytools::bind.tip(new.tree, PUT, edge.length = bind.pos$length,
                                            where = bind.pos$where , position = bind.pos$position )
        }
        if(rand.type=="polytomy"){

            if(is.null(perm.nodes)){
                perm.nodes<- get_permitted_nodes(new.tree, input, MDCC, rank, MDCC.type,
                                                 polyphyly.scheme, use.paraphyletic,
                                                 use.singleton, use.stem)
            }

            if(length(perm.nodes)==1){node <- get_parent_siblings(new.tree, perm.nodes)[[1]]}
            if(length(perm.nodes) >1){node <- ape::getMRCA(new.tree, perm.nodes)}


            bind.pos<- binding_position(new.tree, node,  insertion = "polytomy",  prob)

            new.tree <- phytools::bind.tip(new.tree, PUT, edge.length = bind.pos$length,
                                           where = bind.pos$where , position = bind.pos$position )

        }
    }

    complete.taxa.list <- originalinput$taxon[originalinput$keep.tip=="1"]
    complete.taxa.list.in.tree <- complete.taxa.list[complete.taxa.list %in% new.tree$tip.label]
    not.included <- complete.taxa.list[!(complete.taxa.list %in% complete.taxa.list.in.tree)]
    if(length(not.included) > 0){
        message("The following taxa were not bound to the tree:\n",
                paste0(not.included, "\n"))
    }

    if(isTRUE(prune)){new.tree <- ape::keep.tip(new.tree, complete.taxa.list.in.tree)}
    if(is.null(tree$edge.length)){new.tree$edge.length<-NULL}

    end <- Sys.time()

    if(verbose){
        cat(paste0("\n","\U2713", " PUT binding completed in ",
                  round(as.numeric(difftime(end, start,units = "mins")), 2), " mins\n"))
    }

    return(new.tree)
}


