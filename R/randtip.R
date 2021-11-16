
#' randtip RANDOMIZATION FUNCTIONS
#' @export
rand.tip <- function(input, tree,rand.type = "random",
                    polyphyly.scheme="largest", use.paraphyletic=TRUE,use.singleton=TRUE, use.stem=FALSE,
                    respect.mono=TRUE, respect.para=TRUE, clump.puts = TRUE, prob=TRUE,
                    prune=TRUE, forceultrametric=FALSE, verbose = TRUE){

  if(rand.type == "r"){rand.type <- "random"}
  if(rand.type == "p"){rand.type <- "polytomy"}

  if(polyphyly.scheme == "c"){polyphyly.scheme <- "complete"}
  if(polyphyly.scheme == "f"){polyphyly.scheme <- "frequentist"}
  if(polyphyly.scheme == "l"){polyphyly.scheme <- "largest"}

  if (!inherits(tree, "phylo")) {stop("Backbone tree must be an object of class \"phylo\"")}
  if(!(rand.type %in% c("random", "polytomy"))) {stop("Argument 'rand.type' must be \"random\" or \"polytomy\" ")}
  if(!(polyphyly.scheme %in% c("frequentist", "complete", "largest"))) {stop("Argument 'polyphyly.scheme' must be \"frequentist\", \"complete\" or \"largest\" ")}

  if(length(tree$tip.label[duplicated(tree$tip.label)])==1){stop("Tip ", paste0(tree$tip.label[duplicated(tree$tip.label)], collapse=", "),  " is duplicated in the backbone tree.")}
  if(length(tree$tip.label[duplicated(tree$tip.label)])>1 ) {stop("Tips ",  paste0(tree$tip.label[duplicated(tree$tip.label)], collapse=", "), " are duplicated in the backbone tree.")}

    start<- Sys.time()

    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    input<- randtip::correct.DF(input)
    originalinput<- input
    input<- input[!is.na(input$MDCC),]
    input$taxon <- gsub(" ", "_", input$taxon)

    tree$tip.label <- gsub("_x_", "_x-", tree$tip.label)
    tree$tip.label <- gsub("_X_", "_x-", tree$tip.label)

    new.tree <- tree
    if(is.null(tree$edge.length)){new.tree$edge.length<-rep(1, nrow(new.tree$edge))}

    if(forceultrametric & !ape::is.ultrametric(new.tree)){new.tree<- phytools::force.ultrametric(new.tree)}
    if(isFALSE(forceultrametric) & !ape::is.ultrametric(new.tree)){
      message("The backbone tree is not ultrametric. \nPlease, set the argument 'forceultrametric' to TRUE if the tree is genuinely ultrametric.")}


    input.rand <- NULL
    input.poly <- NULL



    if(prune){
      if(length(input$taxon[input$MDCC=="Tip"])>0){
      trimming.species<- input$taxon[input$MDCC=="Tip"]}else{
        trimming.species<- as.vector(NULL)}

      for(using.mdcc in as.character(unique(input$MDCC))){

        if(using.mdcc=="Tip"){next}

        spp.df<-input[input$MDCC==using.mdcc,]
        using.rank<- as.character(randtip::notNA(unique(spp.df$MDCC.rank)))

        if(!(using.rank%in%names(input))){
          mdcc.genera<-randtip::firstword(spp.df[,c("taxon1","taxon2")]) }else{
          mdcc.genera<-randtip::firstword(input$taxon[input[,using.rank]==using.mdcc])}

        mdcc.species<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)%in%mdcc.genera]
        trimming.species<- c(trimming.species, mdcc.species)
        trimming.species<-trimming.species[trimming.species%in%new.tree$tip.label]
      }
        trimming.species<- randtip::notNA(trimming.species)
        new.tree <- ape::keep.tip(new.tree, trimming.species)}

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



    input.poly    <- input[input$rand.type=="polytomy",]
    input.rand    <- input[input$rand.type=="random",]



    #Phase 1. Random insertions
    if(!is.null(input.rand)){if(nrow(input.rand)>0){

        input.rand.bind<- input.rand[!(input.rand$taxon %in% new.tree$tip.label),]
        input.rand.bind<- input.rand.bind[!is.na(input.rand.bind$MDCC),]

        manual.mdcc.taxa<-input.rand.bind$taxon[!is.na(input.rand.bind$taxon1)|!is.na(input.rand.bind$taxon2)]
        rand.PUTs<- input.rand.bind$taxon
        rand.PUTs<-sample(rand.PUTs, length(rand.PUTs), replace = F)
        rand.PUTs<- c(rand.PUTs[rand.PUTs%in%manual.mdcc.taxa], rand.PUTs[!(rand.PUTs%in%manual.mdcc.taxa)])


        for(i in seq_along(rand.PUTs)){
            PUT <- rand.PUTs[i]

            MDCC  <- randtip::inputfinder(input.rand.bind,PUT, "MDCC")
            rank <- randtip::inputfinder(input.rand.bind,PUT, "MDCC.rank")
            if(rank%in%randtip::randtip_ranks()){
            MDCC.type <- randtip::MDCC.phyleticity(input, new.tree,
                                                   MDCC.info = list(rank=rank,MDCC=MDCC), trim=F)}else{
            MDCC.type <- rank
                                                   }

            use.singleton <- as.logical(randtip::inputfinder(input.rand.bind, PUT, "use.singleton"))
            polyphyly.scheme<- as.character(randtip::inputfinder(input.rand.bind,PUT, "polyphyly.scheme"))

            respect.mono <- as.logical(randtip::inputfinder(input.rand.bind, PUT, "respect.mono"))
            respect.para <- as.logical(randtip::inputfinder(input.rand.bind, PUT, "respect.para"))

            clump.PUT.i<-as.logical(randtip::inputfinder(input.rand.bind, PUT, "clump.puts"))

            prob<-as.logical(randtip::inputfinder(input.rand.bind, PUT, "prob"))




            if(isTRUE(clump.PUT.i)){
              clump<- randtip::bind.clump(new.tree, tree, input, PUT)
              if(!is.null(unlist(clump))){
                MDCC<-clump$MDCC
                rank<-clump$rank
                MDCC.type<-clump$MDCC.type


              }}


            if(verbose){
              cat(paste0("\n",i, "/", length(rand.PUTs),
                           " (",round(i/length(rand.PUTs)*100, 2), "%). ",
                           "Binding ", PUT, " to ", MDCC ,"\r")) }


            #Manual additions
            if(rank=="Sister species"){
              new.tree <- add.to.singleton(tree=new.tree, singleton = MDCC ,new.tips =  PUT)
              next
            }

            if(rank=="Sister genus"){
              sister.genus.tips<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)==MDCC]
              sister.genus.mrca<- ape::getMRCA(new.tree, sister.genus.tips)
              new.tree<-add.over.node(tree= new.tree, new.tip = PUT, node = sister.genus.mrca)
              next
              }

            if(rank=="Manual clade"){
            sp1<- randtip::inputfinder(input.rand.bind, PUT, "taxon1")
            sp2<- randtip::inputfinder(input.rand.bind, PUT, "taxon2")
            clade.mrca<- ape::getMRCA(new.tree, c(sp1, sp2))

            new.tree <- add.into.node(tree = new.tree, node = clade.mrca,
                                        new.tip = PUT, prob = prob )
            next
              }




            #Automatically searched MDCCs additions



            if(rank=="species"){
              new.tree<-add.to.singleton(new.tree, clump$taxa, PUT, use.singleton=T)
              next
            }
            perm.nodes<- get.permitted.nodes(new.tree, input, MDCC, rank, MDCC.type,
                         polyphyly.scheme, use.paraphyletic, use.singleton, use.stem=T)
            forbidden.nodes<- get.forbidden.nodes(new.tree,input, MDCC, rank, perm.nodes,
                                                  respect.mono, respect.para)

              if(all(perm.nodes[-1]%in%forbidden.nodes)){
                perm.nodes <- perm.nodes[1]
              }else{
                if (!use.stem){
                perm.nodes<- perm.nodes[-1]
                perm.nodes<-perm.nodes[!(perm.nodes%in%forbidden.nodes)]}
                if (use.stem){
                  perm.nodes<-perm.nodes[!(perm.nodes%in%forbidden.nodes)]}
              }
            if(is.null(perm.nodes)){
              perm.nodes<- get.permitted.nodes(new.tree, input, MDCC, rank, MDCC.type,
                          polyphyly.scheme, use.paraphyletic, use.singleton, use.stem)}

            adding.DF<- data.frame("parent"=new.tree$edge[,1], "node"=new.tree$edge[,2], "length"= new.tree$edge.length )
            adding.DF<- adding.DF[adding.DF$node %in% perm.nodes,]
            adding.DF$id<- 1:nrow(adding.DF)
            if(nrow(adding.DF)==1){node <- adding.DF$node}
            if(nrow(adding.DF) >1 & prob) {node<-sample(adding.DF$node, 1, prob = adding.DF$length)}
            if(nrow(adding.DF) >1 & !prob){node<-sample(adding.DF$node, 1)}

            bind.pos<- randtip::binding.position(new.tree, node,  insertion = "random",  prob)

            new.tree <- phytools::bind.tip(new.tree, PUT, edge.length = bind.pos$length,
                                           where = bind.pos$where , position = bind.pos$position )

        }
        }}
    #Polytomies
    if(nrow(input.poly)>0){

      input.poly<-input.poly[!(input.poly$taxon %in% new.tree$tip.label),]


      MDCCs<- randtip::notNA(unique(input.poly$MDCC))
      for(MDCCs.i in MDCCs){

        input.poly<- input.poly[!is.na(input.poly$MDCC),]
        MDCCs.i.rank<- unique(input.poly$MDCC.rank[input.poly[,"MDCC"]==MDCCs.i])
        if(!MDCCs.i.rank%in%randtip::randtip_ranks()){
          MDCCs.i.rank<- "MDCC"
          MDCC.type<- "Manual clade"
        }else{

        MDCC.type <- randtip::MDCC.phyleticity(input, new.tree,
                                               MDCC.info = list(rank=MDCCs.i.rank,MDCC=MDCCs.i), trim=F)
        }


        MDCC.taxa.toAdd <- input.poly$taxon[input.poly[,"MDCC"]==MDCCs.i]
        MDCC.taxa.ininput <- input$taxon[input[, MDCCs.i.rank]==MDCCs.i]
        MDCC.genera <- unique(randtip::firstword(MDCC.taxa.ininput))
        MDCC.genera <- randtip::notNA(MDCC.genera)
        MDCC.taxa.inTree<- randtip::sp.genus.in.tree(new.tree, MDCC.genera)



        if(verbose){
          cat(paste0(which(MDCCs==MDCCs.i), "/", length(MDCCs),
                     " (",round(which(MDCCs==MDCCs.i)/length(MDCCs)*100, 2), " %). ",
                     "Adding ", MDCCs.i ," (",
                     length(MDCC.taxa.toAdd)," tips).\n")) }


        if(length(MDCC.taxa.inTree)==1){
          new.tree<- randtip::polytomy.to.singleton(tree = new.tree, singleton = MDCC.taxa.inTree,
                                        new.tip =MDCC.taxa.toAdd, insertion = "long")}

        if(length(MDCC.taxa.inTree)>1) {
          MDCC.mrca<- ape::getMRCA(new.tree, MDCC.taxa.inTree)
          if(MDCC.type=="Polyphyletic"){
            if(length(unique(input.poly[input.poly$taxon%in%MDCC.taxa.toAdd, "polyphyly.scheme"]))==1){
              polyphyly.scheme<- unique(input.poly[input.poly$taxon%in%MDCC.taxa.toAdd, "polyphyly.scheme"])
              if(polyphyly.scheme=="largest"){

                mrca<- ape::getMRCA(new.tree, MDCC.taxa.inTree)
                descs<- phytools::getDescendants(new.tree, mrca,curr = NULL)
                forbidden.nodes<- randtip::get.forbidden.nodes(new.tree, input, MDCCs.i, MDCCs.i.rank, descs, respect.mono, respect.para)
                descs<- descs[!(descs%in%forbidden.nodes)]
                table<-sharingtaxa.descs(tree=new.tree, nodes=descs, MDCC.genera = MDCC.genera)
                table<- table[table$number==max(table$number),]
                table<- table[table$tot.number==min(table$tot.number),]
                node<-table[,"node"]
                MDCC.mrca<-node

              }
            }
          }
          new.tree<- randtip::polytomy.into.node(tree=new.tree,
                            new.tip =MDCC.taxa.toAdd, node =  MDCC.mrca)}
            }



    }




    complete.taxa.list <- originalinput$taxon[originalinput$keep.tip=="1"]
    complete.taxa.list.in.tree <- complete.taxa.list[complete.taxa.list %in% new.tree$tip.label]
    not.included <- complete.taxa.list[!(complete.taxa.list %in% complete.taxa.list.in.tree)]
    if(length(not.included) > 0){
        message("The following taxa were not bound to the tree:\n", paste0(not.included, "\n"))}

    if(isTRUE(prune)){new.tree <- ape::keep.tip(new.tree, complete.taxa.list.in.tree)}
    if(is.null(tree$edge.length)){new.tree$edge.length<-NULL}

    end <- Sys.time()


    if(verbose){
      cat(paste0("\n","\n","\U2713", " PUT binding completed in ",
                 round(as.numeric(difftime(end, start,units = "mins")), 2), " mins\n"))
    }
    return(new.tree)
}
