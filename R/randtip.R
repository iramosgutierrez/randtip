
#' randtip RANDOMIZATION FUNCTIONS
#' @export
rand.list <- function(tree, DF1,
                    rand.type = "random",agg.ssp = FALSE, poly.ins="large",
                    resp.mono=FALSE, resp.para=FALSE, resp.sing=FALSE,
                    prob = TRUE, verbose = FALSE, prune=TRUE, forceultrametric=F){
  if (!inherits(tree, "phylo")) {stop("object \"tree\" is not of class \"phylo\"")}
  if(!(rand.type %in% c("random", "polytomy"))) {stop("rand.type must be \"random\" or \"polytomy\" ")}
  if(!(poly.ins %in% c("freq", "all", "large"))) {stop("poly.ins must be \"freq\", \"all\" or \"large\" ")}

    start<- Sys.time()

    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    DF1<- randtip::correct.DF(DF1)
    originalDF1<- DF1
    DF1<- DF1[!is.na(DF1$using.MDCC),]
    DF1$taxon <- gsub(" ", "_", DF1$taxon)

    tree$tip.label <- gsub("_x_", "_x-", tree$tip.label)
    tree$tip.label <- gsub("_X_", "_x-", tree$tip.label)

    new.tree <- tree

    if(forceultrametric & !ape::is.ultrametric(new.tree)){new.tree<- phytools::force.ultrametric(new.tree)}
    if(isFALSE(forceultrametric) & !ape::is.ultrametric(new.tree)){
      message("Specified tree is not ultrametric. \nTo force the randomization as an ultrametric tree plase set forceultrametric=TRUE")}

    DF1.rand <- NULL
    DF1.poly <- NULL
    DF1.dupl <- NULL


    if(prune){
      if(length(DF1$taxon[DF1$using.MDCC=="Tip"])>0){
      trimming.species<- DF1$taxon[DF1$using.MDCC=="Tip"]}else{
        trimming.species<- as.vector(NULL)}

      for(using.mdcc in as.character(unique(DF1$using.MDCC))){

        if(using.mdcc=="Tip"){next}

        spp.df<-DF1[DF1$using.MDCC==using.mdcc,]
        using.level<- as.character(randtip::notNA(unique(spp.df$using.MDCC.lev)))

        if(!(using.level%in%names(DF1))){
          mdcc.genera<-randtip::firstword(spp.df[,c("taxon1","taxon2")]) }else{
          mdcc.genera<-randtip::firstword(DF1$taxon[DF1[,using.level]==using.mdcc])}

        mdcc.species<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)%in%mdcc.genera]
        trimming.species<- c(trimming.species, mdcc.species)
      }
        trimming.species<- randtip::notNA(trimming.species)
        new.tree <- ape::keep.tip(new.tree, trimming.species)}

    if(isTRUE(agg.ssp)){DF1$agg.ssp[is.na(DF1$agg.ssp)]<-"1"}else{DF1$agg.ssp[is.na(DF1$agg.ssp)]<-"0"}
    if(rand.type=="random"){DF1$rand.type[is.na(DF1$rand.type)]<-"0"}else{DF1$rand.type[is.na(DF1$rand.type)]<-"1"}
    {DF1$poly.ins[is.na(DF1$poly.ins)]<-poly.ins}

    DF1[is.na(DF1$resp.mono), "resp.mono"]<- resp.mono
    DF1[is.na(DF1$resp.para), "resp.para"]<- resp.para
    DF1[is.na(DF1$resp.sing), "resp.sing"]<- resp.sing

    DF1[DF1$resp.mono==0, "resp.mono"]<- FALSE
    DF1[DF1$resp.mono==1, "resp.mono"]<- TRUE

    DF1[DF1$resp.para==0, "resp.para"]<- FALSE
    DF1[DF1$resp.para==1, "resp.para"]<- TRUE

    DF1[DF1$resp.sing==0, "resp.sing"]<- FALSE
    DF1[DF1$resp.sing==1, "resp.sing"]<- TRUE

    DF1.poly    <- DF1[DF1$rand.type=="1",]
    DF1.nonpoly <- DF1[DF1$rand.type=="0",]


    if(nrow(DF1.nonpoly)>0){
    DF1.nonpoly$using.taxa <- get.taxa.to.use(DF1.nonpoly)
    DF1.nonpoly <- DF1.nonpoly[order(DF1.nonpoly$taxon),]
    is.duplicated <- duplicated(DF1.nonpoly$using.taxa)
    DF1.dupl <- DF1.nonpoly[ is.duplicated,]
    DF1.rand <- DF1.nonpoly[!is.duplicated,]
    }


    #Phase 1. Random insertions, non-aggregated
    if(!is.null(DF1.rand)){if(nrow(DF1.rand)>0){

        DF1.rand.bind<- DF1.rand[!(DF1.rand$using.taxa %in% new.tree$tip.label),]
        DF1.rand.bind<- DF1.rand.bind[!is.na(DF1.rand.bind$using.MDCC),]

        manual.mdcc.taxa<-DF1.rand.bind$taxon[!is.na(DF1.rand.bind$taxon1)|!is.na(DF1.rand.bind$taxon2)]
        rand.PUTs<- DF1.rand.bind$taxon
        rand.PUTs<-sample(rand.PUTs, length(rand.PUTs), replace = F)
        rand.PUTs<- c(rand.PUTs[rand.PUTs%in%manual.mdcc.taxa], rand.PUTs[!(rand.PUTs%in%manual.mdcc.taxa)])

        if( verbose){
          cat(paste0("\n", "STARTING RANDOMIZATION \n"))
        }

        for(i in seq_along(rand.PUTs)){
            PUT <- rand.PUTs[i]

            MDCC  <- randtip::DF1finder(DF1.rand.bind,PUT, "using.MDCC")
            level <- randtip::DF1finder(DF1.rand.bind,PUT, "using.MDCC.lev")
            MDCC.type <- randtip::DF1finder(DF1.rand.bind,PUT, "using.MDCC.phylstat")

            resp.mono <- as.logical(randtip::DF1finder(DF1.rand.bind, PUT, "resp.mono"))
            resp.para <- as.logical(randtip::DF1finder(DF1.rand.bind, PUT, "resp.para"))
            resp.sing <- as.logical(randtip::DF1finder(DF1.rand.bind, PUT, "resp.sing"))

            poly.ins<- as.character(randtip::DF1finder(DF1.rand.bind,PUT, "poly.ins"))


            if(verbose){
              cat(paste0(i, "/", length(rand.PUTs),
                           " (",round(i/length(rand.PUTs)*100, 2), "%). ",
                           "Adding ", PUT, " to ", MDCC ," (", MDCC.type, " ", level,")\n")) }


            #Manual additions
            if(level=="Sister species"){
              new.tree <- add.to.singleton(tree=new.tree, singleton = MDCC ,new.tips =  PUT)
              next
            }

            if(level=="Sister genus"){
              sister.genus.tips<- new.tree$tip.label[randtip::firstword(new.tree$tip.label)==MDCC]
              sister.genus.mrca<- ape::getMRCA(new.tree, sister.genus.tips)
              new.tree<-add.over.node(tree= new.tree, new.tip = PUT, node = sister.genus.mrca)
              next
              }

            if(level=="Manual clade"){
            sp1<- randtip::DF1finder(DF1.rand.bind, PUT, "taxon1")
            sp2<- randtip::DF1finder(DF1.rand.bind, PUT, "taxon2")
            clade.mrca<- ape::getMRCA(new.tree, c(sp1, sp2))

            new.tree <- add.into.node(tree = new.tree, node = clade.mrca,
                                        new.tip = PUT, prob = prob )
            next
              }




            #Automatically searched MDCCs additions
            #Add to genus
            if(isTRUE(resp.sing)){
              if(length(sp.genus.in.tree(tree, randtip::firstword(PUT)))==0 &
                 length(sp.genus.in.tree(new.tree, randtip::firstword(PUT)))> 0){
                singleton<-sp.genus.in.tree(new.tree, randtip::firstword(PUT))
                except<- DF1[randtip::firstword(DF1$taxon)==randtip::firstword(PUT)&DF1$resp.sing=="FALSE",]
                singleton<- singleton[!(singleton%in%except)]
                new.tree<- add.to.singleton(new.tree, singleton, PUT, T, T)
                next
              }}

            if(level=="genus"){
              MDCC.type <- randtip::MDCC.phyleticity(DF1, new.tree,
                                                     MDCC.info = list(level=level,MDCC=MDCC), trim=F)
            if(MDCC.type=="Monophyletic"){

            new.tree<-add.to.monophyletic(tree = new.tree, new.tip = PUT, prob=prob)

            }else if(MDCC.type=="Paraphyletic"){
            if(isFALSE(resp.para)){new.tree<-add.to.monophyletic(tree = new.tree, new.tip = PUT, prob=prob)}
            if(isTRUE (resp.para)){new.tree <- add.to.paraphyletic(tree = new.tree,
                                                    new.tip = PUT, prob = prob)}
            }else if(MDCC.type=="Polyphyletic"){

              new.tree<- add.to.polyphyletic(tree = new.tree, new.tip = PUT,
                                               poly.ins , prob, resp.mono=resp.mono)
            }else if(MDCC.type=="Singleton MDCC"){
                # All tips added in one step
                singleton <- new.tree$tip.label[(randtip::firstword(new.tree$tip.label) == MDCC)]

                new.tree <- add.to.singleton(tree = new.tree,
                              singleton = singleton, new.tips = PUT,
                              resp.sing = resp.sing, resp.mono = resp.mono)}
            }
            #Add to other taxonomic MDCC
            if(level%in% randtip::randtip_levels()[-1]){
              MDCC.type <- randtip::MDCC.phyleticity(DF1, new.tree,
                              MDCC.info = list(level=level,MDCC=MDCC), trim=F)

              DF1.mdcc<-  DF1[!is.na(DF1[,level]),]
              MDCC.taxa<- DF1.mdcc$taxon[DF1.mdcc[,level]==MDCC]
              MDCC.genera<- randtip::notNA(unique(randtip::firstword(MDCC.taxa)))
              MDCC.intree<- sp.genus.in.tree(new.tree, MDCC.genera)

              if(MDCC.type=="Singleton MDCC"){
                singleton<-c(MDCC.taxa, MDCC.intree)
                singleton<-unique(singleton[singleton%in%new.tree$tip.label])
                new.tree<- randtip::add.to.singleton(new.tree, singleton , PUT,
                                                     resp.sing, resp.mono)}

              if(MDCC.type=="Monophyletic"){
              MDCC.mrca<- ape::getMRCA(new.tree, MDCC.intree)
              if(isTRUE(resp.mono)){
              permitted.nodes<-randtip::get.permitted.nodes(tree=new.tree, node = MDCC.mrca)
              forbidden.nodes<- randtip::get.forbidden.MDCC.nodes(new.tree, DF1, level, MDCC)
              if(!all(permitted.nodes%in%forbidden.nodes)){
              permitted.nodes<- permitted.nodes[!(permitted.nodes%in%forbidden.nodes)]}
              }
              if(isFALSE(resp.mono)){
                permitted.nodes<- phytools::getDescendants(new.tree, MDCC.mrca,curr = NULL)
              }
              if(length(permitted.nodes)==1){nd<-permitted.nodes}else{
                nd<-sample(permitted.nodes,1)}
              pos<- randtip::binding.position(new.tree, node = nd,df = NULL,
                                              insertion = "random", prob = prob)
              new.tree<- phytools::bind.tip(new.tree, PUT, pos$length, pos$where, pos$position)
              }

              if(MDCC.type=="Paraphyletic"){

                if(isFALSE(resp.para)){
                  MDCC.mrca<- ape::getMRCA(new.tree, MDCC.intree)
                  if(isTRUE(resp.mono)){
                    permitted.nodes<-randtip::get.permitted.nodes(tree=new.tree, node = MDCC.mrca)
                    forbidden.nodes<- randtip::get.forbidden.MDCC.nodes(new.tree, DF1, level, MDCC)
                    if(!all(permitted.nodes%in%forbidden.nodes)){
                      permitted.nodes<- permitted.nodes[!(permitted.nodes%in%forbidden.nodes)]}}
                  if(isTRUE(resp.mono)){
                    permitted.nodes<- phytools::getDescendants(new.tree, MDCC.mrca,curr = NULL)
                  }
                  if(length(permitted.nodes)==1){nd<-permitted.nodes}else{
                    nd<-sample(permitted.nodes,1)}
                  pos<- randtip::binding.position(new.tree, node = nd,df = NULL,
                                                  insertion = "random", prob = prob)
                  new.tree<- phytools::bind.tip(new.tree, PUT, pos$length, pos$where, pos$position)
                }#same as monophyletic
                if(isTRUE (resp.para)){
                  MDCC.mrca<- ape::getMRCA(new.tree, MDCC.intree)

                  permitted.nodes<- get.permitted.nodes(new.tree, MDCC.mrca)
                  forbidden.nodes<- randtip::get.intruder.nodes(new.tree, DF1, level, MDCC)
                  if(!all(permitted.nodes%in%forbidden.nodes)){
                    permitted.nodes<- permitted.nodes[!(permitted.nodes%in%forbidden.nodes)]}
                  if(randtip::findRoot(new.tree)%in%permitted.nodes){
                    permitted.nodes<-permitted.nodes[-which(permitted.nodes==randtip::findRoot(new.tree))]}

                  if(length(permitted.nodes)==0){
                    new.tree<-add.over.node(new.tree, new.tip = unit.taxa, node = MDCC.mrca)
                  }else{
                    if(length(permitted.nodes)==1){nd<-permitted.nodes}else{
                      nd<-sample(permitted.nodes,1)}
                    pos<- randtip::binding.position(new.tree, nd,df = NULL,
                                    "random", prob = prob)
                    new.tree<-phytools::bind.tip(new.tree, PUT, pos$length, pos$where, pos$position)

                  }
                  }
                  }

              if(MDCC.type=="Polyphyletic"){
                if(poly.ins=="all"){
                  MDCC.mrca<- ape::getMRCA(new.tree, MDCC.intree)
                  if(isTRUE(resp.mono)){
                    permitted.nodes<-randtip::get.permitted.nodes(tree=new.tree, node = MDCC.mrca)
                    forbidden.nodes<- randtip::get.forbidden.MDCC.nodes(new.tree, DF1, level, MDCC)
                    if(!all(permitted.nodes%in%forbidden.nodes)){
                      permitted.nodes<- permitted.nodes[!(permitted.nodes%in%forbidden.nodes)]}}
                  if(isFALSE(resp.mono)){
                    permitted.nodes<- phytools::getDescendants(new.tree, MDCC.mrca,curr = NULL)
                  }
                  if(length(permitted.nodes)==1){nd<-permitted.nodes}else{
                    nd<-sample(permitted.nodes,1)}
                  pos<- randtip::binding.position(new.tree, node = nd,df = NULL,
                                                  insertion = "random", prob = prob)
                  new.tree<- phytools::bind.tip(new.tree, PUT, pos$length, pos$where, pos$position)
                }#same as monophyletic
                if(poly.ins=="freq"){
                  mrca<- ape::getMRCA(new.tree, MDCC.intree)
                  descs<- phytools::getDescendants(new.tree, mrca,curr = NULL)
                  forbidden.nodes<- randtip::get.forbidden.MDCC.nodes(new.tree, DF1, level, MDCC)
                  descs<- descs[!(descs%in%forbidden.nodes)]
                  table<-sharingtaxa.descs(tree=new.tree, nodes=descs, MDCC.genera = MDCC.genera)
                  table.row<- sample(1:length(table$node), 1, prob=as.numeric(table$number))
                  node<-table[table.row,"node"]


                  if(isTRUE(resp.mono)){
                    nodes<- randtip::get.permitted.nodes(new.tree, node)
                  }
                  if(isFALSE(resp.mono)){
                    nodes<- phytools::getDescendants(new.tree, node,curr = NULL)

                  }
                  if(length(nodes)==1){nd<-nodes}else{
                    nd<-sample(nodes,1)}
                  pos<- randtip::binding.position(new.tree, nd,df = NULL,
                                  insertion = "random", prob=prob)
                  new.tree<-phytools::bind.tip(new.tree, PUT, pos$length, pos$where, pos$position)


                }
                if(poly.ins=="large"){
                  mrca<- ape::getMRCA(new.tree, MDCC.intree)
                  descs<- phytools::getDescendants(new.tree, mrca,curr = NULL)
                  forbidden.nodes<- randtip::get.forbidden.MDCC.nodes(new.tree, DF1, level, MDCC)
                  descs<- descs[!(descs%in%forbidden.nodes)]
                  table<-sharingtaxa.descs(tree=new.tree, nodes=descs, MDCC.genera = MDCC.genera)
                  table<- table[table$number==max(table$number),]
                  table.row<- sample(1:length(table$node), size = 1)
                  node<-table[table.row,"node"]

                  if(isTRUE(resp.mono)){
                    nodes<- randtip::get.permitted.nodes(new.tree, node)
                  }
                  if(isFALSE(resp.mono)){
                    nodes<- phytools::getDescendants(new.tree, node,curr = NULL)

                  }
                  if(length(nodes)==1){nd<-nodes}else{
                    nd<-sample(nodes,1)}
                  pos<- randtip::binding.position(new.tree, nd,df = NULL,
                                                  insertion = "random", prob=prob)
                  new.tree<-phytools::bind.tip(new.tree, PUT, pos$length, pos$where, pos$position)


                }

              }}}
        }}

    #Phase 2 - Random insertions, aggregated subspecies
    if(!is.null(DF1.dupl)){if(nrow(DF1.dupl)>0){
      DF1.masdupl<- DF1.rand[(DF1.rand$using.taxa%in%tree$tip.label) &
                             !(DF1.rand$taxon%in%tree$tip.label) , ]
      DF1.dupl<- rbind(DF1.dupl, DF1.masdupl)
      new.tree <- randtip.subsp(tree = new.tree, DF1.dupl, verbose)
            }}

    if(verbose){cat("\n")}
    if(nrow(DF1.nonpoly)>0){new.tree <- get.original.names(tree = new.tree, DF1 = DF1.nonpoly, verbose=F )}


    #Polytomies
    if(nrow(DF1.poly)>0){
      if(verbose){
        cat(paste0("\n","Starting polytomies addition \n"))
      }
      DF1.poly<-DF1.poly[!(DF1.poly$taxon %in% new.tree$tip.label),]


      MDCCs<- randtip::notNA(unique(DF1.poly$using.MDCC))
      for(MDCCs.i in MDCCs){

        DF1.poly<- DF1.poly[!is.na(DF1.poly$using.MDCC),]
        MDCCs.i.level<- unique(DF1.poly$using.MDCC.lev[DF1.poly[,"using.MDCC"]==MDCCs.i])


        MDCC.taxa.toAdd <- DF1.poly$taxon[DF1.poly[,"using.MDCC"]==MDCCs.i]
        MDCC.taxa.inDF1 <- DF1$taxon[DF1[, MDCCs.i.level]==MDCCs.i]
        MDCC.genera <- unique(randtip::firstword(MDCC.taxa.inDF1))
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
          new.tree<- randtip::polytomy.into.node(tree=new.tree,
                            new.tip =MDCC.taxa.toAdd, node =  MDCC.mrca)}
            }



    }




    complete.taxa.list <- originalDF1$taxon
    complete.taxa.list.in.tree <- complete.taxa.list[complete.taxa.list %in% new.tree$tip.label]
    not.included <- complete.taxa.list[!(complete.taxa.list %in% complete.taxa.list.in.tree)]
    if(length(not.included) > 0){
        message("The following taxa were not included in the tree:\n", paste0(not.included, "\n"))}

    if(isTRUE(prune)){new.tree <- ape::keep.tip(new.tree, complete.taxa.list.in.tree)}


    end <- Sys.time()


    if(verbose){
      cat(paste0("\n","\n","\U2713", " Randomization completed in ",
                 round(as.numeric(difftime(end, start,units = "mins")), 2), " mins\n"))
    }
    return(new.tree)
}
