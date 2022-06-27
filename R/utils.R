notNA <- function(x){
    vect<- x[!is.na(x)]
    return(vect)
}

############ String manipulation ###################
first.word<- function(string){
    word <- stringr::str_extract(string, "[A-Za-z]+")
    return(word)
}

remove.spaces <- function(species, sep = "_"){
    sp <- stringr::str_trim(gsub("\\s+", " ", species))
    sp <- stringr::str_replace(sp, " ", sep)
    return(sp)
}

############### Dataframe ###################
randtip.ranks<- function(){
    return(as.vector(c("genus","subtribe","tribe",
                      "subfamily","family","superfamily",
                      "order","class")))
}

correct.DF<- function(df){

    df_corrected <- apply(df, 2,
                          function(x){
                                      x[x == ""] <- NA
                                      x <- as.character(x)
    })
    
    return(as.data.frame(df_corrected))
}

isRoot<- function(tree, node){
    tips<- length(tree$tip.label)
    if(node==(tips+1L)){root<-TRUE}else{root<-FALSE}
    return(root)
}

findRoot<- function(tree){
    tips<- length(tree$tip.label)
    return(tips+1L)
}

randtip_ranks<- function(){
    return(as.vector(c("genus","subtribe","tribe","subfamily","family","superfamily","order","class")))
}

is.node<-function(tree, node){
    if(!(node %in% c(tree$edge[,1], tree$edge[,2]))){
        stop("Node number is not in your tree")
    }
    if(length(phytools::getDescendants(tree = tree, node = node, curr=NULL )) > 1){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

is.tip <-function(tree, node){
    if(!(node %in% c(tree$edge[,1], tree$edge[,2]))){
        stop("Node number is not in your tree")
    }
    if(length(phytools::getDescendants(tree = tree, node = node, curr=NULL)) == 1){
        return(TRUE)
    }else{
        return(FALSE)
    }
}


#SPECIFIC FUNCTIONS####
# Return a vector of all species in the tree that match a given genus
sp.genus.in.tree <- function(tree, genus){
    sp <- tree$tip.label
    taxa.vector <- sp[first.word(sp)%in%genus]

    return(taxa.vector)
}

inputfinder<- function(input, taxon, column){
    return(as.character(input[input$taxon==taxon, column]))
}

usingMDCCfinder<- function(input, taxon=NULL, tree, silent = FALSE){

    if(is.null(taxon)){taxon<- input$taxon}
    input<-correct.DF(input)
    MDCC.vect<- vector( mode="character", length = length(taxon))
    MDCC.lev.vect<- vector(mode="character", length = length(taxon))
    #MDCC.phyletictype.vect<- vector(mode="character", length = length(taxon))

    if(!silent){cat(paste0("Searching MDCCs...\n"))}

    #manual MDCC search
    taxa<- input[!is.na(input$taxon1)|!is.na(input$taxon2),]
    if(nrow(taxa)>0){
        for(tx in seq_along(taxa$taxon)){
            if(isTRUE(length(strsplit(taxa$taxon1[tx], "_")[[1]])==1 &
                      length(strsplit(taxa$taxon2[tx], "_")[[1]])==1 &
                      taxa$taxon1[tx]==taxa$taxon2[tx] &
                      length(sp.genus.in.tree(tree, taxa$taxon1[tx]))>0)){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- taxa$taxon1[tx]
                MDCC.lev.vect[pos] <- "Sister genus"
                #MDCC.phyletictype.vect[pos]<-"-"
                next
            }

            if(!(any(tree$tip.label == taxa$taxon1[tx]))){taxa$taxon1[tx]<-NA}
            if(!(any(tree$tip.label == taxa$taxon2[tx]))){taxa$taxon2[tx]<-NA}

            if(isTRUE(any(tree$tip.label == taxa$taxon1[tx]) &
                      any(tree$tip.label == taxa$taxon2[tx]) &
                      taxa$taxon1[tx]==taxa$taxon2[tx])){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- taxa$taxon1[tx]
                MDCC.lev.vect[pos] <- "Sister species"

                next
            }

            if(isTRUE(any(tree$tip.label == taxa$taxon1[tx]) &
                      any(tree$tip.label == taxa$taxon2[tx]) &
                      taxa$taxon1[tx]!=taxa$taxon2[tx])){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- paste0("Clade ", taxa$taxon1[tx], "-", taxa$taxon2[tx])
                MDCC.lev.vect[pos] <- "Manual setting"

                next
            }

            if(isTRUE(any(tree$tip.label == taxa$taxon1[tx]) &
                      is.na(taxa$taxon2[tx]))){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- taxa$taxon1[tx]
                MDCC.lev.vect[pos] <- "Sister species"
                next
            }

            if(isTRUE(any(tree$tip.label == taxa$taxon2[tx]) &
                      is.na(taxa$taxon1[tx]))){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- taxa$taxon2[tx]
                MDCC.lev.vect[pos] <- "Sister species"
                next
            }

        }
    }

    #automatic MDCC search
    ranks<- randtip_ranks()
    taxa<- input[!(!is.na(input$taxon1)|!is.na(input$taxon2)),]

    if(nrow(taxa)>0){
        vect<- which(taxon%in%taxa$taxon)
        for(v in vect){

            if(!silent){

                if(v==vect[1]){
                    cat(paste0("0%       25%       50%       75%       100%", "\n",
                        "|---------|---------|---------|---------|",   "\n"))
                }

                vec<- seq(from=0, to=40, by=40/length(vect))
                vec<-ceiling(vec)
                vec<- diff(vec)
                cat(strrep("*", times=vec[which(vect==v)]))

                if(v ==vect[length(vect)]){cat("*\n")}

            }

            if(any(tree$tip.label == taxon[v])){
              MDCC.vect[v]<- "Tip"
              MDCC.lev.vect[v]<-"Tip"
              next
            }

            i<- which(input$taxon==taxon[v])
            if((MDCC.vect[v])==""){

                MDCC<-as.character(NA)
                MDCC.ranks<-as.character(NA)

                for(rank in ranks){
                    if(is.na(MDCC)){
                        MDCC<-as.character(input[i, rank])
                        if(!is.na(MDCC)){
                            phyleticity<-MDCC.phyleticity(input, tree = tree,
                                                      MDCC.info = list(rank=rank, MDCC= MDCC))
                            if(phyleticity=="Missing"){MDCC<-NA}
                        }

                        lev<-rank
                    }else{next}
                }
                MDCC.vect[v]<-as.character(MDCC)
                MDCC.lev.vect[v]<-as.character(lev)
                if(is.na(MDCC)){MDCC.lev.vect[v]<-NA}

            }
        }
    }

    return(list(MDCC=MDCC.vect,MDCC.ranks=MDCC.lev.vect) )
}

get.parent.siblings <- function(tree, tip){
    tree.sp <- tree$tip.label
    # Direct ancestor
    parent <- tree$edge[tree$edge[,2] == tip, 1]
    # Ancestor's descendants
    parent.desc <- phytools::getDescendants(tree, parent,curr = NULL)
    siblings <- tree.sp[parent.desc]
    siblings <- siblings[!is.na(tree.sp[parent.desc])]

    return(list(parent = parent,
                siblings = siblings))
}

get.groups <- function(tree, genus){
    species <- sp.genus.in.tree(tree, genus)
    sp.mrca<- ape::getMRCA(tree, species)
    mrca.descs <- phytools::getDescendants(tree, sp.mrca,curr = NULL)

    node.descs<- rep(list(NA), times=length(mrca.descs))
    names(node.descs)<- mrca.descs

    node.types<- rep(list(NA), times=length(mrca.descs))
    names(node.types)<- mrca.descs

    for(n in seq_along(mrca.descs)){
        nd<- mrca.descs[n]

        if(is.tip(tree = tree,node = nd)){
            node.descs[[n]]<- tree$tip.label[nd]
            node.types[[n]]<- "tip"
            next
        }

        nd.descs<- phytools::getDescendants(tree, nd,curr = NULL)
        nd.descs <- notNA(tree$tip.label[nd.descs])
        nd.genera<- first.word(nd.descs)

        if(!(genus%in%nd.genera)){
            node.descs[[n]]<- NA
            node.types[[n]]<- NA
            next
        }


        siblings<-get.parent.siblings(tree, tip = nd)$siblings
        siblings<-siblings[-which(siblings%in%nd.descs)]
        siblings.genus<- first.word(siblings)
        if(all(nd.genera==genus)){
            # if(genus%in%siblings.genus){
            #   node.descs[[n]]<- NA
            #   node.types[[n]]<- NA}else{
            node.descs[[n]]<- nd.descs
            node.types[[n]]<- "monophyletic"#}
            next
        }

        intruders<- nd.descs[nd.genera!=genus]

        if(length(intruders)==1){
            node.descs[[n]]<- nd.descs
            node.types[[n]]<- "paraphyletic"
            next
        }

        intruders.mrca<- ape::getMRCA(tree, intruders)

        intruders.descs<- phytools::getDescendants(tree, intruders.mrca,curr = NULL)
        intruders.descs <- notNA(tree$tip.label[intruders.descs])
        intruders.genera<- first.word(intruders.descs)

        if(genus %in% intruders.genera){
            node.descs[[n]]<- nd.descs
            node.types[[n]]<- "polyphyletic"
            next
        }
        if(length(nd.genera[nd.genera==genus])==1){
            node.descs[[n]]<- nd.descs[nd.genera==genus]
            node.types[[n]]<- "singleton"
            next
        }

        group<- nd.descs[nd.genera==genus]
        group.mrca<- ape::getMRCA(tree, group)
        group.descs<- phytools::getDescendants(tree, group.mrca,curr = NULL)
        group.descs<- notNA(tree$tip.label[group.descs])
        group.descs.gen <- first.word(group.descs)
        if(!all(group.descs.gen==genus)){
            node.descs[[n]]<- nd.descs
            node.types[[n]]<- "paraphyletic"
        }else{
            node.descs[[n]]<- NA
            node.types[[n]]<- NA
        }
    }

    node.types<-node.types[as.vector(which(!is.na(node.types)))]
    node.descs<-node.descs[as.vector(which(!is.na(node.descs)))]

    permitted.type<- as.vector(unlist(node.types))
    permitted.type<- which(!(permitted.type%in% c("polyphyletic")))

    permitted.gen<- grep(paste0(genus,"_"), (node.descs))

    permitted<- intersect(permitted.type,permitted.gen)

    node.types<-node.types[permitted]
    node.descs<-node.descs[permitted]

    return(list(species=node.descs, type=node.types))
}

get.position <- function(tree, node){
    if(isRoot(tree, node)){
        position=0
    }else{
        df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                          "length"= tree$edge.length, "id"=1:length(tree$edge[,1]) )
        edge.length <- df[df$node==node,"length"]
        # Bind at a random point of the branch
        position <- edge.length * stats::runif(1, 0, 1)
    }

    return(position)
}

binding.position<- function(tree, node,  insertion,  prob){
    position<-list("length"=NA, "where"=NA, "position"=NA)
    df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                     "length"= tree$edge.length, "id"=1:length(tree$edge[,1]))

    if(ape::is.ultrametric(tree)){
        position$length<-NULL
    }else{
        position$length<-abs(runif(1, 0, max(tree$edge.length)))
    }

    df<- df[df$node==node,]
    position$where<- node

    if(insertion=="polytomy"){
        position$position<- 0
        position$where <- node
    }

    if(insertion=="random"){

        index<- df[df$node==node,"id"]
        pos <- get.position(tree = tree, node = node)

        position$position<- pos
        position$where <- df[df$id==index,"node"]
    }

    return(position)
}

sharingtaxa.descs<-function(tree, nodes, MDCC.genera){
    table<- data.frame("node"=as.numeric(nodes), "number"=rep(NA, length(nodes)),
                       "tot.number"=rep(NA, length(nodes)))
    for(i in seq_along(table$node)){
        node<- table$node[i]
        if(is.tip(tree,node)){
            table$number[i]<-1;table$tot.number[i]<-1
            next
        }
        descs<- phytools::getDescendants(tree, node,curr = NULL)
        table$tot.number[i]<-length(descs)
        descs<- tree$tip.label[descs]
        descs<- notNA(descs)
        descs<- descs[first.word(descs)%in%MDCC.genera]
        table$number[i]<-length(descs)
    }
    table<- table[table$number>0,]

    return(table)
}

get.permitted.nodes <- function (tree, input, MDCC, rank, MDCC.type,
                                 polyphyly.scheme, use.paraphyletic, use.singleton, use.stem){

    input.mdcc <-  input[!is.na(input[,rank]),]
    MDCC.taxa  <- input.mdcc$taxon[input.mdcc[,rank]==MDCC]
    MDCC.genera<- notNA(unique(first.word(MDCC.taxa)))
    MDCC.intree<- sp.genus.in.tree(tree, MDCC.genera)

    if(MDCC.type=="Monophyletic"){
        MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
        nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
        if(use.stem){nodes <- c(MDCC.mrca,nodes)}
    }

    if(MDCC.type=="Paraphyletic"){
        if(!use.paraphyletic){
            MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
            nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
            if(use.stem){nodes <- c(MDCC.mrca,nodes)}
        }
        if(use.paraphyletic){

            MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)

            descendants.nodes<- phytools::getDescendants(tree, node=MDCC.mrca,curr = NULL)
            descendants.tips <- tree$tip.label[descendants.nodes]
            descendants.tips<- notNA(descendants.tips)

            non.mdcc <- input[!is.na(input[,rank]),]
            non.mdcc <- non.mdcc[non.mdcc[,rank]!=MDCC,]
            non.mdcc.genera<- unique(first.word(non.mdcc$taxon))

            intruder.descs<- descendants.tips[!(first.word(descendants.tips)%in%MDCC.genera)]
            if(length(intruder.descs)==1){
                intruder.descs.nodes<-NULL
            }else{
                intruder.mrca<- ape::getMRCA(tree, intruder.descs)
                intruder.descs.nodes<- phytools::getDescendants(tree, intruder.mrca, curr=NULL)
            }

            nodes <- descendants.nodes[!(descendants.nodes%in%intruder.descs.nodes)]
            if(use.stem){nodes <- c(MDCC.mrca,nodes)}
        }

    }

    if(MDCC.type == "Singleton"){
        if(use.singleton){
            nodes<- which(tree$tip.label==MDCC.intree)
            if(length(nodes)>1){
                nodes<- phytools::getDescendants(ape::getMRCA(tree, MDCC.intree))
                if(use.stem){nodes <- c(ape::getMRCA(tree, MDCC.intree),nodes)}
            }
        }
        if(!use.singleton){
            tip<- which(tree$tip.label==MDCC.intree)
            nodes<- phytools::getDescendants(tree, get.parent.siblings(tree, tip)$parent, curr=NULL)
        }
    }

    if(MDCC.type=="Polyphyletic"){
        if(polyphyly.scheme=="complete"){
            MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
            nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
            if(use.stem){nodes <- c(MDCC.mrca,nodes)}
        }

        if(polyphyly.scheme == "frequentist"){
            MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
            nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)

            if(use.stem){nodes <- c(MDCC.mrca,nodes)}

            table<- data.frame("node"=as.numeric(nodes), "descs"=NA, "total.descs"=NA,
                              "sharing.descs"=NA, "eligible"=NA)
            for(i in seq_along(table$node)){
                node<- table$node[i]
                if(is.tip(tree,node)){
                    table$descs[i] <- NA
                    table$total.descs[i]<-1
                    if(tree$tip.label[node] %in% input[input[,rank]==MDCC,"taxon"]){
                        table$sharing.descs[i]<-1;table$eligible[i]<- "TRUE"
                    }else{
                            table$sharing.descs[i]<-0;table$eligible[i]<- "FALSE"
                    }
                }
                if(is.node(tree,node)){
                    descs<- phytools::getDescendants(tree, node, curr=NULL)
                    table$descs[i] <- paste0(descs, collapse = ",")
                    descs<-notNA(tree$tip.label[descs])
                    table$total.descs[i]   <- length(descs)
                    table$sharing.descs[i] <- length(descs[descs %in% input[input[,rank]==MDCC,"taxon"]])
                    if(length(descs[descs %in% input[input[,rank]!=MDCC,"taxon"]])>0){
                        table$eligible[i]<-"FALSE"
                    }else{
                        table$eligible[i]<-"TRUE"
                    }
                }
            }

            table<- table[table$eligible=="TRUE",]
            for(i in seq_along(table$node)){
                node<- table$node[i]
                descs<- unique(unlist(strsplit(table$descs, split=",")))
                if(as.character(node)%in%descs ){table$eligible[i]<- "FALSE"}
            }

            table<- table[table$eligible=="TRUE",]
            if(use.stem){table$descs<- paste0(table$node,",", table$descs)}
            table$descs<- gsub(",NA", "", table$descs)
            if(nrow(table)==1){
                nodes <- table$descs
            }else{
                node<- sample(table$node, size = 1,prob=table$total.descs)
                nodes <- table$descs[table$node==node]
                nodes <- as.numeric(strsplit(nodes, split=",")[[1]])
            }

        }

        if(polyphyly.scheme == "largest"){
            MDCC.mrca<- ape::getMRCA(tree, MDCC.intree)
            nodes <- phytools::getDescendants(tree, MDCC.mrca, curr=NULL)
            if(use.stem){nodes <- c(MDCC.mrca,nodes)}
            table<- data.frame("node"=as.numeric(nodes), "descs"=NA, "total.descs"=NA,
                              "sharing.descs"=NA, "eligible"=NA)
            for(i in seq_along(table$node)){
                node<- table$node[i]
                if(is.tip(tree,node)){
                    table$descs[i]   <- NA
                    table$total.descs[i]<-1
                    if(tree$tip.label[node]%in% input[input[,rank]==MDCC,"taxon"]){
                        table$sharing.descs[i]<-1;table$eligible[i]<- "TRUE"
                    }else{
                        table$sharing.descs[i]<-0;table$eligible[i]<- "FALSE"
                    }
                }
                if(is.node(tree,node)){
                    descs<- phytools::getDescendants(tree, node, curr=NULL)
                    table$descs[i] <- paste0(descs, collapse = ",")
                    descs<-notNA(tree$tip.label[descs])
                    table$total.descs[i]   <- length(descs)
                    table$sharing.descs[i] <- length(descs[descs %in% input[input[,rank]==MDCC,"taxon"]])
                    if(length(descs[descs %in% input[input[,rank]!=MDCC,"taxon"]])>0){
                        table$eligible[i]<-"FALSE"
                    }else{
                        table$eligible[i]<-"TRUE"
                    }
                }
            }
            table<- table[table$eligible=="TRUE",]
            for(i in seq_along(table$node)){
                node<- table$node[i]
                descs<- unique(unlist(strsplit(table$descs, split=",")))
                if(node%in%descs ){table$eligible[i]<- "FALSE"}
            }
            table<- table[table$eligible=="TRUE",]
            table<-table[table$sharing.descs== max(table$sharing.descs),]
            if(use.stem){table$descs<- paste0(table$node,",", table$descs)}
            table$descs<- gsub(",NA", "", table$descs)
            if(nrow(table)>1){table<-table[sample(1:nrow(table), size=1),]}

            nodes <- as.numeric(strsplit(table$descs, split=",")[[1]])
        }
    }

    return(nodes)
}

get.forbidden.nodes <- function(tree,input, MDCC, rank, perm.nodes, respect.mono, respect.para){

    forbidden.nodes<- vector("numeric")
    #groups which must not be forbid
    perm.groups<- input[input[,rank]==MDCC,randtip_ranks()[(which(randtip_ranks()==rank)):8]]
    perm.groups<- notNA(unique(as.vector(unlist(perm.groups))))

    #respect monophyletic clusters of species
    if(respect.mono){
        perm.tips<- notNA(tree$tip.label[perm.nodes])
        perm.species<-paste0(first.word(perm.tips), "_",
                             stringr::word(perm.tips, 2, sep="_"))
        ssps<-perm.species[duplicated(perm.species)]
        if(length(ssps)>0){
            for(ssp in ssps){
                nodes<-which(paste0(first.word(tree$tip.label), "_",
                                    stringr::word(tree$tip.label, 2, sep="_"))==ssp)
                if(length(nodes)==2){
                    forbidden.nodes<-c(forbidden.nodes, nodes)
                    next
                }
                if(length(nodes)> 2){
                    ssp.mrca<- ape::getMRCA(tree, nodes)
                    ssp.descs<- phytools::getDescendants(tree, ssp.mrca, curr=NULL)
                    for(n in ssp.descs){
                        if(n %in%forbidden.nodes){next}
                        node.descs<- tree$tip.label[notNA(phytools::getDescendants(tree, n, curr=NULL))]
                        node.descs<- paste0(first.word(node.descs), "_",
                                            stringr::word(node.descs, 2, sep="_"))
                        if(all(node.descs==ssp)){forbidden.nodes<-c(forbidden.nodes, nodes)}
                    }
                }
            }
        }
    }


    if(respect.mono){
        for(i in seq_along(perm.nodes)){
            nd<- perm.nodes[i]
            if(nd %in% forbidden.nodes){next}
            if(is.tip(tree, nd)){next}
            descs.nd<- phytools::getDescendants(tree, nd)
            descs <- notNA(tree$tip.label[descs.nd])
            genera<- unique(first.word(descs))

            #respect monophyletic genera, disregarding input info
            if(length(genera)==1){
                if(!(genera%in%perm.groups)){
                    forbidden.nodes<- c(forbidden.nodes, descs.nd )
                    next
                }
            }

            #respect monophyletic groups, according to input info
            sub.input <- input[first.word(input$taxon)%in%genera,]

            for(rk in randtip_ranks()){
                rk.vals<-unique(sub.input[,rk])
                if(all(is.na(rk.vals))){next}
                if(length(rk.vals)==1){
                    if(!(rk.vals%in%perm.groups)){
                    forbidden.nodes<- c(forbidden.nodes, descs.nd )
                    next
                    }
                }
            }

        }
    }

    if(respect.para){
        for(i in seq_along(perm.nodes)){

            nd<- perm.nodes[i]
            if(nd %in% forbidden.nodes){next}
            if(is.tip(tree, nd)){next}
            descs.nd<- phytools::getDescendants(tree, nd, curr=NULL)
            descs <- notNA(tree$tip.label[descs.nd])
            genera<- unique(first.word(descs))

            if(length(genera)==1){next}

            #account for singleton tips in otherwise monophyletic clusters
            if(length(genera)==2){
                table<- table(first.word(descs))
                uniq.gen<- names(table)[which(table==1)]
                nest.gen<- names(table)[which(table!=1)]
                if(length(uniq.gen)==1){
                    tip<-descs[first.word(descs)%in%uniq.gen]
                    tip.n<- which(tree$tip.label==tip)

                    tip.par<- get.parent.siblings(tree, tip.n)$parent

                    if(tip.par!=nd){
                        if(!(nest.gen%in%perm.groups)){
                            forbidden.nodes<- c(forbidden.nodes, descs.nd[which(descs.nd!=tip.n)] )
                            next
                        }
                    }
                }
            }

            #respect paraphyletic genera, disregarding input info
            rk.vals.mrca<- vector("numeric", length = length(genera))
            rk.vals.desc<- list( NA)
            for(v in genera){
                ds <- descs[first.word(descs)%in%v]
                if(length(ds)==1){
                    rk.vals.mrca[which(genera==v)]<-which(tree$tip.label==ds)
                    rk.vals.desc[[which(genera==v)]]<- NA
                    next
                }

                ds.mrca<- ape::getMRCA(tree, ds)
                rk.vals.mrca[which(genera==v)]<-ds.mrca
                rk.vals.desc[[which(genera==v)]]<- phytools::getDescendants(tree, ds.mrca, curr=NULL)
            }
            nest<- genera[which(rk.vals.mrca==nd)]
            if(length(nest)==1){
                if(phyleticity(tree, nest)=="Paraphyletic" &
                   !(nest%in%perm.groups)){

                    nested <-rk.vals.desc; nested[which(rk.vals.mrca==nd)]<-NA
                    nested <- tree$tip.label[unlist(nested)]

                    if(!any(first.word(nested)%in%nest)){
                        p<-which(rk.vals.mrca==nd); if(length(p)>1){p<-p[1]}
                        para.nodes<- rk.vals.desc[[p]]
                        intr.nodes <- rk.vals.desc
                        intr.nodes[[p]]<-NA
                        intr.nodes<- notNA(unlist(intr.nodes))
                        para.nodes<- para.nodes[!(para.nodes%in%intr.nodes)]
                        forbidden.nodes<- c(forbidden.nodes, para.nodes )
                        next
                    }
                }
            }

            #respect paraphyletic groups, according to input info
            sub.input <- input[first.word(input$taxon)%in%genera,]

            for( rk in randtip_ranks()){
                rk.vals<-notNA(unique(sub.input[,rk]))

                if(length(rk.vals)>1){
                    rk.vals.mrca<- vector("numeric", length = length(rk.vals))
                    rk.vals.desc<- list( NA)
                    for(v in rk.vals){
                        gen <- unique(first.word(sub.input$taxon[sub.input[,rk]==v]))
                        ds <- descs[first.word(descs)%in%gen]
                        if(length(ds)==1){
                            ds.mrca<-which(tree$tip.label==ds)
                            rk.vals.mrca[which(rk.vals==v)]<- ds.mrca
                            rk.vals.desc[[which(rk.vals==v)]]<- NA
                        }
                        if(length(ds)>1){ds.mrca<- ape::getMRCA(tree, ds)}
                        if(is.null(ds.mrca)){
                            rk.vals.mrca[which(rk.vals==v)]<- which(tree$tip.label==sub.input$taxon[sub.input[,rk]==v])
                            rk.vals.desc[[which(rk.vals==v)]]<- NA
                        }
                        if(!is.null(ds.mrca)){
                            rk.vals.mrca[which(rk.vals==v)]<-ds.mrca
                            rk.vals.desc[[which(rk.vals==v)]]<- phytools::getDescendants(tree, ds.mrca, curr=NULL)
                        }

                    }

                    nest<- rk.vals[rk.vals.mrca==nd]
                    if(length(nest)==1){
                        if(MDCC.phyleticity(input, tree, list(rank=rk, MDCC=nest))=="Paraphyletic"&
                           !(nest %in% perm.groups)){

                            nest.tips<- sub.input$taxon[sub.input[,rk]==nest]
                            nest.gen<- first.word(nest.tips)
                            nest.mrca<- nd

                            nested.nodes<- phytools::getDescendants(tree, node=nest.mrca,curr = NULL)
                            nest.desc.tips <- tree$tip.label[nested.nodes]
                            nest.desc.tips<- notNA(nest.desc.tips)

                            non.nest <- input[!is.na(input[,rank]),]
                            non.nest <- non.nest[non.nest[,rk]!=nest,]
                            non.nest.genera<- unique(first.word(non.nest$taxon))

                            intruders <- nest.desc.tips[first.word(nest.desc.tips)%in%non.nest.genera]
                            if(length(intruders)==1){
                                intruder.mrca<- which(tree$tip.label==intruders)
                            }else{
                                intruder.mrca<- ape::getMRCA(tree, intruders)
                            }

                            if(any(first.word(intruders) %in% nest.gen)){next}


                            int.nodes<-phytools::getDescendants(tree, intruder.mrca, curr=NULL)

                            para.nodes<- nested.nodes[!(nested.nodes%in%c(int.nodes, intruder.mrca))]

                            forbidden.nodes<- c(forbidden.nodes, para.nodes )

                            next
                        }
                    }
                }
            }
        }
    }

    return(unique(forbidden.nodes))
}


bind.clump<- function(new.tree, tree, input, PUT){

    clumplist<- list("MDCC"=NULL, "rank"=NULL,"MDCC.type"=NULL, "taxa"=NULL)
    sp<- paste0(first.word(PUT), "_", stringr::word(PUT, 2,sep = "_"))
    treespp<- paste0(first.word(tree$tip.label), "_",
                     stringr::word(tree$tip.label, 2,sep = "_"))
    new.treespp<- paste0(first.word(new.tree$tip.label), "_",
                         stringr::word(new.tree$tip.label, 2,sep = "_"))

    if(sp%in%new.treespp){
        DFspp<- paste0(first.word(input$taxon), "_",
                       stringr::word(input$taxon, 2,sep = "_"))
        clump<- tree$tip.label[treespp==sp]
        clumpDF<-input[DFspp==sp,]
        clumpDF<-clumpDF[clumpDF$clump.puts==TRUE,"taxon"]
        clump<- unique(c(clump, clumpDF))
        clump<- clump[clump%in%new.tree$tip.label]
        if(length(clump)>1){
            mrca<- ape::getMRCA(new.tree, clump)
            descs<- new.tree$tip.label[phytools::getDescendants(new.tree, mrca, curr=NULL)]
            descs<-notNA(descs)
            if(any(!(descs%in%clump))){clump<- sample(clump, 1)}
        }
        clumplist$MDCC<- sp
        clumplist$rank<- "species"
        clumplist$MDCC.type<- "Singleton"
        clumplist$taxa<- clump
        return(clumplist)
    }

    if(input[input$taxon==PUT, "MDCC.rank"]=="genus"){return(clumplist)}


    tree.genus.tips<-sp.genus.in.tree(tree, first.word(PUT))
    new.tree.genus.tips<-sp.genus.in.tree(new.tree, first.word(PUT))
    if(length(tree.genus.tips)==0 & length(new.tree.genus.tips)>0){
        clumplist$MDCC<-first.word(PUT)
        clumplist$rank<-"genus"
        clumplist$MDCC.type<- phyleticity(new.tree, first.word(PUT))
        clumplist$taxa<-new.tree.genus.tips

        return(clumplist)
    }

    newsearch<- usingMDCCfinder(input, PUT, new.tree, silent=TRUE)
    MDCC<-input[input$taxon==PUT,"MDCC"]
    MDCC.rank<- input[input$taxon==PUT,"MDCC.rank"]

    if(newsearch$MDCC!=MDCC){
        newclump<- input$taxon[input[,newsearch$MDCC.ranks]==newsearch$MDCC]
        newclump<- unique(first.word(newclump))
        newclump <- sp.genus.in.tree(new.tree, newclump)

        clumplist$MDCC<-newsearch$MDCC
        clumplist$rank<-newsearch$MDCC.ranks
        clumplist$taxa<-newclump
        clumplist$MDCC.type<- MDCC.phyleticity(input, new.tree,
                                  MDCC.info = list(rank=newsearch$MDCC.ranks,
                                                   MDCC=newsearch$MDCC, rank), trim = T)

    }

    return(clumplist)
}


add.to.singleton <- function(tree, singleton, new.tips, use.singleton=F, respect.mono=F, respect.para=F){
    singleton<-gsub(" ", "_", singleton)
    singleton<-singleton[singleton%in%tree$tip.label]

    new.tree <- tree

    if(length(singleton)==1){

        node<-which(new.tree$tip.label==singleton)

        if(isTRUE(use.singleton)){
            pos<- binding.position(new.tree, node = node, insertion = "random",prob = T)
            new.tree <- phytools::bind.tip(new.tree,
                                          new.tips,
                                          edge.length = pos$length,
                                          where = pos$where,
                                          position = pos$position)
        }
        if(isFALSE(use.singleton)){
            if(!(isRoot(new.tree, node))){
                parent<- get.parent.siblings(new.tree, node)[[1]]
                if(isFALSE(respect.mono)){nodes<- phytools::getDescendants(new.tree, parent)}
                if(isTRUE(respect.mono)){
                    nodes<- get.permitted.nodes(new.tree, parent)
                    nodes<- nodes[nodes!=parent]
                    if(length(nodes)==0){nodes<-node}
                 }
            }else{
                nodes<-node
            }
            pos<- binding.position(new.tree, node = sample(nodes,1), insertion = "random",prob = T)
            new.tree <- phytools::bind.tip(new.tree,
                                          new.tips,
                                          edge.length = pos$length,
                                          where = pos$where,
                                          position = pos$position)
        }
    }

    if(length(singleton)> 1){
        node<-which(new.tree$tip.label%in%singleton)
        mrca<- ape::getMRCA(new.tree, singleton)

        if(isTRUE(use.singleton)){
            nodes<- c(mrca, node)
            pos<- binding.position(new.tree, node = sample(nodes,1), insertion = "random",prob = T)
            new.tree <- phytools::bind.tip(new.tree,
                                          new.tips,
                                          edge.length = pos$length,
                                          where = pos$where,
                                          position = pos$position)
        }
        if(isFALSE(use.singleton)){
            if(!(isRoot(new.tree, mrca))){
                parent<- get.parent.siblings(new.tree, mrca)[[1]]
                if(isFALSE(respect.mono)){
                    nodes<- phytools::getDescendants(new.tree, parent)
                }
                if(isTRUE(respect.mono)& isFALSE(respect.para)){
                    nodes<-get.permitted.nodes(new.tree, mrca, respect.para = F)
                }
                if(isTRUE(respect.mono)& isTRUE(respect.para)){
                    nodes<-get.permitted.nodes(new.tree, mrca, respect.para = T)
                }

              nodes<- nodes[nodes!=parent]
              if(length(nodes)==0){nodes<-c(node,mrca)}
            }else{
                nodes<-c(node,mrca)
            }

            pos<- binding.position(new.tree, node = sample(nodes,1), insertion = "random",prob = TRUE)
            new.tree <- phytools::bind.tip(new.tree,
                                          new.tips,
                                          edge.length = pos$length,
                                          where = pos$where,
                                          position = pos$position)
        }

    }

    return(new.tree)
}

phyleticity<- function(tree, genus){

    if(length(genus) != 1){stop("Only one genus accepted.") }

    sp <- tree$tip.label
    taxa.vector <- sp[first.word(sp)==genus]

    if(length(taxa.vector) == 0){
        genus.type <- "Missing"
        return(genus.type)
    }

    if(length(taxa.vector) == 1){
        genus.type <- "Singleton"
        return(genus.type)
    }

    mrca <- phytools::findMRCA(tree = tree, tips = taxa.vector)
    descend <- phytools::getDescendants(tree, mrca)
    desc.tips <- notNA(sp[descend])
    desc.genera <- first.word(desc.tips)

    if(length(unique(desc.genera)) == 1){
        genus.type <- "Monophyletic"
        return(genus.type)
    }
    intruder.tips <- desc.tips[desc.genera != genus]
    if(length(intruder.tips) == 1){
        # there is only one intruder: paraphyletic by singleton)
        genus.type <- "Paraphyletic"
    }else{
        intruder.mrca <- phytools::findMRCA(tree = tree, tips = intruder.tips)
        # if there are more than one intruder...
        desc.intruder.mrca <- phytools::getDescendants(tree, intruder.mrca)
        desc.intruder.tips <- sp[desc.intruder.mrca]
        desc.intruder.tips <- notNA(desc.intruder.tips)
        #intruders grouped: paraphyletic group by monophyletic intruder
        suppressWarnings({
            grouped.intruders <- all(sort(desc.intruder.tips) == sort(intruder.tips))
        })
        if(grouped.intruders){
            genus.type <- "Paraphyletic"
        }else{
            genus.type <- "Polyphyletic"
        }
    }

    return(genus.type)
}

MDCC.phyleticity<-function(input, tree, MDCC.info=list("rank"=NA, "MDCC"=NA),
                           trim=TRUE){

    rank<- MDCC.info$rank
    MDCC <- MDCC.info$MDCC
    input<-input[!is.na(input[,rank]),]

    if(rank=="genus"){
        MDCC.type<-phyleticity(tree, MDCC)
        return(MDCC.type)
    }

    tips<- tree$tip.label[first.word(tree$tip.label) %in% first.word(input$taxon)]
    if(isTRUE(trim) & length(tips)>0){
        tree<- ape::keep.tip(phy = tree, tip = tips)
    }

    species<- input[input[,rank]==MDCC,]
    MDCC.genera<- unique(first.word(species$taxon))
    genera.in.tree<-first.word(tree$tip.label)
    genera.in.tree<-MDCC.genera[MDCC.genera %in% genera.in.tree]
    spp.in.tree<- tree$tip.label[first.word(tree$tip.label) %in% genera.in.tree]

    if(length(spp.in.tree)==0){
        MDCC.type<-"Missing"
        return(MDCC.type)
    }

    if(length(spp.in.tree)==1){
        MDCC.type<-"Singleton"
        return(MDCC.type)
    }

    sp.mrca<- phytools::findMRCA(tree, tips = spp.in.tree)
    descs.num<- phytools::getDescendants(tree,sp.mrca)
    descs.name<-notNA(tree$tip.label[descs.num])

    if(all(first.word(descs.name)%in%MDCC.genera)){
        MDCC.type<-"Monophyletic"
        return(MDCC.type)
    }else{
        intruders<-descs.name[!(first.word(descs.name)%in%MDCC.genera)]
        if(length(intruders)==1){
            MDCC.type<-"Paraphyletic"
            return(MDCC.type)
        }

        intruders.mrca<- ape::getMRCA(phy = tree, tip = intruders)
        intruders.descs.num<- phytools::getDescendants(tree,intruders.mrca)
        intruders.descs.name<-notNA(tree$tip.label[intruders.descs.num])

        if(!any(first.word(intruders.descs.name)%in%MDCC.genera)){
            MDCC.type<-"Paraphyletic"
            return(MDCC.type)
        }else{
            MDCC.type<-"Polyphyletic"
            return(MDCC.type)
        }
    }
}
