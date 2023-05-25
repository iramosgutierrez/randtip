#### String manipulation ####
first_word<- function(string){
    return(stringr::str_extract(string, "[A-Za-z]+"))
}

second_word <- function(string){
    return(sapply(stringr::str_extract_all(string, "[A-Za-z\\-\\.]+"), `[`, 2))
}

remove_spaces <- function(species, sep = "_"){
    sp <- stringr::str_trim(gsub("\\s+", " ", species))
    sp <- stringr::str_replace_all(sp, " ", sep)
    return(sp)
}

notNA <- function(x){
  vect<- x[!is.na(x)]
  return(vect)
}

#### Auxiliary functions ####
randtip_ranks<- function(){
    return(as.vector(c("genus","subtribe","tribe",
                      "subfamily","family","superfamily",
                      "order","class")))
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

is_node<-function(tree, node){
    if(!(any(c(tree$edge[,1], tree$edge[,2]) == node))){
        stop("Node number is not in your tree")
    }
    #if(length(phytools::getDescendants(tree = tree, node = node, curr=NULL )) > 1){
    if(node > length(tree$tip.label)){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

is_tip <-function(tree, node){
    if(!(any(c(tree$edge[,1], tree$edge[,2]) == node))){
        stop("Node number is not in your tree")
    }
    #if(length(phytools::getDescendants(tree = tree, node = node, curr=NULL)) == 1){
    if(node <= length(tree$tip.label)){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

name_tree_nodes <- function(tree){
  if(is.null(tree$node.label)){tree$node.label <- rep("", times=tree$Nnode)}
  if(length(tree$node.label)==0){tree$node.label <- rep("", times=tree$Nnode)}
  tree$node.label[is.na(tree$node.label)]<-""
  tree$node.label[which(tree$node.label=="")] <- paste0("ON_",which(tree$node.label==""))
  if(any(duplicated(tree$node.label))){
    tree$node.label[duplicated(tree$node.label)] <- paste0(
      tree$node.label[duplicated(tree$node.label)], "_",
           1:sum(duplicated(tree$node.label)))
    }
  return(tree)
}

listnodes2realnodes <- function(listnodes, tree){
  realnodes <- c(which(tree$tip.label %in% listnodes),
                 (which(tree$node.label %in% listnodes)+length(tree$tip.label)))
  if(is.null(listnodes)){realnodes<-NULL}
  return(realnodes)
}

realnodes2listnodes <- function(realnodes, tree){

tips <- notNA(tree$tip.label[realnodes])
nodes <- realnodes[realnodes>length(tree$tip.label)]
nodes <- nodes-length(tree$tip.label)
nodes <- tree$node.label[nodes]
listnodes <- c(tips,nodes)
return(listnodes)

}

extend2ultrametric <- function(phy){
  if (is.null(phy$edge.length)){stop("the tree has no branch lengths")}
  n <- ape::Ntip(phy)
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  EL <- phy$edge.length
  ## xx: distance from a node or a tip to the root
  xx <- numeric(n + phy$Nnode)
  ## the following must start at the root and follow the
  ## edges contiguously; so the tree must be either in cladewise
  ## order (or in pruningwise but the for loop must start from
  ## the bottom of the edge matrix)
  for (i in seq_len(length(e1))){
    xx[e2[i]] <- xx[e1[i]] + EL[i]
  }
  xx.tip <- xx[1:n]
  xx.max <- max(xx.tip)
  xx.dif <- xx.max - xx.tip
  tip.pos <- which(phy$edge[,2] %in% 1:n)
  phy$edge.length[tip.pos] <- phy$edge.length[tip.pos] + xx.dif
  return(phy)
}

#### Specific functions  ####

correct_DF<- function(df){

  df_corrected <- apply(df, 2,
                        function(x){
                          x[x == ""] <- NA
                          x <- as.character(x)
                        })

  return(as.data.frame(df_corrected))
}

# Return a vector of all species in the tree that match a given genus
sp_genus_in_tree <- function(tree, genus){
    sp <- tree$tip.label
    taxa.vector <- sp[first_word(sp)%in%genus]

    return(taxa.vector)
}

inputfinder<- function(input, taxon, column){
    return(as.character(input[input$taxon==taxon, column]))
}

usingMDCCfinder<- function(input, taxon=NULL, tree, silent = FALSE){

    if(is.null(taxon)){taxon<- input$taxon}
    input<-correct_DF(input)
    MDCC.vect<- vector( mode="character", length = length(taxon))
    MDCC.lev.vect<- vector(mode="character", length = length(taxon))


    if(!silent){cat(paste0("Searching MDCCs\n"))}

    #manual MDCC search
    taxa <- input[!(input$taxon %in% tree$tip.label),]
    taxa<- taxa[!is.na(taxa$taxon1)|!is.na(taxa$taxon2),]
    if(nrow(taxa)>0){
        for(tx in seq_along(taxa$taxon)){
            if(length(strsplit(taxa$taxon1[tx], "_")[[1]])==1 &
               length(strsplit(taxa$taxon2[tx], "_")[[1]])==1 &
               taxa$taxon1[tx]==taxa$taxon2[tx] &
               length(sp_genus_in_tree(tree, taxa$taxon1[tx]))>0){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- taxa$taxon1[tx]
                MDCC.lev.vect[pos] <- "Sister genus"
                #MDCC.phyletictype.vect[pos]<-"-"
                next
            }

            if(!(any(tree$tip.label == taxa$taxon1[tx]))){taxa$taxon1[tx]<-NA}
            if(!(any(tree$tip.label == taxa$taxon2[tx]))){taxa$taxon2[tx]<-NA}

            if(any(tree$tip.label == taxa$taxon1[tx], na.rm = TRUE) &
               any(tree$tip.label == taxa$taxon2[tx], na.rm = TRUE) &
               taxa$taxon1[tx]==taxa$taxon2[tx]){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- taxa$taxon1[tx]
                MDCC.lev.vect[pos] <- "Sister species"

                next
            }

            if(any(tree$tip.label == taxa$taxon1[tx], na.rm = TRUE) &
               any(tree$tip.label == taxa$taxon2[tx], na.rm = TRUE) &
               taxa$taxon1[tx]!=taxa$taxon2[tx]){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- paste0("Clade ", taxa$taxon1[tx], "-", taxa$taxon2[tx])
                MDCC.lev.vect[pos] <- "Manual setting"

                next
            }

            if(any(tree$tip.label == taxa$taxon1[tx], na.rm = TRUE) &
               is.na(taxa$taxon2[tx])){

                pos<- which(taxon==taxa$taxon[tx])
                MDCC.vect[pos] <- taxa$taxon1[tx]
                MDCC.lev.vect[pos] <- "Sister species"
                next
            }

            if(any(tree$tip.label == taxa$taxon2[tx], na.rm = TRUE) &
               is.na(taxa$taxon1[tx])){

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

    if(nrow(taxa) == 0){
        return(list(MDCC=MDCC.vect,MDCC.ranks=MDCC.lev.vect) )
    }
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
                        #  phyleticity<-MDCC_phyleticity(input, tree = tree,
                        #          MDCC.info = list(rank=rank, MDCC= MDCC))
                        # if(phyleticity=="Missing"){MDCC<-NA}
                        #supressed for optimization

                        treegenera <- unique(first_word(tree$tip.label))
                        tree.input <- input[first_word(input$taxon)%in%treegenera,]
                        tree.input <- tree.input[!is.na(tree.input[,rank]),]
                        {if(sum(tree.input[, rank]==MDCC)==0){MDCC<-NA}}

                    }

                    lev<-rank
                }else{next}
            }
            MDCC.vect[v]<-as.character(MDCC)
            MDCC.lev.vect[v]<-as.character(lev)
            if(is.na(MDCC)){MDCC.lev.vect[v]<-NA}

        }
    }

    return(list(MDCC=MDCC.vect,MDCC.ranks=MDCC.lev.vect) )
}


get_parent_siblings <- function(tree, tip){
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

get_groups <- function(tree, genus){
    species <- sp_genus_in_tree(tree, genus)
    sp.mrca<- ape::getMRCA(tree, species)
    mrca.descs <- phytools::getDescendants(tree, sp.mrca,curr = NULL)

    node.descs<- rep(list(NA), times=length(mrca.descs))
    names(node.descs)<- mrca.descs

    node.types<- rep(list(NA), times=length(mrca.descs))
    names(node.types)<- mrca.descs

    for(n in seq_along(mrca.descs)){
        nd<- mrca.descs[n]

        if(is_tip(tree = tree,node = nd)){
            node.descs[[n]]<- tree$tip.label[nd]
            node.types[[n]]<- "tip"
            next
        }

        nd.descs<- phytools::getDescendants(tree, nd,curr = NULL)
        nd.descs <- notNA(tree$tip.label[nd.descs])
        nd.genera<- first_word(nd.descs)

        if(!(genus%in%nd.genera)){
            node.descs[[n]]<- NA
            node.types[[n]]<- NA
            next
        }


        siblings<-get_parent_siblings(tree, tip = nd)$siblings
        siblings<-siblings[-which(siblings%in%nd.descs)]
        siblings.genus<- first_word(siblings)
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
        intruders.genera<- first_word(intruders.descs)

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
        group.descs.gen <- first_word(group.descs)
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

get_position <- function(tree, node){
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

binding_position<- function(tree, node,  insertion,  prob, ultrametric = FALSE){
    position<-list("length"=NA, "where"=NA, "position"=NA)
    df <- data.frame("parent"=tree$edge[,1], "node"=tree$edge[,2],
                     "length"= tree$edge.length, "id"=1:length(tree$edge[,1]))

    if(ultrametric){
        position$length<-NULL
    }else{
      tiplength <- tree$edge.length[which(tree$edge[,2] %in% 1:length(tree$tip.label))] #we use tip lengths
      position$length<-abs(stats::rexp (n = 1, rate = 1 / mean (tiplength))) #negative binomial value
      }

    df<- df[df$node==node,]
    position$where<- node

    if(insertion=="polytomy"){
        position$position<- 0
        position$where <- node
    }

    if(insertion=="random"){

        index<- df[df$node==node,"id"]
        pos <- get_position(tree = tree, node = node)

        position$position<- pos
        position$where <- df[df$id==index,"node"]
    }

    return(position)
}

sharingtaxa_descs<-function(tree, nodes, MDCC.genera){
    table<- data.frame("node"=as.numeric(nodes), "number"=rep(NA, length(nodes)),
                       "tot.number"=rep(NA, length(nodes)))
    for(i in seq_along(table$node)){
        node<- table$node[i]
        if(is_tip(tree,node)){
            table$number[i]<-1;table$tot.number[i]<-1
            next
        }
        descs<- phytools::getDescendants(tree, node,curr = NULL)
        table$tot.number[i]<-length(descs)
        descs<- tree$tip.label[descs]
        descs<- notNA(descs)
        descs<- descs[first_word(descs)%in%MDCC.genera]
        table$number[i]<-length(descs)
    }
    table<- table[table$number>0,]

    return(table)
}

#given the group and specifications, which nodes can be selected to bind a tip over
get_permitted_nodes <- function (tree, input, MDCC, rank, MDCC.type,
                                 polyphyly.scheme, use.paraphyletic, use.singleton, use.stem){

    input.mdcc <-  input[!is.na(input[,rank]),]
    MDCC.taxa  <- input.mdcc$taxon[input.mdcc[,rank]==MDCC]
    MDCC.genera<- notNA(unique(first_word(MDCC.taxa)))
    MDCC.intree<- sp_genus_in_tree(tree, MDCC.genera)

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
            non.mdcc.genera<- unique(first_word(non.mdcc$taxon))

            intruder.descs<- descendants.tips[!(first_word(descendants.tips)%in%MDCC.genera)]
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
            nodes<- phytools::getDescendants(tree, get_parent_siblings(tree, tip)$parent, curr=NULL)
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
                if(is_tip(tree,node)){
                    table$descs[i] <- NA
                    table$total.descs[i]<-1
                    if(tree$tip.label[node] %in% input[input[,rank]==MDCC,"taxon"]){
                        table$sharing.descs[i]<-1;table$eligible[i]<- "TRUE"
                    }else{
                            table$sharing.descs[i]<-0;table$eligible[i]<- "FALSE"
                    }
                }
                if(is_node(tree,node)){
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
                if(is_tip(tree,node)){
                    table$descs[i]   <- NA
                    table$total.descs[i]<-1
                    if(tree$tip.label[node]%in% input[input[,rank]==MDCC,"taxon"]){
                        table$sharing.descs[i]<-1;table$eligible[i]<- "TRUE"
                    }else{
                        table$sharing.descs[i]<-0;table$eligible[i]<- "FALSE"
                    }
                }
                if(is_node(tree,node)){
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

#given a set of nodes with get.permitted.nodes, must some of them should not be considered?
get_forbidden_nodes <- function(tree,input, MDCC, rank, perm.nodes, respect.mono, respect.para){

    forbidden.nodes<- vector("numeric")
    #groups which must not be forbid
    perm.groups<- input[input[,rank]==MDCC,randtip_ranks()[(which(randtip_ranks()==rank)):8]]
    perm.groups<- notNA(unique(as.vector(unlist(perm.groups))))

    #respect monophyletic clusters of species
    if(respect.mono){
        perm.tips<- notNA(tree$tip.label[perm.nodes])
        perm.species<-paste0(first_word(perm.tips), "_",
                             second_word(perm.tips))
        ssps<-perm.species[duplicated(perm.species)]
        if(length(ssps)>0){
            for(ssp in ssps){
                nodes<-which(paste0(first_word(tree$tip.label), "_",
                                    second_word(tree$tip.label))==ssp)
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
                        node.descs<- paste0(first_word(node.descs), "_",
                                            second_word(node.descs))
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
            if(is_tip(tree, nd)){next}
            descs.nd<- phytools::getDescendants(tree, nd)
            descs <- notNA(tree$tip.label[descs.nd])
            genera<- unique(first_word(descs))

            #respect monophyletic genera, disregarding input info
            if(length(genera)==1){
                if(!(genera%in%perm.groups)){
                    forbidden.nodes<- c(forbidden.nodes, descs.nd )
                    next
                }
            }

            #respect monophyletic groups, according to input info
            sub.input <- input[first_word(input$taxon)%in%genera,]

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
            if(is_tip(tree, nd)){next}
            descs.nd<- phytools::getDescendants(tree, nd, curr=NULL)
            descs <- notNA(tree$tip.label[descs.nd])
            genera<- unique(first_word(descs))

            if(length(genera)==1){next}

            #account for singleton tips in otherwise monophyletic clusters
            if(length(genera)==2){
                table<- table(first_word(descs))
                uniq.gen<- names(table)[which(table==1)]
                nest.gen<- names(table)[which(table!=1)]
                if(length(uniq.gen)==1){
                    tip<-descs[first_word(descs)%in%uniq.gen]
                    tip.n<- which(tree$tip.label==tip)

                    tip.par<- get_parent_siblings(tree, tip.n)$parent

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
                ds <- descs[first_word(descs)%in%v]
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

                    if(!any(first_word(nested)%in%nest)){
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
            sub.input <- input[first_word(input$taxon)%in%genera,]

            for( rk in randtip_ranks()){
                rk.vals<-notNA(unique(sub.input[,rk]))

                if(length(rk.vals)>1){
                    rk.vals.mrca<- vector("numeric", length = length(rk.vals))
                    rk.vals.desc<- list( NA)
                    for(v in rk.vals){
                        gen <- unique(first_word(sub.input$taxon[sub.input[,rk]==v]))
                        ds <- descs[first_word(descs)%in%gen]
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
                        if(MDCC_phyleticity(input, tree, list(rank=rk, MDCC=nest))=="Paraphyletic"&
                           !(nest %in% perm.groups)){

                            nest.tips<- sub.input$taxon[sub.input[,rk]==nest]
                            nest.gen<- first_word(nest.tips)
                            nest.mrca<- nd

                            nested.nodes<- phytools::getDescendants(tree, node=nest.mrca,curr = NULL)
                            nest.desc.tips <- tree$tip.label[nested.nodes]
                            nest.desc.tips<- notNA(nest.desc.tips)

                            non.nest <- input[!is.na(input[,rank]),]
                            non.nest <- non.nest[non.nest[,rk]!=nest,]
                            non.nest.genera<- unique(first_word(non.nest$taxon))

                            intruders <- nest.desc.tips[first_word(nest.desc.tips)%in%non.nest.genera]
                            if(length(intruders)==1){
                                intruder.mrca<- which(tree$tip.label==intruders)
                            }else{
                                intruder.mrca<- ape::getMRCA(tree, intruders)
                            }

                            if(any(first_word(intruders) %in% nest.gen)){next}


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


bind_clump<- function(new.tree, tree, input, PUT){

    clumplist<- list("MDCC"=NULL, "rank"=NULL,"MDCC.type"=NULL, "taxa"=NULL)
    sp<- paste0(first_word(PUT), "_", second_word(PUT))
    treespp<- paste0(first_word(tree$tip.label), "_",
                     second_word(tree$tip.label))
    new.treespp<- paste0(first_word(new.tree$tip.label), "_",
                         second_word(new.tree$tip.label))

    if(sp%in%new.treespp){
        DFspp<- paste0(first_word(input$taxon), "_",
                       second_word(input$taxon))
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


    tree.genus.tips<-sp_genus_in_tree(tree, first_word(PUT))
    new.tree.genus.tips<-sp_genus_in_tree(new.tree, first_word(PUT))
    if(length(tree.genus.tips)==0 & length(new.tree.genus.tips)>0){
        clumplist$MDCC<-first_word(PUT)
        clumplist$rank<-"genus"
        clumplist$MDCC.type<- phyleticity(new.tree, first_word(PUT))
        clumplist$taxa<-new.tree.genus.tips

        return(clumplist)
    }

    newsearch<- usingMDCCfinder(input, PUT, new.tree, silent=TRUE)
    MDCC<-input[input$taxon==PUT,"MDCC"]
    MDCC.rank<- input[input$taxon==PUT,"MDCC.rank"]

    if(newsearch$MDCC!=MDCC){
        newclump<- input$taxon[input[,newsearch$MDCC.ranks]==newsearch$MDCC]
        newclump<- unique(first_word(newclump))
        newclump <- sp_genus_in_tree(new.tree, newclump)

        clumplist$MDCC<-newsearch$MDCC
        clumplist$rank<-newsearch$MDCC.ranks
        clumplist$taxa<-newclump
        clumplist$MDCC.type<- MDCC_phyleticity(input, new.tree,
                                  MDCC.info = list(rank=newsearch$MDCC.ranks,
                                                   MDCC=newsearch$MDCC, rank), trim = T)

    }

    return(clumplist)
}


add_to_singleton <- function(tree, singleton, new.tips, use.singleton=F){
    singleton<-gsub(" ", "_", singleton)
    singleton<-singleton[singleton%in%tree$tip.label]

    new.tree <- tree

    if(length(singleton)==1){

        nodes<-which(new.tree$tip.label==singleton)
        if(!(use.singleton)){
          if(!(isRoot(new.tree, nodes))){
            parent<- get_parent_siblings(new.tree, nodes)[[1]]
            nodes <- c(nodes, parent)}
          }

     }

    if(length(singleton)> 1){
        nodes<-which(new.tree$tip.label%in%singleton)
        mrca<- ape::getMRCA(new.tree, singleton)
        nodes<- c(mrca, nodes)
    }
    adding.DF<- data.frame("parent"=new.tree$edge[,1], "node"=new.tree$edge[,2],
                           "length"= new.tree$edge.length )
    adding.DF<- adding.DF[adding.DF$node %in% nodes,]
    adding.DF$id<- 1:nrow(adding.DF)

    if(length(nodes)>1){nodes<-sample(adding.DF$node, 1, prob = adding.DF$length)}

    ultrametric <- ape::is.ultrametric(tree)


      pos<- binding_position(new.tree, node = nodes, insertion = "random",
                             prob = T, ultrametric = ultrametric)
      new.tree <- phytools::bind.tip(new.tree,
                                     new.tips,
                                     edge.length = pos$length,
                                     where = pos$where,
                                     position = pos$position)



    return(new.tree)
}

phyleticity<- function(tree, genus){

    if(length(genus) != 1){stop("Only one genus accepted.") }

    sp <- tree$tip.label
    taxa.vector <- sp[first_word(sp)==genus]

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
    desc.genera <- first_word(desc.tips)

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

#function to obtain the phyletic nature of a group (genus or above)
MDCC_phyleticity<-function(input, tree, MDCC.info=list("rank"=NA, "MDCC"=NA),
                           trim=TRUE){

    rank<- MDCC.info$rank
    MDCC <- MDCC.info$MDCC
    input<-input[!is.na(input[,rank]),] #NAs must be erased to avoid errors

    if(rank=="genus"){
        MDCC.type<-phyleticity(tree, MDCC) #normal genus-level function
        return(MDCC.type)
    }

    tips<- tree$tip.label[first_word(tree$tip.label) %in% first_word(input$taxon)]
    if(trim & length(tips)>0){
        tree<- ape::keep.tip(phy = tree, tip = tips)
    }

    species<- input[input[,rank]==MDCC,]
    MDCC.genera<- unique(first_word(species$taxon))
    genera.in.tree<-first_word(tree$tip.label)
    genera.in.tree<-MDCC.genera[MDCC.genera %in% genera.in.tree]
    spp.in.tree<- tree$tip.label[first_word(tree$tip.label) %in% genera.in.tree]

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

    if(all(first_word(descs.name)%in%MDCC.genera)){
        MDCC.type<-"Monophyletic"
        return(MDCC.type)
    }else{
        intruders<-descs.name[!(first_word(descs.name)%in%MDCC.genera)]
        if(length(intruders)==1){
            MDCC.type<-"Paraphyletic"
            return(MDCC.type)
        }

        intruders.mrca<- ape::getMRCA(phy = tree, tip = intruders)
        intruders.descs.num<- phytools::getDescendants(tree,intruders.mrca)
        intruders.descs.name<-notNA(tree$tip.label[intruders.descs.num])

        if(!any(first_word(intruders.descs.name)%in%MDCC.genera)){
            MDCC.type<-"Paraphyletic"
            return(MDCC.type)
        }else{
            MDCC.type<-"Polyphyletic"
            return(MDCC.type)
        }
    }
}
