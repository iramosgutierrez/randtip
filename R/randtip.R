#PREPARATION####
# Load packages

library(taxize)
library(openxlsx)
library(ape)
library(stringr)
library(phytools)

rm(list = setdiff(ls(), lsf.str()))#eliminar todo menos funciones

#checklist<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/IberoBalearic_Checklist_R1.xlsx")
#checklist$match.name<- gsub(" ","_", checklist$match.name)
#tipos_generos<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/Tipos_generos_Filet.xlsx")
#genus.corrections<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/Corrections.xlsx")
#taxonomy<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/Taxonomic_categories.xlsx")

Tree<- read.tree("data/25tree.tre")
tabla.info<- read.xlsx("data/phylo.table.xlsx")

plot(Tree, label.offset = 2)
nodelabels()
tiplabels(adj = c(-10, 0.5))

#CREACION DE FUNCIONES####


'%!in%' <- function(x,y)!('%in%'(x,y))  #funcion opuesta de %in%


# function phylo.taxonomy returns a table with family, order and class for every genus gien in the species vector.
phylo.taxonomy<- function(species, db="ncbi"){
  message(paste0("Taxonomic information is being looked up in ",db," Taxonomy Database. This may take a while!"))
  taxa.genera<- word(species, 1, sep="_")[!duplicated(word(species, 1, sep="_"))]     #genera are extracted
  taxonomy.df<- data.frame("Genus"=taxa.genera, "Family"=NA, "Order"=NA, "Class"=NA)  #table is created
  for(i in 1:length(taxa.genera)){    #for every genus...
    tryCatch({
      search<-classification(as.character(taxa.genera[i]), db = db)[[1]]  #taxonomic categories are searched and saved in the table
      fam<-search[search$rank=="family","name"]
      ord<-search[search$rank=="order" ,"name"]
      cla<-search[search$rank=="class" ,"name"]
      taxonomy.df$Family[taxonomy.df$Genus==taxa.genera[i]]<-fam
      taxonomy.df$Order[taxonomy.df$Genus==taxa.genera[i]] <-ord
      taxonomy.df$Class[taxonomy.df$Genus==taxa.genera[i]] <-cla
    }, error=function(e){                                               #errors are saved as NA to avoid stopping
      taxonomy.df$Family[taxonomy.df$Genus==taxa.genera[i]]<-NA
      taxonomy.df$Order[taxonomy.df$Genus==taxa.genera[i]] <-NA
      taxonomy.df$Class[taxonomy.df$Genus==taxa.genera[i]] <-NA
    })
    Sys.sleep(0.33)         #taxize allows only 3 searches per second; this avoids collapsing
  }
  return(taxonomy.df)
}
# phylo.taxonomy(species=c("Abies_pinsapo", "Pinus_pinaster", "Achyllea_milleifolium"))


#functionn phyleticity retrieves genus type given a certain tree
phyleticity<- function(tree, genus){
  if(length(genus)!=1){stop("Only one genus accepted")}

  list<- data.frame(species=tree$tip.label)                                 #lists trees species
  taxa.vector<- as.vector(list[word(list$species, 1, sep="_")==genus,])     #lists trees genera
  if(length(taxa.vector)==0){                                               #if the genus is not in the tree
    type<-"NOT INCLUDED"
    #message(paste0("Genus ", genus, " is ", type, " in your tree"))
  }else{
    if(length(taxa.vector)==1){                                             #if this genus appears in the tree once
      type<-"SINGLETON GENUS"
      #message(paste0("Genus ", genus, " is a ", type))
    }else{                                                                 #if this genus appears in the tree more than once..

      MRCA<-findMRCA(tree = tree,tips = taxa.vector)                       #...MRCA is searched...
      Desc.Tips<-tree$tip.label[getDescendants(tree, MRCA)][!is.na(tree$tip.label[getDescendants(tree, MRCA)])] #...descendants are searched...
      Desc.Genera<-word(Desc.Tips, 1, sep="_")[!duplicated(word(Desc.Tips, 1, sep="_"))] #if all MCRA descendants are fron the same genus: monophyletic.
      if(length(Desc.Genera)==1){
        type<-"MONOPHYLETIC"
        #message(paste0("Genus ", genus, " is ", type))
      }else{
        Intruder.Tips<-as.vector(Desc.Tips[word(Desc.Tips, 1, sep = "_")!= genus]) #if there is an intruder...
        Intruder.MRCA<-findMRCA(tree = tree,tips = Intruder.Tips)
        if(length(Intruder.Tips)==1){type<-"PARAPHYLETIC"                          #(and there is only one intruder: paraphyletic by singleton)
        #message(paste0("Genus ", genus, " is ", type))
        }else{                                                                     # if there are more than one intruder...
          Desc.Intruder.Tips<-tree$tip.label[getDescendants(tree, Intruder.MRCA)][!is.na(tree$tip.label[getDescendants(tree, Intruder.MRCA)])] #descendants and intruders MRCA

          if(length(Desc.Intruder.Tips)==length(Intruder.Tips)&all(sort(Desc.Intruder.Tips)==sort(Intruder.Tips))){        #intruders grouped: paraphyletic group by monophyletic intruder
            type<-"PARAPHYLETIC"
            # message(paste0("Genus ", genus, " is ", type))
          } else{                                                                 #if everything else has failed, it has to be polyphyletic
            type<-"POLYPHYLETIC"
            #message(paste0("Genus ", genus, " is ",type))
          }}}}}
  return(type)
}
#phyleticity(Tree, "Inventum")


#function phylo.table returns a table which the main function can work with (editable by hand!).
phylo.table<- function(species, tree=NULL, taxonomy=NULL, genera.phyleticity=NULL){
if(is.null(genera.phyleticity)&(!inherits(tree, "phylo"))){   stop("tree should be an object of class \"phylo\".")}
  phylo.df<- data.frame(matrix(nrow = length(species),ncol=10))
  names(phylo.df)<-c("taxon", "genus", "genus.type", "family", "order", "class", "aggregate.subspecies",
                     "relative.species","synonim.genus","sibling.genus")

  if(length(species[duplicated(species)])){stop("There are duplicated species")}
  phylo.df$taxon<-species                                     #species listed on the species column
  phylo.df$genus<-word(phylo.df$taxon, 1, sep="_")            #genera obtention

  taxa.genera<- phylo.df$genus[!duplicated(phylo.df$genus)]   #genera list

  if(is.null(taxonomy)){  taxonomy<-phylo.taxonomy(species=species) } #if taxonomy is not given, it is searched for.

  for(i in 1:length(taxa.genera)){
    genus<-as.character(taxa.genera[i])
    fam<-taxonomy$Family[taxonomy$Genus==genus ]
    ord<-taxonomy$Order [taxonomy$Genus==genus ]
    cla<-taxonomy$Class [taxonomy$Genus==genus ]
    phylo.df$family[phylo.df$genus==taxa.genera[i]]<-fam
    phylo.df$order[phylo.df$genus==taxa.genera[i]] <-ord
    phylo.df$class[phylo.df$genus==taxa.genera[i]] <-cla

  }



  if(is.null(genera.phyleticity)){                 #if phyleticity is not given, it is searched for.
    if(is.null(tree)){stop("Tree object required")}
    genera.phyleticity<- data.frame("Genus"=taxa.genera, "Type"=NA)
    for(f in 1:nrow(genera.phyleticity)){
      genera.phyleticity$Type[f]<- phyleticity(tree,  genera.phyleticity$Genus[f])
    }}

  for(i in 1:length(taxa.genera)){         #obtained before, or given, data are imputted to our table.
    genus<-as.character(taxa.genera[i])
    genus.type<-genera.phyleticity$Type[genera.phyleticity$Genus==genus]
    phylo.df$genus.type[phylo.df$genus==genus]<-genus.type
  }



  return(phylo.df)
}
#df<-phylo.table(species = checklist$match.name, Tree,  genera.phyleticity = tipos_generos, taxonomy = taxonomy)

#this function reveals if a given node at a given tree is an internal node
is.node<-function(tree, node){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(node %!in% tree$edge){stop("Node number is not in your tree")}
  if (length(getDescendants(tree=tree, node = node))>1){return(TRUE)}else{return(FALSE)}
}

#this function reveals if a given node at a given tree is a tip
is.tip <-function(tree, node){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(node %!in% tree$edge){stop("Node number is not in your tree")}
  if (length(getDescendants(tree=tree, node = node))==1){return(TRUE)}else{return(FALSE)}
}

#modification from phytools "add.random". if prob=TRUE, it remains the same; else a random node is selected without probability vector
add.random<- function (tree, n = NULL, tips = NULL, edge.length = NULL, order = c("random", "input"), prob=TRUE) {
  if (!inherits(tree, "phylo"))
    stop("tree should be an object of class \"phylo\".")

  if(prob==TRUE){randomPosn <- function(tree) {
    cum.edge <- cumsum(tree$edge.length)
    index <- tree$edge[, 2]
    pos <- runif(1) * sum(tree$edge.length)
    edge <- 1
    while (pos > cum.edge[edge]) edge <- edge + 1
    return(list(node = index[edge], posn = cum.edge[edge] - pos))}

  if (is.ultrametric(tree))
    um <- TRUE
  else um <- FALSE
  if (is.null(tips)) {
    if (is.null(n))
      n <- 1
    tips <- paste("t", length(tree$tip) + 1:n, sep = "")
  }
  else n <- length(tips)
  if (is.null(edge.length))
    if (!um)
      edge.length <- runif(n = n, min = min(tree$edge.length),
                           max = max(tree$edge.length))
  if (order[1] == "random") {
    o <- sample(1:n)
    tips <- tips[o]
    if (!is.null(edge.length))
      edge.length <- edge.length[o]
  }
  for (i in 1:n) {
    where <- randomPosn(tree)
    if (is.null(edge.length))
      tree <- bind.tip(tree, tips[i], where = where$node,
                       position = where$posn)
    else tree <- bind.tip(tree, tips[i], where = where$node,
                          position = where$posn, edge.length = edge.length[i])
  }

  }else{
        if(order!="input"){tips<- sample(tips, length(tips), replace = F)}

    for(t in 1: length(tips)){
    new.tip<- tips[t]
    DF <- data.frame(tree$edge,tree$edge.length,1:length(tree$edge.length)) #tree data
    EDGES <- DF[sample(x=1:nrow(DF), size = 1),] # select a random node
    to_index<-EDGES[,4] #indexing position

    WHERE <- tree$edge[to_index,2]   #indexing position
    LENGTH <- tree$edge.length[to_index]   #maximum branch length

    MIN<-0
    MAX<-LENGTH
    POS<-runif(1,MIN,MAX)

    while(POS==MIN | POS==MAX){POS<- runif(1,MIN,MAX)} #intermediate value

    tree<-bind.tip(tree, new.tip, edge.length=NULL, where=WHERE, position=POS)
  }}
  return(tree)
}


# Function to add tips at random from a given node
add_into_node <- function(tree,node,new.tip) {
  new.tip<-gsub(" ", "_",new.tip) #modification in case names are separated with blanks
  tt<-splitTree(tree,split=list(node=node, bp=tree$edge.length[which(tree$edge[,2]==node)])) #tree is splitted at a given node
  if (is.ultrametric(tt[[2]])==TRUE){
    tt[[2]]<-add.random(tt[[2]],tips=new.tip)   #randomly add into the cut tree
    new.tree<-paste.tree(tt[[1]],tt[[2]])       #trees are pasted together
    return(new.tree)
  } else
  {
    warning("Found issue in sticking, forcing ultrametric")
    tt[[2]]<-force.ultrametric(tt[[2]], method="extend")
    tt[[2]]<-add.random(tt[[2]],tips=new.tip)
    new.tree<-paste.tree(tt[[1]],tt[[2]]) }

 return(new.tree)}

# Function to add tips at random from a given node WITH FORBIDDEN BRANCHES
#exception.list<-list(c("Abies_alba","Abies_pinsapo"),c("Achillea_pyrenaica","Achillea_santolinoides","Achillea_chamaemelifolia"))
add_into_node_exceptions <- function(tree,node,new.tip, exception.list) {
  new.tip<-gsub(" ", "_",new.tip) #modification in case names are separated with blanks
  tt<-splitTree(tree,split=list(node=node, bp=tree$edge.length[which(tree$edge[,2]==node)])) #tree is splitted at a given node

  tree1<-tt[[1]]
  tree2<-tt[[2]]
forbiddentaxa<- Reduce(c,exception.list)

if(all(tree2$tip.label %in% forbiddentaxa)==TRUE){ #in case all groups arein a borbidden clade
 randsp<-sample(tree2$tip.label, 1) #one species is selected
 randgr<-exception.list[[grep(randsp,exception.list ) ]]#its group is selected
 new.tree<- add_over_node(tree, new.tip = new.tip, node = findMRCA(tree, tips = randgr))
}else{
  DF <- data.frame(tree2$edge,1:dim(tree2$edge)[1],tree2$edge.length)
  colnames(DF) <- c("Parent Node","Child Node","Edge ID","Edge Length")

  # Workflow

  List_Edge <- vector(mode="list",length=length(exception.list)) # En esta lista voy a almacenar la identidad de las ramas prohibidas

  for (i in 1:length(exception.list)) {

    Group <- exception.list[[i]]  # Saco el grupo correspondiente

    ID_tip <- NA

    for(k in 1:length(Group)) {
      ID_tip <- c(ID_tip,which(tree2$tip.label==Group[k])) # Números que R asigna a las especies
    }

    ID_tip <- ID_tip[!is.na(ID_tip)]

    KK_Pare_tip <- NA

    for(j in 1:length(ID_tip)) {  # Con este bucle voy a identificar las ramas terminales del grupo, así como la identidad de los nodos
      #internos a partir de los cuales vamos a extraer los identificadores de las ramas prohibidas en la aleatorización

      KK_Pare_tip <- c(KK_Pare_tip,getParent(tree2, ID_tip[j]))

    }

    KK_Pare_tip <- unique(KK_Pare_tip[!is.na(KK_Pare_tip)])

    List_Edge[[i]] <- DF[DF[,2] <= max(KK_Pare_tip) & DF[,2] > getMRCA(tree2,ID_tip) | DF[,2] %in% ID_tip,3] # Si borras el "3" de la línea de código, verás que el trozo de DF que incluye las ramas prohibidas del grupo correspondiente

  }


  DF_to_stick <- DF[DF[,3] %in% unlist(List_Edge) == FALSE,] # Subset de DF incluyendo sólo las ramas a aleatorizar

  to_index <- sample(DF_to_stick[,3],1,prob=DF_to_stick[,4]) # selecciono una rama con probabilidad proporcional a la longitud de las mismas

  WHERE = DF[to_index,2] # Nodo (i.e. rama) en el que voy a insertar la especie (puede ser interno o terminal)

  LENGTH <- tree2$edge.length[to_index] # Longitud de la rama seleccionada para la inserción

  MIN=0
  MAX=LENGTH
  runif(1,MIN,MAX)->POS # Seleccionamos una posición dentro de la rama a partir de distribución uniforme

  while(POS==MIN | POS==MAX){ # Bucle para asegurarnos que la posición de inserción no coincide con el valor mínimo o máximo de la longitud de la rama
    runif(1,MIN,MAX)->POS
  }

  tree2 <- bind.tip(tree2, new.tip, edge.length=NULL, where=WHERE, position=POS)

  new.tree<-paste.tree(tree1, tree2)}

  return(new.tree)}


# Function to add tips at random monophyletic clade
add_to_monophyletic <- function(tree,new.tip) {
  new.tip<-gsub(" ", "_",new.tip)       #modification in case names are separated with blanks
  genus<- word(new.tip, 1, sep="_")     #genus obtention (it has to be previously checked to be monophyletic)
  list<- data.frame(species=tree$tip.label)  #species list inside the tree
  taxa.vector<- as.vector(list[word(list$species, 1, sep="_")==genus,]) #species list inside the tree from the genus
  if(length(taxa.vector)==0){stop(paste0("Genus ", genus, " is not included in yout tree."))}
  MRCA<-findMRCA(tree = tree,tips = taxa.vector) #genus MCRA
  new.tree<-add_into_node(tree=tree, node = MRCA, new.tip = new.tip)  #Added into th tree defined by MCRA
  return(new.tree)}
#plot(add_to_monophyletic(Tree, "Invent_inventata"))


# Function to add tips over a node. Valid for "sibling" species or genera.
add_over_node<- function(tree, new.tip, node){
  DF <- data.frame(tree$edge,tree$edge.length,1:length(tree$edge.length)) #tree data
  EDGES <- DF[which.edge(tree, node),] # node especifications
  to_index<-EDGES[,4] #indexing position

  WHERE <- tree$edge[to_index,2]   #indexing position
  LENGTH <- tree$edge.length[to_index]   #maximum branch length

  MIN<-0
  MAX<-LENGTH
  POS<-runif(1,MIN,MAX)

  while(POS==MIN | POS==MAX){POS<- runif(1,MIN,MAX)} #intermediate value

  newtree<-bind.tip(tree, new.tip, edge.length=NULL, where=WHERE, position=POS)
  return(newtree)
}

# A function to stick species at random within a paraphyletic clade
add_into_paraphyletic_node<- function(tree, new.tip, group.node, intern.node){

  intern.tips<- tree$tip.label[getDescendants(tree, intern.node)][!is.na(tree$tip.label[getDescendants(tree, intern.node)])] #tips in intruder node

  group.tree<-splitTree(tree,split=list(node=group.node,bp=tree$edge.length[which(tree$edge[,2]==group.node)])) #tree is cut at the deep node

  new.intern.node<- findMRCA(tree = group.tree[[2]],tips = intern.tips) #node is searched in the new tree

  intern.tree <-splitTree(group.tree[[2]],split=list(node=new.intern.node,bp=group.tree[[2]]$edge.length[which(group.tree[[2]]$edge[,2]==new.intern.node)])) #cut tree is cut again

  intern.tree[[1]]<-add.random(intern.tree[[1]],tips=new.tip) #The new tip is added to the cut tree without the intruder group

  group.tree[[2]]<-paste.tree( intern.tree[[1]], intern.tree[[2]] )
  group.tree[[2]]$root.edge<-0
  new.tree<- paste.tree( group.tree[[1]], group.tree[[2]]) #3 trees are pasted together


  if (is.ultrametric(new.tree)==FALSE){new.tree<-force.ultrametric(new.tree, method="extend")}
  return(new.tree)
}


add_to_paraphyletic <- function(tree,new.tip){
  new.tip<-gsub(" ", "_",new.tip) #Same process than paraphyletic node, bt detecting automatically the deep and shallow nodes
  genus<- word(new.tip, 1, sep="_")
  sp.list<- data.frame(species=tree$tip.label)
  taxa.vector<- as.vector(sp.list[word(sp.list$species, 1, sep="_")==genus,])
  if(length(taxa.vector)==0){stop(paste0("Genus ", genus, " is no included in yout tree."))}
  MRCA<-findMRCA(tree = tree,tips = taxa.vector)
  Desc.Tips<-tree$tip.label[getDescendants(tree, MRCA)][!is.na(tree$tip.label[getDescendants(tree, MRCA)])]
  Intruder.Tips<-as.vector(Desc.Tips[word(Desc.Tips, 1, sep = "_")!= genus])


  genus.tree<-splitTree(tree,split=list(node=MRCA,bp=tree$edge.length[which(tree$edge[,2]==MRCA)]))
if(length(Intruder.Tips)>1){
  Intruder.MRCA<- findMRCA(tree = genus.tree[[2]],tips = Intruder.Tips)

  pure.tree <-splitTree(genus.tree[[2]],split=list(node=Intruder.MRCA,bp=genus.tree[[2]]$edge.length[which(genus.tree[[2]]$edge[,2]==Intruder.MRCA)]))

  pure.tree[[1]]<-add.random(pure.tree[[1]],tips=new.tip)

  genus.tree[[2]]<-paste.tree( pure.tree[[1]], pure.tree[[2]] )
  genus.tree[[2]]$root.edge<-0
  new.tree<- paste.tree( genus.tree[[1]], genus.tree[[2]])}  else{
    Intruder.node<- which(genus.tree[[2]]$tip.label== Intruder.Tips)

    pure.tree <-splitTree(genus.tree[[2]],split=list(node=Intruder.node,bp=0))

    pure.tree[[1]]<-add.random(pure.tree[[1]],tips=new.tip)

    genus.tree[[2]]<-paste.tree( pure.tree[[1]], pure.tree[[2]] )
    genus.tree[[2]]$root.edge<-0
    new.tree<- paste.tree( genus.tree[[1]], genus.tree[[2]])

  }


  if (is.ultrametric(new.tree)==FALSE){new.tree<-force.ultrametric(new.tree, method="extend")}




  return(new.tree)}
#plot(add_to_paraphyletic(Tree, "Achillea_3_inventatus"))



polytomy_into_node<- function(tree, new.tip, node){

  if(!exists("tree")){stop("Tree object not found")}
  newtree<-tree #new tree is created to be modified

  DF <- data.frame(newtree$edge,newtree$edge.length,"id"=1:length(newtree$edge.length)) #DF characteistics
  EDGES <- DF[which.edge(newtree, node),] #given node row
  to_index<-EDGES[,4] #indexing position

  WHERE <- newtree$edge[to_index,2]   #indexing position



  newtree<-bind.tip(newtree, new.tip, edge.length=NULL, where=WHERE, position=0) #indexing position; node length is 0
  return(newtree)
}

polytomy_over_node<- function(tree, node, species, insertion=c("random","long","middle") ){

  newtree<-tree #tree to be modified
  for(i in 1:length(species)){
    if(i==1){   #first species is added to singleton

      DF <- data.frame(newtree$edge,newtree$edge.length,"id"=1:length(newtree$edge.length))
      EDGES <- DF[which.edge(newtree, node),]
      to_index<-EDGES[,4]
      WHERE <- newtree$edge[to_index,2]

      if(insertion=="random"){position<-EDGES[,3]*runif(1, 0, 1)} #bound at a random point of the branch
      if(insertion=="middle"){position<-EDGES[,3]*0.5}  #bound at the middle of the branch
      if(insertion=="long")  {position<-EDGES[,3]}      #bound at the beggining of the branch

      newtree<-bind.tip(newtree, tip.label =  species[i], edge.length=NULL, where=WHERE, position=position)
    }else{
      node<-getParent(tree = newtree, node=which (newtree$tip.label==species[i-1])) #the rest of the species are added to the node
      newtree<- polytomy_into_node(tree = newtree, new.tip  =  species[i], node = node)
    }}

  return(newtree)}

polytomy_to_singleton<- function(tree, singleton, species, insertion=c("long","middle") ){

  newtree<-tree #tree to be modified
  for(i in 1:length(species)){
    if(i==1){   #first species is added to singleton
      sing.node<-which(newtree$tip.label==singleton)
      DF <- data.frame(newtree$edge,newtree$edge.length,"id"=1:length(newtree$edge.length))
      EDGES <- DF[which.edge(newtree, sing.node),]
      to_index<-EDGES[,4]
      WHERE <- newtree$edge[to_index,2]

      if(insertion=="random"){position<-EDGES[,3]*runif(1, 0, 1)} #bound at a random point of the branch
      if(insertion=="middle"){position<-EDGES[,3]*0.5}  #bound at the middle of the branch
      if(insertion=="long")  {position<-EDGES[,3]}      #bound at the beggining of the branch

      newtree<-bind.tip(newtree, tip.label =  species[i], edge.length=NULL, where=WHERE, position=position)
    }else{
      node<-findMRCA(newtree, tips = c(singleton, species[1:(i-1)])) #the rest of the species are added to the node
      newtree<- polytomy_into_node(tree = newtree, new.tip  =  species[i], node = node)
    }}

  return(newtree)}




add_to_singleton<-function(tree, singleton,new.tips){
  newtree<-tree
  new.tips<-gsub(" ", "_",new.tips)
  if(!exists("tree")){stop("Tree object not found")}
  if(length(newtree$tip.label[newtree$tip.label==singleton])==0){stop("Singleton tip not found")}
  To_add <- new.tips
  added <- singleton

  for(i in 1:length(To_add)){

    if(length(added)==1) {    #similar to poltomy, but random branch length

      new.tip<-To_add[i]
      added_R<-which(newtree$tip.label==added[1])

      DF <- data.frame(newtree$edge,newtree$edge.length,1:length(newtree$edge.length))
      EDGES <- DF[which.edge(newtree, added_R),]
      to_index<-EDGES[,4]

      WHERE <- newtree$edge[to_index,2]
      LENGTH <- newtree$edge.length[to_index]

      MIN<-0
      MAX<-LENGTH
      POS<-runif(1,MIN,MAX)

      while(POS==MIN | POS==MAX){POS<- runif(1,MIN,MAX)}

      newtree<-bind.tip(newtree, new.tip, edge.length=NULL, where=WHERE, position=POS)
      added <- c(added,new.tip)

    } else {

      new.tip<-To_add[i]
      species_R <- vector(mode="numeric",length=length(added))
      for(k in 1:length(added)) { species_R[k] <- which(newtree$tip.label==added[k])}

      DF <- data.frame(newtree$edge,newtree$edge.length,1:length(newtree$edge.length))
      EDGES <- DF[which.edge(newtree, species_R),]
      EDGES2<-DF[DF[,1]==getParent(newtree,min(EDGES[,1])) & DF[,2]==min(EDGES[,1]),]
      EDGES<-rbind(EDGES2,EDGES)
      to_index<-sample(EDGES[,4],1,prob=EDGES[,3])
      WHERE <- newtree$edge[to_index,2]
      LENGTH <- newtree$edge.length[to_index]

      MIN<-0
      MAX<-LENGTH
      POS<-runif(1,MIN,MAX)

      while(POS==MIN | POS==MAX){POS<-runif(1,MIN,MAX)}

      newtree <- bind.tip(newtree, new.tip, edge.length=NULL, where=WHERE, position=POS)
      added <- c(added,new.tip)
    }
  }
  return(newtree)
}

#function to add species at specific branches
stick.to.branch <- function(tree,edges,new.tip=new.tip,prob=TRUE){

  DF <- data.frame(tree$edge,tree$edge.length)

  if(prob == TRUE) {

    for(i in 1:dim(edges)[1]) {

      tips1 = edges[i,1:2]
      tips2 = edges[i,3:4]

      tips1.char<- c(as.character(edges[i,1]), as.character(edges[i,2]))
      tips2.char<- c(as.character(edges[i,3]), as.character(edges[i,4]))

      N1<-findMRCA(tree, tips=tips1.char, "node")

      if(duplicated(tips2.char)[2]==TRUE) {N2<-which(tree$tip.label==tips2.char[1])} else {N2<-findMRCA(tree, tips=tips2.char, "node")}

      if(tree$edge[which.edge(tree, N2),1]==N1) {DF[which.edge(tree, N2),3] -> edges[i,5]} else {

        EDGES <- c(NA,NA)
        Child=N2
        getParent(tree,Child) -> Parent
        EDGES <- rbind(c(Parent,Child))

        while(Parent!=N1) {
          getDescendants(tree, Parent, curr=NULL)->Desc_Par
          Desc_Par<-c(Desc_Par[Desc_Par%in%getDescendants(tree, Child, curr=NULL)==FALSE & Desc_Par!=Child],Parent)
          EDGES<-rbind(EDGES,tree$edge[tree$edge[,1]%in%Desc_Par & tree$edge[,2]!=Child,])
          Child<-Parent
          Parent<-getParent(tree,Parent)
          EDGES <- rbind(EDGES,c(Parent,Child))
        }

        sum(DF[DF[,1]%in%EDGES[,1]&DF[,2]%in%EDGES[,2],3]) -> edges[i,5]

      }
    }

    SS <- sample(1:dim(edges)[1],1,replace=FALSE,edges[,5])

    tips1 = edges[SS,1:2]
    tips2 = edges[SS,3:4]


    tips1.char<- c(as.character(edges[SS,1]), as.character(edges[SS,2]))
    tips2.char<- c(as.character(edges[SS,3]), as.character(edges[SS,4]))

    N1<-findMRCA(tree, tips=tips1.char, "node")

    if(duplicated(tips2.char)[2]==TRUE) {N2<-which(tree$tip.label==tips2.char[1])} else{N2<-findMRCA(tree, tips=tips2.char, "node")}

    EDGES <- c(NA,NA)
    Child=N2
    getParent(tree,Child) -> Parent
    EDGES <- rbind(c(Parent,Child))

    while(Parent!=N1) {

      getDescendants(tree, Parent, curr=NULL)->Desc_Par
      Desc_Par<-c(Desc_Par[Desc_Par%in%getDescendants(tree, Child, curr=NULL)==FALSE & Desc_Par!=Child],Parent)
      EDGES<-rbind(EDGES,tree$edge[tree$edge[,1]%in%Desc_Par & tree$edge[,2]!=Child,])
      Child<-Parent
      Parent<-getParent(tree,Parent)
      EDGES <- rbind(EDGES,c(Parent,Child))
    }

    Lengths <- vector(mode="numeric",length=length(EDGES[,1]))

    for(j in 1:dim(EDGES)[1]){
      which.edge(tree, EDGES[j,2]) -> Lengths[j]
    }

    if(length(Lengths)==1){Lengths->to_index} else {sample(Lengths,1,replace=FALSE,DF[Lengths,3])->to_index}

    WHERE = tree$edge[to_index,2]

    tree$edge.length[to_index]->LLL

    MIN=0
    MAX=LLL
    runif(1,MIN,MAX)->POS

    while(POS==MIN | POS==MAX){
      runif(1,MIN,MAX)->POS
    }

    bind.tip(tree, new.tip, edge.length=NULL, where=WHERE, position=POS)->tree

  } else {

    SS <- sample(1:dim(edges)[1],1)

    tips1 = edges[SS,1:2]
    tips2 = edges[SS,3:4]

    tips1.char<- c(as.character(edges[SS,1]), as.character(edges[SS,2]))
    tips2.char<- c(as.character(edges[SS,3]), as.character(edges[SS,4]))

    N1<-findMRCA(tree, tips=tips1.char, "node")

    if(duplicated(tips2.char)[2]==TRUE) {N2<-which(tree$tip.label==tips2.char[1])} else{N2<-findMRCA(tree, tips=tips2.char, "node")}

    EDGES <- c(NA,NA)
    Child=N2
    getParent(tree,Child) -> Parent
    EDGES <- rbind(c(Parent,Child))

    while(Parent!=N1) {

      getDescendants(tree, Parent, curr=NULL)->Desc_Par
      Desc_Par<-c(Desc_Par[Desc_Par%in%getDescendants(tree, Child, curr=NULL)==FALSE & Desc_Par!=Child],Parent)
      EDGES<-rbind(EDGES,tree$edge[tree$edge[,1]%in%Desc_Par & tree$edge[,2]!=Child,])
      Child<-Parent
      Parent<-getParent(tree,Parent)
      EDGES <- rbind(EDGES,c(Parent,Child))
    }

    Lengths <- vector(mode="numeric",length=length(EDGES[,1]))

    for(j in 1:dim(EDGES)[1]){
      which.edge(tree, EDGES[j,2]) -> Lengths[j]
    }

    if(length(Lengths)==1){Lengths->to_index} else {sample(Lengths,1)->to_index}

    WHERE = tree$edge[to_index,2]

    tree$edge.length[to_index]->LLL

    MIN=0
    MAX=LLL
    runif(1,MIN,MAX)->POS

    while(POS==MIN | POS==MAX){
      runif(1,MIN,MAX)->POS
    }

    bind.tip(tree, new.tip, edge.length=NULL, where=WHERE, position=POS)->tree

  }
  return(tree)
}

# A function to stick species at random within a polyphyletic clade
add_to_polyphyletic<-function(tree, species, polyphyletic.insertion=c("freq", "large", "all")){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  species<-gsub(" ", "_",species)
  species.genus<- word(species, 1, sep="_")[!duplicated(word(species, 1, sep="_"))]
  if(length(species.genus)>1){stop("Species from more than 1 genera are included")}
  new.tree<-tree

  genus<- species.genus
  list<- data.frame(species=tree$tip.label)
  taxa.vector<- as.vector(list[word(list$species, 1, sep="_")==genus,]) #genus species within the tree. This vector is used for all species to be added
  if(length(taxa.vector)==0){stop(paste0("Genus ", genus, " is not included in your tree."))}

  groups<- rep(list(NA), times=length(taxa.vector))
  names(groups)<-taxa.vector
  for( t in 1:length(taxa.vector)){
    taxon<- taxa.vector[t]
    taxon.tip<- which(new.tree$tip.label==taxon) #tip value
    parent<-new.tree$edge[new.tree$edge[,2]==taxon.tip,1] #direct ancestor
    siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])] #ancestor's descendants
    if(length(siblings[word(siblings,1,sep="_")==genus])==1){groups[[t]]<- taxon} #its siblings are from a different genus; it is a singleton inside another clade
    if(length(siblings[word(siblings,1,sep="_")==genus]) >1){ #at least one sibling is from the same genus
      while(length(word(siblings,1,sep="_")[!duplicated(word(siblings,1,sep="_"))])==1){ #tip and parent upstream until they are from different genera
        taxon.tip<-parent
        parent<-new.tree$edge[new.tree$edge[,2]==taxon.tip,1]
        siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])]
      }

      gen.MRCA<- findMRCA(tree = new.tree,tips = siblings[word(siblings, 1, sep="_")==genus]) #MCRA form same genus siblings; probably always the same as "sticking tip"; but may change
      grouped<- new.tree$tip.label[getDescendants(new.tree, gen.MRCA)][!is.na(new.tree$tip.label[getDescendants(new.tree, gen.MRCA)])] #non-node descendants
      grouped.gen<- word(grouped,1,sep="_")[!duplicated(word(grouped,1,sep="_"))] #MRCA descendants
      if(length(grouped.gen)==1){groups[[t]]<- grouped}  #monophyletic subgruoup
      if(length(grouped.gen) >1){ #always paraphyletic subgroup, it can't be mono or poly
        groups[[t]]<- grouped[word(grouped, 1, sep="_")==genus]}}
  }

  group.types<-rep(list(NA), times=length(taxa.vector))
  names(group.types)<-taxa.vector
  for( t in 1:length(taxa.vector)){
    taxon<- taxa.vector[t]
    taxon.tip<- which(new.tree$tip.label==taxon) #tip value
    parent<-new.tree$edge[new.tree$edge[,2]==taxon.tip,1] #direct ancestor
    siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])] #ancestor's descendants
    if(length(siblings[word(siblings,1,sep="_")==genus])==1){group.types[[t]]<- "singleton"} #its siblings are from a different genus; it is a singleton inside another clade
    if(length(siblings[word(siblings,1,sep="_")==genus]) >1){ #at least one sibling is from the same genus
      while(length(word(siblings,1,sep="_")[!duplicated(word(siblings,1,sep="_"))])==1){ #tip and parent upstream until they are from different genera
        taxon.tip<-parent
        parent<-new.tree$edge[new.tree$edge[,2]==taxon.tip,1]
        siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])]
      }

      gen.MRCA<- findMRCA(tree = new.tree,tips = siblings[word(siblings, 1, sep="_")==genus]) #MCRA form same genus siblings; probably always the same as "sticking tip"; but may change
      grouped<- new.tree$tip.label[getDescendants(new.tree, gen.MRCA)][!is.na(new.tree$tip.label[getDescendants(new.tree, gen.MRCA)])] #non-node descendants
      grouped.gen<- word(grouped,1,sep="_")[!duplicated(word(grouped,1,sep="_"))] #MRCA descendants
      if(length(grouped.gen)==1){group.types[[t]]<- "monophyletic"}  #monophyletic subgruoup
      if(length(grouped.gen) >1){ #always paraphyletic subgroup, it can't be mono or poly
        group.types[[t]]<- "paraphyletic"} }}

  if(polyphyletic.insertion=="freq"){
  for(p in 1:length(species)){
    sticked.species<- species[p]
    sticking.species<- sample(taxa.vector, 1) #a species from the vector is chosen (this allows a probability based on frequencies inside each group)
    sticking.tip<- which(new.tree$tip.label==sticking.species) #tip value

    parent<-new.tree$edge[new.tree$edge[,2]==sticking.tip,1] #direct ancestor
    siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])] #ancestor's descendants

    if(length(siblings[word(siblings,1,sep="_")==genus])==1){new.tree<- add_to_singleton(new.tree,sticking.species,new.tips=sticked.species)} #its siblings are from a different genus; it is a singleton inside another clade
    if(length(siblings[word(siblings,1,sep="_")==genus]) >1){ #at least one sibling is from the same genus
      while(length(word(siblings,1,sep="_")[!duplicated(word(siblings,1,sep="_"))])==1){ #tip and parent upstream until they are from different genera
        sticking.tip<-parent
        parent<-new.tree$edge[new.tree$edge[,2]==sticking.tip,1]
        siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])]
      }

      gen.MRCA<- findMRCA(tree = new.tree,tips = siblings[word(siblings, 1, sep="_")==genus]) #MCRA form same genus siblings; probably always the same as "sticking tip"; but may change
      grouped<- new.tree$tip.label[getDescendants(new.tree, gen.MRCA)][!is.na(new.tree$tip.label[getDescendants(new.tree, gen.MRCA)])] #non-node descendants
      grouped.gen<- word(grouped,1,sep="_")[!duplicated(word(grouped,1,sep="_"))] #MRCA descendants
      if(length(grouped.gen)==1){new.tree<- add_into_node(new.tree, sticking.tip, sticked.species)} #monophyletic subgruoup
      if(length(grouped.gen) >1){ #always paraphyletic subgroup, it can't be mono or poly
        intruders<-grouped[word(grouped, 1, sep="_")!=genus]
        if(length(intruders)==1){new.tree<- add_into_node(new.tree, sticking.tip, sticked.species)}else{ # singleton intruders; added as monophyletic
          intruders.MRCA<- findMRCA(tree = new.tree,tips = intruders)
          new.tree<- add_into_paraphyletic_node(new.tree, new.tip=sticked.species, group.node = gen.MRCA, intern.node = intruders.MRCA) #monophyletic intruders; added as paraphyletic


        }}}}
  }
  if(polyphyletic.insertion=="large"){
    m<-max(as.numeric(lengths(groups)))
    group.types<- group.types[lengths(groups)==m]
    groups<- groups[lengths(groups)==m]

    slot<-sample(1:length(groups), size = 1)
    group<- groups[[slot]]
    group.type <- group.types[[slot]]

    for(s in 1:length(species)){
      sp<-species[s]
      if(group.type=="monophyletic"){
        MRCA<- findMRCA(new.tree, tips=group)
        new.tree<-add_into_node(tree = new.tree, node = MRCA, new.tip = sp)}
      if(group.type=="paraphyletic"){
        MRCA<- findMRCA(new.tree, tips=group)
        intruders<-new.tree$tip.label[getDescendants(new.tree, MRCA)][!is.na(new.tree$tip.label[getDescendants(new.tree, MRCA)])]
        intruders<- intruders[word(intruders, 1, sep="_")!=genus]
        intruders.MRCA<- findMRCA(new.tree, tips=intruders)

        new.tree<-add_into_paraphyletic_node(tree=new.tree, new.tip = sp, group.node = MRCA, intern.node = intruders.MRCA)}
    }

  }
  if(polyphyletic.insertion=="all"){
    for(s in 1:length(species)){
      sp<- species[s]
      MRCA<- findMRCA(new.tree, taxa.vector)
      new.tree<- add_into_node(new.tree, node = MRCA, new.tip = sp)
    }
  }
  return(new.tree)

}





#randtip "MOTHER FUNCTION" ####


RANDTIP<- function(tree, species.table, type=c("random", "genus.polytomy", "family.polytomy", "order.polytomy", "class.polytomy"),
                   aggregate.subspecies=TRUE, insertion=c("random", "middle","long"), prob=TRUE){

  if(type=="random"){
  species.table$using.taxa<-NA #New column with the name for the first round
  for(i in 1:nrow(species.table)){
    if(is.na(species.table$aggregate.subspecies[i])){    #para los NA (en la tabla no hay especificadas excepciones) se aplica lo marado en la funcion
      if(aggregate.subspecies==TRUE){species.table$using.taxa[i]<-paste0(word(species.table$taxon[i],1, sep="_"), "_",word(species.table$taxon[i],2, sep="_"))}
      else{species.table$using.taxa[i]<-species.table$taxon[i]}
      next}

    if(species.table$aggregate.subspecies[i]==0){species.table$using.taxa[i]<-species.table$taxon[i] #para los especificados como 0 no se agrupa
    next}
    if(species.table$aggregate.subspecies[i]==1){species.table$using.taxa[i]<-paste0(word(species.table$taxon[i],1, sep="_"), "_",word(species.table$taxon[i],2, sep="_"))}} # the ones specified with a 1 are clustered

  species.table.dupl <- species.table[ duplicated(species.table$using.taxa),]
  species.table <-      species.table[!duplicated(species.table$using.taxa),]



  newtree<- tree
  taxa <- species.table$using.taxa[!duplicated(species.table$using.taxa)]
  taxa<- gsub(" ", "_", taxa)
  taxa <- taxa[taxa%!in%tree$tip.label]

  taxa.genera<- word(taxa, 1, sep="_")[!duplicated(word(taxa, 1, sep="_"))]
  taxa.genera<-sample(taxa.genera, length(taxa.genera), replace = F)
  }else{

    species.table$using.taxa <-species.table$taxon  #In polytomy cases, names are not changed


    newtree<- tree
    taxa <- species.table$using.taxa[!duplicated(species.table$using.taxa)]
    taxa<- gsub(" ", "_", taxa)
    taxa <- taxa[taxa%!in%tree$tip.label]

    taxa.genera<- word(taxa, 1, sep="_")[!duplicated(word(taxa, 1, sep="_"))]
    taxa.genera<-sample(taxa.genera, length(taxa.genera), replace = F)
  } #using name preparation

  if(type=="genus.polytomy"){
    genera<- species.table$genus[!duplicated(species.table$genus)]

    for(p in 1:length(genera)){
      genus<- genera[p]
      genus.taxa<- species.table$taxon[species.table$genus==genus]
      genus.taxa<- genus.taxa[genus.taxa%!in%tree$tip.label]
      genus.genera<- species.table$genus[species.table$genus==genus] #this is redundant, but keeps the structure
      union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%genus.genera]   #species (tips) within class IN ORIGINAL TREE
      if(length(genus.taxa)==0){next}
      if(length(union.tips)==0){
        family<- species.table$family[species.table$genus==genus][!duplicated(species.table$family[species.table$genus==genus])]
        family.taxa<- species.table$taxon[species.table$family==family]
        family.taxa<- family.taxa[family.taxa%!in%tree$tip.label]
        family.genera<- species.table$genus[species.table$family==family]
        union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%family.genera]   #species (tips) within class IN ORIGINAL TREE
        if(length(family.taxa)==0){next}
        if(length(union.tips)==0){
          order<-species.table$order[species.table$family==family][!duplicated(species.table$order[species.table$family==family])]
          order.taxa<- species.table$taxon[species.table$order==order]
          order.taxa<- order.taxa[order.taxa%!in%tree$tip.label]
          order.genera<- species.table$genus[species.table$order==order]
          union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%order.genera]   #species (tips) within class IN ORIGINAL TREE
          if(length(order.taxa)==0){next}
          if(length(union.tips)==0){
            class<- species.table$class[species.table$order==order][!duplicated(species.table$class[species.table$order==order])]
            class.taxa<- species.table$taxon[species.table$class==class]
            class.taxa<- class.taxa[class.taxa%!in%tree$tip.label]
            class.genera<- species.table$genus[species.table$class==class]
            union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%class.genera]   #species (tips) within class IN ORIGINAL TREE
            if(length(class.taxa)==0){next}
            if(length(union.tips)==0){message(paste0("ATTENTION: genus ", class.genera, " was not included as no Class coincidences were found"))####JOIN TO UPPER TAXONOMIC LEVEL
              next}
            newtree<- polytomy_over_node(tree = newtree, species = class.taxa, node=findMRCA(newtree, tips=union.tips), insertion = insertion)####JOIN TO UPPER TAXONOMIC LEVEL
            next}
          if(length(union.tips)==1){node<- which(newtree$tip.label==union.tips)}else{node<-findMRCA(newtree, tips=union.tips)}
          newtree<- polytomy_over_node(tree = newtree, species = order.taxa, node=node, insertion = insertion)####JOIN TO UPPER TAXONOMIC LEVEL
          next}

        if(length(union.tips)==1){node<- which(newtree$tip.label==union.tips)}else{node<-findMRCA(newtree, tips=union.tips)}
        newtree<- polytomy_over_node(tree = newtree, species = family.taxa, node=node, insertion = insertion)####JOIN TO UPPER TAXONOMIC LEVEL
        next}

      if(length(union.tips)==1){node<- which(newtree$tip.label==union.tips)}else{node<-findMRCA(newtree, tips=union.tips)}
      newtree<- polytomy_over_node(tree = newtree, species = genus.taxa, node=node, insertion = insertion)
      }}

  if(type=="family.polytomy"){
    families<- species.table$family[!duplicated(species.table$family)]

    for(p in 1:length(families)){
      family<- families[p]
      family.taxa<- species.table$taxon[species.table$family==family]
      family.taxa<- family.taxa[family.taxa%!in%tree$tip.label]
      family.genera<- species.table$genus[species.table$family==family]
      union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%family.genera]   #species (tips) within class IN ORIGINAL TREE
      if(length(family.taxa)==0){next}
      if(length(union.tips)==0){
        order<-species.table$order[species.table$family==family][!duplicated(species.table$order[species.table$family==family])]
        order.taxa<- species.table$taxon[species.table$order==order]
        order.taxa<- order.taxa[order.taxa%!in%tree$tip.label]
        order.genera<- species.table$genus[species.table$order==order]
        union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%order.genera]   #species (tips) within class IN ORIGINAL TREE
        if(length(order.taxa)==0){next}
        if(length(union.tips)==0){
          class<- species.table$class[species.table$order==order][!duplicated(species.table$class[species.table$order==order])]
          class.taxa<- species.table$taxon[species.table$class==class]
          class.taxa<- class.taxa[class.taxa%!in%tree$tip.label]
          class.genera<- species.table$genus[species.table$class==class]
          union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%class.genera]   #species (tips) within class IN ORIGINAL TREE
          if(length(class.taxa)==0){next}
          if(length(union.tips)==0){message(paste0("ATTENTION: genus ", class.genera, " was not included as no Class coincidences were found"))####JOIN TO UPPER TAXONOMIC LEVEL
            next}
          newtree<- polytomy_over_node(tree = newtree, species = class.taxa, node=findMRCA(newtree, tips=union.tips), insertion = insertion)####JOIN TO UPPER TAXONOMIC LEVEL
          next}
        if(length(union.tips)==1){node<- which(newtree$tip.label==union.tips)}else{node<-findMRCA(newtree, tips=union.tips)}
        newtree<- polytomy_over_node(tree = newtree, species = order.taxa, node=node, insertion = insertion)####JOIN TO UPPER TAXONOMIC LEVEL
        next}

      if(length(union.tips)==1){node<- which(newtree$tip.label==union.tips)}else{node<-findMRCA(newtree, tips=union.tips)}
      newtree<- polytomy_over_node(tree = newtree, species = family.taxa, node=node, insertion = insertion)}}

  if(type=="order.polytomy"){
    orders<- species.table$order[!duplicated(species.table$order)]

    for(p in 1:length(orders)){
      order<-orders[p]
      order.taxa<- species.table$taxon[species.table$order==order]
      order.taxa<- order.taxa[order.taxa%!in%tree$tip.label]
      order.genera<- species.table$genus[species.table$order==order]
      union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%order.genera]   #species (tips) within class IN ORIGINAL TREE
      if(length(order.taxa)==0){next}
      if(length(union.tips)==0){
        class<- species.table$class[species.table$order==order][!duplicated(species.table$class[species.table$order==order])]
        class.taxa<- species.table$taxon[species.table$class==class]
        class.taxa<- class.taxa[class.taxa%!in%tree$tip.label]
        class.genera<- species.table$genus[species.table$class==class]
        union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%class.genera]   #species (tips) within class IN ORIGINAL TREE
        if(length(class.taxa)==0){next}
        if(length(union.tips)==0){message(paste0("ATTENTION: genus ", class.genera, " was not included as no Class coincidences were found"))####JOIN TO UPPER TAXONOMIC LEVEL
          next}
        newtree<- polytomy_over_node(tree = newtree, species = class.taxa, node=findMRCA(newtree, tips=union.tips), insertion = insertion)####JOIN TO UPPER TAXONOMIC LEVEL
        next}
      if(length(union.tips)==1){node<- which(newtree$tip.label==union.tips)}else{node<-findMRCA(newtree, tips=union.tips)}
      newtree<- polytomy_over_node(tree = newtree, species = order.taxa, node=node, insertion = insertion)
      }}

  if(type=="class.polytomy"){
    classes<- species.table$class[!duplicated(species.table$class)]

 for(p in 1:length(classes)){
      class<- classes[p]
      class.taxa<- species.table$taxon[species.table$class==class]
      class.taxa<- class.taxa[class.taxa%!in%tree$tip.label]
      class.genera<- species.table$genus[species.table$class==class]
      union.tips<-tree$tip.label[word(tree$tip.label, 1, sep="_")%in%class.genera]   #species (tips) within class IN ORIGINAL TREE
      if(length(class.taxa)==0){next}
      if(length(union.tips)==0){message(paste0("ATTENTION: genus ", class.genera, " was not included as no Class coincidences were found"))####JOIN TO UPPER TAXONOMIC LEVEL
        next}

       if(length(union.tips)==1){node<- which(newtree$tip.label==union.tips)}else{node<-findMRCA(newtree, tips=union.tips)}
       newtree<- polytomy_over_node(tree = newtree, species = class.taxa, node=node, insertion = insertion)}}



  if(type=="random"){
    for(i in 1: length(taxa.genera)){        #loop 1
      genus<- taxa.genera[i]
      genus.type<-species.table$genus.type[species.table$genus==genus][!duplicated(species.table$genus.type[species.table$genus==genus])]

      genus.taxa <- taxa[word(taxa, 1, sep="_")==genus] #genus taxa selection
      genus.taxa <- sample(genus.taxa, length(genus.taxa), replace = F) #order randomized

      grouped.taxa<-species.table$using.taxa[!is.na(species.table$relative.species)]#tips with no "relatives" information selected (they will be bound as monophyletic)

      grouped.taxa<-grouped.taxa[grouped.taxa%in%genus.taxa]
      genus.taxa<-genus.taxa[genus.taxa%!in%grouped.taxa]




      if(length(grouped.taxa)>0){    #if "grouped.taxa" exist
        for( j in 1:length(grouped.taxa)){
          grouping.taxa<- species.table$relative.species[species.table$using.taxa==grouped.taxa[j]]
          grouping.taxa<-gsub(" ","",grouping.taxa)
          grouping.taxa<- strsplit(grouping.taxa, split = ",")[[1]]
          grouping.taxa<-grouping.taxa[grouping.taxa%in%tree$tip.label] #taxa to be grouped hich are in the tree are selected
          if(length(grouping.taxa)==1){
            newtree<- add_to_singleton(newtree, singleton = grouping.taxa, new.tips = grouped.taxa[j])  #first is added as singleton
          } else{
            node<-findMRCA(newtree, tips = grouping.taxa)
            newtree<- add_into_node(newtree, new.tip = grouped.taxa[j],node = node)} #the rest ar added as monophyletic
          rm(j, node, grouping.taxa)}
        next}

      #forbidden.groups
      forb.genera<- species.table[species.table$genus.type=="MONOPHYLETIC"|species.table$genus.type=="PARAPHYLETIC","genus"]
      forb.genera<-forb.genera[!duplicated(forb.genera)]
      forbidden.groups<- rep(list(NA),length(forb.genera))
      for(l in 1:length(forb.genera)){
        forbidden.groups[[l]]<- species.table[species.table$genus==forb.genera[l], "taxon"]
      }


      #Hereon we will work with genus.taxa; i.e., no grouped tips.
      if(genus.type=="MONOPHYLETIC"){
        for( j in 1:length(genus.taxa)){
          newtree<- add_to_monophyletic(newtree, new.tip = genus.taxa[j])#for the given genus, being MP, taxa are added one by one as MP
          rm(j)}}

      if(genus.type=="PARAPHYLETIC"){
        for( j in 1:length(genus.taxa)){
          newtree<- add_to_paraphyletic(newtree, new.tip = genus.taxa[j])#for the given genus, being PaP, taxa are added one by one as PaP
          rm(j)}}

      if(genus.type=="POLYPHYLETIC"){
        for( j in 1:length(genus.taxa)){
          newtree<- add_to_polyphyletic(newtree, new.tip = genus.taxa[j])#for the given genus, being PoP, taxa are added one by one as PoP
          rm(j)}}

      if(genus.type=="SINGLETON GENUS"){   #for the given genus, being singleton, all tips are added in one step
        newtree<- add_to_singleton(newtree, singleton = tree$tip.label[word(tree$tip.label,1,sep="_")== genus], new.tips = genus.taxa)
      }

      if(genus.type=="NOT INCLUDED"){
        genus.synonyms<- species.table$synonim.genus[species.table$genus==genus&!is.na(species.table$synonim.genus)]#looking for synonim information...
        genus.siblings<- species.table$sibling.genus[species.table$genus==genus&!is.na(species.table$sibling.genus)]#... or sibling information

        genus.synonyms.tips<-species.table$using.taxa[species.table$genus==genus&!is.na(species.table$synonim.genus)]
        genus.siblings.tips<-species.table$using.taxa[species.table$genus==genus&!is.na(species.table$sibling.genus)]

        genus.taxa.NO.except<-genus.taxa[genus.taxa%!in%c(genus.synonyms.tips,genus.siblings.tips)]#taxa without exceptions

        if(length(genus.synonyms)>0){#for genera with synonym exceptions
          for(i in 1:length(genus.synonyms.tips)){
            genus.synonyms.tips[i]<-paste0(genus.synonyms[i], "_",genus.synonyms.tips[i]) #synonym genus pasted to the begginig of the name; must be in the tree
            new.genus<- word(genus.synonyms.tips[i], 1, sep="_")
            new.genus.type<-species.table$genus.type[species.table$genus==new.genus][!duplicated(species.table$genus.type[species.table$genus==new.genus])]

            if(new.genus.type=="MONOPHYLETIC"){
              newtree<- add_to_monophyletic(newtree, new.tip = genus.synonyms.tips[i])}
            if(new.genus.type=="PARAPHYLETIC"){
              newtree<- add_to_paraphyletic(newtree, new.tip = genus.synonyms.tips[i])}
            if(new.genus.type=="POLYPHYLETIC"){
              newtree<- add_to_polyphyletic(newtree, species = genus.synonyms.tips[i])}
            if(new.genus.type=="SINGLETON GENUS"){
              newtree<- add_to_singleton(newtree,
                                         singleton = newtree$tip.label[word(newtree$tip.label,1,sep="_")== new.genus],
                                         new.tips = genus.synonyms.tips[i])}

            newtree$tip.label<-gsub(pattern = genus.synonyms.tips[i],
                                    replacement = species.table$using.taxa[species.table$genus==genus&!is.na(species.table$synonim.genus)][i],
                                    x =newtree$tip.label )
          }}
        if(length(genus.siblings)>0){#for genera with sibling exceptions
          if(length(genus.siblings[!duplicated(genus.siblings)])>1){stop(paste0("Multiple genus siblings for genus ",genus))}
          new.genus<- genus.siblings[!duplicated(genus.siblings)]

          genus.siblings.tips<-sample(genus.siblings.tips, length(genus.siblings.tips), replace = F)

          adding.node<- findMRCA(newtree, tips = newtree$tip.label[word(newtree$tip.label, 1, sep="_")==new.genus])
          newtree<-add_over_node(newtree, new.tip = genus.siblings.tips[1], node=adding.node)#it is added over the sibling genus MCRA

          newtree<-add_to_singleton(newtree, singleton = genus.siblings.tips[1],new.tips =  genus.siblings.tips[2:length(genus.siblings.tips)])#te rest are added to singleton
        }
        if(length(genus.taxa.NO.except)>0){#for those with no exceptions, they will be added randomly within family, order or class

          fam<-species.table$family[species.table$using.taxa%in%genus.taxa.NO.except][!duplicated(species.table$family[species.table$using.taxa%in%genus.taxa.NO.except])]
          genera<- species.table$genus[species.table$family==fam]
          if(length(genera)>0){
            family.taxa<- newtree$tip.label[word(newtree$tip.label,1, sep="_")%in%genera]
            newtree<- add_into_node_exceptions(tree=newtree, node=findMRCA(tree = newtree,tips = family.taxa),new.tip = genus.taxa.NO.except[1], exception.list = forbidden.groups)
            if(length(genus.taxa.NO.except)>1){newtree<- add_to_singleton(tree=newtree, singleton = genus.taxa.NO.except[1],new.tips = genus.taxa.NO.except[2:length(genus.taxa.NO.except)])}

          }else{
            ord<-species.table$order[species.table$using.taxa%in%genus.taxa.NO.except][!duplicated(species.table$family[species.table$using.taxa%in%genus.taxa.NO.except])]
            genera<- species.table$genus[species.table$order==ord]
            if(length(genera)>0){
              order.taxa<- newtree$tip.label[word(newtree$tip.label,1, sep="_")%in%genera]
              newtree<- add_into_node_exceptions(tree=newtree, node=findMRCA(tree = newtree,tips = order.taxa),new.tip = genus.taxa.NO.except[1], exception.list = forbidden.groups)
              if(length(genus.taxa.NO.except)>1){newtree<- add_to_singleton(tree=newtree, singleton = genus.taxa.NO.except[1],new.tips = genus.taxa.NO.except[2:length(genus.taxa.NO.except)])}
            }else{
              cla<-species.table$class[species.table$using.taxa%in%genus.taxa.NO.except][!duplicated(species.table$family[species.table$using.taxa%in%genus.taxa.NO.except])]
              genera<- species.table$genus[species.table$class==cla]
              if(length(genera)>0){
                class.taxa<- newtree$tip.label[word(newtree$tip.label,1, sep="_")%in%genera]
                newtree<- add_into_node_exceptions(tree=newtree, node=findMRCA(tree = newtree,tips = class.taxa),new.tip = genus.taxa.NO.except[1], exception.list = forbidden.groups)
                if(length(genus.taxa.NO.except)>1){newtree<- add_to_singleton(tree=newtree, singleton = genus.taxa.NO.except[1],new.tips = genus.taxa.NO.except[2:length(genus.taxa.NO.except)])}
              }else{
                stop(paste0("Genus ", genus, " has no family, order or class relatives"))}}}
        }
      }

      print(paste0(i, " out of ", length(taxa.genera), ". Genus ", genus, " (", genus.type,"): ", round(i/length(taxa.genera)*100, 2), " %"))
    }

    # phase 2. subspecies
    if(nrow(species.table.dupl)>0){ #species grouped
      rep.taxa <-species.table.dupl$taxon #species are found
      rep.taxa.species<- species.table.dupl$using.taxa[!duplicated(species.table.dupl$using.taxa)]#subspecies' species

      for (r in 1:length(rep.taxa.species)){
        ssps <- rep.taxa[paste0(word(rep.taxa, 1, sep="_"),"_",word(rep.taxa, 2, sep="_"))==rep.taxa.species[r]] #subspecies within species are selected
        ssps <- sample(ssps, length(ssps), replace = F) #order randomized
        sp.to.add<- newtree$tip.label[paste0(word(newtree$tip.label,1, sep="_"), "_", word(newtree$tip.label,2, sep="_"))== rep.taxa.species[r] ]

        newtree<- add_to_singleton(newtree, singleton = sp.to.add, new.tips = ssps) #subspecies are added to their sister as singleton

      }

    }

    for(n in 1:nrow(species.table)){  #original names are returned
      if(species.table$using.taxa[n]== species.table$taxon[n]){next}
      if(species.table$using.taxa[n]%!in% newtree$tip.label){next}
      if(species.table$using.taxa[n]%in% newtree$tip.label){
        newtree$tip.label[newtree$tip.label==species.table$using.taxa[n]]<-species.table$taxon[n]
      }}


  }
  return(newtree)
}

plot(RANDTIP(tree=Tree, species.table=tabla.info, aggregate.subspecies=TRUE, type = "genus.polytomy", insertion = "middle")) #, insertion="long"
plot(RANDTIP(tree=Tree, species.table=tabla.info, aggregate.subspecies=TRUE, type = "genus.polytomy", insertion = "long")) #, insertion="long"
plot(RANDTIP(tree=Tree, species.table=tabla.info, aggregate.subspecies=TRUE, type = "genus.polytomy")) #, insertion="long"
plot(RANDTIP(tree=Tree, species.table=tabla.info, aggregate.subspecies=TRUE, type = "random"))

png(paste0("C:/Users/Alumno/Desktop/mierdas/Arbol_prueba3.png"), height = 1500, width = 2500, units = "px")
plot(prueba)
dev.off()

sub.species<- checklist[checklist$match.name%!in%GBOTB.extended$tip.label,]
sub.species<- sub.species[word(sub.species$match.name,1, sep="_")%in%word(Tree$tip.label,1, sep="_"),]

{
  png(paste0("C:/Users/Alumno/Desktop/mierdas/Arbol_prueba2.png"), height = 1500, width = 2500, units = "px")
  cut<- Tree$tip.label[word(Tree$tip.label, 1, sep="_")%in%word(sub.species$match.name[50:75], 1, sep="_")]
  par(mfrow=c(1,2))
  plot( drop.tip(Tree, Tree$tip.label[-match(cut, Tree$tip.label)]), label.offset = 2)
  #nodelabels()
  #tiplabels(adj = c(-10, 0.5))
  plot(RANDTIP(tree = Tree, species = sub.species[50:75, "match.name"],genera.phyleticity=tipos_generos,aggregate.subspecies = TRUE))
  #nodelabels()
  #tiplabels(adj = c(-10, 0.5))
  dev.off()}
