#install.packages("taxize")
library(taxize)
library(openxlsx)
library(ape)
library(stringr)
library(phytools)

rm(list = setdiff(ls(), lsf.str()))#eliminar todo menos funciones
#comentario de prueba

#checklist<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/IberoBalearic_Checklist_R1.xlsx")
#checklist$match.name<- gsub(" ","_", checklist$match.name)
#tipos_generos<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/Tipos_generos_Filet.xlsx")
#genus.corrections<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/Corrections.xlsx")
#taxonomy<- read.xlsx("E:/UNI/4. DOCTORADO/4. Regionalization/LISTADOS/Taxonomic_categories.xlsx")

Tree<- read.tree("E:/UNI/4. DOCTORADO/3. Randtip/25tree.tre")
tabla.info<- read.xlsx("E:/UNI/4. DOCTORADO/3. Randtip/phylo.table.xlsx")

plot(Tree, label.offset = 2)
nodelabels()
tiplabels(adj = c(-10, 0.5))

#CREACIÓN DE FUNCIONES####


'%!in%' <- function(x,y)!('%in%'(x,y))  #función opuesta de %in%


#la función phylo.taxonomy da una tabla de familia, orden y clase de cada género dado en el vector de especies.
phylo.taxonomy<- function(species, db="ncbi"){
  message(paste0("Taxonomic information is being looked up in ",db," Taxonomy Database. This may take a while!"))
  taxa.genera<- word(species, 1, sep="_")[!duplicated(word(species, 1, sep="_"))]     #se extraen los géneros
  taxonomy.df<- data.frame("Genus"=taxa.genera, "Family"=NA, "Order"=NA, "Class"=NA)  #se crea la tabla
  for(i in 1:length(taxa.genera)){    #por cada uno de los géneros...
    tryCatch({
      search<-classification(as.character(taxa.genera[i]), db = db)[[1]]  #se buscan las categorías taxonómicas y se guardan en la tabla
      fam<-search[search$rank=="family","name"]
      ord<-search[search$rank=="order" ,"name"]
      cla<-search[search$rank=="class" ,"name"]
      taxonomy.df$Family[taxonomy.df$Genus==taxa.genera[i]]<-fam
      taxonomy.df$Order[taxonomy.df$Genus==taxa.genera[i]] <-ord
      taxonomy.df$Class[taxonomy.df$Genus==taxa.genera[i]] <-cla
    }, error=function(e){                                               #si hay un error se guarda como NA para evitar parones
      taxonomy.df$Family[taxonomy.df$Genus==taxa.genera[i]]<-NA
      taxonomy.df$Order[taxonomy.df$Genus==taxa.genera[i]] <-NA
      taxonomy.df$Class[taxonomy.df$Genus==taxa.genera[i]] <-NA
    })
    Sys.sleep(0.33)         #taxize solo permite 3 búsquedas por segundo para no colapsar; esto evita que se pare
  }
  return(taxonomy.df)
}
# phylo.taxonomy(species=c("Abies_pinsapo", "Pinus_pinaster", "Achyllea_milleifolium"))


#la función phyleticity busca qué tipo de género es el dado, dentro de un árbol dado.
phyleticity<- function(tree, genus){
  if(length(genus)!=1){stop("Only one genus accepted")}

  list<- data.frame(species=tree$tip.label)                                 #lista las especies del árbol
  taxa.vector<- as.vector(list[word(list$species, 1, sep="_")==genus,])     #lista los géneros  del árbol
  if(length(taxa.vector)==0){                                               #si no aparece el género en el árbol
    type<-"NOT INCLUDED"
    #message(paste0("Genus ", genus, " is ", type, " in your tree"))
  }else{
    if(length(taxa.vector)==1){                                             #si aparece el género en el árbol una sola vez
      type<-"SINGLETON GENUS"
      #message(paste0("Genus ", genus, " is a ", type))
    }else{                                                                 #si aparece el género en el árbol más de una vez...

      MRCA<-findMRCA(tree = tree,tips = taxa.vector)                       #... se busca el MRCA...
      Desc.Tips<-tree$tip.label[getDescendants(tree, MRCA)][!is.na(tree$tip.label[getDescendants(tree, MRCA)])] #...se buscan sus descendientes...
      Desc.Genera<-word(Desc.Tips, 1, sep="_")[!duplicated(word(Desc.Tips, 1, sep="_"))] # si todos los descendientes son del mismo género: monof.
      if(length(Desc.Genera)==1){
        type<-"MONOPHYLETIC"
        #message(paste0("Genus ", genus, " is ", type))
      }else{
        Intruder.Tips<-as.vector(Desc.Tips[word(Desc.Tips, 1, sep = "_")!= genus]) #si hay algún "intruso" en el grupo...
        Intruder.MRCA<-findMRCA(tree = tree,tips = Intruder.Tips)
        if(length(Intruder.Tips)==1){type<-"PARAPHYLETIC"                          #(si solo hay un intruso: parafilético por singleton)
        #message(paste0("Genus ", genus, " is ", type))
        }else{                                                                     # si hay más de un intruso...
          Desc.Intruder.Tips<-tree$tip.label[getDescendants(tree, Intruder.MRCA)][!is.na(tree$tip.label[getDescendants(tree, Intruder.MRCA)])] #MRCA de los intrusos y descendientes

          if(length(Desc.Intruder.Tips)==length(Intruder.Tips)&all(sort(Desc.Intruder.Tips)==sort(Intruder.Tips))){        #si estan agrupados: parafilético por intruso monofilético
            type<-"PARAPHYLETIC"
            # message(paste0("Genus ", genus, " is ", type))
          } else{                                                                 #si nada de lo anterior ha resultado, solo queda polifilético
            type<-"POLYPHYLETIC"
            #message(paste0("Genus ", genus, " is ",type))
          }}}}}
  return(type)
}
#phyleticity(Tree, "Inventum")


#la función phylo.table devuelve una tabla que nutrirá a la función randtip.
phylo.table<- function(species, tree=NULL, taxonomy=NULL, genera.phyleticity=NULL){

  phylo.df<- data.frame(matrix(nrow = length(species),ncol=10))
  names(phylo.df)<-c("taxon", "genus", "genus.type", "family", "order", "class", "aggregate.subspecies",
                     "relative.species","synonim.genus","sibling.genus")

  if(length(species[duplicated(species)])){stop("There are duplicated species")}
  phylo.df$taxon<-species                                     #se listan las especies en la columna de especies
  phylo.df$genus<-word(phylo.df$taxon, 1, sep="_")            #se obtienen los géneros

  taxa.genera<- phylo.df$genus[!duplicated(phylo.df$genus)]   #se listan los géneros

  if(is.null(taxonomy)){  taxonomy<-phylo.taxonomy(species=species) } #si no se aporta una tabla con los datos de taxonomía, se obtienen

  for(i in 1:length(taxa.genera)){
    genus<-as.character(taxa.genera[i])
    fam<-taxonomy$Family[taxonomy$Genus==genus ]
    ord<-taxonomy$Order [taxonomy$Genus==genus ]
    cla<-taxonomy$Class [taxonomy$Genus==genus ]
    phylo.df$family[phylo.df$genus==taxa.genera[i]]<-fam
    phylo.df$order[phylo.df$genus==taxa.genera[i]] <-ord
    phylo.df$class[phylo.df$genus==taxa.genera[i]] <-cla

  }    #esto es por si le das una tabla (p. ej. obtenida con la función phylo.taxonomy)



  if(is.null(genera.phyleticity)){                 #si no se aporta una tabla con los datos de "fileticidad", se obtiene...
    if(is.null(tree)){stop("Tree object required")}
    genera.phyleticity<- data.frame("Genus"=taxa.genera, "Type"=NA)
    for(f in 1:nrow(genera.phyleticity)){
      genera.phyleticity$Type[f]<- phyleticity(tree,  genera.phyleticity$Genus[f])
    }}

  for(i in 1:length(taxa.genera)){         #obtenida antes, o dada desde el principio, se imputan los datos a nuestra tabla.
    genus<-as.character(taxa.genera[i])
    genus.type<-genera.phyleticity$Type[genera.phyleticity$Genus==genus]
    phylo.df$genus.type[phylo.df$genus==genus]<-genus.type
  }



  return(phylo.df)
}
#df<-phylo.table(species = checklist$match.name, Tree,  genera.phyleticity = tipos_generos, taxonomy = taxonomy)

#esta función sirve para ver si un "label" es un nodo interno
is.node<-function(tree, node){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(node %!in% tree$edge){stop("Node number is not in your tree")}
  if (length(getDescendants(tree=tree, node = node))>1){return(TRUE)}else{return(FALSE)}
}

#esta función sirve para ver si un "label" es un tip externo
is.tip <-function(tree, node){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(node %!in% tree$edge){stop("Node number is not in your tree")}
  if (length(getDescendants(tree=tree, node = node))==1){return(TRUE)}else{return(FALSE)}
}

# Function to add tips at random from a given node
add_into_node <- function(tree,node,new.tip) {
  new.tip<-gsub(" ", "_",new.tip) #modificación por si se ha puesto los nombres co espacios...
  tt<-splitTree(tree,split=list(node=node, bp=tree$edge.length[which(tree$edge[,2]==node)])) #se parte el árbol por el nodo especificado
  if (is.ultrametric(tt[[2]])==TRUE){
    tt[[2]]<-add.random(tt[[2]],tips=new.tip)   #se añade aleatoriamente dentro del arbol cortado
    new.tree<-paste.tree(tt[[1]],tt[[2]])       #se pegan los árboles
    return(new.tree)
  } else
  {
    warning("Found issue in sticking, forcing ultrametric") #warning por si no es ultramétrico
    tt[[2]]<-force.ultrametric(tt[[2]], method="extend")
    tt[[2]]<-add.random(tt[[2]],tips=new.tip)
    new.tree<-paste.tree(tt[[1]],tt[[2]])
    return(new.tree)
  }}


# Function to add tips at random monophyletic clade
add_to_monophyletic <- function(tree,new.tip) {
  new.tip<-gsub(" ", "_",new.tip)       #errores de nombres
  genus<- word(new.tip, 1, sep="_")     #se otiene el género (previamente tiene que haber sido identificado como MonoF)
  list<- data.frame(species=tree$tip.label)  #listado de especies dentro del árbol
  taxa.vector<- as.vector(list[word(list$species, 1, sep="_")==genus,]) #especies del árbol del género
  if(length(taxa.vector)==0){stop(paste0("Genus ", genus, " is not included in yout tree."))}
  MRCA<-findMRCA(tree = tree,tips = taxa.vector) #MCRA de las especies
  new.tree<-add_into_node(tree=tree, node = MRCA, new.tip = new.tip)  #se añade desde MRCA
  return(new.tree)}
#plot(add_to_monophyletic(Tree, "Invent_inventata"))


# Function to add tips over a node. Vale para especies o géneros "hermanos".
add_over_node<- function(tree, new.tip, node){
  DF <- data.frame(tree$edge,tree$edge.length,1:length(tree$edge.length)) #datos del arbol
  EDGES <- DF[which.edge(tree, node),] #especificaciones del nodo dado
  to_index<-EDGES[,4] #posición a indexar

  WHERE <- tree$edge[to_index,2]   #posición a indexar
  LENGTH <- tree$edge.length[to_index]   #longitud de rama máxima

  MIN<-0
  MAX<-LENGTH
  POS<-runif(1,MIN,MAX)

  while(POS==MIN | POS==MAX){POS<- runif(1,MIN,MAX)} #encontrar un valor interno

  newtree<-bind.tip(tree, new.tip, edge.length=NULL, where=WHERE, position=POS)
  return(newtree)
}

### A function to stick species at random within a paraphyletic clade ###
add_into_paraphyletic_node<- function(tree, new.tip, group.node, intern.node){

  intern.tips<- tree$tip.label[getDescendants(tree, intern.node)][!is.na(tree$tip.label[getDescendants(tree, intern.node)])] #tips del nodo de "intruso"

  group.tree<-splitTree(tree,split=list(node=group.node,bp=tree$edge.length[which(tree$edge[,2]==group.node)])) #se rompe el árbol por el nodo profundo

  new.intern.node<- findMRCA(tree = group.tree[[2]],tips = intern.tips) #se encuentra en nodo interno en el nuevo arbol

  intern.tree <-splitTree(group.tree[[2]],split=list(node=new.intern.node,bp=group.tree[[2]]$edge.length[which(group.tree[[2]]$edge[,2]==new.intern.node)])) #el arbol cortado se vuelve a romper

  intern.tree[[1]]<-add.random(intern.tree[[1]],tips=new.tip) #se añade el new tip al árbol cortado sin el grupo intruso

  group.tree[[2]]<-paste.tree( intern.tree[[1]], intern.tree[[2]] )
  group.tree[[2]]$root.edge<-0
  new.tree<- paste.tree( group.tree[[1]], group.tree[[2]])#se pegan los 3 árboles


  if (is.ultrametric(new.tree)==FALSE){new.tree<-force.ultrametric(new.tree, method="extend")} #se fuerza el ultramétrico
  return(new.tree)
}


add_to_paraphyletic <- function(tree,new.tip){
  new.tip<-gsub(" ", "_",new.tip) #lo mismo que con el nodo, pero detectando de manera automática el nodo profundo y el intruso
  genus<- word(new.tip, 1, sep="_")
  sp.list<- data.frame(species=tree$tip.label)
  taxa.vector<- as.vector(sp.list[word(sp.list$species, 1, sep="_")==genus,])
  if(length(taxa.vector)==0){stop(paste0("Genus ", genus, " is no included in yout tree."))}
  MRCA<-findMRCA(tree = tree,tips = taxa.vector)
  Desc.Tips<-tree$tip.label[getDescendants(tree, MRCA)][!is.na(tree$tip.label[getDescendants(tree, MRCA)])]
  Intruder.Tips<-as.vector(Desc.Tips[word(Desc.Tips, 1, sep = "_")!= genus])


  genus.tree<-splitTree(tree,split=list(node=MRCA,bp=tree$edge.length[which(tree$edge[,2]==MRCA)]))

  Intruder.MRCA<- findMRCA(tree = genus.tree[[2]],tips = Intruder.Tips)

  pure.tree <-splitTree(genus.tree[[2]],split=list(node=Intruder.MRCA,bp=genus.tree[[2]]$edge.length[which(genus.tree[[2]]$edge[,2]==Intruder.MRCA)]))

  pure.tree[[1]]<-add.random(pure.tree[[1]],tips=new.tip)

  genus.tree[[2]]<-paste.tree( pure.tree[[1]], pure.tree[[2]] )
  genus.tree[[2]]$root.edge<-0
  new.tree<- paste.tree( genus.tree[[1]], genus.tree[[2]])


  if (is.ultrametric(new.tree)==FALSE){new.tree<-force.ultrametric(new.tree, method="extend")}




  return(new.tree)}
#plot(add_to_paraphyletic(Tree, "Acer_inventatus"))



polytomy_into_node<- function(tree, new.tip, node){

  if(!exists("tree")){stop("Tree object not found")}
  newtree<-tree #se crea un árbol nuevo para ser modificado

  DF <- data.frame(newtree$edge,newtree$edge.length,"id"=1:length(newtree$edge.length)) #características del DF
  EDGES <- DF[which.edge(newtree, node),] #la fila del nodo dado
  to_index<-EDGES[,4] #posicion a indexar

  WHERE <- newtree$edge[to_index,2]   #posicion a indexar



  newtree<-bind.tip(newtree, new.tip, edge.length=NULL, where=WHERE, position=0) #la posicion a indexar; la longitud del nodo es 0
  return(newtree)
}

polytomy_to_singleton<- function(tree, singleton, species, politomy.insertion=c("long","middle") ){

  newtree<-tree #arbol para ir editando
  for(i in 1:length(species)){
    if(i==1){   #la primera especie del listado se une al singleton
      sing.node<-which(newtree$tip.label==singleton)
      DF <- data.frame(newtree$edge,newtree$edge.length,"id"=1:length(newtree$edge.length))
      EDGES <- DF[which.edge(newtree, sing.node),]
      to_index<-EDGES[,4]
      WHERE <- newtree$edge[to_index,2]

      if(politomy.insertion=="middle"){position<-EDGES[,3]*0.5}  #se une a la mitad de la rama
      if(politomy.insertion=="long")  {position<-EDGES[,3]}      # se une al final de la rama

      newtree<-bind.tip(newtree, tip.label =  species[i], edge.length=NULL, where=WHERE, position=position)
    }else{
      node<-findMRCA(newtree, tips = c(singleton, species[1:(i-1)])) #el resto de especies se unen al nodo unido de singleton y añadidas
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

    if(length(added)==1) {    #similar a la politomía, pero con longitud de rama aleatoria

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




#para grupos polifiléticos
add_to_polyphyletic<-function(tree, species){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  species<-gsub(" ", "_",species)
  species.genus<- word(species, 1, sep="_")[!duplicated(word(species, 1, sep="_"))]
  if(length(species.genus)>1){stop("Species from more than 1 genera are included")}
  new.tree<-tree

  genus<- species.genus
  list<- data.frame(species=tree$tip.label)
  taxa.vector<- as.vector(list[word(list$species, 1, sep="_")==genus,]) #especies del género dentr del arbol. Este vector vale para todas las especies a añadir
  if(length(taxa.vector)==0){stop(paste0("Genus ", genus, " is no included in your tree."))}


  for(p in 1:length(species)){
    sticked.species<- species[p]
    sticking.species<- sample(taxa.vector, 1) #se selecciona una especie del vector de incluidas (esto conlleva una probabilidad basada en frecuencias dentro de nodos)
    sticking.tip<- which(new.tree$tip.label==sticking.species) #su valor de tip

    parent<-new.tree$edge[new.tree$edge[,2]==sticking.tip,1] #cuál es su ancestro
    siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])] #descendientes del ancestro

    if(length(siblings[word(siblings,1,sep="_")==genus])==1){new.tree<- add_to_singleton(new.tree,sticking.species,new.tips=sticked.species)} #los hemanos son de otro género; es un singleton en otro clado
    if(length(siblings[word(siblings,1,sep="_")==genus]) >1){#alguno de sus hermanos es congénere
      while(length(word(siblings,1,sep="_")[!duplicated(word(siblings,1,sep="_"))])==1){ #va subiendo de nivel el tip y el parent hasta que no sean del mismo género
        sticking.tip<-parent
        parent<-new.tree$edge[new.tree$edge[,2]==sticking.tip,1]
        siblings<- new.tree$tip.label[getDescendants(new.tree, parent)][!is.na(new.tree$tip.label[getDescendants(new.tree, parent)])]
      }

      gen.MRCA<- findMRCA(tree = new.tree,tips = siblings[word(siblings, 1, sep="_")==genus]) #nodo MRCA de los hermanos congéneres; probablemente siempre sea el "sticking tip"; pero podrá ser que no
      grouped<- new.tree$tip.label[getDescendants(new.tree, gen.MRCA)][!is.na(new.tree$tip.label[getDescendants(new.tree, gen.MRCA)])] #descendientes (no nodos; por eso en !is.na)
      grouped.gen<- word(grouped,1,sep="_")[!duplicated(word(grouped,1,sep="_"))] #género de los descendientes del MRCA
      if(length(grouped.gen)==1){new.tree<- add_into_node(new.tree, sticking.tip, sticked.species)} #subgrupo monofilético
      if(length(grouped.gen) >1){ #subgrupo parafilético siempre, por la metodología de la búsqueda no puede ser ni polifilético ni singleton
        intruders<-grouped[word(grouped, 1, sep="_")!=genus]
        if(length(intruders)==1){new.tree<- add_into_node(new.tree, sticking.tip, sticked.species)}else{ #intrusos singletons; a efectos prácticos se añade como si fuera monofilético
          intruders.MRCA<- findMRCA(tree = new.tree,tips = intruders)
          new.tree<- add_into_paraphyletic_node(new.tree, new.tip=sticked.species, group.node = gen.MRCA, intern.node = intruders.MRCA) #intruso monofiletico; se añade como parafilético


        }}}}

  return(new.tree)

}





#PROTOCOLO COMPLETO - FUNCIÓN MADRE randtip ####

#tabla<- read.xlsx("E:/UNI/4. DOCTORADO/3. Randtip/phylo.table.xlsx")


RANDTIP<- function(tree, species.table, type=c("random", "genus.polytomy", "family.polytomy", "order.polytomy", "class.polytomy"), aggregate.subspecies=TRUE, politomy.insertion="long"){

  if(type=="random"){
  species.table$using.taxa<-NA #se crea una columna nueva en la tabla con el nombre a utilizar para la primera ronda
  for(i in 1:nrow(species.table)){
    if(is.na(species.table$aggregate.subspecies[i])){    #para los NA (en la tabla no hay especificadas excepciones) se aplica lo marado en la funcion
      if(aggregate.subspecies==TRUE){species.table$using.taxa[i]<-paste0(word(species.table$taxon[i],1, sep="_"), "_",word(species.table$taxon[i],2, sep="_"))}
      else{species.table$using.taxa[i]<-species.table$taxon[i]}
      next}

    if(species.table$aggregate.subspecies[i]==0){species.table$using.taxa[i]<-species.table$taxon[i] #para los especificados como 0 no se agrupa
    next}
    if(species.table$aggregate.subspecies[i]==1){species.table$using.taxa[i]<-paste0(word(species.table$taxon[i],1, sep="_"), "_",word(species.table$taxon[i],2, sep="_"))}} #para los especificados como 1 sí se agrupa

  species.table.dupl <- species.table[ duplicated(species.table$using.taxa),]
  species.table <-      species.table[!duplicated(species.table$using.taxa),]



  newtree<- tree
  taxa <- species.table$using.taxa[!duplicated(species.table$using.taxa)]
  taxa<- gsub(" ", "_", taxa)
  taxa <- taxa[taxa%!in%tree$tip.label]

  taxa.genera<- word(taxa, 1, sep="_")[!duplicated(word(taxa, 1, sep="_"))]
  taxa.genera<-sample(taxa.genera, length(taxa.genera), replace = F)
  }else{

    species.table$using.taxa <-species.table$taxon  #En casos de politomía se utilizan los nombres tal cual


    newtree<- tree
    taxa <- species.table$using.taxa[!duplicated(species.table$using.taxa)]
    taxa<- gsub(" ", "_", taxa)
    taxa <- taxa[taxa%!in%tree$tip.label]

    taxa.genera<- word(taxa, 1, sep="_")[!duplicated(word(taxa, 1, sep="_"))]
    taxa.genera<-sample(taxa.genera, length(taxa.genera), replace = F)
  } #preparación de listado<-using.name

  if(type=="class.polytomy"){
    for(p in 1:length(taxa.genera)){
      genus<- taxa.genera[p]#se selecciona un genero

      genus.taxa <- taxa[word(taxa, 1, sep="_")==genus] #seleccionamos los taxones del género
      genus.taxa <- sample(genus.taxa, length(genus.taxa), replace = F) #aleatorizamos el orden

      genus.class<- species.table$class[species.table$genus==genus][!duplicated(species.table$class[species.table$genus==genus])]#identificamos la clase
      class.genera<- species.table$genus[species.table$class==genus.class]   #encontramos los géneros de esa clase
      union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%class.genera]   #encontramos los tips dentro de la clase
      if(length(union.tips)==0){message(paste0("ATTENTION: genus ", genus, " was not included as no Class coincidences were found")) #no se encuentra donde juntar
        next}else{
          if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, union.tips, genus.taxa, politomy.insertion = politomy.insertion)}else{
            for(x in 1:length(genus.taxa)){
              tip<-genus.taxa[x]
              newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=union.tips))
            }}}}}

  if(type=="order.polytomy"){
    for(p in 1:length(taxa.genera)){
      genus<- taxa.genera[p]

      genus.taxa <- taxa[word(taxa, 1, sep="_")==genus] #seleccionamos los taxones del género
      genus.taxa <- sample(genus.taxa, length(genus.taxa), replace = F) #aleatorizamos el orden

      genus.order<- species.table$order[species.table$genus==genus][!duplicated(species.table$order[species.table$genus==genus])]
      order.genera<- species.table$genus[species.table$order==genus.order]
      union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%order.genera]
      if(length(union.tips)==0){
        genus.class<- species.table$class[species.table$genus==genus][!duplicated(species.table$class[species.table$genus==genus])]
        class.genera<- species.table$genus[species.table$class==genus.class]
        union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%class.genera]
        if(length(union.tips)==0){message(paste0("ATTENTION: genus ", genus, " was not included as no Class coincidences were found"))}else{
          for(x in 1:length(genus.taxa)){
            tip<-genus.taxa[x]
            newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=))
          }}
      }else{
        if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, union.tips, genus.taxa,politomy.insertion = politomy.insertion)}else{
          for(x in 1:length(genus.taxa)){
            tip<-genus.taxa[x]
            newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=union.tips))
          }}}}}

  if(type=="family.polytomy"){
    for(p in 1:length(taxa.genera)){
      genus<- taxa.genera[p]

      genus.taxa <- taxa[word(taxa, 1, sep="_")==genus] #seleccionamos los taxones del género
      genus.taxa <- sample(genus.taxa, length(genus.taxa), replace = F) #aleatorizamos el orden

      genus.family<- species.table$family[species.table$genus==genus][!duplicated(species.table$family[species.table$genus==genus])]
      family.genera<- species.table$genus[species.table$family==genus.family]

      union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%family.genera]
      if(length(union.tips)==0){
        genus.order<- species.table$order[species.table$genus==genus][!duplicated(species.table$order[species.table$genus==genus])]
        order.genera<- species.table$genus[species.table$order==genus.order]
        union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%order.genera]
        if(length(union.tips)==0){
          genus.class<- species.table$class[species.table$genus==genus][!duplicated(species.table$class[species.table$genus==genus])]
          class.genera<- species.table$genus[species.table$class==genus.class]
          union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%class.genera]
          if(length(union.tips)==0){message(paste0("ATTENTION: genus ", genus, " was not included as no Class coincidences were found"))}else{
            if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, union.tips, genus.taxa,politomy.insertion = politomy.insertion)}else{
              for(x in 1:length(genus.taxa)){
                tip<-genus.taxa[x]
                newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=))
              }}}
        }else{

          if(length(union.tips)==1){newtree<-add_to_singleton(newtree, union.tips, genus.taxa)}else{
            for(x in 1:length(genus.taxa)){
              tip<-genus.taxa[x]
              newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=union.tips))
            }}}

      }else{
        if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, union.tips, genus.taxa, politomy.insertion = politomy.insertion)}else{
          for(x in 1:length(genus.taxa)){
            tip<-genus.taxa[x]
            newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=union.tips))
          }}}}}

  if(type=="genus.polytomy"){
    for(p in 1:length(taxa.genera)){
      genus<- taxa.genera[p]

      genus.taxa <- taxa[word(taxa, 1, sep="_")==genus] #seleccionamos los taxones del género
      genus.taxa <- sample(genus.taxa, length(genus.taxa), replace = F) #aleatorizamos el orden

      union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")==genus]
      if(length(union.tips)==0){

        genus.family<- species.table$family[species.table$genus==genus][!duplicated(species.table$family[species.table$genus==genus])]
        family.genera<- species.table$genus[species.table$family==genus.family]

        union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%family.genera]
        if(length(union.tips)==0){
          genus.order<- species.table$order[species.table$genus==genus][!duplicated(species.table$order[species.table$genus==genus])]
          order.genera<- species.table$genus[species.table$order==genus.order]
          union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%order.genera]
          if(length(union.tips)==0){
            genus.class<- species.table$class[species.table$genus==genus][!duplicated(species.table$class[species.table$genus==genus])]
            class.genera<- species.table$genus[species.table$class==genus.class]
            union.tips<-newtree$tip.label[word(newtree$tip.label, 1, sep="_")%in%class.genera]
            if(length(union.tips)==0){message(paste0("ATTENTION: genus ", genus, " was not included as no Class coincidences were found"))}else{
              if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, union.tips, genus.taxa, politomy.insertion = politomy.insertion)}else{
                for(x in 1:length(genus.taxa)){
                  tip<-genus.taxa[x]
                  newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=))
                }}}
          }else{

            if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, union.tips, genus.taxa, politomy.insertion = politomy.insertion)}else{
              for(x in 1:length(genus.taxa)){
                tip<-genus.taxa[x]
                newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=union.tips))
              }}}

        }else{
          if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, union.tips, genus.taxa, politomy.insertion = politomy.insertion)}else{
            for(x in 1:length(genus.taxa)){
              tip<-genus.taxa[x]
              newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=union.tips))
            }}}}else{
              if(length(union.tips)==1){newtree<-polytomy_to_singleton(newtree, singleton = union.tips,species =  genus.taxa, politomy.insertion = politomy.insertion)}else{
                for(x in 1:length(genus.taxa)){
                  tip<-genus.taxa[x]
                  newtree<- polytomy_into_node(newtree, tip, node=findMRCA(newtree, tips=union.tips))
                }}}}}

  if(type=="random"){
    for(i in 1: length(taxa.genera)){        #loop 1
      genus<- taxa.genera[i]
      genus.type<-species.table$genus.type[species.table$genus==genus][!duplicated(species.table$genus.type[species.table$genus==genus])]

      genus.taxa <- taxa[word(taxa, 1, sep="_")==genus] #seleccionamos los taxones del género
      genus.taxa <- sample(genus.taxa, length(genus.taxa), replace = F) #aleatorizamos el orden

      grouped.taxa<-species.table$using.taxa[!is.na(species.table$relative.species)]#se seleccionan los que no tengan información de relatives (se juntarán como grupo MF->monofilético)

      grouped.taxa<-grouped.taxa[grouped.taxa%in%genus.taxa]
      genus.taxa<-genus.taxa[genus.taxa%!in%grouped.taxa]

      if(length(grouped.taxa)>0){    #si existen "grouped.taxa"
        for( j in 1:length(grouped.taxa)){
          grouping.taxa<- species.table$relative.species[species.table$using.taxa==grouped.taxa[j]]
          grouping.taxa<-gsub(" ","",grouping.taxa)
          grouping.taxa<- strsplit(grouping.taxa, split = ",")[[1]]
          grouping.taxa<-grouping.taxa[grouping.taxa%in%tree$tip.label]#se seleccionan los taxones a agrupar que estan en el árbol
          if(length(grouping.taxa)==1){
            newtree<- add_to_singleton(newtree, singleton = grouping.taxa, new.tips = grouped.taxa[j])  #se añade el primero como singleton
          } else{
            node<-findMRCA(newtree, tips = grouping.taxa)
            newtree<- add_into_node(newtree, new.tip = grouped.taxa[j],node = node)} #se añade el resto como grupo MF
          rm(j, node, grouping.taxa)}
        next}

      #a partir de aqui se trabaja con genus.taxa; es decir, los no agrupados
      if(genus.type=="MONOPHYLETIC"){
        for( j in 1:length(genus.taxa)){
          newtree<- add_to_monophyletic(newtree, new.tip = genus.taxa[j])#para el género dado, si es MF, se añaden los taxones uno a uno como MF
          rm(j)}}

      if(genus.type=="PARAPHYLETIC"){
        for( j in 1:length(genus.taxa)){
          newtree<- add_to_paraphyletic(newtree, new.tip = genus.taxa[j])#para el género dado, si es PaF, se añaden los taxones uno a uno como PaF
          rm(j)}}

      if(genus.type=="POLYPHYLETIC"){
        for( j in 1:length(genus.taxa)){
          newtree<- add_to_polyphyletic(newtree, new.tip = genus.taxa[j])#para el género dado, si es PoF, se añaden los taxones uno a uno como PoF
          rm(j)}}

      if(genus.type=="SINGLETON GENUS"){   #para el género dado, si es singleton, se añaden los taxones de una vez
        newtree<- add_to_singleton(newtree, singleton = tree$tip.label[word(tree$tip.label,1,sep="_")== genus], new.tips = genus.taxa)
      }

      if(genus.type=="NOT INCLUDED"){
        genus.synonyms<- species.table$synonim.genus[species.table$genus==genus&!is.na(species.table$synonim.genus)]#primero busca si es sinónimo...
        genus.siblings<- species.table$sibling.genus[species.table$genus==genus&!is.na(species.table$sibling.genus)]#... o hermano

        genus.synonyms.tips<-species.table$using.taxa[species.table$genus==genus&!is.na(species.table$synonim.genus)]
        genus.siblings.tips<-species.table$using.taxa[species.table$genus==genus&!is.na(species.table$sibling.genus)]

        genus.taxa.NO.except<-genus.taxa[genus.taxa%!in%c(genus.synonyms.tips,genus.siblings.tips)]#taxones si estas excepciones

        if(length(genus.synonyms)>0){#para g géneros con excepción de sinónimo
          for(i in 1:length(genus.synonyms.tips)){
            genus.synonyms.tips[i]<-paste0(genus.synonyms[i], "_",genus.synonyms.tips[i]) #añade el género sinónimo al principio del nombre; es necesario que esté en el árbol
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
        if(length(genus.siblings)>0){#si las excepciones incluyen hermanos
          if(length(genus.siblings[!duplicated(genus.siblings)])>1){stop(paste0("Multiple genus siblings for genus ",genus))}
          new.genus<- genus.siblings[!duplicated(genus.siblings)]

          genus.siblings.tips<-sample(genus.siblings.tips, length(genus.siblings.tips), replace = F)

          adding.node<- findMRCA(newtree, tips = newtree$tip.label[word(newtree$tip.label, 1, sep="_")==new.genus])
          newtree<-add_over_node(newtree, new.tip = genus.siblings.tips[1], node=adding.node)#se añade sobre el MRCA del género hermano el primer taxon

          newtree<-add_to_singleton(newtree, singleton = genus.siblings.tips[1],new.tips =  genus.siblings.tips[2:length(genus.siblings.tips)])#el resto se añaden al primero como singleton
        }
        if(length(genus.taxa.NO.except)>0){#para los que no hay excepciones se añaden aleatoriamente dentro de la familia, orden o clase

          fam<-species.table$family[species.table$using.taxa%in%genus.taxa.NO.except][!duplicated(species.table$family[species.table$using.taxa%in%genus.taxa.NO.except])]
          genera<- species.table$genus[species.table$family==fam]
          if(length(genera)>0){
            family.taxa<- newtree$tip.label[word(newtree$tip.label,1, sep="_")%in%genera]
            newtree<- add_into_node(tree=newtree, node=findMRCA(tree = newtree,tips = family.taxa),new.tip = genus.taxa.NO.except[1])
            if(length(genus.taxa.NO.except)>1){newtree<- add_to_singleton(tree=newtree, singleton = genus.taxa.NO.except[1],new.tips = genus.taxa.NO.except[2:length(genus.taxa.NO.except)])}

          }else{
            ord<-species.table$order[species.table$using.taxa%in%genus.taxa.NO.except][!duplicated(species.table$family[species.table$using.taxa%in%genus.taxa.NO.except])]
            genera<- species.table$genus[species.table$order==ord]
            if(length(genera)>0){
              order.taxa<- newtree$tip.label[word(newtree$tip.label,1, sep="_")%in%genera]
              newtree<- add_into_node(tree=newtree, node=findMRCA(tree = newtree,tips = order.taxa),new.tip = genus.taxa.NO.except[1])
              if(length(genus.taxa.NO.except)>1){newtree<- add_to_singleton(tree=newtree, singleton = genus.taxa.NO.except[1],new.tips = genus.taxa.NO.except[2:length(genus.taxa.NO.except)])}
            }else{
              cla<-species.table$class[species.table$using.taxa%in%genus.taxa.NO.except][!duplicated(species.table$family[species.table$using.taxa%in%genus.taxa.NO.except])]
              genera<- species.table$genus[species.table$class==cla]
              if(length(genera)>0){
                class.taxa<- newtree$tip.label[word(newtree$tip.label,1, sep="_")%in%genera]
                newtree<- add_into_node(tree=newtree, node=findMRCA(tree = newtree,tips = class.taxa),new.tip = genus.taxa.NO.except[1])
                if(length(genus.taxa.NO.except)>1){newtree<- add_to_singleton(tree=newtree, singleton = genus.taxa.NO.except[1],new.tips = genus.taxa.NO.except[2:length(genus.taxa.NO.except)])}
              }else{
                stop(paste0("Genus ", genus, " has no family, order or class relatives"))}}}
        }
      }

      print(paste0(i, " out of ", length(taxa.genera), ". Genus ", genus, " (", genus.type,"): ", round(i/length(taxa.genera)*100, 2), " %"))
    }

    #fase 2. las subespecies
    if(nrow(species.table.dupl)>0){#las especies agrupadas
      rep.taxa <-species.table.dupl$taxon#se encuentra la especie
      rep.taxa.species<- species.table.dupl$using.taxa[!duplicated(species.table.dupl$using.taxa)]#especies de las subespecies

      for (r in 1:length(rep.taxa.species)){
        ssps <- rep.taxa[paste0(word(rep.taxa, 1, sep="_"),"_",word(rep.taxa, 2, sep="_"))==rep.taxa.species[r]] #seleccionamos las subespecies de cada especie
        ssps <- sample(ssps, length(ssps), replace = F) #aleatorizamos el orden
        sp.to.add<- newtree$tip.label[paste0(word(newtree$tip.label,1, sep="_"), "_", word(newtree$tip.label,2, sep="_"))== rep.taxa.species[r] ]

        newtree<- add_to_singleton(newtree, singleton = sp.to.add, new.tips = ssps)#se añaden todas a su subespecie hermana como singleton

      }

    }

    for(n in 1:nrow(species.table)){  #se devuelven los nombres originales
      if(species.table$using.taxa[n]== species.table$taxon[n]){next}
      if(species.table$using.taxa[n]%!in% newtree$tip.label){next}
      if(species.table$using.taxa[n]%in% newtree$tip.label){
        newtree$tip.label[newtree$tip.label==species.table$using.taxa[n]]<-species.table$taxon[n]
      }}


  }
  return(newtree)
}

plot(RANDTIP(tree=Tree, species.table=tabla.info, aggregate.subspecies=TRUE, type = "class.polytomy", insertion="long"))

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
