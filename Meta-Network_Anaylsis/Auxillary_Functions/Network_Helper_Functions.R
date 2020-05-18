 

###########
require(igraph)
require(ggnetwork)
network_whole_precision<- function(x,physeq){
  ## Make network   
  corr1<- x
  cormat<- corr1[[5]]
  base::diag(cormat)<-0
  g1<- graph_from_adjacency_matrix(cormat,
                                   mode="undirected",weighted="correlation")
  #Attach taxonomy
  V(g1)$Phylum <- tax_table(physeq)[,3]
  V(g1)$Class <- tax_table(physeq)[,4]
  if(length(V(g1))==0){
    return(NA)
  }else{
  return(g1)}
  }
  


network_whole_graph<- function(x,physeq){
  ## Make network   
  corr1<- x
  cormat<- corr1[[2]]
  base::diag(cormat)<-0
  g1<- graph_from_adjacency_matrix(cormat,
                                   mode="undirected",weighted="correlation")
  #Attach taxonomy
  V(g1)$Phylum <- tax_table(physeq)[,3]
  V(g1)$Class <- tax_table(physeq)[,4]
  #V(g1)$Prevelance <- tax_table(physeq)[,1]
  ## Match edge names with taxa names after remove non interacting taxa
  if(length(V(g1))==0){
    return(NA)
  }else{
  cor_edge_list<- igraph::as_data_frame(g1,"edges")
  
  E(g1)$color <- ifelse(cor_edge_list[,3]>0,"Co-occur Positive","Co-occur Negative")
  
  E(g1)$weight <- sqrt(cor_edge_list[,3]^2)
  
  return(g1)}
}

ggplotter<- function(x,which.pal=prevcolors,nodecolor=TRUE, node_color="Prevelance", weights=TRUE){
  ggplot(x, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(size=weight,color=color, linetype="solid"),ncp=2,curvature = 0.01,alpha=0.3)+
 scale_color_manual(values=which.pal)+
    theme_blank()+
  geom_nodes(aes_string(color=node_color,size=1.2,alpha=0.5),alpha=0.9)+
 geom_nodes(aes_string(size=1.2,alpha=0.5),alpha=0.9)+
 theme_light()+xlab("")+ylab("")+
    theme(axis.text = element_blank(),legend.position = "none")+
    guides(size=FALSE,linetype=FALSE)+
    scale_y_continuous(limits=c(-0.05,1.05))+scale_x_continuous(limits=c(-0.05,1.05))
}




ggplotter<- function(x,which.pal=prevcolors,nodecolor=FALSE, node_color="Prevelance", weights=TRUE){
  ggplot(x, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size=weight,color=color, linetype="solid"),ncp=2,curvature = 0.01,alpha=0.3)+
    scale_color_manual(values=which.pal)+
    theme_blank()+
    geom_nodes(aes_string(color=node_color,size=1.2,alpha=0.5),alpha=0.9)+
    geom_nodes(aes_string(size=1.2,alpha=0.5),alpha=0.9)+
    theme_light()+xlab("")+ylab("")+
    theme(axis.text = element_blank(),legend.position = "none")+
    guides(size=FALSE,linetype=FALSE)+
    scale_y_continuous(limits=c(-0.05,1.05))+scale_x_continuous(limits=c(-0.05,1.05))
}

networkspieceasi<- function(x,physeq,whole=TRUE){
  ## Make network   
  g1<- x
  #Attach taxonomy
  V(g1)$name<- taxa_names(physeq)
  V(g1)$Phylum <- tax_table(physeq)[,3]
  V(g1)$Class <- tax_table(physeq)[,4]
  #V(g1)$Prevelance <- tax_table(physeq)[,1]
  if(whole==TRUE){
    return(g1)
  }
  if(whole==FALSE){cor_edge_list<- igraph::as_data_frame(g1,"edges")
  g2<- graph_from_data_frame(cor_edge_list,directed=FALSE)
  ## Match edge names with taxa names after remove non interacting taxa
  if(length(V(g2))==0){
    return(NA)
  }
  connections <- sort(table(c(cor_edge_list$from, cor_edge_list$to)), decreasing = TRUE)
  tokeep<- match(rownames(tax_table(physeq)),names(V(g2)))
  tokeep<- match(names(V(g2)),rownames(tax_table(physeq)))# new

  V(g2)$names1 <- rownames(tax_table(physeq))[na.omit(tokeep)]
  V(g2)$Phylum <- tax_table(physeq)[na.omit(tokeep),3]
  V(g2)$Class<- tax_table(physeq)[na.omit(tokeep),4]
  #V(g2)$Prevelance<- tax_table(physeq)[na.omit(tokeep),1]
  
  cor_edge_list2<- igraph::as_data_frame(g2,"edges")
  E(g2)$color <- ifelse(cor_edge_list2[,3]>0, "Co-occur Positive","Co-occur Negative")
  E(g2)$weight <- sqrt(cor_edge_list2[,3]^2)
  return(g2)}}
  
renamer<- function(physeq){
  newNames<- paste(tax_table(physeq)[,2],tax_table(physeq)[,3],tax_table(physeq)[,4],tax_table(physeq)[,5])
  if(anyDuplicated(newNames)!=0){
    newNames<- make.unique(newNames)
  }
  taxa_names(physeq)<- newNames
  return(physeq)
}


taxwording<- function(physeq){
  tax_table(physeq)[tax_table(physeq)=="unclassified"]<-"Unassigned"
  tax_table(physeq)[tax_table(physeq)=="unassigned"]<-"Unassigned"
  tax_table(physeq)[grepl("incertae_sedis",
                          tax_table(physeq))] <- "Unassigned"
  tax_table(physeq)[tax_table(physeq)=="NA"]<-"Unassigned"
  
  tax_table(physeq)<- tax_table(physeq)[,1:5]
  
  return(physeq)
}

taxwordingGG<- function(physeq){
  tax_table(physeq)[tax_table(physeq)=="k__"]<-"__Unassigned"
  tax_table(physeq)[tax_table(physeq)=="p__"]<-"__Unassigned"
  tax_table(physeq)[tax_table(physeq)=="o__"]<-"__Unassigned"
  tax_table(physeq)[tax_table(physeq)=="c__"]<-"__Unassigned"
  tax_table(physeq)[tax_table(physeq)=="f__"]<-"__Unassigned"
  tax_table(physeq)[tax_table(physeq)=="NA"]<-"__Unassigned"
  tax_table(physeq)[is.na(tax_table(physeq))]<-"__Unassigned"
  
  tax_table(physeq)<- tax_table(physeq)[,1:5]
  
  return(physeq)
}


brassica_renamer<- function(physeq){
for(i in 1:nrow(tax_table(physeq))){
  for(j in 1:ncol(tax_table(physeq))){
   if(is.na(pmatch("Root;",tax_table(physeq)[i,j]))==FALSE){
    tax_table(physeq)[i,j]<- "Unassigned"
    print(paste("altered",i,j))
  }
  }
}
  return(physeq)
}



pnas_renamer<- function(physeq){
  
  physeq<- pnasM2
  foralteration<- grep(names(table(tax_table(physeq)[,"Family"])),pattern="e_")
        tax_table(physeq)[foralteration,"Family"]<- sapply(strsplit(tax_table(physeq)[foralteration,"Family"],"e_"),"[",1)
  return(physeq)
}


taxformat<- function(physeq){
  tax_table(physeq)[tax_table(physeq)=="unclassified"]<-"__Unassigned"
  
  tax_table(physeq)<-  apply(tax_table(physeq),2,function(y) sapply(strsplit(y,"__"),"[",2))
  return(physeq)}




################## Delete below? nah I dont think so


renameNetworks <- function(x,physeq,type="Precision"){
  
  
  fulltaxonomy<- taxa_names(physeq)
  
  if(type=="Precision"){
    cor1<- x[[5]]
    rownames(cor1)<- fulltaxonomy
    colnames(cor1)<- fulltaxonomy
    x[[5]]<-cor1}
  if(type=="Residual"){
    cor1<- x[[2]]
    rownames(cor1)<- fulltaxonomy
    colnames(cor1)<- fulltaxonomy
    x[[2]]<-cor1
  }
  if(type=="mb"){
    cor1<- x[[4]]
    rownames(cor1)<- fulltaxonomy
    colnames(cor1)<- fulltaxonomy
    x[[4]]<-cor1
  }
  return(x)}



graph_to_data_frame <- function(x,physeq,type="Precision"){
  
  fulltaxonomy<- taxa_names(physeq)
  
  if(type=="Precision"){
    
    cor1<- x[[5]]}
  
  if(type=="Residual"){
    cor1<- x[[2]]}
  
  
  rownames(cor1)<- fulltaxonomy
  colnames(cor1)<- fulltaxonomy
  cor2<-cor1[sort(rownames(cor1),decreasing = TRUE), ]
  cor3<- t(cor2)
  cor4<- cor3[sort(rownames(cor3),decreasing = TRUE), ]
  
  mtx<- cor4
  mtx[lower.tri(mtx, diag=TRUE)] <- NA
  dx <- as.data.frame(as.table(mtx))
  dx <- subset(dx, !is.na(Freq))
  
  taxaNames<- paste0(dx$Var2,"--" ,dx$Var1)
  dx2<- data.frame(dx,taxaNames)
  colnames(dx2)<-c("taxa1_name","taxa2_name","PrecisionValue","edgeID")
  return(dx2)
}



renameSpiec <- function(x,physeq){
  fulltaxonomy<- taxa_names(physeq)
  colanames(x)<- fulltaxonomy
  rownames(x)<-fulltaxonomy
  return(x)}

spiecConvert<- function(x,physeq){
  secor  <-    symBeta(getOptBeta(x), mode='maxabs')
  fulltaxonomy<- taxa_names(physeq)
  rownames(secor)<- fulltaxonomy
  colnames(secor)<- fulltaxonomy
  return(secor)
}
  

network_emp<- function(x,physeq){
  ## Make network   
  corr1<- x
  cormat<- corr1[[2]]
  diag(cormat)<-0
  g1<- graph_from_adjacency_matrix(cormat,
                                   mode="undirected",weighted="correlation")
  #Attach taxonomy
  cor_edge_list<- igraph::as_data_frame(g1,"edges")
  g2<- graph_from_data_frame(cor_edge_list,directed=FALSE)
  ## Match edge names with taxa names after remove non interacting taxa
  if(length(V(g2))==0){
    return(NA)
  }
  return(g2)}

network_mb<- function(x,physeq){
  ## Make network   
  corr1<- x
  cormat<- corr1[[4]]
  diag(cormat)<-0
  g1<- graph_from_adjacency_matrix(cormat,
                                   mode="undirected",weighted="correlation")
  #Attach taxonomy
 
  cor_edge_list<- igraph::as_data_frame(g1,"edges")
  g2<- graph_from_data_frame(cor_edge_list,directed=FALSE)
  if(length(V(g2))==0){
    return(NA)}
 return(g2)}
  
network_precision<- function(x,physeq){
  ## Make network   
  corr1<- x
  cormat<- corr1[[5]]
  diag(cormat)<-0
  g1<- graph_from_adjacency_matrix(cormat,
                                   mode="undirected",weighted="correlation")
  #Attach taxonomy
 
  cor_edge_list<- igraph::as_data_frame(g1,"edges")
  g2<- graph_from_data_frame(cor_edge_list,directed=FALSE)
  if(length(V(g2))==0){
    return(NA)}
  return(g2)} 



network_precision2<- function(x,physeq){
  ## Make network   
  cormat<- x

  diag(cormat)<-0
  g1<- graph_from_adjacency_matrix(cormat,
                                   mode="undirected",weighted="correlation")
  #Attach taxonomy
  V(g1)$Phylum <- tax_table(physeq)[,2]
  V(g1)$Class <- tax_table(physeq)[,3]
  cor_edge_list<- igraph::as_data_frame(g1,"edges")
  g2<- graph_from_data_frame(cor_edge_list,directed=FALSE)
  if(length(V(g2))==0){
    return(NA)}
  
  
  ## Match edge names with taxa names after remove non interacting taxa
  connections <- sort(table(c(cor_edge_list$from, cor_edge_list$to)), decreasing = TRUE)
  tokeep<- match(rownames(tax_table(physeq)),names(V(g2)))
  tokeep<- match(names(V(g2)),rownames(tax_table(physeq)))# new
  
  
  V(g2)$names1 <- rownames(tax_table(physeq))[na.omit(tokeep)]
  V(g2)$Phylum <- tax_table(physeq)[na.omit(tokeep),3]
  V(g2)$Class<- tax_table(physeq)[na.omit(tokeep),4]
  V(g2)$Prevelance<- tax_table(physeq)[na.omit(tokeep),1]
  
  
  cor_edge_list2<- igraph::as_data_frame(g2,"edges")
  
  E(g2)$color <- ifelse(cor_edge_list[,3]>0, "Co-occur Positive","Co-occur Negative")
  
  
  E(g2)$weight <- sqrt(cor_edge_list[,3]^2)
  
  
  return(g2)} 

dmatrix<- function(x,w1=0.45,w2=0.45,w3=0.1){
  require(Matrix)
  if(sum(is.na(x))>0){
 remove<- which(is.na(x))
 x<- x[-remove]}
  matres<- matrix(ncol=length(x),nrow=length(x))
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      matres[i,j]<- D(x[[i]],x[[j]],w1,w2,w3)
    }
  }
  rownames(matres)<- names(x)
  colnames(matres)<- names(x)
  return(matres)}

dmatrix2<- function(x,w1=0.45,w2=0.45,w3=0.1){
  require(Matrix)
  if(sum(is.na(x))>0){
 remove<- which(is.na(x))
 x<- x[-remove]}
  matres<- matrix(ncol=length(x),nrow=length(x))
  for(i in 1:(length(x)-1)){
    for(j in 2:length(x)){
      matres[i,j]<- D(x[[i]],x[[j]],w1,w2,w3)
    }
  }
  rownames(matres)<- names(x)
  colnames(matres)<- names(x)
  return(matres)}
  

#Could also weight by correlation strength





entropia<-function(a)
{
  a<-a[which(a>0)];
  -sum(a*log(a));
}



#returns the node distance matrix

node_distance<-function(g){
  
  n<-length(V(g))
  
  if(n==1){
    
    retorno=1
    
  }
  
  if(n>1){
    
    a<-Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    
    m<-shortest.paths(g,algorithm=c("unweighted"))
    
    m[which(m=="Inf")]<-n
    
    quem<-setdiff(intersect(m,m),0)
    
    for(j in (1:length(quem))){
      
      l<-which(m==quem[j])/n
      
      linhas<-floor(l)+1
      
      posicoesm1<-which(l==floor(l))
      
      if(length(posicoesm1)>0){
        
        linhas[posicoesm1]<-linhas[posicoesm1]-1
        
      }
      
      
      a[1:n,quem[j]]<-hist(linhas,plot=FALSE,breaks=(0:n))$counts
      
    }
    
    #m<-c()
    
    retorno=(a/(n-1))
    
  }
  
  return(retorno)
  
}


# nnd

nnd<-function(g){
  
  N<-length(V(g))
  
  nd<-node_distance(g)
  
  pdfm<-colMeans(nd)
  
  norm<-log(max(c(2,length(which(pdfm[1:(N-1)]>0))+1)))
  
  return(c(pdfm,max(c(0,entropia(pdfm)-entropia(nd)/N))/norm))
  
  
}

#function

alpha<-function(g){
  
  N<-length(V(g))
  
  r<-sort(alpha.centrality(g,exo=igraph::degree(g)/(N-1),alpha=1/N))/((N^2))
  
  return(c(r,max(c(0,1-sum(r)))))
  
}

#function

D<-function(g,h,w1,w2,w3){
  
  first<-0
  
  second<-0
  
  third<-0
  g<- g
  h<-h
  
  
  
  N<-length(V(g))
  
  M<-length(V(h))
  
  PM<-matrix(0,ncol=max(c(M,N)))
  
  if(w1+w2>0){
    
    pg=nnd(g)
    
    PM[1:(N-1)]=pg[1:(N-1)]
    
    PM[length(PM)]<-pg[N]
    
    ph=nnd(h)
    
    PM[1:(M-1)]=PM[1:(M-1)]+ph[1:(M-1)]
    
    PM[length(PM)]<-PM[length(PM)]+ph[M]
    
    PM<-PM/2
    
    first<-sqrt(max(c((entropia(PM)-(entropia(pg[1:N])+entropia(ph[1:M]))/2)/log(2),0)))
    
    second<-abs(sqrt(pg[N+1])-sqrt(ph[M+1]))
    
    
  }
  
  if(w3>0){
    
    pg<-alpha(g)
    
    ph<-alpha(h)
    
    m<-max(c(length(pg),length(ph)))
    
    Pg<-matrix(0,ncol=m)
    
    Ph<-matrix(0,ncol=m)
    
    Pg[(m-length(pg)+1):m]<-pg
    
    Ph[(m-length(ph)+1):m]<-ph
    
    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2
    
    g<-graph.complementer(g)
    
    h<-graph.complementer(h)
    
    
    pg<-alpha(g)
    
    ph<-alpha(h)
    
    m<-max(c(length(pg),length(ph)))
    
    Pg<-matrix(0,ncol=m)
    
    Ph<-matrix(0,ncol=m)
    
    Pg[(m-length(pg)+1):m]<-pg
    
    Ph[(m-length(ph)+1):m]<-ph
    
    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2
    
    
  }
  
  return(w1*first+w2*second+w3*third)
  
  
}

#function

d<-function(g,h,w1,w2,w3){
  
  first<-0
  
  second<-0
  
  third<-0
  
  g<-read.graph(g,format=c("edgelist"),directed=FALSE)
  
  h<-read.graph(h,format=c("edgelist"),directed=FALSE)
  
  N<-length(V(g))
  
  M<-length(V(h))
  
  PM<-matrix(0,ncol=max(c(M,N)))
  
  if(w1+w2>0){
    
    pg=nnd(g)
    
    PM[1:(N-1)]=pg[1:(N-1)]
    
    PM[length(PM)]<-pg[N]
    
    ph=nnd(h)
    
    PM[1:(M-1)]=PM[1:(M-1)]+ph[1:(M-1)]
    
    PM[length(PM)]<-PM[length(PM)]+ph[M]
    
    PM<-PM/2
    
    first<-sqrt(max(c((entropia(PM)-(entropia(pg[1:N])+entropia(ph[1:M]))/2)/log(2),0)))
    
    second<-abs(sqrt(pg[N+1])-sqrt(ph[M+1]))
    
    
  }
  
  if(w3>0){
    
    pg<-alpha(g)
    
    ph<-alpha(h)
    
    m<-max(c(length(pg),length(ph)))
    
    Pg<-matrix(0,ncol=m)
    
    Ph<-matrix(0,ncol=m)
    
    Pg[(m-length(pg)+1):m]<-pg
    
    Ph[(m-length(ph)+1):m]<-ph
    
    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))
    
    
  }
  
  return(w1*first+w2*second+w3*third)
  
  
}









  
  
  ggNetworker<- function(x,labels=TRUE,which.pal=prevcolors){
    if(labels==TRUE){
      g1<- ggplot(fortify(x), aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(size=1, linetype="solid"),ncp=1,curvature = 0.01,alpha=0.4) + 
        theme_blank()+geom_nodes(aes(size=1.2,alpha=0.5),alpha=0.9)+theme_light()+xlab("")+ylab("")+
        theme( axis.text = element_blank(),legend.position = "none")+ geom_nodelabel(aes(label = vertex.names ), fontface = "bold")+ scale_y_continuous(limits=c(-0.05,1.05))+
        scale_x_continuous(limits=c(-0.05,1.05))}
    if(labels==FALSE){
      g1<- ggplot(fortify(x), aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(color=color,size=1, linetype="solid"),ncp=1,curvature = 0.01,alpha=0.4) + 
        theme_blank()+geom_nodes(aes(size=1.2,alpha=0.5),alpha=0.9)+
        scale_color_manual(values=which.pal)+theme_light()+xlab("")+ylab("")+
        theme( axis.text = element_blank(),legend.position = "none")+  scale_y_continuous(limits=c(-0.05,1.05))+
        scale_x_continuous(limits=c(-0.05,1.05))}
    print(g1)
  }
  
  
  ggSpiec<- function(x,labels=TRUE){
    if(labels==TRUE){
      g1<- ggplot(fortify(x), aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(size=1, linetype="solid"),ncp=1,curvature = 0.01,alpha=0.4) + 
        theme_blank()+
        geom_nodes(aes(size=1.2,alpha=0.5),alpha=0.9)+
        theme_light()+xlab("")+ylab("")+
        theme( axis.text = element_blank(),legend.position = "none")+ geom_nodelabel(aes(label = vertex.names ), fontface = "bold")+ scale_y_continuous(limits=c(-0.05,1.05))+
        scale_x_continuous(limits=c(-0.05,1.05))}
    if(labels==FALSE){
      g1<- ggplot(x, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(size=1, linetype="solid"),ncp=1,curvature = 0.01,alpha=0.4) + 
        theme_blank()+geom_nodes(aes(size=1.2,alpha=0.5),alpha=0.9)+
        theme_light()+xlab("")+ylab("")+
        theme( axis.text = element_blank(),legend.position = "none")+  scale_y_continuous(limits=c(-0.05,1.05))+
        scale_x_continuous(limits=c(-0.05,1.05))}
    print(g1)
  }
  
  
  networkspieceasi2<- function(x,physeq,whole=TRUE){
    ## Make network   
    
    g1<- igraph::graph_from_adjacency_matrix(x,mode="undirected",weighted = "correlation")
    #Attach taxonomy
    V(g1)$name<- taxa_names(physeq)
    V(g1)$Phylum <- tax_table(physeq)[,3]
    V(g1)$Class <- tax_table(physeq)
    
    #V(g1)$Prevelance <- tax_table(physeq)[,1]
    if(whole==TRUE){
      return(g1)
    }
    if(whole==FALSE){cor_edge_list<- igraph::as_data_frame(g1,"edges")
    g2<- graph_from_data_frame(cor_edge_list,directed=FALSE)
    ## Match edge names with taxa names after remove non interacting taxa
    if(length(V(g2))==0){
      return(NA)
    }
    connections <- sort(table(c(cor_edge_list$from, cor_edge_list$to)), decreasing = TRUE)
    tokeep<- match(rownames(tax_table(physeq)),names(V(g2)))
    tokeep<- match(names(V(g2)),rownames(tax_table(physeq)))# new
    
    V(g2)$names1 <- rownames(tax_table(physeq))[na.omit(tokeep)]
    V(g2)$Phylum <- tax_table(physeq)[na.omit(tokeep),3]
    V(g2)$Class<- tax_table(physeq)[na.omit(tokeep),4]
    #V(g2)$Prevelance<- tax_table(physeq)[na.omit(tokeep),1]
    
    cor_edge_list2<- igraph::as_data_frame(g2,"edges")
    E(g2)$color <- ifelse(cor_edge_list2[,3]>0, "Co-occur Positive","Co-occur Negative")
    E(g2)$weight <- sqrt(cor_edge_list2[,3]^2)
    return(g2)}}
  

  
