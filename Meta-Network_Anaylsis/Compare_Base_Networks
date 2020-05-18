library(sna);library(network);library(centiserve);library(intergraph)
library(CINNA)
library(igraph);library(phyloseq)
pathP<- "P:/Surrogate/Network_Comparison/Data/"
pathOld<- "H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/"
pathS<-"P:/Surrogate/Network_Comparison/Scripts/"
#source(paste0(pathS,"network_functions13_05.R"))
source(paste0(pathS,"d_statistic.R"))
source(paste0(pathS,"network_ps.R"))
allphyslv2<- readRDS(paste0(pathP,"renamed_phylo_list.RDS"))

pathPlots<-"P:/Surrogate/Network_Comparison/Plots/"

baseGraphs<- readRDS(file=paste0(pathP,"baseGraphs.RDS"))
baseAdjs<- readRDS(file=paste0(pathP,"baseAdjs.RDS"))

#theEnsembleAdjs<-readRDS(paste0(pathP,"ensemble_adjs.RDS"))
#saveRDS(notstringent_adj,paste0(pathP,"notstringent_adjs.RDS"))
notstng<- readRDS(paste0(pathP,"notstringent_ntks.RDS"))

dbase<- lapply(baseGraphs,dmatrix)
saveRDS(dbase,file=paste0(pathP,"baseDstats.RDS"))


```

#Here, I compare the base networks to one another, to see if they correlate. 
```{r}

edge1<- lapply(notstng,function(x)length(E(x)))
edge2<- lapply(theEnsembleAdjs$Less_Stringent,function(x)length(E(x)))
edge3<- lapply(theEnsembleAdjs$More_Stringent,function(x)length(E(x)))
#Remove unconnected vertices

# Extract edges and vertices
allVEs<- lapply(baseGraphs2,function(x){
  lapply(x,function(y){
    
    ifelse(anyNA(y)== FALSE & is.null(y)==FALSE,
           return(cbind(length(E(y)),length(V(y)))),
           return(NA))
  })
})
# Compare edges/nodes between networks

allVEs3<- lapply(allVEs,function(x)plyr::ldply(x,as.data.frame))

allVEs2<- plyr::ldply(lapply(allVEs,function(x)plyr::ldply(x,as.data.frame)),as.data.frame)


Graphtype<- c(rep(names(allVEs)[1],108),rep(names(allVEs)[2],108),rep(names(allVEs)[3],108),
              rep(names(allVEs)[4],108))#,rep(names(allVEs)[5],108),rep(names(allVEs)[6],108))

allVEs2$type_sample<- paste0(Graphtype,"__",allVEs2$.id)
allVEs2$type<- Graphtype

ggNode<- ggplot(allVEs2,aes(x=type,y=V1))+geom_boxplot(size=1.5)+
  theme_light()+scale_y_log10()+labs(y="log10 Edges")

ggEdge<- ggplot(allVEs2,aes(x=type,y=V2))+geom_boxplot(size=1.5)+
  theme_light()+scale_y_log10()+labs(y="log10 Nodes")

ggsave(arrangeGrob(ggNode,ggEdge,ncol=2),
       file=paste0(pathPlots,"edgenodcompare.png"),dpi=300)

# Make correlation plots to compare 


cormat<- matrix(ncol=4,nrow=4)
rownames(cormat)<-colnames(cormat)<- names(allVEs3)[1:4]
nodecor<- cormat
for(i in  1:4){
  for(j in 1:4){
    cormat[i,j]<- cor(log(allVEs3[[i]]$V1+1),log(allVEs3[[j]]$V1+1),use="complete.obs")
    nodecor[i,j]<- cor(log(allVEs3[[i]]$V2+1),log(allVEs3[[j]]$V2+1),use="complete.obs")
    
  }
}

ggedgecor<- ggcorrplot::ggcorrplot(cormat,type="lower",outline.color="white")
ggnodecor<-  ggcorrplot::ggcorrplot(nodecor,type="lower",outline.color="white")
grid.arrange(ggedgecor,ggnodecor,ncol=2)
ggsave(arrangeGrob(ggedgecor,ggnodecor,ncol=2),
       file=paste0(pathPlots,"base_ntk_corrplot.png"),dpi=300) 

name_mat<- matrix(ncol=4,nrow=4)
for(i in  1:4){
  for(j in 1:4){
    name_mat[i,j]<- paste0(names(allVEs3)[i],"-",names(allVEs3)[j])
    
  }
}

name_vector<- as.vector(name_mat[lower.tri(name_mat)])
bray_df<- matrix(ncol=6,nrow=108)

# Calculate the similarity between degree dist between each graph nb: 1- BC to make it similarity 
for(i in 1:108){ 
  rm(bry1)
  require(vegan)
  print(i)
  gl<- list(baseGraphs2[[1]][[i]], baseGraphs2[[2]][[i]] ,
            baseGraphs2[[3]][[i]] , baseGraphs2[[4]][[i]])
  tryCatch({
    nodeNames<- lapply(gl,function(x)names(V(x)))
    derepn<- unique(unlist(nodeNames))
    node_mat<- matrix(0,nrow=4,ncol=length(derepn))
    
    for(j in seq_along(gl)){
      thenodes<-  na.omit(match(names(V(gl[[j]])),derepn))
      node_mat[j,thenodes ]<-  degree(gl[[j]])}
  },error=function(e)NA)
  
  bry1<- as.vector(vegdist(log(node_mat+1)))
  bray_df[i,]<- 1-bry1
}


colnames(bray_df)<- name_vector
bray_long<- reshape2::melt(bray_df)
ggbraycompare<- ggplot(bray_long,aes(x=Var2,y=value))+
  geom_boxplot(size=1.2)+theme_light()+labs(x="Graph Type",y="1-Bray-Curtis")

ggsave(ggbraycompare,
       file=paste0(pathPlots,"base_ntk_bray.png"),dpi=300) 

dbase<- lapply(baseGraphs,dmatrix)



basenames<- names(baseGraphs)
baseGraphs2<- baseGraphs
network_types<- vector("list",length=4)
for(i in 1:4){
  names(baseGraphs2[[i]])<- paste0(basenames[i], names(baseGraphs2[[i]]))
  network_types[[i]]<- rep(names(baseDstats[i]),nrow(baseDstats[[i]])) 
  
}
allgraphs<- c(baseGraphs2[[1]],baseGraphs2[[2]],
              baseGraphs2[[3]],baseGraphs2[[4]])
all_ds<- dmatrix(allgraphs)
ntkmeta<- ldply(network_types,as.data.frame)
allcmds<- cmdscale(as.dist(all_ds),eig=TRUE)
allcmds<- data.frame(allcmds$points,ntkmeta)

ggplot(data=allcmds,aes(x=X1,y=X2,color=X..i..))+
  geom_point(size=2)+
  theme_light()+
  labs(x="Axis 1 [ 32.8 %] ",y="Axis 2 [ 26.3 % ]" )


```


Comparing the base networks to the ensemble entworks
```{r}
### D stats for networks 
dstats<- lapply(ensembleGraphs,dmatrix)
saveRDS(dstats,file=paste0(pathP,"ensembleDstats.RDS"))
dstats<- lapply(baseGraphs,dmatrix)
saveRDS(dtats,file=paste0(pathP,"baseDstats.RDS"))

# Read in ensemble networks, D stats 


## Compare mantel of the network structure. 

m_stat_results<- matrix(ncol=4,nrow=4)
for(k in 1:4){
  keep<- na.omit(match(rownames(baseDstats[[k]]),rownames(ensembleDstats[[2]])))
  m_stat<-  mantel(as.dist(ensembleDstats[[2]][keep,keep]),as.dist(baseDstats[[k]]))
  m_stat_results[k,1:2]<- c(m_stat$statistic,m_stat$signif)
  
  keep<- na.omit(match(rownames(baseDstats[[k]]),rownames(ensembleDstats[[1]])))
  keep2<-match(rownames(ensembleDstats[[2]])[7],rownames(baseDstats[[k]]))
  m_stat<-  mantel(ensembleDstats[[2]][keep,keep],baseDstats[[k]][-keep2,-keep2])
  m_stat_results[k,3:4]<- c(m_stat$statistic,m_stat$signif)
  
}
colnames(m_stat_results)<- c("NS r","NS Pr.","S r","S Pr.")
rownames(m_stat_results) <- names(baseGraphs)

round(m_stat_results,3) %>% 
  kable()%>%
  kable_styling(full_width = F) %>%
  as_image(width=0.1,"png")

#An index of samples in each set of networks for later
nameList<- vector("list",length=4)
for(i in 1:4){
  nameList[[i]]<- match(rownames(baseDstats[[i]]),rownames(ensembleDstats[[2]]))
}

base_cmds<- lapply(baseDstats,function(x)cmdscale(as.dist(x)))

# Combined sample data to each df 

for(i in 1:4){
  base_cmds[[i]]<-   data.frame(base_cmds[[i]],Author[nameList[[i]]],Organ[nameList[[i]]])
}

allcmds<- plyr::ldply(base_cmds,as.data.frame)


# Plot each dfs D-stat 
ggdstatbase<- ggplot(allcmds,aes(X1,X2,color=Author.nameList..i...))+
  geom_point(size=1.5,alpha=0.6)+
  facet_wrap(~.id,scale="free")+theme_light()+
  theme(legend.title = element_blank())

ggsave(ggdstatbase,file=paste0(pathPlots,"ggdstat_base.png"),dpi=300)

```
