library(limma);library(boot);library(igraph);library(dplyr)


## function for splitting ps objects
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
  
  # require(phyloseq)
  # require(plyr)
  
  ## Check the input
  if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
    stop("Sample data is missing in the phyloseq-object.\n")
  }
  
  ## Extract sample meta-data
  mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  
  if(!variable %in% colnames(mtd)){
    stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
  }
  
  if(class(mtd[, variable]) %in% c("integer", "numeric") ){
    if( length( unique(mtd[, variable]) ) > 5){
      stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
    } else {
      warning("Groupping variable is numeric and it was coerced to factor.\n")
      mtd[, variable] <- factor(mtd[, variable])
    }
  }
  
  if(length(table(mtd[, variable])) == 1){
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }
  
  ## Add sample IDs to the meta-data
  smp <- data.frame(
    SID = phyloseq::sample_names(physeq),
    mtd,
    stringsAsFactors = F)
  
  ## Extract sample names by the specified variable
  svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })
  
  ## Extract samples by groupping variable
  res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })
  
  ## Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }
  
  return(res)
}

##########################################################################\
library(limma);library(boot);library(igraph);library(dplyr)

# First need to generate weights, 3 options
# 1: weightd by network sample size NSS

# 2: Weight studies equally  EW
# 3: weight by prevelance of the taxa in the original communites - 

networkMeans<- function(physeq,type="NSS"){
### Split data by author, calculate the weights 

# log transform data for a reasonable scale 
physeq2<- transform_sample_counts(physeq,function(x)log(x+1))
physeq3<- as.data.frame(otu_table(physeq2))
# Convert data from s4 vector, and add the weights
physeq3<- data.frame(physeq3,weight=sample_data(physeq2)$N,study=sample_data(physeq2)$Author)

#Predefined results matrix 
n_taxa<- ntaxa(physeq)
#Boostrap sample statistic 
samplewmean <- function(data, whicht, i) {
  w <- data[i,n_taxa+1]
  d <- data[i,whicht]
  wm<- weighted.mean(d, w)
  return(wm)
  
}

studywmean <- function(data, whicht, i) {
  w <- data[i,c(whicht,c(n_taxa+1,n_taxa+2))]
  ds <- split(w,w[,3])
  
  ds2<- lapply(ds,function(x){
    return(mean(x[,1]))
    })
  
  ds3<- unlist(ds2)
  dwm<- sum(ds3)/6
  return(dwm)
  
}
cat("Starting Bootstrap")
if(type=="NSS" ){
## Calculate boostrap and P value (where the null is mu 0,10000 bootstraps)
  bs_results<- matrix(nrow=n_taxa,ncol=5)
  
for(k in 1:n_taxa){ 
  bs <- boot(data= physeq3, whicht=k,
             statistic = samplewmean, 
             R=10000) 
  h0<-  bs$t - mean(bs$t)
  p_val<- mean(abs(h0) > abs(bs$t0))
  bs_results[k,]<- c(p_val,bs$t0, quantile(bs$t,c(0.025,0.975)),var(bs$t))
}
bs_results<- as.data.frame(bs_results)  
colnames(bs_results)<- c("P_value","t0","CI[0.025%]","CI[0.975%]","Variance")

}else if(type=="SW"){
  
######## Study weighted ############
bs_results<- matrix(nrow=n_taxa,ncol=6)
for(k in 1:n_taxa){ 
    bs <- boot(data= physeq3, whicht=k,
             statistic = studywmean, 
             R=10000) 
  torm<- sum(is.na(bs$t))
  bst<- na.omit(bs$t)
  h0<-  bst - mean(bst)
  p_val<- sum(abs(h0) > abs(bs$t0))/(10000-torm)
  bs_results[k,]<- c(p_val,bs$t0, quantile(bst,c(0.025,0.975)),var(bst), torm)
}
bs_results<- as.data.frame(bs_results)  

colnames(bs_results)<- c("P_value","t0","CI[0.025%]","CI[0.975%]","Variance","Nremoved")
}
cat("Finished Bootstrap")
bs_results$P_fdr_adj<- p.adjust(bs_results[,1],"fdr")
bs_results$taxa <- taxa_names(physeq2)

return(bs_results)
}


############################ End #######################################################



## Start here 
## Define functions 

## Zero p.adjusted means, to incoporate multiple testing 
merge_network_estimates<- function(x){

allmeasures<- plyr::ldply(x,by=.id) # Merge lists 

allmeasures<- allmeasures %>% 
  mutate(P_adj= if_else(P_fdr_adj== 0,0.000000001, P_fdr_adj)) # Replace completes zeros with a pseudocount for fishers test

# Remove taxa only detected in 1 network, group by taxa and calc fishers, grand mean pooled var and Cis. 

res<- allmeasures %>% 
  
  add_count(taxa,name="taxacount")%>%
  
  filter(taxacount> 1)%>%
  
  group_by(taxa)%>% 
  
  summarise(se=sqrt(poolVar(Variance,n=107)$var)/sqrt(4),
            grand_mean=sum(t0)/4,
            fisher_P=metap::sumlog(P_adj)$p,
            LCI=grand_mean-se*1.96,
            UCI=grand_mean+se*1.96)

return(res)
}

############## Edge base



edgeMeans<- function(graphList,physeq=ls_ntkdf,type="unweighted"){

  
  
  
samplew_edgemean <- function(data, i) {
  w <- data[i,6]
  d <- data[i,4]
  wm<- weighted.mean(d, w)
  return(wm)
  
}


sample_edgemean <- function(data, i) {
  return(mean(data[i,4]))
}


# Convert data from s4 vector, and add the weights
theNs<- data.frame(weight=physeq$N,Community=rownames(physeq),Study=physeq$Author)

index1<- match(graphList$.id, theNs$Community)
graphList$N<- theNs$weight[index1]
graphList$Study<- theNs$Study[index1]

#Predefined results matrix 
#Boostrap sample statistic 
graphList2<- graphList %>% 
  add_count(edgeID,name="edgefreq")%>% 
  filter(edgefreq > 5) %>% 
  group_by(edgeID) %>%
  mutate(Nstudies= length(unique(Study)))%>%
  filter(Nstudies > 1  )%>%
  ungroup()

 

n_edges<- length(unique(graphList2$edgeID))

graphList2$edgeID<-  factor(graphList2$edgeID)
graphList3<- split(graphList2,graphList2$edgeID)


bs_results<- base::data.frame(matrix(nrow=n_edges,ncol=6))

if(type=="unweighted"){
for(k in 1:n_edges){ 
  bs <- boot(data= data.frame(graphList3[[k]]),
             statistic = sample_edgemean, 
             R=10000)
  
  h0<-  na.omit(bs$t) - mean(na.omit(bs$t))
  len<- sum(is.na(bs$t))
  p_val<- sum(abs(h0) > abs(bs$t0))/(length(bs$t)-len)
  
  bs_results[k,]<- c(p_val,bs$t0, quantile(na.omit(bs$t),c(0.025,0.975)),var(na.omit(bs$t)),
                    nrow(graphList3[[k]]))
}
}else{
  for(k in 1:n_edges){ 
    bs <- boot(data= data.frame(graphList3[[k]]),
               statistic = samplew_edgemean, 
               R=10000)
    
    h0<-  na.omit(bs$t) - mean(na.omit(bs$t))
    len<- sum(is.na(bs$t))
    p_val<- sum(abs(h0) > abs(bs$t0))/(length(bs$t)-len)
    
    bs_results[k,]<- c(p_val,bs$t0, quantile(na.omit(bs$t),c(0.025,0.975)),var(na.omit(bs$t)),
                       nrow(graphList3[[k]]))
  }
}

colnames(bs_results)<- c("P_value","t0","CI[0.025%]","CI[0.975%]","Variance","N")
bs_results$P_fdr_adj<- p.adjust(bs_results[,1],"fdr")

alldfs<- vector("list",length=n_edges)
for(k in 1:n_edges){
 d_f<- data.frame( to=graphList3[[k]][1,"to"],from=graphList3[[k]][1,"from"],edgeID=graphList3[[k]][1,"edgeID"])
 alldfs[[k]]<- d_f}

bs_results<- cbind(bs_results,plyr::ldply(alldfs))


return(bs_results)
}
#############
merge_edge_estimates<- function(x){
  
  allmeasures<- plyr::ldply(x,by=.id) # Merge lists 
  allmeasures$edge_ID<- paste0(allmeasures$to,"--",allmeasures$from)
  
  allmeasures<- allmeasures %>% 
    mutate(P_adj= if_else(P_fdr_adj== 0,0.000000001, P_fdr_adj),
           scaled_vals=scale(t0)) # Replace completes zeros with a pseudocount for fishers test
  
  # Remove taxa only detected in 1 network, group by taxa and calc fishers, grand mean pooled var and Cis. 
    res<- allmeasures %>% 
      
      add_count(edgeID,name="taxacount")%>%
      
      filter(taxacount> 1)%>%
      
      group_by(edge_ID)%>% 
      
      summarise(se=sqrt(poolVar(Variance,n=length(N))$var)/sqrt(length(N)),
                grand_mean=sum(t0)/4,
                fisher_P=metap::sumlog(P_adj)$p,
                LCI=grand_mean-se*1.96,
                UCI=grand_mean+se*1.96,
                grand_mean_scaled=sum(scaled_vals)/4,
                se_scaled=sd(scaled_vals)/sqrt(length(N)),
                LCI_scaled=grand_mean_scaled-se_scaled*1.96,
                UCI_scaled=grand_mean_scaled+ se_scaled*1.96)
    
  return(res)
}


##############################


edgeFreqs<- function(graphList,physeq=ls_ntkdf){

  
  
  # Convert data from s4 vector, and add the weights
  theNs<- data.frame(weight=physeq$N,Community=rownames(physeq),Study=physeq$Author)
  
  index1<- match(graphList$.id, theNs$Community)
  graphList$N<- theNs$weight[index1]
  graphList$Study<- theNs$Study[index1]
  
  #Predefined results matrix 
  #Boostrap sample statistic 
  graphList2<- graphList %>% 
    add_count(edgeID,name="edgefreq")%>% 
    filter(edgefreq > 5) %>% 
    group_by(edgeID,to,from) %>%
    mutate(Nstudies= length(unique(Study)))%>%
    filter(Nstudies > 1  )  %>%
    summarise(Pos_int = sum(PrecisionValue > 0 ), 
              Neg_int = sum(PrecisionValue < 0),
              ratio=  (Pos_int-Neg_int)/length(PrecisionValue),
              Total = length(PrecisionValue)  )

  return(graphList2)
}
  
 
#############
edge_freq_merge<- function(edges){
  
  
  
  ratio_edgemean <- function(data, i) {
    return(weighted.mean(data[i,7],data[i,8]))
  }
  
  
 edges<-  edges %>%
    add_count(edgeID,name="edgeCount")%>%
    filter(edgeCount> 1)
    edges$edgeID<- droplevels(edges$edgeID)
    edgeSplit<- split(edges,edges$edgeID)
    n_edges<- length(edgeSplit)
    bs_results<- base::data.frame(matrix(nrow=n_edges,ncol=5))
    
    for(k in 1:n_edges){ 
      
      bs <- boot(data= data.frame(edgeSplit[[k]]),
                 statistic = ratio_edgemean, 
                 R=10000)
      
      h0<-  na.omit(bs$t) - mean(na.omit(bs$t))
      len<- sum(is.na(bs$t))
      p_val<- sum(abs(h0) > abs(bs$t0))/(length(bs$t)-len)
      
      bs_results[k,]<- c(p_val,bs$t0, 
                         quantile(na.omit(bs$t),c(0.025,0.975)),
                         var(na.omit(bs$t)))
    }
    
    colnames(bs_results)<- c("P_value","t0","CI[0.025%]","CI[0.975%]","Variance")
    
    bs_results$P_fdr_adj<- p.adjust(bs_results[,1],"fdr")
    
    alldfs<- vector("list",length=n_edges)
    for(k in 1:n_edges){
      d_f<- data.frame( to=edgeSplit[[k]][1,"to"],from=edgeSplit[[k]][1,"from"],edgeID=edgeSplit[[k]][1,"edgeID"])
      alldfs[[k]]<- d_f}
    
    bs_results<- cbind(bs_results,plyr::ldply(alldfs))
    
    
    return(bs_results) 
  
}

#################


edgePrevalance<- function(graphList,physeq=ls_ntkdf){
  
  
  # Convert data from s4 vector, and add the weights
  theNs<- data.frame(weight=physeq$N,Community=rownames(physeq),Study=physeq$Author)
  
  index1<- match(graphList$.id, theNs$Community)
  graphList$N<- theNs$weight[index1]
  graphList$Study<- theNs$Study[index1]
  
  #Predefined results matrix 
  #Boostrap sample statistic 
  graphList2<- graphList %>% 
    add_count(edgeID,name="edgefreq")%>% 
    filter(edgefreq > 5) %>% 
    group_by(edgeID,to,from) %>%
    mutate(Nstudies= length(unique(Study)))%>%
    filter(Nstudies > 1  )%>%
    ungroup()
  
  graphList2$edgeID<-droplevels(graphList2$edgeID)
  
bystudy<-  split(graphList2,graphList2$Study)
theNs<- table(ls_ntkdf$Author)
bystudy2<-lapply(1:6,function(x){
  s<- theNs[x]
  i<- bystudy[[x]]
  i %>% group_by(edgeID,to,from)%>%
    add_count(edgeID,name="edgeprev")%>%
 summarise(Prev=sum(edgeprev)/length(edgeprev)/s*.1666667 )
    })

allstudies<- plyr::ldply(bystudy2,by=.id)

edgeprev<- allstudies %>% 
  group_by(edgeID,to,from)%>%
  summarise(total_prev=sum(Prev))

return(edgeprev)
}
  
 

prev_edge_mean<- function(x){

prev_edgemean <- function(data, i) {
  return(sum(data[i,5])/4)
}


edgeSplit<- split(x,x$edgeID)
n_edges<- length(edgeSplit)
bs_results<- base::data.frame(matrix(nrow=n_edges,ncol=5))

for(k in 1:n_edges){ 
  
  bs <- boot(data= data.frame(edgeSplit[[k]]),
             statistic = prev_edgemean, 
             R=10000)
  
  h0<-  na.omit(bs$t) - mean(na.omit(bs$t))
  len<- sum(is.na(bs$t))
  p_val<- sum(abs(h0) > abs(bs$t0))/(length(bs$t)-len)
  
  bs_results[k,]<- c(p_val,bs$t0, 
                     quantile(na.omit(bs$t),c(0.025,0.975)),
                     var(na.omit(bs$t))
                     )
}

colnames(bs_results)<- c("P_value","t0","CI[0.025%]","CI[0.975%]","Variance")

bs_results$P_fdr_adj<- p.adjust(bs_results[,1],"fdr")

alldfs<- vector("list",length=n_edges)
for(k in 1:n_edges){
  d_f<- data.frame( to=as.character(edgeSplit[[k]][1,"to"]),from=as.character(edgeSplit[[k]][1,"from"]))
  alldfs[[k]]<- d_f}

bs_results<- cbind(bs_results,plyr::ldply(alldfs))


return(bs_results) 

}
