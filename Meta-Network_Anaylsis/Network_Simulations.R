
library(mvtnorm)



species_number<- c(100)
sample_number<- c(25,50)
latent_number <- c(2,3,5)
latent_number <- c(2)

enviro_number <- c(1,3,5,7)


param_df<- expand.grid(species_number,sample_number,latent_number,enviro_number)
colnames(param_df)<- c("NSpecies","Nsamples","Nlatent","Ncoefs")



## Function to check for overdispersion (i.e., var > mean)
nb.overdispersion <- function(eta, phi, returnVals=FALSE){ 
  nb.mean <- apply(eta, 2, mean)
  nb.var <- nb.mean + phi*nb.mean^2 
  bool <- nb.var > nb.mean
  if(returnVals==TRUE){
    return(rbind(mean=nb.mean, var=nb.var))
  }else{
    return(bool)
  }
}

list1<- vector("list",length=3)
names(list1)<- c("True_counts","Enviro_counts","Enviro_cov")
list2<- vector("list",length=100)
for(k in 1:100){
  list2[[k]]<- list1
}
names(list2)<- 1:100
list3<- vector("list",length=nrow(param_df))
names(list3)<- paste(param_df[1,],param_df[2,],param_df[,3],param_df[,4])
for(k in 1:nrow(param_df)){
  list3[[k]]<- list2
}
allseeds<- 1:1e+6

for( a in  1:length(list3)){ 
  print(a)
  
  for( l in 1:length(list3[[1]])){
    theSeed<- sample(allseeds,1)  
    set.seed(theSeed)
    
  n.species<- param_df[a,1]
  n.samples<- param_df[a,2] 
  n.coefs<- param_df[a,4]
  n.latent<- param_df[a,3]

  eta_true <- z_true <- y_true <- matrix(0, nrow = n.samples, ncol = n.species)
  eta_enviro<- y_enviro<- z_envrio<-  matrix(0, nrow = n.samples, ncol = n.species)
  gamma_true <- rnorm(n.species,-1,1) 

  LV_true <- matrix(0,nrow=n.samples,ncol=n.latent)
  for(k in 1:n.samples) {LV_true[k,] <- mvtnorm::rmvnorm(1, mean = rep(0,n.latent))}
  Loading_true <- matrix(0,nrow=n.species,ncol=n.latent)
  for(k in 1:n.species) {Loading_true[k,] <- mvtnorm::rmvnorm(1, mean = rep(0,n.latent))}

  epsilon_true <- as.matrix(LV_true)%*%t(as.matrix(Loading_true))
 phi_true <- runif( n.species,0.0001, 1) # overdispersion parameter

 
 
 spp_response<- matrix(runif((n.species*n.coefs),-1,1),nrow=n.species,ncol=n.coefs)
 
 #This generates what the sample environment is 
 if(n.coefs==1){
   coef1<- c(runif(n.samples,0,1))#Continous
   the_coefs<- matrix(coef1,nrow=1)
   
 }else if(n.coefs==3){
   coef1<- runif(n.samples,0,1) #Continous
   coef2<- rbinom(n.samples,1,0.5)
   coef3<- runif(n.samples,0,1)#Continous
   the_coefs<- rbind(coef1,coef2,coef3)
   
 }else if (n.coefs==5){
   coef1<- c(runif(n.samples,0,1))#Continous
   coef2<- rbinom(n.samples,1,0.5)
   coef3<- c(runif(n.samples,0,1))#Continous 
   coef4<- rbinom(n.samples,1,0.5)
   coef5<- c(runif(n.samples,0,1))#Continous     
   the_coefs<- rbind(coef1,coef2,coef3,coef4,coef5)
   
 }
 else if (n.coefs==7){
   coef1<- c(runif(n.samples,0,1))#Continous
   coef2<- rbinom(n.samples,1,0.5)
   coef3<- c(runif(n.samples,0,1))#Continous 
   coef4<- rbinom(n.samples,1,0.5)
   coef5<- c(runif(n.samples,0,1))#Continous   
   coef6<- rbinom(n.samples,1,0.5)
   coef7<- c(runif(n.samples,0,1))
   the_coefs<- rbind(coef1,coef2,coef3,coef4,coef5,coef6,coef7)
   
 }
 response<-t(spp_response%*%the_coefs)
  
 enviro<- crossprod(response)
  for(i in 1:n.samples) {
    for(j in 1:n.species) {
      eta_true[i,j] <- alpha_true[i] + gamma_true[j] + epsilon_true[i,j] 
      y_true[i,j] <- rnbinom(n = 1, mu = exp(eta_true[i,j]), size = 1/phi_true[j]) # alternative parametrization via mean
      eta_enviro[i,j] <- alpha_true[i] + gamma_true[j] + epsilon_true[i,j] + response[i,j]
      y_enviro[i,j] <- rnbinom(n = 1, mu = exp(eta_enviro[i,j]), size = 1/phi_true[j]) # alternative parametrization via mean
      
    } 
  }
  

  y_true[is.na(y_true)] <- 0 # if NA creeps in; sometimes happens when the mean gets too large
  y_enviro[is.na(y_enviro)] <- 0 # if NA creeps in; sometimes happens when the mean gets too large
  
  
  list3[[a]][[l]][[1]]<- y_true
  list3[[a]][[l]][[2]]<- y_enviro
  list3[[a]][[l]][[3]]<- enviro
  }
}
  
  #nb_overdisp1 <- nb.overdispersion(eta_true,phi_true,returnVals=T) # check overdispersion

  
  mber<- function(cov1){ 
   require(pulsar)
    require(huge)
      lmax <- getMaxCov(cov1)
      lams <- getLamPath(lmax, lmax*.001, len=100)
      hugeargs <- list(lambda=lams, verbose=FALSE,lambda.min.ratio=1e-5,scr=TRUE,scr.num=(ncol(cov1)-1))
      out.mb    <- pulsar(cov1, fun=huge, fargs=hugeargs, rep.num=30,
                          criterion='stars',lb.stars = TRUE,ub.stars = TRUE)
      refit.mb<- refit(out.mb)
      optimal<- out.mb$stars$opt.index
      adj<- SpiecEasi::symBeta(refit.mb$est$beta[[optimal]])
      
      return(adj)
  }
 
  results<- vector("list",length=length(list3))
  
for(j in 1:length(list3)){ 
   
  df_1<- matrix(nrow=100,ncol=6)
  colnames(df_1)<- c("Interaction Edge Overlap","Enviromental Edge Overlap","Interaction edges","LV enviro Edge","Enviro Edges","Int-Env Overlap")
  
  
  for(k in 1:length(list3[[1]])){
  #for(k in 1:10){
    
    require(igraph)
    
  tryCatch({
  a1<- SpiecEasi::getRefit(SpiecEasi::spiec.easi(as.matrix(list3[[j]][[k]][[1]]), method="mb",nlambda=100,lambda.min.ratio=1e-5 ))
  g1<- graph_from_adjacency_matrix(a1,mode="undirected",diag=FALSE,weighted="Precision")
  
  a2<- SpiecEasi::getRefit(SpiecEasi::spiec.easi(as.matrix(list3[[j]][[k]][[2]]),method="mb",nlambda=100,lambda.min.ratio=1e-5   ))
  g2<- graph_from_adjacency_matrix(a2,mode="undirected",diag=FALSE,weighted="Precision")
  
  a3<- mber(as.matrix(list3[[j]][[k]][[3]]))
  g3<- graph_from_adjacency_matrix(a3,mode="undirected",diag=FALSE,weighted="Precision")
  
  int_g1_g2<- igraph::intersection(g1,g2)# between models
  int_g1_g3<- igraph::intersection(g1,g3)# true LV and enviro
  int_g2_g3<- igraph::intersection(g2,g3) # Enviro 
  
  interactions_edges <-  length(E(int_g1_g2))
  enviro_edges <-  length(E(int_g2_g3))
  int_enviro_true<- length(E(int_g1_g3))
  df_1[k,1]<- interactions_edges
  df_1[k,2]<- enviro_edges
  df_1[k,3]<- length(E(g1))
  df_1[k,4]<- length(E(g2))
  df_1[k,5]<- length(E(g3))
  df_1[k,6]<- int_enviro_true
  rm(g1);rm(g2);rm(g3);rm(a1);rm(a2);rm(a3)},error=function(e)NA)
  }
  
  results[[j]]<- df_1
  
 }
  names(results)<- paste0(param_df[,4], "--",param_df[,2])
#saveRDS(results,file=paste0(pathP,"_nb_environtk_Simulations.RDS" ))
  pathP<- "P:/Surrogate/Network_Comparison/Data/"
  
results2<- plyr::ldply(results,as.data.frame)

library(ggplot2)
library(gridExtra)
results2$enviro_frac<- results2$`Enviromental Edge Overlap`/results2$`LV enviro Edge`*100
results2$lv_frac<- results2$`Interaction Edge Overlap`/results2$`LV enviro Edge` *100
results2$lv_env_frac <- results2$`Int-Env Overlap`/results2$`Interaction edges` *100

ggenv<- ggplot(results2,aes(x=.id,y= enviro_frac))+geom_boxplot(size=1.1)+theme_light()+labs(x="",y="Environmental Correlations")
gglv<- ggplot(results2,aes(x=.id,y= lv_frac))+geom_boxplot(size=1.1)+theme_light()+labs(x="",y="True Interactions")
gglvenv<- ggplot(results2,aes(x=.id,y= lv_env_frac))+geom_boxplot(size=1.1)+theme_light()+labs(x="",y="True Interactions")
pathPlots<- "P:/Surrogate/Network_Comparison/Plots/"

ggsave(arrangeGrob(gglv,ggenv,ncol=2),file=paste0(pathP,"env_ntk_simulation_results.png"),dpi=300)
ggN<- ggplot(results2,aes(x=.id,y= results2$`LV enviro Edge`))+geom_boxplot(size=1.1)+theme_light()+labs(x="",y="True Interactions")
ggsave(arrangeGrob(gglv,ggenv,gglvenv,ncol=2),file=paste0(pathPlots,"env_ntk_simulation_results.png"),dpi=300)

