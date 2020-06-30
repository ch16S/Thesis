
ppes<- readRDS(file=paste0("P:/Surrogate/ecological modelling/",
                                                  "gompertz_hs_tauhalfnormal.RDS"))

sig_alpha<- get_sigs(ppes,"alpha")
g1<- graph_from_adjacency_matrix(sig_alpha,mode="directed",weighted="Density",
                                 diag=FALSE)
V(g1)$Family<- (ps_data$tax[,"Family"])
V(g1)$Genus<-  (ps_data$tax[,"Genus"])
V(g1)$Class<-  (ps_data$tax[,"Class"])
g1.2<- g1
g1.2 <- set.vertex.attribute(g1.2, "name", value=as.character(ps_data$tax[,"Genus"]))
edgedf<- as_data_frame(g1.2)
edgedf
library(phyloseq)
g2.2<- delete_vertices(g1.2,degree(g1.2)==0)
ps_data<- readRDS(file=paste0("P:/Surrogate/ecological modelling/","contam_model_data_by_rajsd.RDS"))

g2<- delete_vertices(g1,degree(g1)==0)

V(g2)$color <- ifelse(V(g2)$Class == "c__[Saprospirae]","lightblue",
                      ifelse(V(g2)$Class == "c__Gammaproteobacteria","lightgreen",
                             ifelse(V(g2)$Class ==  "c__Alphaproteobacteria","orange",
                                    ifelse(V(g2)$Class ==  "c__Betaproteobacteria","yellow",
                                           ifelse(V(g2)$Class == "c__Cytophagia" ,"pink" ,"green" )))))

V(g2)$Membership<- ifelse(infoclust$membership==1,"lightblue","orange")

E(g2)$lty<- ifelse(E(g2)$Density >0,2,1)



png("DM_gompertz_density_class.png", width = 465, height = 225, units='mm', res = 300)

plot(g2, arrow.size=0.2, vertex.frame.color=NA,vertex.label.dist=0.0,vertex.label=NA, vertex.color=V(g2)$color,
     vertex.label.cex=0.7,vertex.label.font=.02,edge.arrow.size=0.4,vertex.size= 6, layout=layout_nicely(g2))
# Add a legend
legend("bottomleft", legend=unique(V(g2)$Class)  , col = unique(V(g2)$color) , bty = "n", pch=20 , pt.cex = 3, cex = .75*2, , horiz = FALSE, inset = c(0.1, 0.1))
dev.off()

dm_plv<- readRDS(file=paste0("P:/Surrogate/ecological modelling/",
                             "DM_PLV_GLLVM.RDS"))


sort_axis<- function(axes){
  axes_ordered<- vector("list",length= ncol(axes)*2)
  for(i in 1:ncol(axes)){
    axes_ordered[[i]]<-  order(axes[,i],decreasing=TRUE)
    axes_ordered[[(i+ncol(axes))]] <- order(axes[,i], decreasing=FALSE)
  }
  results<- as.data.frame(axes_ordered)
  return(results)
}

# Run the same models with random partitioned
lvs<- sort_axis(dm_plv$PPE$lvs )
lo<- c(36:39)

ps<- data.frame(ps_data$otus)[lvs[,5],]

n_obs1<- nrow(ps)
forcast_len<- 4
test<-  ps[(n_obs1-forcast_len+1):n_obs1,]
y<- ps[1:(n_obs1-forcast_len),]
n_sp<- ncol(y)
n_obs<- nrow(y)

samps<- nrow(ppes$Draws[[1]])
mod_cor<- pred_cor<- mod_mae<-mod_rmse<- pred_mae<- pred_rmse<- vector(length=samps)
preds<- array(dim=c(forcast_len,n_sp,samps))


for( i in 1:samps){
  mc_lv <- rmvnorm(1* lv.mc, mean = rep(0, n_lv))
  mc_lv <- array(c(mc_lv), dim = c(lv.mc, 1, n_lv))
  mc_lv1 <- rmvnorm(forcast_len* lv.mc, mean = rep(0, n_lv))
  mc_lv1 <- array(c(mc_lv1), dim = c(lv.mc, forcast_len, n_lv))
  linpred<- array(dim=c(n_sp,1,lv.mc))
  iexpeta<- matrix(ppes$Draws$expeta[i,],ncol=n_sp,nrow=n_obs)
  mpred<- iexpeta/rowSums(iexpeta)*rowSums(y)
  mod_rmse[i]<- RMSE(as.vector(mpred),as.vector(as.matrix(y)))
  mod_cor[i]<- cor(as.vector(mpred),as.vector(as.matrix(y)))^2
  mod_mae[i]<-MAE(as.vector(mpred),as.vector(as.matrix(y)))
  ieta<- log(iexpeta)
  A<- sweep(matrix(ppes$Draws$alpha[i,],ncol=n_sp,nrow=n_sp),1,ppes$Draws$k[i,],"/")
  
  for (b in 1:lv.mc) {
    
    random_test<- matrix(mc_lv[b,,],ncol=1) %*% t(matrix(ppes$Draws$loadings[i,],ncol=1))
    
    if(exists("X")){
      X_test<-model.matrix(~ps_data$sdata$Treatment)[lo,-1]
      X_coef<- X_test %*% matrix(ppes$Draws$coef[i,],ncol=n_sp,nrow=nrow(ppes$PPE$coef))
      eta_test<- X_coef+ random_test
    }else{
      eta_test<- random_test
    }
    
    est<-  ieta[35, ]+ppes$Draws$R_diag[i,] *
      (1 - A %*% ieta[35, ]) + 
      t(eta_test)[,1] #
    linpred[,,b]<-  est
  }
  
  prederror<- apply(linpred,c(1,2),mean)
  predmore1<- predmore<- vector("list",length=(forcast_len))
  predmore[[1]]<- prederror
  
  for(k in 2:length(predmore)){
    predmore1[[k]]<- array(dim=c(1,35,lv.mc))
  }
  
  
  for(k in 2:(forcast_len)){
    for(b in 1:lv.mc){
      
      random <-  matrix(mc_lv1[b,,],ncol=1) %*% t(matrix(ppes$Draws$loadings[i,],ncol=1))
      
      predmore1[[k]][,,b]<-  predmore[[(k-1)]] + ppes$Draws$R_diag[i,]* (1-A %*% predmore[[(k-1)]]) +
        random[(k-1),]
      
    }
    predmore[[k]]<- t(apply(predmore1[[k]],c(1,2),mean))
    
  }
  allpreds<- t(do.call(cbind,predmore))
  pred<- exp(allpreds)/rowSums(exp(allpreds))*rowSums(test)
  
  
  expest<- pred
  preds[,,i]<- expest
  pred_rmse[i]<- RMSE(as.vector(expest),as.vector(as.matrix(test)))
  pred_mae[i]<- MAE(as.vector(expest),as.vector(as.matrix(test)))
  pred_cor[i]<- cor(as.vector(expest),as.vector(as.matrix(test)))^2
  rm(mc_lv)
}

onecor<- onermse<-onemae<- vector("list",length=30000)

for( i in 1:30000){
  onecor[i]<- cor(as.vector(preds[1,,i]),as.vector(as.matrix(test[1,])))^2
  onemae[i]<- MAE(as.vector(preds[1,,i]),as.vector(as.matrix(test[1,])))
  onermse[i]<- RMSE(as.vector(preds[1,,i]),as.vector(as.matrix(test[1,])))
}
dfp<- data.frame(RMSE= pred_rmse,MAE=pred_mae,R2=pred_cor,Type=rep("Out of fit all",length(pred_cor)))
dfm<- data.frame(RMSE= mod_rmse,MAE=mod_mae,R2=mod_cor,Type=rep("In fit",length(mod_cor)))
dfo<- data.frame(RMSE= unlist(onermse),MAE=unlist(onemae),R2=unlist(onecor),Type=rep("Out of fit one",length(mod_cor)))

dfs<- rbind(dfp,dfm,dfo)

gg_mae<-    ggplot(data=dfs,aes(x=Type,y=MAE))+
  geom_violin(size=1)+
  stat_summary(fun.y=median, geom="point",shape=20, size=2, color="black")+
  scale_y_log10()+
  theme_light()+labs(y="MAE",x="")
gg_rmse<-    ggplot(data=dfs,aes(x=Type,y=RMSE))+
  geom_violin(size=1)+
  stat_summary(fun.y=median, geom="point",shape=20, size=2, color="black")+
  scale_y_log10()+
  theme_light()+labs(y="RMSE",x="")
gg_cor<-    ggplot(data=dfs,aes(x=Type,y=R2))+
  geom_violin(size=1)+
  stat_summary(fun.y=mean, geom="point",shape=20, size=2, color="black")+
  theme_light()+labs(y="R2",x="")
error_grid<- (arrangeGrob(grobs=list(gg_mae,gg_cor,gg_rmse),nrow=1))

ggsave(error_grid, file=paste0("P:/Surrogate/ecological modelling/Plots/",
                               "Gompertz_Error_plot_HS.png"),dpi=300)
