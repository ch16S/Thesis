
```{r}
library(phyloseq)
library(parallel)
library(smooth)
library(pracma)
library(ggplot2)
library(gridExtra)

pathP<- "P:/Surrogate/Network_Comparison/Data/"
pathS<- "P:/Surrogate/Universality/"
pathPlots<- "P:/Surrogate/Network_Comparison/Plots/"
source(paste0(pathS,"doc_calc15_06.R"))
source(paste0(pathS,"DOC_bootstrap04_07.R"))
source(paste0(pathS,"detect_negative_slope_Lowess_04.R"))
source(paste0(pathS,"make_null_com11_06.R"))



mainra <- readRDS(file=paste0(pathP,"main_raps_uni.RDS"))
#grapenew<- readRDS(paste0(pathP,"new_ps_grape.RDS"))
#rootra<-transform_sample_counts(grapenew$roots,function(x)x/sum(x))
#mainra$grape_grape_grape_roots<- rootra
#saveRDS(mainra,file=paste0(pathP,"main_raps_uni.RDS"))

```

```{r}

mainra2<- lapply(mainra,function(physeq){
  if(phyloseq::taxa_are_rows(physeq)==TRUE){
    physeq<- otu_table(physeq)
    physeq<- t(physeq)
    m1<- as.matrix(physeq)
  }else{
    physeq<- otu_table(physeq)
    m1<-as.matrix(physeq)
  }
  return(m1)
}
 )



```

```{r}

overlap1_5<- lapply(mainra2[1:5],function(x){print("done");func_Cal_Overlap_rJSD_from_relative_abundance(x)})
overlap6_8<- lapply(mainra2[6:8],function(x){print("done");func_Cal_Overlap_rJSD_from_relative_abundance(x)})
overlap11_15<- lapply(mainra2[11:15],function(x){print("done");func_Cal_Overlap_rJSD_from_relative_abundance(x)})
overlap16_20<- lapply(mainra2[16:20],function(x){print("done");func_Cal_Overlap_rJSD_from_relative_abundance(x)})
overlap9_10<- lapply(mainra2[9:10],function(x){print("done");func_Cal_Overlap_rJSD_from_relative_abundance(x)})
gprtdoc<- func_Cal_Overlap_rJSD_from_relative_abundance(graprt)

saveRDS(overlap1_5,file=paste0(pathP,"overlap1-5.RDS"))
saveRDS(overlap6_8,file=paste0(pathP,"overlap6_8.RDS"))
saveRDS(overlap11_15,file=paste0(pathP,"overlap11-15.RDS"))
saveRDS(overlap16_20,file=paste0(pathP,"overlap16-20.RDS"))
saveRDS(overlap9_10,file=paste0(pathP,"overlap9_10.RDS"))


overlap1_5<- readRDS(file=paste0(pathP,"overlap1-5.RDS"))
overlap6_8<- readRDS(file=paste0(pathP,"overlap6_8.RDS"))
overlap11_15<- readRDS(file=paste0(pathP,"overlap11-15.RDS"))
overlap16_20<- readRDS(file=paste0(pathP,"overlap16-20.RDS"))
overlap9_10<- readRDS(file=paste0(pathP,"overlap9_10.RDS"))
overlaps<- c(overlap1_5,overlap6_8 ,overlap9_10,overlap11_15,overlap16_20)

saveRDS(overlaps,file=paste0(pathP,"all_real_overlaps.RDS"))

#bootstrap1_5<- lapply(overlaps1_5,function(x){print("done"); func_cal_rlowess_bootstrap(x)})
#bootstraps6_8<- lapply(overlaps6_8,function(x){print("done"); func_cal_rlowess_bootstrap(x)})
#bootstraps11_15<- lapply(overlap11_15,function(x){print("done"); func_cal_rlowess_bootstrap(x)})
#bootstraps16_20<- lapply(overlap16_20,function(x){print("done"); func_cal_rlowess_bootstrap(x)})
#bootstraps9_10<- lapply(overlap9_10,function(x){print("done"); func_cal_rlowess_bootstrap(x)})

 all_bs2<- lapply(overlaps,function(x){print("done");func_cal_rlowess_bootstrapv2(x)}) ## Modified the bootstrap function 
saveRDS(all_bs2,file=paste0(pathP,"all_bs2.RDS"))
saveRDS(all_bs3,file=paste0(pathP,"all_bs3.RDS"))


saveRDS(bootstraps1_5,file=paste0(pathP,"bootstrap1_5.RDS"))
saveRDS(bootstraps6_8,file=paste0(pathP,"bootstrap6_8.RDS"))
saveRDS(bootstraps11_15,file=paste0(pathP,"bootstrap11_15.RDS"))
saveRDS(bootstraps15_20,file=paste0(pathP,"bootstrap16_20.RDS"))
saveRDS(bootstraps9_10,file=paste0(pathP,"bootstrap9_10.RDS"))


bs1_5<- readRDS(file=paste0(pathP,"bootstrap1_5.RDS"))
bs6_8<- readRDS(file=paste0(pathP,"bootstrap6_8.RDS"))
bs11_15<-readRDS(file=paste0(pathP,"bootstrap11_15.RDS"))
bs16_20<- readRDS(file=paste0(pathP,"bootstrap16_20.RDS"))
bs9_10<- readRDS(file=paste0(pathP,"bootstrap9_10.RDS"))

all_bs<- c(bs1_5,bs6_8,bs9_10, bs11_15,bs16_20)

mainnull<- readRDS(file=paste0(pathP,"null_main.RDS"))
null_bs<- readRDS(file=paste0(pathP,"nullmain_bs.RDS"))
nulloverlap<-  readRDS(file=paste0(pathP,"null_main_doc.RDS"))

saveRDS(mainnull,file=paste0(pathP,"null_main.RDS"))
saveRDS(null_bs,file=paste0(pathP,"nullmain_bs.RDS"))
saveRDS(nulloverlap,file=paste0(pathP,"null_main_doc.RDS"))
theones<- c(1:3,5:20)
all_bs3<- vector("list",length=20)
for(i in seq_along(theones)){
  x<- theones[i]
  newbs<-  list(all_bs2[[x]][[1]], all_bs2[[x]][[2]],all_bs[[x]]$neg_slopes)
  names(newbs)<- names(all_bs[[x]])
  all_bs3[[x]]<- newbs
}


  

```

I fixed the bootstraps and added the grape root data to bs and docs. 
```{r}

all_bs3 <- readRDS(file=paste0(pathP,"all_bs_final.RDS"))
all_docs<- readRDS(file=paste0(pathP,"all_real_overlaps.RDS"))
null_bs<-  readRDS(file=paste0(pathP,"nullmain_bs.RDS"))

Fns<- Ys<- xcs<- vector("list",length(all_bs3))
for(i in 1:length(all_bs3)){
 print(i)
ys<- myrlowess(cbind(as.vector(as.dist(all_docs[[i]]$Overlap)),
                     as.vector(as.dist(all_docs[[i]]$jsd))),
               all_bs3[[i]]$xs,span=0.3) # Using span 0.3 because some was too noisy at one point 
xc <- detect_negative_slope_Lowess_03(all_bs3[[i]]$xs, ys)
xcs[[i]]<- xc
Ys[[i]]<- ys
Fns[[i]] <- sum(as.vector(as.dist(all_docs[[i]]$Overlap)) > xc)/length(as.vector(as.dist(all_docs[[i]]$Overlap)))
}



# all_bs3 1 has neg slopes combine
neg_ps<- lapply(all_bs3,function(x)sum(x$neg_slopes>0)/200)
universality<- data.frame(pval=neg_ps,Fns=unlist(Fns),slope_coef= unlist(lapply(all_bs3,function(x)mean(x$neg_slopes))))

Study<- c(rep("Zarraonaindia",5),
          rep("Wagner",3),
          rep("Fitzpatrick",5),
          rep("Edwards",4),
             rep("Santos-Medellín",3))

Compartment<- c("Grape","Leaves","Rhizosphere","Roots","Soil","Leaves","Roots","Soil",
                "Root Endosphere","Endosphere Toothpick ","Rhizosphere","Rhizosphere Toothpick ", "Soil",
                "Soil","Root Endosphere","Rhizoplane","Rhizosphere",
                "Soil","Rhizosphere","Root Endosphere")

###

sc<- paste0(Compartment,"-",Study)
cs<- paste0(Study,"-",Compartment)


slopes<- lapply(all_bs,function(x) x$neg_slopes)
 names(slopes)<- cs
slopesld<- plyr::ldply(slopes,as.data.frame)

ggneg<- ggplot()+geom_boxplot(data=slopesld,aes(x=slopesld$.id,y=slopesld$`X[[i]]`),size=1.5)+coord_flip()+geom_hline(yintercept = 0)+theme_light()+labs(y="Bootstrap Slope Coefficient",x="")

names(Fns)<- cs
fns2<- plyr::ldply(Fns,as.data.frame)
ggfns<- ggplot()+geom_col(data=fns2, aes(x=.id,y=`X[[i]]`,fill=Study))+coord_flip()+geom_hline(yintercept = c(0,1))+theme_light()+labs(y="Fns",x="")+theme(axis.text.y = element_blank())

gmatrix<- rbind(c(1,1,2,2))

ggsave(arrangeGrob(ggneg,ggfns,layout_matrix  =gmatrix),file=paste0(pathPlots,"slope_nfs.png"),dpi=300)

alldocs_df<- lapply(1:20,function(x){
   dfreal<- data.frame(mean= apply(all_bs3[[x]]$`Y-Bootstrap`,2,median),
   lower=apply(all_bs3[[x]]$`Y-Bootstrap`,2,function(x)quantile(x,c(0.03))),
   upper= apply(all_bs3[[x]]$`Y-Bootstrap`,2,function(x)quantile(x,c(0.97))),
   xs=all_bs3[[x]]$xs,
   type=rep("Real",ncol(all_bs3[[x]]$`Y-Bootstrap`)))
     dfnull<- data.frame(mean= apply(null_bs[[x]]$`Y-Bootstrap`,2,median),
   lower=apply(null_bs[[x]]$`Y-Bootstrap`,2,function(x)quantile(x,c(0.03))),
   upper= apply(null_bs[[x]]$`Y-Bootstrap`,2,function(x)quantile(x,c(0.97))),
   xs=null_bs[[x]]$xs,
   type=rep("Null",ncol(null_bs[[x]]$`Y-Bootstrap`)))
     dfall<- rbind(dfreal,dfnull)
   return(dfall)}
   )

names(alldocs_df)<- sc
alldocunlist<- plyr::ldply(alldocs_df,data.frame)

lapply(alldocs_df,function(x){
  ggplot(data=alldocunlist)+
    geom_line(data=x,aes(x=xs,y=mean,color=type),size=1.1)+
    geom_ribbon(data=x,aes(ymin=lower,ymax=upper,fill=type),alpha=0.1)+
    theme_light()+facet_wrap(~.id)
}

ggdoc<- ggplot(data=alldocunlist)+
  geom_line(data=alldocunlist,aes(x=xs,y=mean,color=type),size=1.2)+
        geom_ribbon(data=alldocunlist ,aes(x=xs,ymin=lower,ymax=upper,fill=type),alpha=0.4)+
        theme_light()+
  facet_wrap(~.id,scale= "free") +
  labs(x="Overlap",y="rJSD")

ggsave(ggdoc,file=paste0(pathPlots,"All_DOCs_wNulls.png"),dpi=400)

### 

library(kableExtra)

tab<- read.csv(paste0(pathP,"author_info.csv"),header=TRUE,row.names = NULL)

 colnames(tab)[1]<- "Study"
colnames(tab)[3]<- "Interesting Aspects of the Study"

tab %>%  kable(align = "l")%>% 
    kable_styling(full_width = T,bootstrap_options = "bordered")%>% 
  column_spec(c(1:3),color="black", bold = T )%>%
  as_image(file="Study-info.png",width=0.1) 
  

```

# For Seqdepths 
```{r}
saveRDS(seqdepths,file=paste0(pathP,"uni_seqdepths.RDS"))

seqdepths<-readRDS(file=paste0(pathP,"uni_seqdepths.RDS"))


seqdfs<-lapply(1:20,function(x){tryCatch({
    seq_mat<-abs(outer(seqdepths[[x]],seqdepths[[x]],"-"))
    seq_vec<- as.vector(as.dist(seq_mat))
    ov<- as.vector(as.dist(overlaps[[x]]$Overlap))
    j<- as.vector(as.dist(overlaps[[x]]$jsd))
    df1<- data.frame(jsd=j,overlap=ov,seq_depth=log(seq_vec+1))
    return(df1)},error=function(e)NA) })


gg_ovseqs<- lapply(seqdfs,function(x){
  ggplot()+geom_smooth(data=x,aes(x=seq_depth,y=overlap), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+
    theme(legend.title = element_blank())+labs(x="Change in Sequncing depth",y="Overlap")
})

seq_sum<- lapply(seqdfs,function(x){
 summary(x$seq_depth)
 qs<- quantile(x$seq_depth,c(0.1,0.5,0.9) )
 
 return(qs)
})

ovmedian<- lapply(seqdfs,function(x){
 median(x$overlap)
})

gg_ovseqs<- lapply(1:20,function(x){
  ggplot()+geom_smooth(data=seqdfs[[x]],aes(x=seq_depth,y=overlap), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+
    geom_vline(aes_(xintercept=seq_sum[[x]][c(1,2,3)]))+
        geom_hline(aes_(yintercept=ovmedian[[x]][1],color="red"),size=1)+
    theme(legend.position = "none")+labs(x="Change in Sequncing depth",y="Overlap")
})

 grid.arrange(arrangeGrob(grobs=gg_ovseqs[1:5],ncol=3))
 grid.arrange(arrangeGrob(grobs=gg_ovseqs[c(9,11,13)],ncol=3))
grid.arrange(arrangeGrob(grobs=gg_ovseqs[14:17],ncol=2))
grid.arrange(arrangeGrob(grobs=gg_ovseqs[18:20],ncol=3))
grid.arrange(arrangeGrob(grobs=gg_ovseqs[6:8],ncol=3))

ggsave(arrangeGrob(grobs=gg_ovseqs[1:5],ncol=3) ,file=paste0(pathPlots,"grape_doc_seq_.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_ovseqs[c(9,11,13)],ncol=3),file=paste0(pathPlots,"pnas_doc_seq_.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_ovseqs[14:17],ncol=2) ,file=paste0(pathPlots,"ricet_doc_seq_.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_ovseqs[18:20],ncol=3) ,file=paste0(pathPlots,"riced_doc_seq_.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_ovseqs[6:8],ncol=3),file=paste0(pathPlots,"bras_doc_seq_.png"),dpi=300)

upper<- lapply(seqdfs,function(x){
  which(x$seq_depth > quantile(x$seq_depth,0.5) )
})

lower<- lapply(seqdfs,function(x){
  which(x$seq_depth < quantile(x$seq_depth,0.5) )
})

gg_upper<- lapply(1:20,function(x){
  ggplot()+geom_smooth(data=seqdfs[[x]][,],aes(x=overlap,y=jsd,color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=seqdfs[[x]][-upper[[x]],],aes(x=overlap,y=jsd,color="Upper"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=seqdfs[[x]][-lower[[x]],],aes(x=overlap,y=jsd,color="Lower"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+
    theme()+labs(x="Overlap",y="rJSD")
})




# Lower = red, orig= green, upper= blue
ggsave(arrangeGrob(grobs=gg_upper[1:5],ncol=3) ,file=paste0(pathPlots,"grape_doc_seqrmv.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_upper[c(9,11,13)],ncol=3),file=paste0(pathPlots,"pnas_doc_seqrmv.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_upper[14:17],ncol=2) ,file=paste0(pathPlots,"ricet_doc_seqrmv.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_upper[18:20],ncol=3) ,file=paste0(pathPlots,"riced_doc_seqrmv.png"),dpi=300)
ggsave(arrangeGrob(grobs=gg_upper[6:8],ncol=3),file=paste0(pathPlots,"bras_doc_seqrmv.png"),dpi=300)


# To remove fabaceae 
fab<- sample_data(mainra$pnas_pnas_pnas_E)$Family %in%c("Fabaceae",""))
seq_mat<-abs(outer(seqdepths[[x]][-fab],seqdepths[[x]][-fab],"-"))
seq_vec<- as.vector(as.dist(seq_mat))
ov<- as.vector(as.dist(overlaps$pnas_pnas_pnas_E$Overlap[-fab,-fab]))
j<- as.vector(as.dist(overlaps$pnas_pnas_pnas_E$jsd[-fab,-fab]))
df1<- data.frame(overlap=ov,jsd=j ,seq_depth=log(seq_vec+1))
uppernofab<- which(df1$seq_depth > quantile(df1$seq_depth,0.5)) 
lowernofab<- which(df1$seq_depth < quantile(df1$seq_depth,0.5)) 

  ggplot()+geom_smooth(data=df1,aes(x=overlap,y=jsd,color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=df1[- uppernofab,],aes(x=overlap,y=jsd,color="Upper"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=df1[-  lowernofab,],aes(x=overlap,y=jsd,color="Lower"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+
    theme()+labs(x="Overlap",y="rJSD")

 ggsave(file=paste0(pathP,"nofab_seq_depth.png"),dpi=300)



```


```{r}
brasRdata<- sample_data(mainra$bras_bras_bras_root)

site_deets<- read.csv(paste0(pathP,"soildata_bras.csv"),header=TRUE)
site_deets2<- scale(site_deets[,-1])
pca1<- princomp(site_deets2)
pcadf<- data.frame(pca1$scores[,1:2],Site= site_deets$Site)
ggplot(data=pcadf,aes(x=Comp.1,Comp.2,color=Site))+geom_point()
test<-HCPC(pcadf[,c(1,2)])

table(test$data.clust$clust,site_deets$Site)

axis1<- tapply(pcadf$Comp.1,pcadf$Site,mean)
axis2<- tapply(pcadf$Comp.2,pcadf$Site,mean)

field<- which(brasRdata$Treatment=="field")
green<- which(brasRdata$Treatment!="field")

fieldsdata<- subset(brasRdata,brasRdata$Treatment=="field")
greensdata<- subset(brasRdata,brasRdata$Treatment!="field")


coords<- read.csv(paste0(pathP,"sites_bras.csv"),header=TRUE,row.names=1)

sdata1<- fieldsdata
soilaxis1<- vector(length=nrow(sdata1))
soilaxis2<-  vector(length=nrow(sdata1))
for(i in 1:nrow(sdata1)){
  one<-   match(sdata1$Site[i],names(axis1))
    two<-   match(sdata1$Site[i],names(axis2))

soilaxis1[i]<- axis1[one]
soilaxis2[i]<- axis2[two]
   }
ages<- list(which(is.na(brasRdata$Age)),
            which(brasRdata$Age==2),
            which(brasRdata$Age==3),
which(brasRdata$Age==4))

thecohort<- ifelse(is.na(brasRdata$Cohort)==TRUE,"3",brasRdata$Cohort)
brasRdata$SiteN<- as.numeric(brasRdata$Site)

brasRdata$SiteN<- as.numeric(brasRdata$Site)
brasRdata$GenotypeN<- as.numeric(brasRdata$Genotype)
sg<- as.numeric(as.factor(soil_group))

field_site<-outer(brasRdata$SiteN[field],brasRdata$SiteN[field],"-") 
field_Geno<-outer(brasRdata$GenotypeN[field],brasRdata$GenotypeN[field],"-") 
field_age<-outer(brasRdata$Age[field],brasRdata$Age[field],"-") 
field_lat<- abs(outer(lat,lat,"-"))
field_long<- abs(outer(long,long,"-"))
field_sg<- outer(sg,sg,"-")
field_axis1<- abs(outer(soilaxis1,soilaxis1,"-"))
field_axis2<- abs(outer(soilaxis2,soilaxis2,"-"))


brasf_rootdf<- data.frame(overlap=as.vector(as.dist(overlaps$bras_bras_bras_root$Overlap[field,field])),
                         jsd=as.vector(as.dist(overlaps$bras_bras_bras_root$jsd[field,field])),
                         site=as.vector(as.dist(field_site)),
                         geno=as.vector(as.dist(field_Geno)),
                         age=abs(as.vector(as.dist(field_age))),
                         soilaxis1=as.vector(as.dist(field_axis1)),
                         soilaxis2=as.vector(as.dist(field_axis2)))#,
                             #    lat=as.vector(as.dist(field_lat)),
                        # long=as.vector(as.dist(field_long)),
                         #  soilgroup=as.vector(as.dist(field_sg)))



fsitermv<- which(brasf_rootdf$site==0)
fgenormv<- which(brasf_rootdf$geno==0)
fsitegenormv<- unique(c(fsitermv,fgenormv))
fagermv<- which(brasf_rootdf$age==0)
fsiteagermv<- unique(c(fsitermv,fagermv))
fgenoagermv<- unique(c(fgenormv,fagermv))
rmlist<-list(fsitermv,fgenormv,fagermv,fgenoagermv,fsitegenormv,fsiteagermv)
fsitermv<- which(brasf_rootdf$soilaxis1 < 2)


g1<- ggplot(data=brasf_rootdf,aes(x=soilaxis1,y=overlap))+geom_point(alpha=0.1)+geom_smooth(size=1.5)+theme_light()
g2<- ggplot(data=brasf_rootdf,aes(x=soilgroup2,y=overlap))+geom_point(alpha=0.1)+geom_smooth(size=1.5)+theme_light()

g3<- ggplot(data=brasf_rootdf,aes(x=long,y=overlap))+geom_point(alpha=0.1)+geom_smooth(size=1.5)+theme_light()
g4<- ggplot(data=brasf_rootdf,aes(x=lat,y=overlap))+geom_point(alpha=0.1)+geom_smooth(size=1.5)+theme_light()

ggsave(arrangeGrob(g1,g2,g3,g4,ncol=2),file=paste0(pathPlots,"bras_enviro_geo.png"),dpi=300)

brasagef<- brasf_rootdf
brasagef<- subset(brasagef,brasagef$age!="NA")
brasagef$age<-as.factor(brasagef$age)

age_ov<- ggplot()+geom_boxplot(size=1.5,data=brasagef,aes(x=age,y=overlap))+
  theme_light()+labs(y="Overlap",x="Age")
age_jsd<- ggplot()+geom_boxplot(size=1.5,data=brasagef,aes(x=age,y=jsd))+
  theme_light()+labs(y="rJSD",x="Age")

ggsave(arrangeGrob(age_ov,age_jsd,ncol=2),file=paste0(pathPlots,"bras-rt_fieldage.png"),dpi=300)

field_gs<-ggplot()+
  geom_smooth(data=brasf_rootdf,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_rootdf[-rmlist[[1]],],aes(x=overlap,y=jsd,fill="Site",color="Site"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_rootdf[-rmlist[[2]],],aes(x=overlap,y=jsd,fill="Genotype",color="Genotype"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_rootdf[-rmlist[[3]],],aes(x=overlap,y=jsd,fill="Age",color="Age"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_rootdf[-rmlist[[4]],],aes(x=overlap,y=jsd,fill="Age Genotype",color="Age Genotype"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_rootdf[-rmlist[[5]],],aes(x=overlap,y=jsd,fill="Site Genotype",color="Site Genotype"), 
                formula=y~s(x),method="gam", size=1.5)+
      geom_smooth(data=brasf_rootdf[-rmlist[[6]],],aes(x=overlap,y=jsd,fill="Age Site",color="Age Site"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")

fsoilrmv2<- which(brasf_rootdf$soilaxis2 < 1)
fsoilrmv1<- which(brasf_rootdf$soilaxis1 < 2)
ggplot()+
  geom_smooth(data=brasf_rootdf,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_rootdf[-c(fsoilrmv1,fsoilrmv2,fagermv,fgenormv),],aes(x=overlap,y=jsd,fill="Age Site",color="Age Site"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")


ggsave(field_gs,
       file=paste0(pathPlots,"bras_fd_rts_site_geno.png"),dpi=300)

ages<- list(which(is.na(brasRdata$Age)),
            which(brasRdata$Age==2),
            which(brasRdata$Age==3),
which(brasRdata$Age==4))





brasf_rootdagesf<- lapply(1:4,function(x){ data.frame(overlap=as.vector(as.dist(overlaps$bras_bras_bras_root$Overlap[ages[[x]],ages[[x]]])),
                   
})

 
ggage_root<- ggplot()+
  geom_smooth(data=brasf_rootdagesf[[1]],aes(x=overlap,y=jsd,fill="<1",color="<1"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_rootdagesf[[2]],aes(x=overlap,y=jsd,fill="2",color="2"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_rootdagesf[[3]],aes(x=overlap,y=jsd,fill="3",color="3"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_rootdagesf[[4]],aes(x=overlap,y=jsd,fill="4",color="4"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 



greensdata$SiteN<- as.numeric(greensdata$Treatment)
greensdata$GenotypeN<- as.numeric(greensdata$Genotype)


gsoil_group<- vector(length=nrow(greensdata))
for(i in 1:nrow(greensdata)){
  if(greensdata$Treatment[i] %in% c("JAMsoil","PARsoil")){
    gsoil_group[i]<-"one"
  }else if(greensdata$Treatment[i] %in% c("MAHsoil","SILsoil")){
     gsoil_group[i]<-"two"
  }
  else if(greensdata$Treatment[i] %in% c("MILsoil")){
    gsoil_group[i]<-"three" 
  }else{
    gsoil_group[i]<-"four" 
  }
}



gsg<- as.numeric(as.factor(gsoil_group))

green_site<-outer(greensdata$SiteN,greensdata$SiteN,"-") 
green_Geno<-outer(greensdata$GenotypeN,greensdata$GenotypeN,"-") 
green_sg<-outer(gsg,gsg,"-") 

brasg_rootdf<- data.frame(overlap=as.vector(as.dist(overlaps$bras_bras_bras_root$Overlap[green,green])),
                         jsd=as.vector(as.dist(overlaps$bras_bras_bras_root$jsd[green,green])),
                         site=as.vector(as.dist(green_site)),
                         geno=as.vector(as.dist(green_Geno)))#,
                        # soilgroup=as.vector(as.dist(green_sg)))
                          



gsitermv<- which(brasg_rootdf$site==0)
ggenormv<- which(brasg_rootdf$geno==0)
gsitegenormv<- unique(c(gsitermv,ggenormv))
grmlist<-list(gsitermv,ggenormv,gsitegenormv)

gh_gs<-ggplot()+
  geom_smooth(data=brasg_rootdf,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasg_rootdf[-grmlist[[1]],],aes(x=overlap,y=jsd,fill="Site",color="Site"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasg_rootdf[-grmlist[[2]],],aes(x=overlap,y=jsd,fill="Genotype",color="Genotype"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasg_rootdf[-grmlist[[3]],],aes(x=overlap,y=jsd,fill="Site Genotype",color="Site Genotype"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")


grid.arrange(field_gs,gh_gs,ncol=2 )

grid.arrange(brasg_root_orig,
             otherg_brasrootgs[[1]],
             otherg_brasrootgs[[2]],
             otherg_brasrootgs[[3]],ncol=2 )

ggsave(gh_gs,
       file=paste0(pathPlots,"bras_gh_rts_site_geno.png"),dpi=300)
######################################### Leaf ############################

brasLdata<- sample_data(mainra$bras_bras_bras_leaf)
lfield<- which(brasLdata$Treatment=="field")


sdata1<- brasLdata

llong<- vector(length=nrow(sdata1))
llat<- vector(length=nrow(sdata1))
for(i in 1:nrow(sdata1)){
  one<-   match(sdata1$Site[i],rownames(coords))
llong[i]<- coords[one,2]
llat[i]<- coords[one,1]
}


soil_group<- vector(length=nrow(lfieldsdata))



soil_group<- vector(length=nrow( brasLdata))
for(i in 1:nrow( brasLdata)){
    if(brasLdata$Site[i] %in% c("Jam","Par")){
        lsoil_group[i]<-"one"
    }else if( brasLdata$Site[i] %in% c("Mah","Sil")){
        lsoil_group[i]<-"two"
    }
    else{ lsoil_group[i]<-"three" }
}


lfieldsdata<- subset(brasLdata,brasLdata$Treatment=="field")

lfieldsdata$SiteN<- as.numeric(lfieldsdata$Site)
lfieldsdata$GenotypeN<- as.numeric(lfieldsdata$Genotype)
lsg<- as.numeric(as.factor(lsoil_group))

lfield_Geno<-outer(lfieldsdata$GenotypeN,lfieldsdata$GenotypeN,"-") 
lfield_Site<-outer(lfieldsdata$SiteN,lfieldsdata$SiteN,"-") 
lfield_age<-outer(lfieldsdata$Age,lfieldsdata$Age,"-") 
agesitegenoe<- c()


lfield_lat<- abs(outer(llat,llat,"-"))
lfield_long<- abs(outer(llong,llong,"-"))
lfield_sg<- outer(lsg,lsg,"-")
lfield_sg1<- abs(outer(soilgroup1,soilgroup1,"-"))
lfield_sg2<- abs(outer(soilgroup1,soilgroup2,"-"))

brasf_leafdf<- data.frame(overlap=as.vector(as.dist(overlaps$bras_bras_bras_leaf$Overlap)),
                         jsd=as.vector(as.dist(overlaps$bras_bras_bras_leaf$jsd)),
                         site=as.vector(as.dist(lfield_Site)),
                         geno=as.vector(as.dist(lfield_Geno)),
                         age=abs(as.vector(as.dist(lfield_age))))#,
                         #long=as.vector(as.dist(lfield_long)),
                         #lat=as.vector(as.dist(lfield_lat)),
                         #soilgroup=as.vector(as.dist(lfield_sg))))
       


lfsitermv<- which(brasf_leafdf$site==0)
lfgenormv<- which(brasf_leafdf$geno==0)
lfagermv<- which(brasf_leafdf$age==0)

lfsiteage<- unique(c(lfsitermv,lfagermv))
lfsitegenormv<- unique(c(lfsitermv,lfgenormv))
lfagegenormv<- unique(c(lfgenormv,lfagermv))
lrmlist<- list(lfsitermv,lfgenormv,lfagermv,lfsitegenormv,lfsiteage,lfagegenormv)

lf_gs<-ggplot()+
  geom_smooth(data=brasf_leafdf,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_leafdf[-lfsitermv,],aes(x=overlap,y=jsd,fill="Site",color="Site"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_leafdf[-lfgenormv,],aes(x=overlap,y=jsd,fill="Genotype",color="Genotype"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_leafdf[-lfagermv,],aes(x=overlap,y=jsd,fill="Age",color="Age"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_leafdf[-lfsitegenormv,],aes(x=overlap,y=jsd,fill="Site Genotype",color="Site Genotype"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_leafdf[-lfsiteage,],aes(x=overlap,y=jsd,fill="Age Site",color="Age Site"), 
                formula=y~s(x),method="gam", size=1.5)+
          geom_smooth(data=brasf_leafdf[-lfagegenormv,],aes(x=overlap,y=jsd,fill="Age Genotype",color="Age Genotyoe"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")

ggsave(lf_gs,file=paste0(pathP,"brs_lf_genoagesite.png"),dpi=300)

jg1<- ggplot(data=brasf_leafdf,aes(x=long,y=overlap))+geom_point()+geom_smooth(size=1.5)+theme_light()
jg2<- ggplot(data=brasf_leafdf,aes(x=lat,y=overlap))+geom_point()+geom_smooth(size=1.5)+theme_light()
jg3<- ggplot(data=brasf_leafdf,aes(x=soilgroup1,y=overlap))+geom_point()+geom_smooth(size=1.5)+theme_light()
jg4<- ggplot(data=brasf_leafdf,aes(x=soilgroup2,y=overlap))+geom_point()+geom_smooth(size=1.5)+theme_light()

ggsave(arrangeGrob(jg1,jg2,jg3,jg4,ncol=2),file=paste0(pathP,"brs_lfov_latlngsgs.png"),dpi=300)

g1<- ggplot(data=brasf_leafdf,aes(x=long,y=jsd))+geom_point()+geom_smooth(size=1.5)+theme_light()
g2<- ggplot(data=brasf_leafdf,aes(x=lat,y=jsd))+geom_point()+geom_smooth(size=1.5)+theme_light()
g3<- ggplot(data=brasf_leafdf,aes(x=soilgroup1,y=jsd))+geom_point()+geom_smooth(size=1.5)+theme_light()
g4<- ggplot(data=brasf_leafdf,aes(x=soilgroup2,y=jsd))+geom_point()+geom_smooth(size=1.5)+theme_light()

ggsave(arrangeGrob(g1,g2,g3,g4,ncol=2),file=paste0(pathP,"brs_lfjsd_latlngsgs.png"),dpi=300)




ggsave(arrangeGrob(field_gs,gh_gs,lf_gs,ncol=2 ),file=paste0(pathP,"brs_rtfl_rtgh_lffl.png"),dpi=300)


brasSdata<- sample_data(mainra$bras_bras_bras_soil)
fieldsdata<- brasSdata$Type

To fix still 

sdata1<- brasSdata

llong<- vector(length=nrow(brasSdata))
llat<- vector(length=nrow(brasSdata))
for(i in 1:nrow(sdata1)){
  one<-   match(sdata1$Site[i],rownames(coords))
llong[i]<- coords[one,2]
llat[i]<- coords[one,1]
}


soil_group<- vector(length=nrow(fieldsdata))

first<-  match(fieldsdata$Site, names(axis1) )
second<-  match(fieldsdata$Site, names(axis2) )

soilgroup1<- axis1[first]
soilgroup2<- axis2[second]


lsoil_group<- vector(length=nrow( brasSdata))
for(i in 1:nrow( brasSdata)){
    if(brasSdata$Site[i] %in% c("Jam","Par")){
        lsoil_group[i]<-"one"
    }else if( brasSdata$Site[i] %in% c("Mah","Sil")){
        lsoil_group[i]<-"two"
    }
    else{ lsoil_group[i]<-"three" }
}


fieldsdata<- brasSdata

fieldsdata$SiteN<- as.numeric(fieldsdata$Site)
fieldsdata$GenotypeN<- as.numeric(fieldsdata$Genotype)
fs<- as.numeric(as.factor(paste0(fieldsdata$Treatment,fieldsdata$Site)))
fs<- as.numeric(as.factor(paste0(fieldsdata$Treatment,fieldsdata$Site)))


field_Site<-outer(fieldsdata$SiteN,fieldsdata$SiteN,"-") 
field_treat<-outer(fieldsdata$Treatment,fieldsdata$Treatment,"-") 
field_lat<- abs(outer(llat,llat,"-"))
field_long<- abs(outer(llong,llong,"-"))
field_sg<- outer(lsg,lsg,"-")
field_sg<- outer(fs,fs,"-")

field_sg1<- abs(outer(soilgroup1,soilgroup1,"-"))
field_sg2<- abs(outer(soilgroup2,soilgroup2,"-"))

brasf_soildf<- data.frame(overlap=as.vector(as.dist(overlaps$bras_bras_bras_soil$Overlap)),
                         jsd=as.vector(as.dist(overlaps$bras_bras_bras_soil$jsd)),
                         site=as.vector(as.dist(field_Site)),
                         age=abs(as.vector(as.dist(field_age))),
                         long=as.vector(as.dist(field_long)),
                         lat=as.vector(as.dist(field_lat)),
                         soilgroup=as.vector(as.dist(field_sg)),
                         soilgroup1=as.vector(as.dist(field_sg1)),
                         soilsg=as.vector(as.dist(field_sg)))
                    
)

lfsitermv<- which(brasf_soildf$site==0)
lflat<- which(brasf_soildf$lat==0)
lflong<- which(brasf_soildf$long==0)
lfsg1rmv<- which(brasf_soildf$soilgroup1==0)
sltrmv<- which(brasf_soildf$soilsg==0)

lfagermv<- which(brasf_soildf$age==0)
lfgrmv<- which(brasf_soildf$soilgroup==0)
lfsiteage<- unique(c(lfsitermv,lfagermv))
lfsitegenormv<- unique(c(lfsitermv,lfgenormv))
lrmlist<- list(lfsitermv,lfgenormv,lfagermv,lfgrmv,lfsiteage)

lf_gs<-ggplot()+
  geom_smooth(data=brasf_soildf,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_soildf[-lfsitermv,],aes(x=overlap,y=jsd,fill="Site",color="Site"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_soildf[-lfagermv,],aes(x=overlap,y=jsd,fill="Age",color="Age"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_soildf[-lflat,],aes(x=overlap,y=jsd,fill="Age",color="Age"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_soildf[-lflong ],aes(x=overlap,y=jsd,fill="Age",color="Age"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")


g1<- ggplot(data=brasf_soildf,aes(x=long,y=jsd))+geom_point()+geom_smooth(size=1.5)+theme_light()
g2<- ggplot(data=brasf_soildf,aes(x=lat,y=jsd))+geom_point()+geom_smooth(size=1.5)+theme_light()
g3<- ggplot(data=brasf_soildf,aes(x=soilgroup1,y=jsd))+geom_point(method="lm")+geom_smooth(size=1.5)+theme_light()
g4<- ggplot(data=brasf_soildf,aes(x=soilgroup2,y=jsd))+geom_point()+geom_smooth(size=1.5)+theme_light()


agesL<- list(
            which(brasLdata$Age==2),
            which(brasLdata$Age==3),
which(brasLdata$Age==4))





brasf_leafdagesf<- lapply(1:4,function(x){ data.frame(overlap=as.vector(as.dist(overlaps$bras_bras_bras_leaf$Overlap[agesL[[x]],agesL[[x]]])),
               jsd=as.vector(as.dist(overlaps$bras_bras_bras_leaf$jsd[agesL[[x]],agesL[[x]]])))

                   
})

 
ggage_leaf<- ggplot()+
    geom_smooth(data=brasf_leafdagesf[[2]],aes(x=overlap,y=jsd,fill="2",color="2"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=brasf_leafdagesf[[3]],aes(x=overlap,y=jsd,fill="3",color="3"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=brasf_leafdagesf[[4]],aes(x=overlap,y=jsd,fill="4",color="4"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 

                                        ggsave(arrangeGrob(ggage_root,ggage_leaf,ncol=2),file=paste0(pathPlots,"bras_age_doc.png"),dpi=300)

```


```{r}

#PNAS
plantTraits<- read.csv(file=paste0(pathP,"Host_plant_phenotypic_distance.csv"),header=TRUE,row.names = 1)
pnaseID2<- matrix(rep(1:nrow(overlaps$pnas_pnas_pnas_E$Overlap),nrow(overlaps$pnas_pnas_pnas_E$Overlap)),ncol=nrow(overlaps$pnas_pnas_pnas_E$Overlap))
pnaseID1<- t(pnaseID2)



func_plant_dist<- function(x,y){
sdata1<- x
plant_dist<- matrix(ncol=nrow(sdata1),nrow=nrow(sdata1))
for(i in 1:nrow(sdata1)){
  for(j in 1:nrow(sdata1)){
  one<-   match(sdata1$species[i],rownames(plantTraits))
   two<-  match(sdata1$species[j],rownames(plantTraits))
   if(any(is.na(c(one,two)))==TRUE){
     plant_dist[i,j]<- NA
   }else{
   plant_dist[i,j]<- plantTraits[one,two]
   }
  }
}
distrmv<- which( is.na(as.vector(as.dist(plant_dist))))
pnasoverlapvec<- as.vector(as.dist(y$Overlap))[-distrmv]

plntdistvec<- as.vector(as.dist(plant_dist))[-distrmv]
pnasjsdvec<- as.vector(as.dist(y$jsd))[-distrmv]
df<- data.frame(plant_dist=plntdistvec,overlap=pnasoverlapvec,jsd=pnasjsdvec)
return(df)
}


func_plant_distwd<- function(x,y,drought){
sdata1<- x
plant_dist<- matrix(ncol=nrow(sdata1),nrow=nrow(sdata1))
for(i in 1:nrow(sdata1)){
  for(j in 1:nrow(sdata1)){
  one<-   match(sdata1$species[i],rownames(plantTraits))
   two<-  match(sdata1$species[j],rownames(plantTraits))
   if(any(is.na(c(one,two)))==TRUE){
     plant_dist[i,j]<- NA
   }else{
   plant_dist[i,j]<- plantTraits[one,two]
   }
  }
}
distrmv<- which( is.na(as.vector(as.dist(plant_dist))))
pnasoverlapvec<- as.vector(as.dist(y$Overlap))[-distrmv]
plntdistvec<- as.vector(as.dist(plant_dist))[-distrmv]
pnasjsdvec<- as.vector(as.dist(y$jsd))[-distrmv]
drghtvec<- as.vector(as.dist(drought))[-distrmv]

df<- data.frame(plant_dist=plntdistvec,overlap=pnasoverlapvec,jsd=pnasjsdvec,drought=drghtvec)
return(df)
}


#### Endo
allpss<- readRDS(file=paste0(pathP,"pnasE_counts_wFam.RDS"))
fab<- which(sample_data(allpss)$Family=="Fabaceae")
pnasdocnofab<- lapply(overlaps$pnas_pnas_pnas_E,function(x)x[-fab,-fab])
drought<- which(sample_data(mainra$pnas_pnas_pnas_E)$treatment=="d")
drought1<- as.numeric(sample_data(allpss)$treatment)
drought_m<- outer(drought1,drought1,"-")
fam1<- as.numeric(sample_data(allpss)$Family[-fab])
famdist<- outer(fam1,fam1,"-")
pnasdoconlyfab<- lapply(overlaps$pnas_pnas_pnas_E,function(x)x[fab,fab])


wfabd<- func_plant_distwd(sample_data(allpss),overlaps$pnas_pnas_pnas_E,drought=drought_m)
drought<- as.numeric(sample_data(mainra$pnas_pnas_pnas_E)$treatment[-fab])
droughtdist<- as.vector(as.dist(outer(drought,drought,"-")))
nofabd<- func_plant_distwd(sample_data(allpss)[-fab],pnasdocnofab,drought=drought_m[-fab,-fab])


onlyfabd<- func_plant_distwd(sample_data(allpss)[fab,],pnasdoconlyfab,drought=drought_m[fab,fab])

wfabdrought<- 
wfab<- func_plant_dist(sample_data(allpss),overlaps$pnas_pnas_pnas_E)
nofab<- func_plant_dist(sample_data(allpss)[-fab,],pnasdocnofab)
nofab<- func_plant_dist(sample_data(allpss)[-fab,],pnasdocnofab,famdist)

nfrmv<-which(nofab$plant_dist==0)
wfrmv<- which(wfab$plant_dist==0)
nffamrmv<-which(nofab$family==0)

endodrought<- func_plant_dist(sample_data(allpss)[drought,],
                               lapply(overlaps$pnas_pnas_pnas_E,function(x)x[drought,drought]) )
endowater<- func_plant_dist(sample_data(allpss)[-drought],
                              lapply(overlaps$pnas_pnas_pnas_E,function(x)x[-drought,-drought]))
nfdrmv<-which(nofab$plant_dist==0)


datanf<- sample_data(allpss)[-fab,]
nfdrought<- which(datanf$treatment=="d")

nfendodrought<- func_plant_dist(datanf[nfdrought,],
                              lapply(pnasdocnofab,function(x)x[nfdrought,nfdrought]))
nfendowater<- func_plant_dist(datanf[-nfdrought,],
                              lapply(pnasdocnofab,function(x)x[-nfdrought,-nfdrought]))




wtnf<- which(nfendowater$plant_dist==0)
dtnf<- which(nfendodrought$plant_dist==0)
drmv<- which(wfabd$drought==0)  
####### Compare similarity of dorught 
#ggplot()+
#   geom_smooth(data=wfabd[,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
#                formula=y~s(x),method="gam", size=1.5)+
#    geom_smooth(data=wfabd[-drmv,],aes(x=overlap,y=jsd,fill="Same host",color="Same host"), 
#                formula=y~s(x),method="gam", size=1.5)+
#    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 
#######
ggfab<- ggplot()+
   geom_smooth(data=wfab[,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=wfab[-wfrmv,],aes(x=overlap,y=jsd,fill="Same host",color="Same host"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 
 

 
 ggnofab<-ggplot()+
   geom_smooth(data=nofab[,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=nofab[-nfrmv,],aes(x=overlap,y=jsd,fill="Same host",color="Same host"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")
  ggnofab<-ggplot()+
   geom_smooth(data=nofab[,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=nofab[-nfrmv,],aes(x=overlap,y=jsd,fill="Same host",color="Same host"), 
                formula=y~s(x),method="gam", size=1.5)+
        #geom_smooth(data=nofab[-nffamrmv,],aes(x=overlap,y=jsd,fill="Same host",color="Same host"), 
                #formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")
 
## Removing the same host species makes no difference    

ggendodrought<- ggplot()+
   geom_smooth(data=wfab,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=endowater,aes(x=overlap,y=jsd,fill="Watered",color="Watered"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=endodrought,aes(x=overlap,y=jsd,fill="Drought",color="Drought"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 


ggnfendodrought<- ggplot()+
    geom_smooth(data=nofab,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=nfendowater,aes(x=overlap,y=jsd,fill="Watered",color="Watered"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=nfendodrought,aes(x=overlap,y=jsd,fill="Drought",color="Drought"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 

ggnofab_d<- ggplot()+
    geom_smooth(data=nofabd[,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=nofabd[-nfrmv,],aes(x=overlap,y=jsd,fill="Host",color="Host"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=nofabd[-drght1,],aes(x=overlap,y=jsd,fill="Drought",color="Drought"), 
                formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=nofabd[-c(drght1,nfrmv),],aes(x=overlap,y=jsd,fill="Host Drought",color="Host Drought"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")


ggsave(arrangeGrob(gwf,gwfrmv,gnf,gnfrmv,ncol=2),
       file=paste0(pathPlots,"pnas_re_fabornot_samehostornot.png"),dpi=300)

ggsave(arrangeGrob(ggfab,ggnofab,ggnofab_d,ggnfendodrought,ncol=2),
       file=paste0(pathPlots,"pnas_root_fbsh_nfsh_fd_fnd.png"),dpi=300)

fabspp<- which(onlyfabd$plant_dist==0)

ggonlyfab<- ggplot()+
    geom_smooth(data=onlyfabd[,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")

######### Rhizo 

drought<- which(sample_data(mainra$pnas_pnas_pnas_R)$treatment=="d")

fab<- which(sample_data(mainra$pnas_pnas_pnas_R)$Family=="Fabaceae")
rhznofab<- lapply(overlaps$pnas_R,function(x)x[-fab,-fab])


pnasrhizo<- func_plant_dist(sample_data(mainra$pnas_pnas_pnas_R),overlaps$pnas_R )

pnasnofab_rhizo<- func_plant_dist(sample_data(mainra$pnas_pnas_pnas_R)[-fab],rhznofab )

rhizodrought<- func_plant_dist(sample_data(mainra$pnas_pnas_pnas_R)[drought,],
                               lapply(overlaps$pnas_R,function(x)x[drought,drought]) )
rhizowater<- func_plant_dist(sample_data(mainra$pnas_pnas_pnas_R)[-drought],
                              lapply(overlaps$pnas_R,function(x)x[-drought,-drought]))


fbrmv<-which(pnasrhizo$plant_dist==0)
nfbrmv<-which(pnasnofab_rhizo$plant_dist==0)
wrhz<-which(rhizowater$plant_dist==0)
drhz<-which(rhizodrought$plant_dist==0)

 ggrhz<-ggplot()+
       geom_smooth(data=pnasrhizo,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=rhizodrought,aes(x=overlap,y=jsd,fill="Drought",color="Drought"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=rhizowater[,],aes(x=overlap,y=jsd,fill="Watered",color="Watered"), 
                formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=pnasrhizo[,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=rhizodrought[-drhz,],aes(x=overlap,y=jsd,fill="Drought Species",color="Drought Species"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=rhizowater[-wrhz,],aes(x=overlap,y=jsd,fill="Watered Species",color="Watered Species"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")

 
 
rhznofab<- ggplot()+
   geom_smooth(data=pnasrhizo,aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=pnasnofab_rhizo[,],aes(x=overlap,y=jsd,fill="No Fabaceae",color="No Fabaceae"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")


wfab
pnasrhizo

ggdistrhz<- ggplot(data=pnasrhizo)+geom_point(aes(x=plant_dist,y=overlap),alpha=0.1)+
   geom_smooth(aes(x=plant_dist,y=overlap,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
         theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Phenotypic Distance",y="Overlap")
ggdistendo<- ggplot(data=wfab)+geom_point(aes(x=plant_dist,y=overlap),alpha=0.1)+
   geom_smooth(aes(x=plant_dist,y=overlap,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
         theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Phenotypic Distance",y="Overlap")
ggdistrhzj<- ggplot(data=pnasrhizo)+geom_point(aes(x=plant_dist,y=jsd),alpha=0.1)+
   geom_smooth(aes(x=plant_dist,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
         theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Phenotypic Distance",y="rJSD")
ggdistendoj<- ggplot(data=wfab)+geom_point(aes(x=plant_dist,y=jsd),alpha=0.1)+
   geom_smooth(aes(x=plant_dist,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
         theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Phenotypic Distance",y="rJSD")



ggsave(arrangeGrob(ggdistendo,ggdistendoj,ggdistrhz,ggdistrhzj,ncol=2),
       file=paste0(pathPlots,"edov_edj_rhzov_rhzj.png"),dpi=300)

ggsave(arrangeGrob(ggrhz,rhznofab,ncol=2),
       file=paste0(pathPlots,"rhzdrought_rhznf.png"),dpi=300)

```

Rice temporal Data 

```{r}




rice_times<- overlaps[14:17]
rice_ps<- mainra[14:17]

rice_sdata<- lapply(rice_ps,function(x)sample_data(x))


agediff<- lapply(rice_sdata,function(x){ 
  abs(outer(x$Age,x$Age,"-"))})

age_dynamic<- lapply(rice_sdata,function(x){ 
  age_state<-ifelse(x$Age > 42,"Older","Younger")
  age_stateN<- as.numeric(as.factor(age_state))
  agestate_mat<- outer( age_stateN, age_stateN,"-")
  return(agestate_mat)
})

age_dyn<- lapply(rice_sdata,function(x){ 
  age_state<-ifelse(x$Age > 42,"Older","Younger")
  oldage<- which(age_state=="Older")
  return(oldage)
})



sitediff<-  lapply(rice_sdata,function(x){ 
  site<-as.numeric(x$Site)
  site_mat<- outer(site,site,"-")
  return(site_mat)})

cultdiff<-  lapply(rice_sdata,function(x){ 
  Cultivar<-as.numeric(x$Cultivar)
  Cultivar_mat<- outer(Cultivar,Cultivar,"-")
  return(site_mat)})

agevector<- lapply(agediff,function(x)as.vector(as.dist(x)))
overlapvector<- lapply(rice_times,function(x)as.vector(as.dist(x$Overlap)))
jsdvector<- lapply(rice_times,function(x)as.vector(as.dist(x$jsd)))
sitevector<- lapply(sitediff,function(x)as.vector(as.dist(x)))
agestvector<- lapply(age_dynamic,function(x)as.vector(as.dist(x)))
cultvector<- lapply(cultdiff,function(x)as.vector(as.dist(x)))

allcoms<- lapply(1:4,function(x){
  data.frame(age=agevector[[x]],jsd=jsdvector[[x]],overlap=overlapvector[[x]],
             site=sitevector[[x]],agestate=agestvector[[x]],
             cultivar=cultvector[[x]])
})
ages<- age_dyn[[2]]

endotimeold<-  data.frame(jsd=as.vector(as.dist(rice_times[[2]]$jsd[ages,ages])),
                          overlap=as.vector(as.dist(rice_times[[2]]$Overlap[ages,ages])),
                          site=as.vector(as.dist(sitediff[[2]][ages,ages])),
                           cultivar=as.vector(as.dist(cultdiff[[2]][ages,ages])))
endotimeyng<-  data.frame(jsd=as.vector(as.dist(rice_times[[2]]$jsd[-ages,-ages])),
                          overlap=as.vector(as.dist(rice_times[[2]]$Overlap[-ages,-ages])),
                          site=as.vector(as.dist(sitediff[[2]][-ages,-ages])),
                          cultivar=as.vector(as.dist(cultdiff[[2]][-ages,-ages])))

names(allcoms)<- c("sl","eo","rp","rs")

ageList<- lapply(allcoms,function(x) which(x$age==0))
siteList<- lapply(allcoms,function(x) which(x$site==0))
agesiteList<- lapply(1:4,function(x)unique(c(ageList[[x]],siteList[[x]])))

cultList<- lapply(allcoms,function(x) which(x$cultivar==0))
sitecultList<- lapply(1:4,function(x)unique(c(siteList[[x]],cultList[[x]])))

allgs<- lapply(c(1,2,4),function(x){
ggplot()+geom_smooth(data=allcoms[[x]][,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                         formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=allcoms[[x]][-ageList[[x]],],aes(x=overlap,y=jsd,fill="Age",color="Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=allcoms[[x]][-siteList[[x]],],aes(x=overlap,y=jsd,fill="Site",color="Site"), 
                         formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=allcoms[[x]][-agesiteList[[x]],],aes(x=overlap,y=jsd,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 
})



ggsave(arrangeGrob(allgs[[1]],allgs[[3]],allgs[[2]],ncol=2),
       file=paste0(pathPlots,"rice_sl_rhz_endo.png"),dpi=300)


cultsiteo<- unique(c(which(endotimeold$site==0),endotimeold$cultivar==0))
cultsitey<- unique(c(which(endotimeyng$site==0),endotimeyng$cultivar==0))

gg_time1<- ggplot()+geom_smooth(data=endotimeold,aes(x=overlap,y=jsd,fill="Older",color="Older"), 
                         formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=endotimeold[-which(endotimeold$site==0) ,]
              ,aes(x=overlap,y=jsd,fill="Older Site",color="Older Site"), 
                  formula=y~s(x),method="gam", size=1.5)+
    #geom_smooth(data=endotimeold[-cultsiteo ,],
        #      aes(x=overlap,y=jsd,fill="Older Site2",color="Older Site2"), 
         #         formula=y~s(x),method="gam", size=1.5)+
     theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 

gg_time2<- ggplot()+
    geom_smooth(data=endotimeyng[-which(endotimeyng$site==0) ,],
aes(x=overlap,y=jsd,fill="Younger Site",color="Younger Site"), 
                         formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=endotimeyng,aes(x=overlap,y=jsd,fill="Younger",color="Younger"), 
                         formula=y~s(x),method="gam", size=1.5)+geom_vline(xintercept=0.6,alpha=0.3,size=1)+
      #geom_smooth(data=endotimeyng[-cultsitey,],aes(x=overlap,y=jsd,fill="Younger2",color="Younger2"), 
                      #   formula=y~s(x),method="gam", size=1.5)+geom_vline(xintercept=0.6,alpha=0.3,size=1)+
       theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 

  
ggsave(arrangeGrob(gg_time1,gg_time2,ncol=2),file=paste0(pathP,"ricendo_time.png"),dpi=300)

ages<- age_dyn[[3]]

endotimeold2<-  data.frame(jsd=as.vector(as.dist(rice_times[[3]]$jsd[ages,ages])),
                          overlap=as.vector(as.dist(rice_times[[3]]$Overlap[ages,ages])),
                          cult=as.vector(as.dist(cultdiff[[3]][ages,ages])))
endotimeyng2<-  data.frame(jsd=as.vector(as.dist(rice_times[[3]]$jsd[-ages,-ages])),
                          overlap=as.vector(as.dist(rice_times[[3]]$Overlap[-ages,-ages])),
                          cult=as.vector(as.dist(cultdiff[[3]][-ages,-ages])))

gg_time12<- ggplot()+geom_smooth(data=endotimeold2,aes(x=overlap,y=jsd,fill="Older",color="Older"), 
                         formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=endotimeold2[-which(endotimeold2$cult==0) ,],
              aes(x=overlap,y=jsd,fill="Older Cultivar",color="Older Cultivar"), 
                  formula=y~s(x),method="gam", size=1.5)+
     theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 

gg_time22<- ggplot()+
    geom_smooth(data=endotimeyng2[-which(endotimeyng2$cult==0) ,],
                aes(x=overlap,y=jsd,fill="Younger Site",color="Younger Cultivar"), 
                         formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=endotimeyng2,aes(x=overlap,y=jsd,fill="Younger",color="Younger"), 
                         formula=y~s(x),method="gam", size=1.5)+
       theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD") 

ggsave(arrangeGrob(gg_time12,gg_time22,ncol=2),file=paste0(pathPlots,"ricerp_time_cult.png"),dpi=300)

so<- ggplot()+ geom_point(data=allcoms[[1]],aes(x=age,y=overlap),alpha=0.1)+
  geom_smooth(data=allcoms[[1]],aes(x=age,y=overlap,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="Overlap") 

sj<- ggplot()+ geom_point(data=allcoms[[1]],aes(x=age,y=jsd),alpha=0.1)+
  geom_smooth(data=allcoms[[1]],aes(x=age,y=jsd,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="rJSD")

eo<- ggplot()+geom_point(data=allcoms[[1]],aes(x=age,y=overlap),alpha=0.1)+
  geom_smooth(data=allcoms[[2]],aes(x=age,y=overlap,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="Overlap")

ej<- ggplot()+geom_point(data=allcoms[[1]],aes(x=age,y=jsd),alpha=0.1)+
geom_smooth(data=allcoms[[2]],aes(x=age,y=jsd,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="rJSD")

rpo<- ggplot()+  geom_point(data=allcoms[[1]],aes(x=age,y=overlap),alpha=0.1)+
  geom_smooth(data=allcoms[[3]],aes(x=age,y=overlap,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="Overlap") 

rpj<- ggplot()+geom_point(data=allcoms[[1]],aes(x=age,y=jsd),alpha=0.1)+
  geom_smooth(data=allcoms[[3]],aes(x=age,y=jsd,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="rJSD")
rso<- ggplot()+geom_point(data=allcoms[[1]],aes(x=age,y=overlap),alpha=0.1)+
  geom_smooth(data=allcoms[[4]],aes(x=age,y=overlap,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="Overlap") 
rsj<- ggplot()+geom_point(data=allcoms[[1]],aes(x=age,y=jsd),alpha=0.1)+
   geom_smooth(data=allcoms[[4]],aes(x=age,y=jsd,fill="Site Age",color="Site Age"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.position = "none")+labs(x="Change in Age",y="rJSD")




grid.arrange(eo,ej,rpo,rpj,rso,rsj,so,sj,ncol=2)

ggsave(arrangeGrob(eo,ej,rpo,rpj,rso,rsj,so,sj,ncol=2),file=paste0(pathPlots,"rice_doc_by_time.png"),dpi=300)



```

rice drought 
```{r}


 
rice_drght<- overlaps[18:20]
riced_ps<- mainra[18:20]

rice_sdata<- lapply(riced_ps,function(x)sample_data(x))


drought<- lapply(rice_sdata,function(x){ 
  drought<- as.numeric(x$Treatment)
  outer(drought,drought,"-")})
Soil<- lapply(rice_sdata,function(x){ 
  Soil<- as.numeric(x$Soil)
  outer(Soil,Soil,"-")})
cultivar<- lapply(rice_sdata,function(x){ 
 cultivar<- as.numeric(x$Cultivar)
  outer(cultivar,cultivar,"-")})
tub<- lapply(rice_sdata,function(x){ 
 Tub<- as.numeric(x$Tub)
  outer(Tub,Tub,"-")})

whichdrought<- lapply(rice_sdata,function(x) which(x$Treatment=="DS"))


overlapdrought<- lapply(1:3,function(x)
  as.vector(as.dist(rice_drght[[x]]$Overlap[-whichdrought[[x]],-whichdrought[[x]]])))
jsddrought<- lapply(1:3,function(x)
  as.vector(as.dist(rice_drght[[x]]$jsd[-whichdrought[[x]],-whichdrought[[x]]])))
soildrought<- lapply(1:3,function(x)
  as.vector(as.dist(Soil[[x]][-whichdrought[[x]],-whichdrought[[x]]])))
cultivardrought<- lapply(1:3,function(x)
  as.vector(as.dist(cultivar[[x]][-whichdrought[[x]],-whichdrought[[x]]])))

overlapnodrought<- lapply(1:3,function(x)as.vector(as.dist(rice_drght[[x]]$Overlap[whichdrought[[x]],whichdrought[[x]]])))
jsdnodrought<- lapply(1:3,function(x)as.vector(as.dist(rice_drght[[x]]$jsd[whichdrought[[x]],whichdrought[[x]]])))
soilnodrought<- lapply(1:3,function(x)as.vector(as.dist(Soil[[x]][whichdrought[[x]],whichdrought[[x]]])))
cultivarnodrought<- lapply(1:3,function(x)as.vector(as.dist(cultivar[[x]][whichdrought[[x]],whichdrought[[x]]])))

nd_df<- lapply(1:3, function(x){ data.frame(jsd=jsdnodrought[[x]],
                                            overlap=overlapnodrought[[x]],
                                            cultivar=cultivarnodrought[[x]],
                                            soil=soilnodrought[[x]])})
  
d_df<- lapply(1:3, function(x){ data.frame(jsd=jsddrought[[x]],
                                           overlap=overlapdrought[[x]],
                                           cultivar=cultivardrought[[x]],
                                           soil=soildrought[[x]])})

ds<- lapply(1:3,function(x)which(d_df[[x]]$soil==0))
ws<- lapply(1:3,function(x)which(d_df[[x]]$soil==0))

drought_gs<- lapply(1:3,function(x){
ggplot()+geom_smooth(data=d_df[[x]][,],aes(x=overlap,y=jsd,fill="Drought",color="Drought"), 
                         formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=nd_df[[x]][,],aes(x=overlap,y=jsd,fill="Watered",color="Watered"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")
})

ggsave(arrangeGrob(drought_gs[[1]],drought_gs[[2]],drought_gs[[3]],ncol=1),file=paste0(pathPlots,"sl_rhz_en_droughtvsnorm_rice.png"),dpi=300)

drought_gs<- lapply(1:3,function(x){
ggplot()+geom_smooth(data=d_df[[x]][,],aes(x=overlap,y=jsd,fill="Drought",color="Drought"), 
                         formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=nd_df[[x]][,],aes(x=overlap,y=jsd,fill="Watered",color="Watered"), 
                         formula=y~s(x),method="gam", size=1.5)+
   geom_smooth(data=d_df[[x]][-ds[[x]],],aes(x=overlap,y=jsd,fill="Drought Site",color="Drought Site"), 
                         formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=nd_df[[x]][-ws[[x]],],aes(x=overlap,y=jsd,fill="Watered Site",color="Watered Site"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")
})


wtfovers<- nd_df[[3]][-ws[[3]],]
dtfovers<- d_df[[3]][-ds[[3]],]

m1<- which(wtfovers$overlap> median(wtfovers$overlap))
m2<- which(dtfovers$overlap> median(dtfovers$overlap))

summary(lm(wtfovers$jsd[ m1]~wtfovers$jsd[m1] ))

droughtvector<- lapply(drought,function(x)as.vector(as.dist(x)))
overlapvector<- lapply(rice_drght,function(x)as.vector(as.dist(x$Overlap)))
jsdvector<- lapply(rice_drght,function(x)as.vector(as.dist(x$jsd)))
Soilvector<- lapply(Soil,function(x)as.vector(as.dist(x)))
cultivarvector<- lapply(cultivar,function(x)as.vector(as.dist(x)))


allcoms<- lapply(1:3,function(x){
  data.frame(drought=droughtvector[[x]],jsd=jsdvector[[x]],overlap=overlapvector[[x]],
             soil=Soilvector[[x]],cultivar=cultivarvector[[x]])
})

soilrmv<- lapply(allcoms,function(x)which(x$soil==0))
cultrmv<- lapply(allcoms,function(x)which(x$cultivar==0))
drghtrmv<- lapply(allcoms,function(x)which(x$drought==0))
soildrght<- lapply(1:3,function(x) c(soilrmv[[x]],drghtrmv[[x]]))
tubrmv<- lapply(allcoms,function(x)which(x$tub==0))

ggplot()+geom_smooth(data=allcoms[[3]][-soilrmv[[3]],],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                         formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")

all_gs<- lapply(1:3,function(x){
  ggplot()+
    geom_smooth(data=allcoms[[x]][,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                         formula=y~s(x),method="gam", size=1.5)+
  geom_smooth(data=allcoms[[x]][-drghtrmv[[x]],],aes(x=overlap,y=jsd,fill="Drought",color="Drought"), 
                    formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=allcoms[[x]][-soilrmv[[x]],],aes(x=overlap,y=jsd,fill="Soil",color="Soil"), 
                    formula=y~s(x),method="gam", size=1.5)+
         # geom_smooth(data=allcoms[[x]][-soildrght[[x]],],aes(x=overlap,y=jsd,fill="Soil Drought",color="Soil Drought"), 
                  #formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")})

ggsave(arrangeGrob(all_gs[[1]],all_gs[[2]],all_gs[[3]]),file=patse0(pathPlots,"sl_rhz,en_drought_rice.png"),dpi=300)
ggsave(arrangeGrob(all_gs[[1]],drought_gs[[1]],
                   all_gs[[2]],drought_gs[[2]],
                   all_gs[[3]],drought_gs[[3]] ),
       file=paste0(pathPlots,"sl_rhz_en_all_drought_rice.png"),dpi=300)




```

Grape
```{r}


 ## Need to add year 
grape<- overlaps[2:5]
grape_ps<- mainra[2:5]

grape_date<- readRDS(file=paste0(pathP,"grape_date.RDS"))

grape_sdata<- lapply(grape_date[2:5],function(x)sample_data(x))
 area<- lapply(grape_sdata,function(x){ 
  area<- as.numeric(x$area)
  outer( area, area,"-")})

 area2<- lapply(grape_sdata,function(x){ 
  area2<- as.numeric(x$area_2)
  outer( area2, area2,"-")})

dates<-  lapply(grape_sdata,function(x){ 
  Date<- as.numeric(x$Date)
  outer(Date, Date,"-")})
  outer( area2, area2,"-")})

years<-  lapply(grape_sdata,function(x){ 
  year<- sapply(strsplit(as.character(x$Date),"/"),"[",3)
  years2<- as.numeric(as.factor(year))
  outer(years2,years2,"-")

})
stage<-  lapply(grape_sdata,function(x){ 
 month<- sapply(strsplit(as.character(x$Date),"/"),"[",2)
   month2<- as.numeric(as.factor( month))
  return(outer( month2, month2,"-"))

})
        
   
theyear<-  lapply(grape_sdata,function(x){ 
  year<- sapply(strsplit(as.character(x$Date),"/"),"[",3) })                    
 
 areavector<- lapply(area,function(x)as.vector(as.dist(x)))
overlapvector<- lapply(grape,function(x)as.vector(as.dist(x$Overlap)))
jsdvector<- lapply(grape,function(x)as.vector(as.dist(x$jsd)))
area2vector<- lapply(area2,function(x)as.vector(as.dist(x)))
datevector<- lapply(dates,function(x)as.vector(as.dist(x)))
stagevector<- lapply(stage,function(x)as.vector(as.dist(x)))
yearvector<- lapply(years,function(x)as.vector(as.dist(x)))

allcoms<- lapply(1:4,function(x){
  data.frame(area=areavector[[x]],jsd=jsdvector[[x]],overlap=overlapvector[[x]],
             area2=area2vector[[x]] ,date=datevector[[x]],stage=stagevector[[x]],
             year=yearvector[[x]]
                                                                            
  )
})

arearmv<- lapply(allcoms,function(x)which(x$area==0))
area2rmv<- lapply(allcoms,function(x)which(x$area2==0))
yearrmv<- lapply(allcoms,function(x)which(x$year==0))
datermv<- lapply(allcoms,function(x)which(x$date==0))
stagermv<- lapply(allcoms,function(x)which(x$stage==0))
yeararea<- lapply(1:4,function(x) c(yearrmv[[x]],arearmv[[x]]))
stagearea<- c


all_gs<- lapply(1:4,function(x){
ggplot()+
    geom_smooth(data=allcoms[[x]][,],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                formula=y~s(x),method="gam", size=1.5)+
    geom_smooth(data=allcoms[[x]][-stagermv[[x]],],aes(x=overlap,y=jsd,fill="Stage",color="Stage"), 
                formula=y~s(x),method="gam", size=1.5)+
      geom_smooth(data=allcoms[[x]][-yearrmv[[x]],],aes(x=overlap,y=jsd,fill="Year",color="Year"), 
                formula=y~s(x),method="gam", size=1.5)+
        geom_smooth(data=allcoms[[x]][-arearmv[[x]],],aes(x=overlap,y=jsd,fill="Area",color="Area"), 
                formula=y~s(x),method="gam", size=1.5)+
        geom_smooth(data=allcoms[[x]][-area2rmv[[x]],],aes(x=overlap,y=jsd,fill="Area 2",color="Area 2"), 
                formula=y~s(x),method="gam", size=1.5)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")})

ggsave(arrangeGrob(grobs=all_gs[1:4],ncol=2),file=paste(pathPlots,"grape_lfrhzrtsl_com_area.png"),dpi=300)

grdata<- grape_sdata$roots

month<- sapply(strsplit(as.character(grdata$Date ),"/"),"[",2)

 ## Need to add year 
gr_ovs<- grape$grape_grape_grape_roots

preharvest<- which(month== c("09","10"))
dormancy<- which(month == "04")
flowering<- which(month== "06")

tokeep<- list(preharvest,dormancy,flowering)

rootsdf<- lapply(tokeep,function(x){
ov<- overlapvector<- as.vector(as.dist(grape[[3]]$Overlap[x,x]))
j<- jsdvector<- as.vector(as.dist(grape[[3]]$jsd[x,x]))
area2vector<- as.vector(as.dist(area2[[3]][x,x] ))
 areavector<- as.vector(as.dist(area[[3]][x,x]))
 yearvector<- as.vector()
 return(data.frame(overlap=ov,jsd=j,area=areavector,area2=area2vector))
})



 areasrmv<- lapply(rootsdf,function(x) which(x$area==0 ))
  areas2rmv<- lapply(rootsdf,function(x) which(x$area2==0 ))

 
 ggrootyear<- lapply(1:3,function(x){
   ggplot()+
    geom_smooth(data=rootsdf[[x]],aes(x=overlap,y=jsd,fill="Original",color="Original"), 
                         formula=y~s(x),method="gam", size=1.5,alpha=0.2)+
       geom_smooth(data=rootsdf[[x]][-areasrmv[[x]],],aes(x=overlap,y=jsd,fill="area",color="area"), 
                    formula=y~s(x),method="gam", size=1.5,alpha=0.2)+
     geom_smooth(data=rootsdf[[x]][-areas2rmv[[x]],],aes(x=overlap,y=jsd,fill="area2",color="area2"), 
                    formula=y~s(x),method="gam", size=1.5,alpha=0.2)+
    theme_light()+guides(fill=FALSE)+theme(legend.title = element_blank())+labs(x="Overlap",y="rJSD")
 })
 ggsave(arrangeGrob(grobs=ggrootyear[1:3],ncol=2),file=paste(pathPlots,"grape_rt_life_stg.png"),dpi=300)


 ggsave(arrangeGrob(ggrootyear,ggrootyeararea,ncol=2),
        file=paste0(pathPlots,"grpe_rt_year_yeararea.png"),dpi=300)
 
 rootrich<- data.frame(richness= apply(otu_table(mainra$grape_grape_grape_roots),1,function(x)sum(x>0)),
                       year=grdata$Year)
 
 rootrichdf
#### For soil 
test3<- test2[,c(7,8,11,17,19,20,21)]

testdist<- lapply(1:7,function(x){
  as.vector(as.dist(abs(outer(test3[[x]],test3[[x]],"-"))))
})

names(testdist)<- colnames(test3)
soil_df <- data.frame(matrix(unlist(testdist), ncol=length((testdist), byrow=F))
colnames(soil_df) <- colnames(test3)  
dist1<- geosphere::distm(test2[,c("longitude","latitude")])

soil_df2<- data.frame(soil_df,overlap=as.vector(as.dist(overlap1_5$grape_grape_grape_soil$Overlap[-removal,-removal])),
                      jsd=as.vector(as.dist(overlap1_5$grape_grape_grape_soil$jsd[-removal,-removal])),
                      geo=as.vector(abs(as.dist(dist1))))

so_cn<- ggplot(data=soil_df2,aes(x=c_n_ratio,y=overlap))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                      
sj_cn<- ggplot(data=soil_df2,aes(x=c_n_ratio,y=jsd))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                      
so_ph<- ggplot(data=soil_df2,aes(x=ph,y=overlap))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                      
sj_ph<- ggplot(data=soil_df2,aes(x=ph,y=jsd))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                      
so_mst<- ggplot(data=soil_df2,aes(x=perc_moisture,y=overlap))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                      
sj_mst<- ggplot(data=soil_df2,aes(x=perc_moisture,y=jsd))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                     
so_tmp<- ggplot(data=soil_df2,aes(x=temperature,y=overlap))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                     
sj_tmp<- ggplot(data=soil_df2,aes(x=temperature,y=jsd))+geom_point(alpha=0.3)+geom_smooth()+theme_light() 
so_geo<- ggplot(data=soil_df2,aes(x=geolog,y=overlap))+geom_point(alpha=0.3)+geom_smooth()+theme_light()                     
sj_geo<- ggplot(data=soil_df2,aes(x=geolog,y=jsd))+geom_point(alpha=0.3)+geom_smooth()+theme_light()


```


make Null communities
```{r}

mainnull<- lapply(mainra2,func_make_null)
nulloverlap<- lapply(mainnull,function(x){print("done");func_Cal_Overlap_rJSD_from_relative_abundance(x)})
saveRDS(mainnull,file=paste0(pathP,"null_main.RDS"))
saveRDS(nulloverlap,file=paste0(pathP,"null_main_doc.RDS"))
nullmain_bs<- lapply(nulloverlap,function(x){print("done");func_cal_rlowess_bootstrap(x)})
saveRDS(nullmain_bs,file=paste0(pathP,"nullmain_bs.RDS"))


```



```{r}

theMeans<- lapply(all_bs,function(x)apply(x[[1]],2,mean))


Prctile_BS<- quantile(bootStraps[[1]], c(.3,  .97)) 
 
mean_BS <- apply(bootstraps1_5[[1]],2,function(x)mean(x))
Prctile_BS <- apply(bootStraps[[1]],2,function(x)quantile(x, c(.03,  .97)))

theMeans<- lapply(doc_bs,function(x)apply(x[[1]],2,mean))


null_docs_bs<- nulls
theNullMeans<- lapply(null_docs_bs,function(x)apply(x[[1]],2,mean))

      


all_null_slopes<- vector("list",length(null_docs_bs))

for(i in 1:length(all_null_slopes)){
  print(names(all_null_slopes[i]))
  nullslopes<- NULL
  for (k in 1:200){
    tryCatch({
lm1<- lm(null_docs_bs[[i]]$`Y-Bootstrap`[k,null_docs_bs[[i]]$xs > xcs[[i]]]~
     poly(null_docs_bs[[i]]$xs[null_docs_bs[[i]]$xs> xcs[[i]] ] , 1))
nullslopes<- c(nullslopes,lm1$coefficients[2])
  },error=function(e)NA) }
  all_null_slopes[[i]]<- nullslopes
  
} 

set.seed(123)
all_slopes2<- vector("list",length(all_bs))
theboots<- sample(1:1000,200)
for(i in 1:length(all_slopes2)){
  slopes<- NULL
  for (j in 1:200){
    k<- theboots[j]
    tomodel<- which(all_bs[[i]]$xs > median(overlaps[[i]]$Overlap))
    lm1<- lm(all_bs[[i]]$`Y-Bootstrap`[k,tomodel ]~
               all_bs[[i]]$xs[tomodel] )
    slopes<- c(slopes,lm1$coefficients[2])
  }
  all_slopes2[[i]]<- slopes
  
} 

means<- lapply(all_slopes2,function(x)mean(x))

pvals<- lapply(all_slopes2,function(x)sum(x>0)/200)
plot(all_bs[[i]]$`Y-Bootstrap`[k,tomodel ],
               all_bs[[i]]$xs[tomodel])
tomodel<- which(all_bs[[i]]$xs > median(overlaps[[i]]$Overlap)

tomodel<- which(all_bs[[i]]$`Y-Bootstrap`[k > median(overlaps[[i]]$Overlap)])

for(i in 1:length(all_slopes)){
  slopes<- NULL
  for (k in 1:200){
    lm1<- lm(all_bs[[i]]$`Y-Bootstrap`[k,all_bs[[i]]$xs > xcs[[i]] ]~
               poly(all_bs[[i]]$xs[all_bs[[i]]$xs> xcs[[i]] ] , 1))
    slopes<- c(slopes,lm1$coefficients[2])
  }
  all_slopes[[i]]<- slopes
  
} 
testres<- vector("list",length=length(all_null_slopes))
for(i in 1:length(all_null_slopes)){
  ttest1 <- t.test(all_slopes[[i]],all_null_slopes[[i]],"less")
  testres[[i]]<- ttest1
}

t.tes

pvals<- lapply(testres,function(x){
  x1<- x$p.value
  return(x1)
})


```


```{r}



ho_cmds<- lapply(1:20,function(x){
  print(x)
ov<- overlaps[[x]]$Overlap
diag(ov)<- 0
seventy5<- quantile(as.vector(ov),0.9)
if(seventy5 >= 0.95){
  retain<- which(apply(ov,1,function(x)sum(x  >= 0.95))> 0)
}else{
  
retain<- which(apply(ov,1,function(x)sum(x  >= seventy5))> 0)}
if(length(retain) < 4){ df1<- NA
}else{cmd1<- cmdscale(as.dist(overlaps[[x]]$Overlap[retain,retain]))
df1<- data.frame(cmd1,sample_data(mainra[[x]])[retain,])}
return(df1)})
  
ho_cmds2<- lapply(ho_cmds,)

cmds<- lapply(overlaps,function(x){cmdscale(as.dist(x$jsd))})
ords<- lapply(1:20,function(x){
  sdata<- sample_data(mainra[[x]])
  df<- data.frame(cmds[[x]],sdata)
  return(df)
})

# thevals<- as.list(c(rep("area_2",5),rep("Site",3),
#                      rep("species",5),rep("Age",2),"Age","Age",rep("Site",3)))


thevals<- as.list(c(rep("area_2",5),rep("Site",3),
                    rep("treatment",5),"Site","Age","Age",
                     "Site",rep("Treatment",3)))


ggs<- lapply(1:20,function(x){
  ggplot(data=ords[[x]],aes_string(x="X1",y="X2",color=thevals[[x]]))+
    geom_point(size=2,alpha=0.7)+
    theme_light()+labs(x="PCoA 1",y="PCoA 2")+
    theme(legend.position="none",axis.text.y=element_text(size=8),                          axis.text.x=element_text(size=8))+
    theme(title = element_text(size=6))+ggtitle(sc[x])
})


ho_ggs<- lapply(1:20,function(x){
  ggplot(data=ho_cmds[[x]],aes_string(x="X1",y="X2",color=thevals[[x]]))+
    geom_point(size=2,alpha=0.7)+
      theme_light()+labs(x="PCoA 1",y="PCoA 2")+
      theme(legend.position="none",axis.text.y=element_text(size=8),                          axis.text.x=element_text(size=8))+
      ggtitle(sc[x])+theme(title = element_text(size=6)) 
})

hlay<- rbind(c(1,2,3,4,5),
      c(6,7,8,NA,NA),
      c(9,10,11,12,13),
      c(14,15,16,17,NA),
      c(18,19,20,NA,NA))



ggsave(arrangeGrob(grobs=ggs[1:20],layout_matrix =hlay),
       file=paste0(pathPlots,"PCoA_jsd_by_compartment.png"),dpi=300)

ggsave(arrangeGrob(grobs=ho_ggs[1:20],layout_matrix=hlay),
       file=paste0(pathPlots,"PCoA_high_overlap_by_compartment.png"),dpi=300)

ords<- lapply(mainra,function(x){plot_ordination(x,                                                 ordinate(x,"MDS","jsd"))+
    geom_point(alpha=0.4)+
    theme_light()+
    xlab("")+ylab("")})

thevals<-c(rep("area",5),rep("Site",3),rep("species",5),rep("Site",5),rep("Soil",5))

theords_coloured<- lapply(c(1:20),function(x){
  p1<- plot_ordination(mainra[[x]],ordinate(mainra[[x]],"MDS","jsd"),color=thevals[x])+
    theme_light()+
    xlab("")+ylab("")
  return(p1)})

theords2<- lapply(theord_coloured,function(x){
  p1<- x+theme(legend.position="none")
  return(p1)
})

grapes<- lapply(c(1:5),function(x){
  sample_data(mainra[[x]])$thearea <- as.factor(paste0(sample_data(mainra[[x]])$area,
                                                       sample_data(mainra[[x]])$area_2))
  return(mainra[[x]])})

mainra[1:5]<- grapes[1:5]

```
Split age in 2 for 16

age > 
For rice > 42 split 

## Growing season -> 
##


```{r}
ag1<-list(ords[[5]],ords[[4]],ords[[3]],ords[[2]],
             ords[[8]],ords[[7]],ords[[8]],
             ords[[13]],ords[[11]],ords[[9]],
             ords[[14]],ords[[17]],ords[[16]],ords[[15]],
             ords[[18]],ords[[19]],ords[[20]])

hlay<- rbind(c(1,2,NA,3,4),
      c(5,NA,NA,6,7),
      c(8,9,NA,10,NA),
      c(11,12,13,14,NA),
      c(15,16,NA,17,NA))

hlay<- rbind(c(1,2,NA,3,4,5),
      c(6,NA,NA,7,NA,8),
      c(9,10,NA,11,NA,NA),
      c(12,13,14,15,NA,NA),
      c(16,17,NA,18,NA,NA))

colTitles<- c("Soil","Rhizosphere","Rhizoplane","Root","Leaf")
rowTitles<- c("Zarraonaindia et al 2015","Wagner et al 2016",
              "Fitzpatrick et al 2018", "Edwards et al 2018" ,"Santos-Medellin et al 2017")


col.titles = paste("C_Title", 1:3)
row.titles = paste("R_Title", 4:6)

# Add row titles
ag1[c(1,4,8,11,15)] = lapply(c(1,4,8,11,15), function(i) arrangeGrob(ag1[[i]], left=rowTitles[i]))

# Add column titles and lay out plots
grid.arrange(grobs=lapply(c(1,4,7), function(i) {
  arrangeGrob(grobs=pl[i:(i+2)], top=col.titles[i/3 + 1], ncol=1)
}), ncol=3)

grid.arrange(ords[[5]],ords[[4]],ords[[3]],ords[[2]],
             ords[[8]],ords[[7]],ords[[8]],
             ords[[13]],ords[[11]],ords[[9]],
             ords[[14]],ords[[17]],ords[[16]],ords[[15]],
             ords[[18]],ords[[19]],ords[[20]], layout_matrix=hlay)
grid.arrange(ag1,,)


lapply(c(1:20),function(x){
  fpc::pamk(as.dist())
})

```
Split time b


```{r}



otherEffects<- readRDS(file=paste0(pathP,"other_ps_uni.RDS"))

othera<- lapply(otherEffects[[1]],function(x)transform_sample_counts(x,function(y)y/sum(y)))
othera2<- lapply(otherEffects[[2]],function(x)transform_sample_counts(x,function(y)y/sum(y)))
othera3<- lapply(otherEffects[[3]],function(x)transform_sample_counts(x,function(y)y/sum(y)))


other2<- lapply(othera,function(physeq){
  if(phyloseq::taxa_are_rows(physeq)==TRUE){
    physeq<- otu_table(physeq)
    physeq<- t(physeq)
    m1<- as.matrix(physeq)
  }else{
    physeq<- otu_table(physeq)
    m1<-as.matrix(physeq)
  }
  return(m1)
}
 )


other2_2<- lapply(othera2,function(physeq){
  if(phyloseq::taxa_are_rows(physeq)==TRUE){
    physeq<- otu_table(physeq)
    physeq<- t(physeq)
    m1<- as.matrix(physeq)
  }else{
    physeq<- otu_table(physeq)
    m1<-as.matrix(physeq)
  }
  return(m1)
}
 )


other3_2<- lapply(othera3,function(physeq){
  if(phyloseq::taxa_are_rows(physeq)==TRUE){
    physeq<- otu_table(physeq)
    physeq<- t(physeq)
    m1<- as.matrix(physeq)
  }else{
    physeq<- otu_table(physeq)
    m1<-as.matrix(physeq)
  }
  return(m1)
}
 )

other_doc<- lapply(other2,function(x)func_Cal_Overlap_rJSD_from_relative_abundance(x))
other_null<- lapply(other2,function(x)func_make_null(x))
other_null_doc<- lapply(other_null,function(x)func_Cal_Overlap_rJSD_from_relative_abundance(x))
saveRDS(other_doc,paste0(pathP,"other_doc.RDS"))
saveRDS(other_null_doc,paste0(pathP,"other_null_doc.RDS"))

other_doc2<- lapply(other2_2,function(x)func_Cal_Overlap_rJSD_from_relative_abundance(x))
other_null2<- lapply(other2_2,function(x)func_make_null(x))
other_null_doc2<- lapply(other_null2,function(x)func_Cal_Overlap_rJSD_from_relative_abundance(x))
saveRDS(other_doc2,paste0(pathP,"other_doc2.RDS"))
saveRDS(other_null_doc2,paste0(pathP,"other_null_doc2.RDS"))

other_doc3<- lapply(other3_2,function(x)func_Cal_Overlap_rJSD_from_relative_abundance(x))
other_null3<- lapply(other3_2,function(x)func_make_null(x))
other_null_doc3<- lapply(other_null3,function(x)func_Cal_Overlap_rJSD_from_relative_abundance(x))
saveRDS(other_doc3,paste0(pathP,"other_null_doc3.RDS"))
saveRDS(other_null_doc3,paste0(pathP,"other_null_doc3.RDS"))

bs_other_null_doc3<- lapply(other_null_doc3,function(x)func_cal_rlowess_bootstrap(x))
bs_other_doc3<- lapply(other_doc3,function(x)func_cal_rlowess_bootstrap(x))
bs_other_null_doc2<- lapply(other_null_doc2,function(x)func_cal_rlowess_bootstrap(x))
bs_other_doc2<- lapply(other_doc2,function(x)func_cal_rlowess_bootstrap(x))
bs_other_null_doc<- lapply(other_null_doc,function(x)func_cal_rlowess_bootstrap(x))
bs_other_doc<- lapply(other_doc,function(x)func_cal_rlowess_bootstrap(x))

saveRDS(bs_other_null_doc,paste0(pathP,"bs_other_null_doc.RDS"))
saveRDS(bs_other_doc,paste0(pathP,"bs_other_doc.RDS"))

saveRDS(bs_other_doc2,paste0(pathP,"bs_other_doc2.RDS"))
saveRDS(bs_other_null_doc2,paste0(pathP,"bs_other_null_doc2.RDS"))

saveRDS(bs_other_doc3,paste0(pathP,"bs_other_doc3.RDS"))
saveRDS(bs_other_null_doc3,paste0(pathP,"bs_other_null_doc3.RDS"))

```


```{r}


sitera<- lapply(sitesEffects$bras,function(x) transform_sample_counts(x,function(u)u/sum(u)))

sites<- read.csv(paste0(pathP,"sites_bras.csv"),header=TRUE)
chems<- read.csv(paste0(pathP,"soildata_bras.csv"),header=TRUE)

chemlist<- split(chems,chems$Site)
chemlist2<- lapply(chemlist,function(x) apply(x,2,function(y)mean(y)))


test<- vector(length=nsamples(mainra$bras_bras_root))

for(i in 1:nsamples(mainra$bras_bras_root)){
if(bdata$Treatment[i] =="field"){
    bdata$Treatment[i]<- bdata$Site[i]
}
  
  else{
    test[i]<- sapply(strsplit(sample_data(mainra$bras_bras_root)$Treatment,"soil"),"[",1)
}
}



 sapply(strsplit(ntk2$.id,"___"),"[",2)
```

```{r}

 Overlap_Vec=nan(NumSamples*(NumSamples-1)/2,1);
    RootJSD_Vec=nan(NumSamples*(NumSamples-1)/2,1);
    Subj1_Vec=nan(NumSamples*(NumSamples-1)/2,1);
    Subj2_Vec=nan(NumSamples*(NumSamples-1)/2,1);

Subj1=repmat((1:ncol(x)),1,ncol(x))
Subj2=repmat(1:NumSamples,NumSamples,1)

X<- mainra2$grape_grape_grape_grapes
mx = dim(X)[1]
nx = dim(X)[2]
matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
n<- 1
sbj1<-sbj2<- matrix(ncol=NumSamples,nrow=NumSamples)
for (i in 1:(NumSamples)){
    for (j in 1:NumSamples){
      sbj1[i,j]<- i
       sbj2[i,j]<- j
    }}
      

```

