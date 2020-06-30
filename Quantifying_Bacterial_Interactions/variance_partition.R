
library(phyloseq)
ps_data<- readRDS(file=paste0("P:/Surrogate/ecological modelling/","contam_model_data_w_intial.RDS"))



y<- data.frame(ps_data$otus)[-1,]
seqDepth<- rowSums(y)



contrasts(ps_data$sample_data$Treatment)=contr.sum
X<- model.matrix(~ps_data$sample_data$Treatment)

library(gllvm)



gllvm1<- gllvm(y= y,X=data.frame(X), num.lv=2,family="negative.binomial",method="VA",row.eff = TRUE,n.init=10)
gllvm3<- gllvm(y= y,X=data.frame(X), num.lv=3,family="negative.binomial",method="VA",row.eff = TRUE,n.init=10)
gllvm5<- gllvm(y= y,X=data.frame(X), num.lv=5,family="negative.binomial",method="VA",row.eff = TRUE,n.init=10)

# 3 IS best based on AIC

library(boral)
mcmcs= list(n.burnin=50000, n.iteration=200000,n.thin=10,seed=123)

priors = list(type = c("normal","normal","cauchy","uniform"),hypparams = c(10, 10, 6.25, 30),ssvs.index = -1)
brl1<- boral(y, X=X[,-1],
             family="negative.binomial",lv.control=list(num.lv=3) ,row.eff="fixed",
             row.ids=matrix(1:nrow(y),ncol=1),control = priors, save.model=TRUE)
saveRDS(brl1,file=paste0("P:/Surrogate/ecological modelling/results/","contam_CR_brl_3lv.RDS"))


vars<- calc.varpart(brl1)

env<- data.frame(type=rep("Environmental",length(vars$varpart.X)),Var= vars$varpart.X)

rand<- data.frame(type=rep("Latent",length(vars$varpart.X)),Var= vars$varpart.lv)

all<- rbind(env,rand)

ggplot(data=all)+geom_density(aes(Var,color=type))

ggplot(data=all)+geom_density(aes(Var,fill=type),alpha=0.3)+theme_light()+labs(x="Variance Explained")
ggsave(file=paste0("P:/Surrogate/ecological modelling/results/","contam_brl_var.png"),dpi=300)
