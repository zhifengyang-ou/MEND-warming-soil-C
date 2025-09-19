
##treatment effects by linear mixed model
env=read.table("environmental variables.txt",header=T,sep="\t",row.names = 1)
env[is.na(env)]<-0

library(lme4)
library(lmerTest)
library(glmulti)
library(MuMIn)
cand.models<-list()
cand.models[[1]]<-lmer(scale(env$variable)~warm*precip*clip+(1|block)+(1|year),data=env)
cand.models[[2]]<-lmer(scale(env$variable)~warm*precip*clip*scale(env$year)+(1|block),data=env)
cand.models[[3]]<-lmer(scale(env$variable)~warm*precip*clip+scale(env$year)+(1|block),data=env)
lm.ave<-model.avg(cand.models[[1]],cand.models[[2]],cand.models[[3]])
summary(lm.ave)

##paired t test
env.warm=read.table("environmental variables under warming.txt",header=T,sep="\t",row.names = 1)
env.unwarming=read.table("environmental variables under unwarming.txt",header=T,sep="\t",row.names = 1)

t.test(env.warm$variable, env.unwarming$variable, paired= TRUE, var.equal = TRUE)

