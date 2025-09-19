library(vegan)
library(picante)
library(lme4)
library(lmerTest)
library(MuMIn)

env=read.table("environmental variable.txt",header=T,sep="\t",row.names = 1)
env[is.na(env)]<-0
fm1=lmer(scale(env$Cvariable)~scale(env$variable)+(1|block)+(1|year),data=env)
summary(fm1)
presult<-car::Anova(fm1,type=2)
coefs<-coef(summary(fm1))[ , "Estimate"]
pvalue=presult[,3]
r2<-r.squaredGLMM(fm1)
r=ifelse(coefs>0,(r2[1,1])^0.5,-(r2[1,1])^0.5)
result=list(r=r,pvalue=pvalue)