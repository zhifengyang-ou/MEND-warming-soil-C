library(vegan)
library(picante)
library(ieggr)
library(lme4)

#permutation test 
otu=read.table("community.txt",header=T,sep="\t",row.names = 1)
T.otu=t(otu[,5:556])
T.otu[is.na(T.otu)]<-0
treat=read.table("treatment_list.txt",header=T,sep="\t",row.names = 1)

Bray=vegdist(T.otu,method="bray")
adonis(Bray~Warm*Precipitation*Clip*year+block,data=treat,perm=999)#overall adonis

#######test for gene abundance
library(vegan)
All=read.table("community.txt", header=T, sep="\t", row.names=1) ##fill0 already
head(All)

All = cbind(All[,1:4],All[,5:280],All[,281:556])
All[1,5:556]
delete = which(rowSums(All[,5:556])==556)
if(length(delete!=0)){
  All = All[-delete,]
}
a=as.matrix(All[,5:556])
mean(a)
gene.name = All[,1]
gene.list=unique(gene.name)
n = length(gene.list)

Gene.name = as.vector(gene.list)
Probe.Num = c(rep(9999,n))
Gene.category = c(rep(9999,n))
Gene.subcategory1 = c(rep(9999,n))
Gene.subcategory2 = c(rep(9999,n))
warming.p = c(rep(9999,n))
warming.pn = c(rep(9999,n))
bonferroni = c(rep(9999,n))
fdr = c(rep(9999,n))
difference = c(rep(9999,n))
std.err = c(rep(9999,n))

for(i in 1:n){
  geneselected = which(gene.name == gene.list[i]) ##for cycle
  genetable = All[geneselected,]
  genetable.d = genetable[,5:556] ## change the row and column
  Y = t(genetable.d)
  lable = t(genetable[,1:4]) ## change please see your All file
  y1= unlist(genetable.d[1:276], use.names = FALSE)
  y2= unlist(genetable.d[277:552], use.names = FALSE)
  
  if(length(y1) > 1){
    
    ttest = t.test(y1,y2,paired = TRUE)
    
    Gene.category[i] = lable[2,1]
    Gene.subcategory1[i] = lable[3,1]
    Gene.subcategory2[i] = lable[4,1]
    Probe.Num[i] = length(genetable[,1])
    warming.p[i] = ttest[[3]]
    warming.pn[i] = ttest[[3]]*length(y1)
    difference[i] = (sum(y2)-sum(y1))/length(y1)
    #difference[i] = sum(y1)/sum(y2)
    std.err[i] = sd(y2 - y1)/sqrt(length(y1))
  }
  else{
    Gene.category[i] = lable[2,1]
    Gene.subcategory1[i] = lable[3,1]
    Gene.subcategory2[i] = lable[4,1]
    Probe.Num[i] = length(genetable[,1])
    warming.p[i] = 9999
    warming.pn[i] = 9999
    difference[i] = 9999
    std.err[i] = 9999
  }
}
bonferroni = p.adjust(warming.p, method = 'bonferroni')
fdr = p.adjust(warming.p, method = 'fdr')
GeneAOVtable = cbind(Gene.name,Probe.Num,Gene.category,Gene.subcategory1,Gene.subcategory2,warming.p,warming.pn,bonferroni,fdr,difference,std.err)
write.csv(GeneAOVtable, file="E:/ten year of OK new warming/OK.NWS.GeoChip.2009-2020/Gene responseANOVA/Gene_warming_responseANOVA_Across all years.csv", row.names=F)


#Link environmental variables to community
otu.env=read.table("environmental variable.txt",header=T,sep="\t",row.names = 1)
#mantel test 
otu.dist=vegdist(T.otu)
variable.dist=vegdist(scale(otu.env$variable),"euclid")
mantel(otu.dist,variable.dist)

##correlation matrix 

library(corrplot)
dat<-read.csv("env factors for correlation matrix.csv",header = T,row.names = 1)

M<-cor(dat,use = "pairwise.complete.obs")
corrplot(M, method= "square",type="upper",diag=FALSE,tl.col="black",tl.srt=45)

write.csv(M,"env correlation matrix.csv")

########correlation of genes and functioning
library(vegan)
library(picante)
library(lme4)
library(lmerTest)
library(MuMIn)
gene=read.table("functional_gene.txt",header=T,sep="\t",row.names = 1)
env=read.table("environmental variable.txt",header=T,sep="\t",row.names = 1)
gene[is.na(gene)]<-0
env[is.na(env)]<-0
fm1=lmer(gene$G1~scale(env$variable)+(1|block)+(1|year),data=env)
summary(fm1)
presult<-car::Anova(fm1,type=2)
coefs<-coef(summary(fm1))[ , "Estimate"]
pvalue=presult[,3]
r2<-r.squaredGLMM(fm1)
r=ifelse(coefs>0,(r2[1,1])^0.5,-(r2[1,1])^0.5)
result=list(r=r,pvalue=pvalue)

#piecewise SEM
library(nlme)
library(lme4)
library(piecewiseSEM)
library(QuantPsyc)
library(piecewiseSEM) 
library(semTools)
library(lavaan)

data1 <- read.csv("summary of variables.csv",header = T)

normalized_data <- data1
normalized_data[, c("warm","precip","clip","half","double","NO3", "NH4", "TN", "TC", "pH", "FlTotl","FlC4","FlC3","plantrichness","annual_moisture","temperature_annual","PC1_GeoChip","metabolic_quotient_RhDNAcon")] <- scale(data1[, c("warm","precip","clip","half","double", "NO3", "NH4", "TN", "TC", "pH", "FlTotl","FlC4","FlC3","plantrichness","annual_moisture","temperature_annual","PC1_GeoChip","metabolic_quotient_RhDNAcon")])
data1 <- normalized_data
data1$block <- factor(data1$block)
data1$year <- factor(data1$year)

###############

model <- psem(
  lmer(temperature_annual ~ warm*half+warm*double +clip+ (1|block)+(1|year), data1),
  lmer(annual_moisture ~ warm*half+warm*double +clip+temperature_annual+ (1|block)+(1|year), data1),
  lmer(pH ~ warm*half+warm*double +clip+temperature_annual+annual_moisture+ (1|block)+(1|year), data1),
  lmer(NO3 ~ warm*half+warm*double +clip+temperature_annual+annual_moisture+ (1|block)+(1|year), data1),
  lmer(FlTotl ~ warm*half+warm*double +clip+temperature_annual+annual_moisture+NO3+plantrichness+ (1|block)+(1|year), data1),
  lmer(plantrichness ~ warm*half+warm*double +clip+temperature_annual+annual_moisture+ (1|block)+(1|year), data1),
  lmer(PC1_GeoChip ~ warm*half+warm*double +clip+temperature_annual+annual_moisture+pH+NO3+FlTotl+plantrichness+ (1|block)+(1|year), data1),
  lmer(metabolic_quotient_RhDNAcon ~ warm*half+warm*double +clip+temperature_annual+annual_moisture+pH+FlTotl+PC1_GeoChip+ (1|block)+(1|year), data1),
  lmer(TC ~ warm*half+warm*double +clip+temperature_annual+annual_moisture+NO3+FlTotl+PC1_GeoChip+metabolic_quotient_RhDNAcon+ (1|block)+(1|year), data1),
  
  data = data1
)

summary(model, .progressBar = F)

