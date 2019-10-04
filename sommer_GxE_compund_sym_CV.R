
library(MASS)
library(rrBLUP)
library(sommer)
library(INLA)


geno=read.table("rrblup_geno.txt",header=T)

dat=read.table("pheno_gxe.txt",header=T)


names(dat)[1] <- "animal"

dat=dat[order(dat$animal),]
dat$animal2=dat$animal


marker=geno
marker$Marker=NULL
marker$Chr=NULL
marker$Pos=NULL
marker=marker[,match(unique(dat$animal),colnames(marker))]

M=t(marker)
colnames(M)=geno$Marker

A=A.mat(M,shrink=T)

dat$Env=as.factor(dat$Env)

co_ARK=vector()


co_KAK=vector()


co_KTI=vector()

reco=dim(dat)[1]/3

rep=10

for ( i in 1:rep)
{

CV=dat

set.seed(i)

val=sample(reco,80)

del=unique(dat$animal)[val]

CV[CV$animal %in% del,3]=NA

mix1=mmer2(Phe~1+Env,random=~g(animal)+diag(Env):g(animal),rcov=~units, G=list(animal=A),data=CV,draw=TRUE)

res1=data.frame(mix1$u.hat$"1:g(animal)"+mix1$u.hat$"g(animal)")

res2=mix1$u.hat$"2:g(animal)"+mix1$u.hat$"g(animal)"

res3=mix1$u.hat$"3:g(animal)"+mix1$u.hat$"g(animal)"

dat1=dat[dat$Env %in% 1,]
dat2=dat[dat$Env %in% 2,]
dat3=dat[dat$Env %in% 3,]


co_ARK[i]=cor(dat1$Phe[del],res1[,1][del])

co_KAK[i]=cor(dat2$Phe[del],res2[,1][del])

co_KTI[i]=cor(dat$Phe[del],res3[,1][del])

}
