
library(rrBLUP)

library(sommer)

dat=read.table("pheno_rrblup.txt",header=T)
geno=read.table("rrblup_geno.txt",header=T)

names(dat)[1] <- "animal"

#dat=dat[order(dat$animal),]

marker=geno
marker$Marker=NULL
marker$Chr=NULL
marker$Pos=NULL
marker=marker[,match(dat$animal,colnames(marker))]

M=t(marker)
colnames(M)=geno$Marker

A=A.mat(M,shrink=T)

#dat$animal=as.factor(dat$animal)
reco=dim(dat)[1]

#WAS2(fixed, random, rcov, data, weights, G=NULL, M=NULL, grouping=NULL, 

###MCMC

#mix2<- mmer2(cbind(EBU,KAK,KTI) ~1, random = ~ us(trait):g(animal), rcov= ~ us(trait):units, G=list(animal=A),data =dat,draw=FALSE)
qtl=GWAS2(cbind(ARK,Abr,Fad) ~1, random = ~ us(trait):g(animal), G=list(animal=A),M=M)

co_EBU_EBU_multi_IDV=vector()
#trait FAD

co_KAK_multi_IDV=vector()

#trait Abr

co_KTI_multi_IDV=vector()

rep=10

for ( i in 1:rep)
{

CV=dat

set.seed(i)

val=sample(reco,80)

CV[val,2]=NA

#single trait

mix3 <- mmer2(cbind(ARK,Abr,Fad) ~1, random = ~ us(trait):g(animal),  rcov= ~ us(trait):units, G=list(animal=A),data =CV,draw=FALSE)

dat_test=data.frame(Id=names(mix3$u.hat[[1]][,1]),val1=mix3$u.hat[[1]][,1],val2=mix3$u.hat[[1]][,2],val3=mix3$u.hat[[1]][,3])

test=dat_test[match(CV$animal,dat_test$Id),]

co_EBU_EBU_multi_IDV[i]=cor(test$val1[val],dat$ARK[val])

co_KAK_multi_IDV[i]=cor(test$val2[val],dat$Abr[val])

co_KTI_multi_IDV[i]=cor(test$val3[val],dat$Fad[val])

}


