
library(MASS)
library(rrBLUP)
library(synbreed)
library(MCMCglmm)
library(INLA)


dat=read.table("pheno_rrblup.txt",header=T)
geno=read.table("rrblup_geno.txt",header=T)

marker=geno
marker$Marker=NULL
marker$Chr=NULL
marker$Pos=NULL
marker=marker[,match(dat$Gid,colnames(marker))]

M=t(marker)
colnames(M)=geno$Marker

A=A.mat(M,shrink=T)

A_inv=ginv(A)
A_spar=as(A_inv,"dgCMatrix")

names(dat)[1] <- "animal"
dat$animal=as.factor(dat$animal)
	dat$ARK=as.numeric(dat$ARK)
	dat$Fad=as.numeric(dat$Fad)
	dat$Abr=as.numeric(dat$Abr)
	

rownames(A_spar)=dat$animal

Phe.var<-matrix(c(var(dat$ARK,na.rm=TRUE),0,0,0,var(dat$Fad,na.rm=TRUE),0,0,0,var(dat$Abr,na.rm=TRUE)),3,3)

#3##########################
#single trait
#3##########################

#############
#MCMCglmm priors
#############




prior2.2 <- list(G = list(G1 = list(V =1, nu = 1)), R = list(V =1, nu = 1))

###########
# first trait
##############

###MCMC

prior_ARK <- list(G = list(G1 = list(V =Phe.var[1][1]/3, n=1)), R = list(V =Phe.var[1][1]/3, n=1))
set.seed(100)
mod_mcmc_ark=MCMCglmm(ARK~1,random=~animal,ginverse=list(animal=A_spar),pr=T,data=dat,nitt=10000,thin=50,burnin=3000,prior=prior_ARK,verbose=TRUE)

val_ark=colMeans (mod_mcmc_ark$Sol)

breed_ark=val_ark[2:(dim(A)[1]+1)]
co_ark=cor(dat$ARK,breed_ark)
co_ark
posterior.mode(mod_mcmc$VCV)
#RRBLUP
ans_ark=kin.blup(dat,K=A,geno="animal",pheno="ARK")

ans_ark$Vg
ans_ark$Ve

cor(ans_ark$g,dat$ARK)

#3##########################


###########
# second trait
##############

###MCMC

prior_fad <- list(G = list(G1 = list(V =Phe.var[2,2]/3, n=1)), R = list(V =Phe.var[2,2]/3, n=1))
set.seed(100)
mod_mcmc_fad=MCMCglmm(Fad~1,random=~animal,ginverse=list(animal=A_spar),pr=T,data=dat,nitt=10000,thin=50,burnin=3000,prior=prior_fad,verbose=TRUE)

val_fad=colMeans (mod_mcmc_fad$Sol)

breed_fad=val_fad[2:(dim(A)[1]+1)]
co_fad=cor(dat$Fad,breed_fad)

posterior.mode(mod_mcmc_fad$VCV)
#RRBLUP
ans_fad=kin.blup(dat,K=A,geno="animal",pheno="Fad")

ans_fad$Vg
ans_fad$Ve

cor(ans_fad$g,dat$Fad)
#3##########################

###########
# Third trait
##############

prior_Abr <- list(G = list(G1 = list(V =Phe.var[3,3]/3, n=1)), R = list(V =Phe.var[3,3]/3, n=1))
set.seed(100)
mod_mcmc_Abr=MCMCglmm(Abr~1,random=~animal,ginverse=list(animal=A_spar),pr=T,data=dat,nitt=10000,thin=50,burnin=3000,prior=prior_Abr,verbose=TRUE)

val_Abr=colMeans (mod_mcmc_Abr$Sol)

breed_Abr=val_Abr[2:(dim(A)[1]+1)]
co_Abr=cor(dat$Abr,breed_Abr)

posterior.mode(mod_mcmc_Abr$VCV)
#RRBLUP
ans_Abr=kin.blup(dat,K=A,geno="animal",pheno="Abr")

ans_Abr$Vg
ans_Abr$Ve

cor(ans_Abr$g,dat$Abr)

#3###########################3###########################3##########################

###########
# Multi trait
###############3###########################3###########################3##########################


Phe.var<-matrix(c(var(dat$ARK,na.rm=TRUE),0,0,0,var(dat$Fad,na.rm=TRUE),0,0,0,var(dat$Abr,na.rm=TRUE)),3,3)

prior2.2 <- list(G = list(G1 = list(V = Phe.var/3, n = 3)), R = list(V = Phe.var/3, n = 3))

#nu=0

#prior2.J <- list(G = list(G1 = list(V = Phe.var/3, n = 3)), R = list(V = Phe.var/3, n = 3,nu=0.002))

prior2.2 <- list(G = list(G1 = list(V = Phe.var/3, n = 3)), R = list(V = Phe.var/3, n = 3))

#prior2.1 <- list(G = list(G1 = list(V = diag(3), n = 1.002)), R = list(V = diag(3), n = 1.002))



prior2.2 <- list(G = list(G1 = list(V = Phe.var/3, n = 3)), R = list(V = Phe.var/3, n = 3))
set.seed(100)
model2.0<-MCMCglmm(cbind(ARK,Fad,Abr) ~trait-1,random=~us(trait):animal,rcov=~us(trait):units,family=c("gaussian","gaussian","gaussian"),pr=T,data=dat,ginverse=list(animal=A_spar),
	nitt=10000,thin=50,burnin=3000,prior=prior2.2,verbose=TRUE)
val_multi=colMeans (model2.0$Sol)
posterior.mode(model2.0$VCV)

#independent error
###################################################

 model2.id<-MCMCglmm(cbind(ARK,Fad,Abr) ~trait-1,random=~us(trait):animal,rcov=~idh(trait):units,family=c("gaussian","gaussian","gaussian"),pr=T,data=CV,ginverse=list(animal=A_spar),
+ nitt=10000,thin=50,burnin=3000,prior=prior2.2,verbose=TRUE)


#breeding values 


val_multi=colMeans (model2.0$Sol)
#trait correlations
#trait correlations
cor(val_multi[4:374],dat$ARK)
cor(val_multi[375:(374+reco)],dat$Fad)
cor(val_multi[746:(745+reco)],dat$Abr)

#Ante dependence mmodel first order 

prior2.2 <- list(G = list(G1 = list(V = Phe.var/3, n = 3)), R = list(V = Phe.var/3, n = 3))

#prior for ante depe fix the regression of trait 3 on trait 2 to be zero

prior2.ant <- list(G = list(G1 = list(V = diag(3), n = 3)), R=list(V=diag(3), nu=0, beta.mu=rep(0,3), beta.V=diag(3)*100))

prior2.ant$R$beta.V[3,3]<-0.00000001

plot(posterior.ante(m1$VCV[,1:9], k=2))


set.seed(100)

model2.ant1<-MCMCglmm(cbind(ARK,Fad,Abr) ~trait-1,random=~us(trait):animal,rcov=~ante1(trait):units,family=c("gaussian","gaussian","gaussian"),pr=T,data=dat,ginverse=list(animal=A_spar),
	nitt=10000,thin=50,burnin=3000,prior=prior2.ant,verbose=TRUE)
posterior.mode(model2.ant1$VCV)

val_ant1=colMeans (model2.ant1$Sol)
#trait correlations
cor(val_ant1[4:416],dat$ARK)
cor(val_ant1[417:(416+413)],dat$Fad)
cor(val_ant1[830:(829+413)],dat$Abr)

#Ante dependence mmodel second order 

prior2.2 <- list(G = list(G1 = list(V = Phe.var/3, n = 3)), R = list(V = Phe.var/3, n = 3))
set.seed(1100)

model2.ant2<-MCMCglmm(cbind(ARK,Fad,Abr) ~trait-1,random=~us(trait):animal,rcov=~ante2(trait):units,family=c("gaussian","gaussian","gaussian"),pr=T,data=dat,ginverse=list(animal=A_spar),
	nitt=10000,thin=50,burnin=3000,prior=prior2.2,verbose=TRUE)

posterior.mode(model2.ant2$VCV)

val_ant2=colMeans (model2.ant2$Sol)
#trait correlations
cor(val_ant2[4:416],dat$ARK)
cor(val_ant2[417:(416+413)],dat$Fad)
cor(val_ant2[830:(829+413)],dat$Abr)



