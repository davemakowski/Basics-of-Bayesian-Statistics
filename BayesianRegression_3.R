# David Makowski
# 2023

#Data and their summaries
TAB<-read.table("Data_Zhao_Wheat.txt", sep="\t", header=T)
summary(TAB)

###Data###

Y=TAB$Sensitivity
Temp=TAB$TGS
REF=as.numeric(as.factor(as.character(TAB$Site_name)))

DATA=data.frame(REF, Y, Temp)

############
#MCMCglmm#
############

library(MCMCglmm)

#prior1<-list(B=list(mu=c(0,-10),V=diag(c(10^8,1))), R=list(V=1,nu=1),G=list(G1=list(V=1,nu=1)))

Mod_mcmc<-MCMCglmm(Y~1+Temp,random=~REF, data=DATA,verbose=F, nitt=50000, 
                   thin=10, burnin=10000, pr=TRUE)

summary(Mod_mcmc)

#Graphique de la chaine de valeur pour le paramètre mu
plot(Mod_mcmc)

plot(Temp,Y, xlab="Average temperature (°C)", ylab="Yield sensitivity (%)", pch=19)

Temp_vec=5:15
Response=Mod_mcmc$Sol[,1]+as.matrix(Mod_mcmc$Sol[,2])%*%t(Temp_vec)

med_rep=apply(Response, 2, median)
q2.5_rep=apply(Response, 2, quantile, 0.025)
q97.5_rep=apply(Response, 2, quantile, 0.975)

lines(Temp_vec, med_rep, col="red", lwd=3)
lines(Temp_vec, q2.5_rep, col="red", lty=2)
lines(Temp_vec, q97.5_rep, col="red", lty=2)
abline(h=0, lty=2)

Mod_mcmc_1<-MCMCglmm(Y~1+Temp,random=~REF, data=DATA,verbose=F, nitt=50000, thin=10, burnin=10000)
Mod_mcmc_2<-MCMCglmm(Y~1+Temp,random=~REF, data=DATA,verbose=F, nitt=50000, thin=10, burnin=10000)
Mod_mcmc_3<-MCMCglmm(Y~1+Temp,random=~REF, data=DATA,verbose=F, nitt=50000, thin=10, burnin=10000)

ChainList<-mcmc.list(Mod_mcmc_1$Sol,Mod_mcmc_2$Sol,Mod_mcmc_3$Sol)
gelman.plot(ChainList)
