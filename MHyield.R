####Estimation using a Metropolis-Hastings algorithm####

###Seed used for random sampling###
set.seed(2)

###Number of simulations###
N<-100

###Prior distribution###
mu<-1.4
tau<-0.2

###Data###
Y<-1.9
sigma<-0.2

###True posterior###
B<-tau^2/(tau^2+sigma^2)
Mu_post<-(1-B)*mu + B*Y
V_post<-(1-B)*tau^2

###Initialisation of the parameter chains###
Theta<-1:N
Theta[1]<-mu
Reject<-0

###Graphic of likelihood

plot(seq(0.5,3,by=0.01),dnorm(Y,seq(0.5,3,by=0.01),sigma), xlab="Yield (t ha-1)", ylab="Density",type="l",lwd=3,lty=3,ylim=c(0,5))
points(Y,0,pch=19,cex=2)
lines(seq(0.5,3,by=0.01),dnorm(seq(0.5,3,by=0.01),mu,tau),col="darkgreen",lwd=3,lty=3)
lines(seq(0.5,3,by=0.01),dnorm(seq(0.5,3,by=0.01),Mu_post,sqrt(V_post)),col="blue",lwd=3)

###Computation of the likelihood for the first set of parameter values###
Likeli<-1:N
Likeli[1]<-dnorm(Y,Theta[1],sigma)

###Tuning parameters###
Lambda<-0.1

###Generation of candidate parameter values and test###

for (i in 2:N) {

	Theta.c<-rnorm(1,Theta[i-1],Lambda)

	Likeli[i]<-dnorm(Y,Theta.c,sigma)
	
	Likeli[i]<-max(Likeli[i], 10^(-10))

	Test.i<-log(Likeli[i]) + log(dnorm(Theta.c,mu,tau)) - log(Likeli[i-1]) - log(dnorm(Theta[i-1],mu,tau))

	u<-runif(1,min=0,max=1)
	
	if (Test.i>log(u)) {Theta[i]<-Theta.c}

	if (Test.i<=log(u)) {Theta[i]<-Theta[i-1]
			Reject<-Reject+1
				}			
		}

####Prior####
par(mfrow=c(2,2))

plot(seq(0,max(Theta),by=max(Theta)/100),dnorm(seq(0,max(Theta),by=max(Theta)/100),mu,tau),col="darkgreen",type="l",lwd=3, xlab="Yield", ylab="Density")

###Chains of parameter values
plot(1:N, Theta, xlab="Iteration", ylab="Yield (t ha-1)", type="l",cex=2)

###Posterior distribution###
hist(Theta[(N/2):N],labels=F,xlab="Yield (t ha-1)",main=" ",freq=F,cex=2,xlim=c(0.5, 2.5), ylim=c(0,3),breaks=seq(0,max(Theta[(N/2):N]),by=max(Theta[(N/2):N])/100))
points(Y,0,pch=20,cex=2)
lines(density(Theta[(N/2):N]),col="red", lwd=2)
lines(seq(min(Theta),max(Theta),by=(max(Theta)-min(Theta))/100),dnorm(seq(min(Theta),max(Theta),by=(max(Theta)-min(Theta))/100),Mu_post,sqrt(V_post)),col="blue",lwd=2)
lines(seq(0,max(Theta),by=max(Theta)/100),dnorm(seq(0,max(Theta),by=max(Theta)/100),mu,tau),col="darkgreen",lwd=2)

print("Posterior mean:")	
print(mean(Theta[1:(N/2)]))
print("Posterior variance:")
print(var(Theta[1:(N/2)]))
print("Acceptance rate (%):")
print(100*(N-Reject)/N)		




