
# Appendix S4R Code to create population dataFunction to create population data with constant vitalrates (i.e.  C data)
sim.dm <- function(nsites=15, nyears=30, lambda=100, phi=0.6, r=0.4, p=0.8){

  # Authors : Edwige Bellier, Marc Kery and Michael Schaub# Date    : March 2015# nsites  : number of sampling sites
  # nyears  : number of years in the time series# lambda  : mean density during t=1
  # phi     : apparent survival (survival)
  # r       : recruitment rate# p       : detection probability
  # y       : count at site 1 at time t

  N <- array(NA,dim=c(nsites, nyears))
  N[,1] <- rpois(nsites, lambda)
  S <- array(NA,dim=c(nsites, nyears-1))
  G <- array(NA,dim=c(nsites, nyears-1))
  phi <- array(phi,dim=c(nsites, nyears-1))
  r <- array(r,dim=c(nsites, nyears-1))

  # Create population data from Dail and Madsen model (2011)

  for (i in 1:nsites){
    for (t in 1:(nyears-1)){
      S[i,t] <- rbinom(1,N[i,t], phi[i,t])
      G[i,t] <- rpois(1, N[i,t]*r[i,t])
      N[i,t+1] <- S[i,t] + G[i,t]
      }
    }
  # Sampling to obtain counts
  y <- matrix(NA,ncol=nyears, nrow=nsites)
  y[] <- rbinom(nsites*nyears, N, p)
  # Make a graph counts
  ypar(mfrow=c(1,1),bty="n",las=1)

  matplot(seq(1,nyears,1),t(y),type="l",lty=1,main="Counts at site i")

  return(list(nsites=nsites, nyears=nyears, lambda=lambda, p=p,phi=phi, r=r,y=y, N=N, Nm=mean(N)))
  }

# Function to create population data with density-dependencein recruitment (i.e.  DDR data)

sim.ddr <- function(nsites=15, nyears=30, lambda=100, phi=0.4,mr=0.8, p=0.8, beta=-0.005, sigma=0.2){
  # Authors : Edwige Bellier, Marc Kery and Michael Schaub
  # Date    : March 2015
  # nsites  : number of sampling sites
  # nyears  : number of years in the time series
  # lambda  : mean density during t=1
  # phi     : apparent survival
  # mr      : recruitment rate on log scale (mean recruitment)
  # p       : detection probability
  # y       : count at site 1 at time t
  # beta    : strength of density-dependence
  # sigma   : environmental stochasticity (SD)

  N <- array(NA,dim=c(nsites, nyears))N[,1] <- rpois(nsites, lambda)r <- array(NA,dim=c(nsites, nyears-1))log.r <- array(NA,dim=c(nsites, nyears-1))phi <- array(phi,dim=c(nsites, nyears-1))S <- array(NA,dim=c(nsites, nyears-1))G <- array(NA,dim=c(nsites, nyears-1))

  # Create population data with density-dependence# and environmental stochasticity in recruitment as in Eq.6
  for (i in 1:nsites){
    for (t in 2:(nyears)){
      eps <- rnorm(1,0,sigma)
      log.r[i,t-1] <- log(mr)+beta*N[i,t-1]+eps
      r[i,t-1] <- exp(log.r[i,t-1])
      S[i,t-1] <- rbinom(1,N[i,t-1], phi[i,t-1])
      G[i,t-1] <- rpois(1, N[i,t-1]*r[i,t-1])
      N[i,t] <- S[i,t-1] + G[i,t-1]}}

  # Sampling to obtain countsy <- matrix(NA,ncol=nyears, nrow=nsites)y[] <- rbinom(nsites*nyears, N, p)# Make a graph of counts ypar(mfrow=c(1,1),bty="n",las=1)matplot(seq(1,nyears,1),t(y),type="l",lty=1,main="Counts at site i")return(list(nsites=nsites, nyears=nyears, lambda=lambda,mr=mr, p=p, beta=beta, sigma=sigma,phi=phi, r=r,y=y, N=N, G=G, S=S, Nm=mean(N)))}

#  Function to create population data with density-dependencein survival (i.e.  DDS data)

  sim.dds <- function(nsites=15, nyears=30, lambda=100, mphi=0.55,r=0.6, p=0.8, beta=-0.003, sigma=0.2){# Authors : Edwige Bellier, Marc Kery and Michael Schaub# Date    : March 2015# nsites  : number of sampling sites# nyears  : number of years in the time series# lambda  : mean abundance during t=1 (mean initial population size)# mphi    : apparent survival on logit scale (mean survival)# r       : recruitment rate# p       : detection probability on logit scales# y       : count at site 1 at time t# beta    : strength of density-dependence# sigma   : environmental stochasticity (SD)
    N <- array(NA, dim = c(nsites, nyears))N[,1] <- rpois(nsites, lambda)r <- array(mr, dim = c(nsites, nyears-1))ls.phi <- array(NA, dim = c(nsites, nyears-1))phi <- array(NA, dim = c(nsites, nyears-1))S <- array(NA, dim = c(nsites, nyears-1))G <- array(NA, dim = c(nsites, nyears-1))

    # Create population data with density-dependence# and environmental stochasticity in survival as in Eq. 8
    for (i in 1:nsites)
      {for (t in 2:(nyears)){
        eps <- rnorm(1, 0, sigma)

    ls.phi[i,t-1] <- mphi+beta*N[i,t-1]+ eps
    phi[i,t-1] <- 1/(1+exp(-ls.phi[i,t-1]))
    S[i,t-1] <- rbinom(1, N[i,t-1], phi[i,t-1])
    G[i,t-1] <- rpois(1, N[i,t-1]*r[i,t-1])
    N[i,t] <- S[i,t-1] + G[i,t-1]4
  }}# Sampling to obtain countsy<-matrix(NA,ncol=nyears,nrow=nsites)y[]<-rbinom(nsites*nyears,N,p)# Make a graph of counts ypar(mfrow=c(1,1),bty="n",las=1)matplot(seq(1,nyears,1),t(y),type="l",lty=1,main="Counts at site i")return(list(nsites=nsites, nyears=nyears, lambda=lambda,p=p,beta=beta, sigma=sigma,phi=phi, r=r,y=y, N=N, G=G, S=S, Nm=mean(N)))}Function to create population data with density-dependencein survival and recruitment (i.e.  DDSR data)sim.ddsr <- function(nsites=10, nyears=15, lambda=100, mphi=0.4, mr=0.8,p=0.8, beta1=-0.003, beta2=-0.005, sigma1=0.2, sigma2=0.1){# Authors   : Edwige Bellier, Marc Kery and Michael Schaub# Date      : March 2015# nsites    : number of sampling sites# nyears    : number of years in the time series# lambda    : mean abundance during t=1 (mean initial population size)# mphi      : apparent survival on logit scale (mean survival)# mr        : recruitment rate on log scale (mean recruitment)# p         : detection probability# y         : count at site 1 at time t# beta1     : strength of density-dependence in recruitment5
# beta2     : strength of density-dependence in survival# sigma1    : environmental stochasticity in recruitment# sigma2    : environmental stochasticity in survivalN <- array(NA,dim=c(nsites, nyears))N[,1] <- rpois(nsites, lambda)r <- array(NA,dim=c(nsites, nyears-1))log.r <- array(NA,dim=c(nsites, nyears-1))ls.phi <- array(NA,dim=c(nsites, nyears-1))phi <- array(NA,dim=c(nsites, nyears-1))S <- array(NA,dim=c(nsites, nyears-1))G <- array(NA,dim=c(nsites, nyears-1))# Create population data with density-dependence# and environmental stochasticity in recruitment and survival as in Eq. 9for (i in 1:nsites){for (t in 2:(nyears)){eps1 <- rnorm(1, 0, sigma1)eps2 <- rnorm(1, 0, sigma2)# Density-dependence and environmental stochasticity# in recruitmentlog.r[i,t-1] <- log(mr)+beta1*N[i,t-1]+eps1r[i,t-1] <- exp(log.r[i,t-1])# Density-dependence and environmental stochasticity# in survivalls.phi[i,t-1] <- mphi+beta2*N[i,t-1]+eps2phi[i,t-1] <- 1/(1+exp(-ls.phi[i,t-1]))S[i,t-1] <- rbinom(1,N[i,t-1], phi[i,t-1])G[i,t-1] <- rpois(1,N[i,t-1]*r[i,t-1])N[i,t] <- S[i,t-1] + G[i,t-1]}}6
# Sampling to obtain countsy <- matrix(NA, ncol=nyears, nrow=nsites)y[] <- rbinom(nsites*nyears, N, p)# Make a graph of counts ypar(mfrow=c(1,1),bty="n",las=1)matplot(seq(1,nyears,1),t(y),type="l",lty=1,main="yit Counts")return(list(nsites=nsites, nyears=nyears, lambda=lambda, mphi=mphi, mr=mr,p=p, beta1=beta1, beta2=beta2, sigma1=sigma1, sigma2=sigma2,phi=phi, r=r,y=y, N=N, G=G, S=S, Nm=mean(N)))}JAGS code of the different modelsConstant model in JAGSsink("DM.jags")cat("model{# Priors and constraintsphi ~ dunif(0, 1)r ~ dunif(0, 2)p ~ dunif(0, 1)lambda ~ dunif(0.1, 200)# Likelihood# Initial population sizefor (g in 1:G){N[g,1]~ dpois(lambda)# State process7
for (t in 1:(T-1)){S[g,t+1] ~ dbin(phi,N[g,t])R[g,t+1] ~ dpois(N[g,t] * r)N[g,t+1] < -S[g,t+1] + R[g,t+1]}# Observation processfor (t in 1:T) {y[g,t]~bin(p, N[g,t])}}# g}",fill = TRUE)sink()DDR model in JAGSsink("DDR.jags")cat("model{# Priors and constraintss ~ dunif(0, 1)mean.r ~ dunif(0, 3)lmr <- log(mean.r)beta ~ dnorm(0, 0.001)p ~ dunif(0, 1)lambda ~ dunif(0.1, 200)sigma ~ dunif(0, 2)tau <- 1/ pow(sigma, 2)# Likelihood# Initial population sizefor (g in 1:G){N[g,1] ~ dpois(lambda)8
# State processfor (t in 1:(T-1)){S[g,t] ~ dbin(s, N[g,t])R[g,t] ~ dpois(N[g,t] * r[g,t])N[g,t+1] <- S[g,t] + R[g,t]

# Density-dependence and environmental stochasticity
# in recruitment
lr[g,t] <- lmr + beta*(N[g,t]-Nm) + eps[g,t]r[g,t] <- exp(lr[g,t])eps[g,t] ~ dnorm(0, tau)}# Observation processfor (t in 1:T) {y[g,t] ~ dbin(p, N[g,t])}} # g}",fill = TRUE)sink()DDS model in JAGSsink("DDS.jags")cat("model{# Priors and constraintsr ~ dunif(0, 2)mean.phi ~ dunif(0, 1)mphi<- log(mean.phi/(1-mean.phi))beta ~ dunif(-1, 1)p ~ dunif(0, 1)lambda ~ dunif(0.1, 200)sigma ~ dunif(0, 0.5)9
tau <- 1/pow(sigma, 2)# Likelihood# Initial population sizefor (g in 1:G){N[g,1] ~ dpois(lambda)# State processfor (t in 1:(T-1)){S[g,t] ~ dbin(phi[g,t], N[g,t])R[g,t] ~ dpois(N[g,t] * r)N[g,t+1] <- S[g,t] + R[g,t]# Density-dependence and environmental stochasticity# in survivallphi[g,t] <- mphi + beta*(N[g,t]-Nm) + eps[g,t]phi[g,t] <- 1/(1+exp(-lphi[g,t]))eps[g,t] ~ dnorm(0, tau)}# Observation processfor (t in 1:T) {y[g,t] ~ dbin(p, N[g,t])}} # g}",fill = TRUE)sink()DDSR model in JAGSThe model can be fitted in R by using the library jagsUIsink("DDSR.jags")cat("model{10
# Priors and constraintsmean.r ~ dunif(0, 3)lmr <- log(mean.r)beta1 ~ dunif(-1, 1)sigma1 ~ dunif(0, 2)tau1 <- 1/pow(sigma1, 2)mean.phi ~ dunif(0, 1)mphi <- log(mean.phi/(1-mean.phi))beta2 ~ dunif(-1, 1)sigma2 ~ dunif(0, 2)tau2 <- 1/pow(sigma2, 2)p ~ dunif(0, 1)lambda ~ dunif(0.1, 200)# Likelihood# Initial population sizefor (g in 1:G){N[g,1] ~ dpois(lambda)# State processfor (t in 1:(T-1)){S[g,t] ~ dbin(phi[g,t], N[g,t])R[g,t] ~ dpois(N[g,t] * r[g,t])N[g,t+1] <- S[g,t] + R[g,t]# Density-dependence and environmental stochasticity in recruitmentlr[g,t] <- lmr + beta1*(N[g,t]-Nm)+ eps1[g,t]r[g,t] <- exp(lr[g,t])eps1[g,t] ~ dnorm(0, tau1)# Density-dependence and environmental stochasticity in survivallphi[g,t] <- mphi + beta2*(N[g,t]-Nm) + eps2[g,t]11

