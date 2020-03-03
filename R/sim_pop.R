#' Dail-Madsen Simulator
#'
#' Appendix S4R Code to create population dataFunction to create population data with constant vitalrates (i.e.  C data)
#' @param nsites  : Integer number of sites. Defaults to 15. Likely need 100 or even 200+ sites for estimation with real data.
#' @param nyears  : number of years in the time series
#' @param lambda  : mean density during t=1
#' @param phi     : apparent survival (survival)
#' @param r       : recruitment rate
#' @param p       : detection probability
#' @keywords abundance
#' @export
#' @examples
#' sim_dm()

sim_dm <- function(nsites=15, nyears=30, lambda=100, phi=0.6, r=0.4, p=0.8, plot = FALSE) {

  # Authors : Edwige Bellier, Marc Kery and Michael Schaub# Date    : March 2015
  # y       : count at site 1 at time t

  N <- array(NA,dim=c(nsites, nyears))
  pres <- array(NA, dim=c(nsites, nyears))
  N[,1] <- rpois(nsites, lambda)
  S <- array(NA,dim=c(nsites, nyears-1))
  G <- array(NA,dim=c(nsites, nyears-1))
  phi <- array(phi,dim=c(nsites, nyears-1))
  r <- array(r,dim=c(nsites, nyears-1))

  pres[ , 1] <- ifelse(N[ , 1] > 0, 1, 0)

  # Create population data from Dail and Madsen model (2011)

  for (i in 1:nsites){
    for (t in 1:(nyears-1)){
      S[i,t] <- rbinom(1,N[i,t], phi[i,t])
      G[i,t] <- rpois(1, N[i,t]*r[i,t])
      N[i,t+1] <- S[i,t] + G[i,t]
      pres[i, t+1] <- ifelse(N[i, t+1] > 0, 1, 0)
    }
  }

  # Sampling to obtain counts
  y <- matrix(NA,ncol=nyears, nrow=nsites)
  y[] <- rbinom(nsites*nyears, N, p)
  # Make a graph counts
if(plot) {
  matplot(seq(1,nyears,1),t(y),type="l",lty=1,main="Counts at site i")
}

  return(list(nsites=nsites, nyears=nyears, lambda=lambda, p=p,phi=phi, r=r,y=y, N=N, Nm=mean(N), pres=pres))
}





#' Open Population Simulation
#'
#' Simulate population over time not separating survival and recruitment
#' @param nsites  : Integer number of sites. Defaults to 15. Likely need 100 or even 200+ sites for estimation with real data.
#' @param nyears  : number of years in the time series
#' @param lambda  : mean density during t=1
#' @param r     : per capita pop birth rate
#' @param K       : carrying capacity
#' @param p       : detection probability
#' @keywords abundance
#' @export
#' @examples
#' sim_dm()

sim_dd <- function(nsites=15, nyears=30, lambda=100, phi=0.6, r=0.4, p=0.8, plot = FALSE) {

  # Authors : Edwige Bellier, Marc Kery and Michael Schaub# Date    : March 2015
  # y       : count at site 1 at time t

  N <- array(NA,dim=c(nsites, nyears))
  pres <- array(NA, dim=c(nsites, nyears))
  S <- array(NA,dim=c(nsites, nyears-1))
  G <- array(NA,dim=c(nsites, nyears-1))
  phi <- array(phi,dim=c(nsites, nyears-1))
  r <- array(r,dim=c(nsites, nyears-1))

  N[,1] <- rpois(nsites, lambda)
  pres[ , 1] <- ifelse(N[ , 1] > 0, 1, 0)

  # Create population data from Dail and Madsen model (2011)

  for (i in 1:nsites){
    for (t in 1:(nyears-1)){
      S[i,t] <- rbinom(1,N[i,t], phi[i,t])
      G[i,t] <- rpois(1, N[i,t]*r[i,t])
      lam[i,t+1] <- N[i, t] * exp(r[i] * (1 - (N[i, t] / K[i])))
      pres[i, t+1] <- ifelse(N[i, t+1] > 0, 1, 0)
    }
  }

  # Sampling to obtain counts
  y <- matrix(NA,ncol=nyears, nrow=nsites)
  y[] <- rbinom(nsites*nyears, N, p)
  # Make a graph counts
  if(plot) {
    matplot(seq(1,nyears,1),t(y),type="l",lty=1,main="Counts at site i")
  }

  return(list(nsites=nsites, nyears=nyears, lambda=lambda, p=p,phi=phi, r=r,y=y, N=N, Nm=mean(N), pres=pres))
}



sim_trend <- function(nsites=100, nyears=30, lambda=5, trend = -0.01, p=0.2, plot = FALSE) {

  N <- array(NA_integer_, dim=c(nsites, nyears))
  pres <- array(NA_integer_, dim=c(nsites, nyears))

  N[,1] <- rpois(nsites, lambda)
  pres[ , 1] <- ifelse(N[ , 1] > 0, 1, 0)

  for (i in 1:nsites){
    for (t in 1:(nyears-1)){
      lam <- N[i, t] + trend * N[i, t]
      N[i, t+1] <- rpois(1, lam)
      pres[i, t+1] <- ifelse(N[i, t+1] > 0, 1, 0)
    }
  }

  # Sampling to obtain counts
  y <- matrix(NA, ncol=nyears, nrow=nsites)
  y[] <- rbinom(nsites*nyears, N, p)
  # Make a graph counts
  if(plot) {
    matplot(seq(1,nyears,1),t(y),type="l",lty=1,main="Counts at site i")
  }

  Nt <- colSums(N)
  occ <- colSums(pres) / nsites

  return(list(nsites=nsites, nyears=nyears, lambda=lambda, p=p, y=y, N=N, Nt = Nt, Nm=mean(N), pres=pres, occ = occ))
}


sim_list <- sim_trend(lambda = 2, p = 0.5, nyears = 20)
library(ggplot2)

g1 <- ggplot(data = data.frame(year = 1:sim_list$nyears, Nt = sim_list$Nt), aes(year, Nt)) + geom_line()
g1

library(unmarked)
umf <- unmarkedFramePCO(y = sim_list$y, numPrimary = 20)
summary(umf)
summary(sim_list$y)

# Fit model and backtransform
m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=30, dynamics = "autoreg", se = FALSE))

# Finite sample inference. Abundance at site i, year t
re <- ranef(m1)
(N_hat1 <- colSums(bup(re)))

# lam <- exp(coef(m1, type="lambda"))
# gam <- exp(coef(m1, type="gamma"))
# om <- plogis(coef(m1, type="omega"))
# p <- plogis(coef(m1, type="det"))
#
# # Expected values of N[i,t]
# N_hat2 <- matrix(NA, 100, 30)
# N_hat2[,1] <- lam
# for(t in 2:30) {
#   N_hat2[,t] <- om*N_hat2[ , t-1] + gam
# }

df <- data.frame(year = rep(1:20, times = 2), Estimate = c(rep("Nt", 20), rep("N_hat", 20)), N = c(sim_list$Nt, N_hat1))
ggplot(df, aes(year, N, color = Estimate)) + geom_line()

