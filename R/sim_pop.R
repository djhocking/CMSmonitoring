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

sim_dm <- function(nsites=15, nyears=30, lambda=100, phi=0.6, r=0.4, p=0.8) {

  # Authors : Edwige Bellier, Marc Kery and Michael Schaub# Date    : March 2015
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
