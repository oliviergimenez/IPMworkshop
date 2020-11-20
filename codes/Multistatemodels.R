################################################
#
# Multistate capture-recapture models
#
# IPM workshop, virtual Montpellier November 2020
#
################################################
# Based on the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" 
# by Marc Kéry & Michael Schaub (2012, Academic Press)
#########################################################################

# Michael Schaub
# Olivier Gimenez 19.11.2020
#
################################################

#########################################################################

# Make sure that you have installed JAGS version 3.2.0
library(jagsUI)  # Load package

#########################################################################

# Specify path
setwd('...')

#########################################################################

# Load multistate capture-histories
data <- read.csv("dat/multistate.csv", sep = ";")
ch <- as.matrix(data)
ch
nrow(ch)
ncol(ch)

###############################################

# A simple multistate model is shown (corresponding to chapter 9.2. in BPA)
# There are 2 sites, at each site individuals are marked and recaptured
# Individuals may move between these sites. The goal is to estimate site-specific survival and recapture probabilities as well as movement rates
# Two approaches are presented: the state-space likelihood and the multinomial likelihood; we will focus on the former, the code is provided for the latter

###############################################



# 1. Approach: state-space likelihood

# 1.1. Specify model in BUGS language
cat(file = "state-space-ms.jags", "
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# psiAB: movement probability from site A to site B
# psiBA: movement probability from site B to site A
# pA: recapture probability at site A
# pB: recapture probability at site B
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 dead
# Observations (O):  
# 1 seen at A 
# 2 seen at B
# 3 not seen
# -------------------------------------------------

# Priors and constraints
   phiA ~ dunif(0, 1)    # Priors for state-spec. survival
   phiB ~ dunif(0, 1)    # Priors for state-spec. survival
   psiAB ~ dunif(0, 1)    # Priors for transitions
   psiBA ~ dunif(0, 1)    # Priors for transitions
   pA ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   pB ~ dunif(0, 1)      # Priors for mean state-spec. recapture

# Define state-transition and observation matrices
   # Define probabilities of state S(t+1) given S(t)
      ps[1,1] <- phiA * (1 - psiAB)
      ps[1,2] <- phiA * psiAB
      ps[1,3] <- 1 - phiA
      ps[2,1] <- phiB * psiBA
      ps[2,2] <- phiB * (1-psiBA)
      ps[2,3] <- 1 - phiB
      ps[3,1] <- 0
      ps[3,2] <- 0
      ps[3,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,1] <- pA
      po[1,2] <- 0
      po[1,3] <- 1 - pA
      po[2,1] <- 0
      po[2,2] <- pB
      po[2,3] <- 1 - pB
      po[3,1] <- 0
      po[3,2] <- 0
      po[3,3] <- 1

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], 1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], 1:3])
      } #t
   } #i
}
")


# 1.2. Data preparation
# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(ch, 1, get.first)

# Recode ch matrix: note, a 0 is not allowed!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rch <- ch          # Recode ch
rch[rch==0] <- 3

# Bundle data 
jags.data <- list(y = rch, f = f, n.occasions = ncol(rch), nind = nrow(rch))


# 1.3. Initial values
# Function to create initial values for unknown z
# Input: ch: multistate capture-recapture data (where "not seen" is not 0); f: vector with the occasion of marking
ms.init.z <- function(ch, f){
   states <- max(ch, na.rm = TRUE)   # identify the state that corresponds to "not seen"
   known.states <- 1:(states-1)
   v <- which(ch==states)
   ch[v] <- sample(known.states, length(v), replace = TRUE)   
   for (i in 1:nrow(ch)){ch[i,1:f[i]] <- NA}
   return(ch)
   }

inits <- function(){list(phiA = runif(1, 0, 1), phiB = runif(1, 0, 1), 
                         psiAB = runif(1, 0, 1), psiBA = runif(1, 0, 1), 
                         pA = runif(1, 0, 1), pB = runif(1, 0, 1),
                         z = ms.init.z(rch, f))}  

# 1.4. Parameters monitored
parameters <- c("phiA", "phiB", "psiAB", "psiBA", "pA", "pB")

# 1.5. MCMC settings
ni <- 400
nt <- 1
nb <- 200
nc <- 3

# 1.6. Call JAGS from R
ms <- jags(jags.data, inits, parameters, "state-space-ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# 1.7. Inspect results
print(ms, digits = 3)

par(mfrow = c(3,3))
traceplot(ms)




# If ran longer...

# 1.5. MCMC settings
ni <- 10000
nt <- 1
nb <- 2000
nc <- 2

# 1.6. Call JAGS from R
ms <- jags(jags.data, inits, parameters, "state-space-ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# 1.7. Inspect results
par(mfrow = c(3,3))
traceplot(ms)

print(ms, digits = 3)


# JAGS output for model 'state-space-ms.jags', generated by jagsUI.
# Estimates based on 2 chains of 10000 iterations,
# adaptation = 100 iterations (sufficient),
# burn-in = 2000 iterations and thin rate = 1,
# yielding 16000 total samples from the joint posterior. 
# MCMC ran for 4.659 minutes at time 2020-11-19 22:33:49.
# 
# mean      sd    2.5%      50%    97.5% overlap0 f  Rhat n.eff
# phiA        0.774   0.023   0.732    0.773    0.822    FALSE 1 1.000 16000
# phiB        0.677   0.033   0.608    0.679    0.738    FALSE 1 1.006   249
# psiAB       0.305   0.090   0.141    0.304    0.477    FALSE 1 1.009   291
# psiBA       0.578   0.096   0.424    0.565    0.789    FALSE 1 1.007   408
# pA          0.761   0.099   0.599    0.753    0.972    FALSE 1 1.012   182
# pB          0.385   0.143   0.218    0.344    0.782    FALSE 1 1.004  2802
# deviance 1435.831 211.512 871.615 1494.200 1702.800    FALSE 1 1.013   129
# 
# Successful convergence based on Rhat values (all < 1.1). 
# Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
# For each parameter, n.eff is a crude measure of effective sample size. 
# 
# overlap0 checks if 0 falls in the parameter's 95% credible interval.
# f is the proportion of the posterior with the same sign as the mean;
# i.e., our confidence that the parameter is positive or negative.
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 22195.9 and DIC = 23631.71 
# DIC is an estimate of expected predictive error (lower is better).


###########################
# Exercise 1:
# Write the model in such a way that survival in A is time-dependent (fixed time effect) 
# and survival at site B it is time-dependent with a random temporal effect.
###########################


