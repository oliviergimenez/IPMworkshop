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

inits <- function(){list(phiA = runif(1, 0, 1), 
                         phiB = runif(1, 0, 1), 
                         psiAB = runif(1, 0, 1), 
                         psiBA = runif(1, 0, 1), 
                         pA = runif(1, 0, 1), 
                         pB = runif(1, 0, 1),
                         z = ms.init.z(rch, f))}  

# 1.4. Parameters monitored
parameters <- c("phiA", "phiB", "psiAB", "psiBA", "pA", "pB")

# 1.5. MCMC settings
ni <- 400
nt <- 1
nb <- 200
nc <- 3

# 1.6. Call JAGS from R
ms <- jags(jags.data, 
           inits, 
           parameters, 
           "state-space-ms.jags", 
           n.chains = nc, 
           n.thin = nt, 
           n.iter = ni, 
           n.burnin = nb,
           parallel = TRUE)

# 1.7. Inspect results
print(ms, digits = 3)

par(mfrow = c(3,3))
traceplot(ms)


# If ran longer...

# 1.5. MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 2

# 1.6. Call JAGS from R
ms <- jags(jags.data, 
           inits, 
           parameters, 
           "state-space-ms.jags", 
           n.chains = nc, 
           n.thin = nt, 
           n.iter = ni, 
           n.burnin = nb,
           parallel = TRUE)

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

#
# You may find the solution code below
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# Specify model in BUGS language
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
for (t in 1:(n.occasions-1)){
   phiA[t] <- alpha[t] # fixed time effect 
   logit(phiB[t]) <- mu + epsilon[t] # random time effect
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   }
for (u in 1:2){
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   }
for (t in 1:(n.occasions-1)){
   alpha[t] ~ dunif(0, 1)    # Priors for survival A
   epsilon[t] ~ dnorm(0, taueps)      # Priors for variance RE survival B
   }
mu ~ dnorm(0, 1)      # Priors for mean survival B
taueps <- 1/(sdeps*sdeps)
sdeps ~ dunif(0, 10)

# Define state-transition and observation matrices
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,t,1] <- phiA[t] * (1-psiAB[t])
      ps[1,t,2] <- phiA[t] * psiAB[t]
      ps[1,t,3] <- 1-phiA[t]
      ps[2,t,1] <- phiB[t] * psiBA[t]
      ps[2,t,2] <- phiB[t] * (1-psiBA[t])
      ps[2,t,3] <- 1-phiB[t]
      ps[3,t,1] <- 0
      ps[3,t,2] <- 0
      ps[3,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,t,1] <- pA[t]
      po[1,t,2] <- 0
      po[1,t,3] <- 1-pA[t]
      po[2,t,1] <- 0
      po[2,t,2] <- pB[t]
      po[2,t,3] <- 1-pB[t]
      po[3,t,1] <- 0
      po[3,t,2] <- 0
      po[3,t,3] <- 1
      } #t

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], t-1, 1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], t-1, 1:3])
      } #t
   } #i
}
")


# Compute vector with occasion of first capture
f <- apply(ch, 1, get.first)

# Recode ch matrix: note, a 0 is not allowed!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rch <- ch          # Recoded ch
rch[rch==0] <- 3

# Bundle data 
jags.data <- list(y = rch, f = f, n.occasions = ncol(rch), nind = nrow(rch))


# Initial values
inits <- function(){list(alpha = runif(ncol(rch) - 1, 0, 1),
                         mu = rnorm(1, 0, 1), 
                         sdeps = runif(1, 1, 5),
                         mean.psi = runif(2, 0, 1), 
                         mean.p = runif(2, 0, 1), 
                         z = ms.init.z(rch, f))}  

# Parameters monitored
parameters <- c("alpha", "mean.psi", "mean.p", "mu", "sdeps")

# MCMC settings
ni <- 2500
nt <- 1
nb <- 500
nc <- 2

# Call JAGS from R
ms <- jags(jags.data, 
           inits, 
           parameters, 
           "state-space-ms.jags", 
           n.chains = nc, 
           n.thin = nt, 
           n.iter = ni, 
           n.burnin = nb,
           parallel = TRUE)

# Inspect results
print(ms, digits = 3)

par(mfrow = c(4,3))
traceplot(ms)



###########################
# Exercise 2:
# The purpose of the exercise is to build models incorporating uncertainty in the assignment of states. 
# We use (historical) data collected between 1940 and 1957 by Lance Richdale on Sooty shearwaters Puffinus griseus (titis hereafter). 
# These data were reanalyzed using multistate capture-recapture models by Scofield et al. (2001) who kindly provided us with the data. 
# Following the way the data were collected, four states were originally considered:
# (1) Breeder (SB);
# (2) Accompanied by another bird in a burrow;
# (3) Alone in a burrow;
# (4) On the surface.
# Because of numerical issues, we pooled 2-3-4 together in a Non-Breeder state (NSB) that includes 
# failed breeders (birds that had bred previously – skip reproduction or divorce) and pre-breeders (birds that had yet to breed). 
# Note that because burrows were not checked before hatching, some birds in the category NSB might have already failed. 
# We therefore regard those birds in the SB state as successful breeders, and those in the NSB state as nonbreeders + prebreeders and failed breeders.
# We artificially generated uncertainty on the assignment of the two states Non-Successful-Breeder and Successful-Breeder as follows. 
# For each individual at each detection occasion: 1) a Non-Successful-Breeder was ascertained as Non-Successful-Breeder with probability 0.2 
# and assigned to obs 1, and not ascertained with the complementary probability 0.8 and assigned to obs 3, while 2) a Successful-Breeder 
# was ascertained as Successful-Breeder with probability 0.7 and assigned to obs 2, and not ascertained with the complementary probability 0.3 
# and assigned to obs 3.
# This procedure was implemented in R, see the script maketitisuncertain.R in the codes/ directory
# Try and fit a model with uncertainty in the state assignment. Consider 2 probabilities to ascertain the breeding status of an individual 
# encountered as Non-Successful-Breeder (deltaNB) or Successful-Breeder (deltaB).
###########################

#
# You may find the solution code below
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

# Read in the data: 
y <- read.table('dat/titis2.txt')
head(y)
dim(y)
nind <- dim(y)[1]
n.occasions <- dim(y)[2]

# Compute the date of first capture for each individual
f <- apply(y, 1, get.first)

# Specify model in BUGS language
cat(file = "multievent.jags", "
model {

   # OBSERVATIONS (+1)
   # 0 = non-detected
   # 1 = seen and ascertained as non-breeder
   # 2 = seen and ascertained as breeder
   # 3 = not ascertained
   
   # STATES
   # 1 = alive non-breeder
   # 2 = alive breeder
   # 3 = dead
   
   # PARAMETERS
   # phiNB  survival prob. of non-breeders
   # phiB  survival prob. of breeders
   # pNB  detection prob. of non-breeders
   # pB  detection prob. of breeders
   # psiNBB transition prob. from non-breeder to breeder
   # psiBNB transition prob. from breeder to non-breeder
   # piNB prob. of being in initial state non-breeder
   # deltaNB prob to ascertain the breeding status of an individual encountered as non-breeder
   # deltaB prob to ascertain the breeding status of an individual encountered as breeder
   
   # DEFINE PARAMETERS	
   # probabilities for each INITIAL STATES
   px0[1] <- piNB # prob. of being in initial state NB
   px0[2] <- 1 - piNB # prob. of being in initial state B
   px0[3] <- 0 # prob. of being in initial state dead
   
   # OBSERVATION PROCESS: probabilities of observations (columns) at a given occasion given states (rows) at this occasion
   # step 1: detection
   po1[1,1] <- 1 - pNB
   po1[1,2] <- pNB
   po1[1,3] <- 0
   po1[2,1] <- 1 - pB
   po1[2,2] <- 0
   po1[2,3] <- pB
   po1[3,1] <- 1
   po1[3,2] <- 0
   po1[3,3] <- 0
   
   po1.init[1,1] <- 0
   po1.init[1,2] <- 1
   po1.init[1,3] <- 0
   po1.init[2,1] <- 0
   po1.init[2,2] <- 0
   po1.init[2,3] <- 1
   po1.init[3,1] <- 1
   po1.init[3,2] <- 0
   po1.init[3,3] <- 0
   # step 2: assignement
   po2[1,1] <- 1
   po2[1,2] <- 0
   po2[1,3] <- 0
   po2[1,4] <- 0
   po2[2,1] <- 0
   po2[2,2] <- deltaNB
   po2[2,3] <- 0
   po2[2,4] <- 1 - deltaNB
   po2[3,1] <- 0
   po2[3,2] <- 0
   po2[3,3] <- deltaB
   po2[3,4] <- 1 - deltaB
   # form the matrix product
   po <- po1 %*% po2
   po.init <- po1.init %*% po2
   
   # STATE PROCESS: probabilities of states at t+1 (columns) given states at t (rows)
   # step 1: survival
   px1[1,1] <- phiNB
   px1[1,2] <- 0
   px1[1,3] <- 1 - phiNB
   px1[2,1] <- 0
   px1[2,2] <- phiB
   px1[2,3] <- 1 - phiB
   px1[3,1] <- 0
   px1[3,2] <- 0
   px1[3,3] <- 1
   # step 2: transition
   px2[1,1] <- 1 - psiNBB
   px2[1,2] <- psiNBB
   px2[1,3] <- 0
   px2[2,1] <- psiBNB
   px2[2,2] <- 1 - psiBNB
   px2[2,3] <- 0
   px2[3,1] <- 0
   px2[3,2] <- 0
   px2[3,3] <- 1
   # form the matrix product
   px <- px1 %*% px2
   
   for (i in 1:nind){ # for each ind
      # estimated probabilities of initial states are the proportions in each state at first capture occasion
      z[i,f[i]] ~ dcat(px0[1:3])
      y[i,f[i]] ~ dcat(po.init[z[i,f[i]],1:4])
      for (j in (f[i]+1):n.occasions){  # loop over time
         ## STATE EQUATIONS ##
         # draw states at j given states at j-1
         z[i,j] ~ dcat(px[z[i,j-1],1:3])
         ## OBSERVATION EQUATIONS ##
         # draw observations at j given states at j
         y[i,j] ~ dcat(po[z[i,j],1:4])
      }
   }
   # PRIORS 
   phiNB ~ dunif(0, 1)
   phiB ~ dunif(0, 1)
   pNB ~ dunif(0, 1)
   pB ~ dunif(0, 1)
   psiNBB ~ dunif(0, 1)
   psiBNB ~ dunif(0, 1)
   piNB ~ dunif(0, 1)
   deltaNB ~ dunif(0, 1)
   deltaB ~ dunif(0, 1)
}
")

# Form the list of data; remember we need to recode the data so that 
# 1 = for not detected, 
# 2 = for seen and ascertained as non-breeder
# 3 = seen and ascertained as breeder
# 4 = not ascertained
# to do that, we just add 1 to y

jags.data <- list(nind = nind, n.occasions = n.occasions, y = as.matrix(y + 1), f = f)

# Create initial values for unknown z
z <- y
# assign 2s to 3s
z[z == 3] <- 2
for (i in 1:nind) {
   for (j in 1:n.occasions) {
      if (j > f[i] & y[i,j]==0) {z[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
      if (j < f[i]) {z[i,j] <- NA}
   }
}
z1 <- as.matrix(z)
z <- y
# assign 1s to 3s
z [z==3] <- 1
for (i in 1:nind) {
   for (j in 1:n.occasions) {
      if (j > f[i] & y[i,j]==0) {z[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
      if (j < f[i]) {z[i,j] <- NA}
   }
}
z2 <- as.matrix(z)

# Now form the list of initial values:
init1 <- list(pB = 0.5, phiNB = 0.3, z = z1)
# second list of inits
init2 <- list(pB = 0.5, phiNB = 0.6, z = z2)
# concatenate list of initial values
inits <- list(init1, init2)

# Specify the parameters to be monitored
parameters <- c("phiB", "phiNB", "psiNBB", "psiBNB", "piNB", "pB", "pNB", "deltaNB", "deltaB")

# MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 2

# Tadaaaaaaan, call JAGS from R:
ms <- jags(jags.data, 
           inits, 
           parameters, 
           "multievent.jags", 
           n.chains = nc, 
           n.thin = nt, 
           n.iter = ni, 
           n.burnin = nb, 
           parallel = TRUE)

# Inspect results
par(mfrow = c(3,4))
traceplot(ms)

# Visualize results
library(lattice)
densityplot(ms)

# Print results
print(ms, digits = 3)



# These results may be compared with the results obtained using E-SURGE 
# (Table 1 in Gimenez et al. 2012):
# deltaB | 0.737576452 0.616347581 0.831002071 0.055236538 
# deltaNB | 0.187797018 0.162175744 0.216420470 0.013831853 
# pB | 0.597718828 0.532890645 0.659302061 0.032413822 
# pNB | 0.564643587 0.508403711 0.619267746 0.028396321 
# phiB | 0.837489708 0.796415345 0.871612997 0.019139418 
# phiNB | 0.814037056 0.779702384 0.844090508 0.016414429 
# piNB | 0.704217686 0.646248176 0.756270393 0.028149150 
# psiBNB | 0.226471935 0.144866984 0.335984429 0.048899140 
# psiNBB | 0.219402142 0.173907703 0.272866821 0.025255272 
