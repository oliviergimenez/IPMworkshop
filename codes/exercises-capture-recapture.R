###############
### EXERCISE 1: Fit a model with constant survival and detection probabilities to data simulated with constant parameters. 
# Fit a model with time-dependent survival and detection probabilities (CJS) to the same data.
###############

###############
### EXERCISE 2: We are going to analyse real capture-recapture data.
# The data concern the European Dipper (Cinclus cinclus). 
# Captures were carried out for 7 years (1981-1987) in eastern France by G. Marzolin who kindly provided us with the data. 
# They consist of initial markings and recaptures of 293 breeding adults each year during the March-June period. 
# Birds were at least 1 year old when initially banded. 
# Estimate sex-specific survival, assuming constant capture probability.
# Add temporal variance on top of survival in the model.
# To get you started, we read in the data and clean up the data. 
###########################

dipper <- read.csv("dat/dipper.csv", sep = ";")
head(dipper)
# Create a sex variable with 1 = female and 2 = male
dipper$sex <- ifelse(dipper$males == 1, 2, 1)
head(dipper)
# Delete the redundant columns
dipper <- dipper[, -c(8,9)]
dipper <- as.matrix(dipper)
head(dipper)

# Specify model in BUGS language
cat(file = "cjs-sex-temp.txt", "
model {
    
    # Constraints
    for (i in 1:nind){
       for (t in 1:(n.occasions-1)){
          logit(phi[i,t]) <- mu[sex[i]] + epsilon[t]
          p[i,t] <- mean.p
          } #t
       } #i

    for (t in 1:(n.occasions-1)){
       epsilon[t] ~ dnorm(0, tau)
       }

    # Priors
    mean.p ~ dunif(0, 1)                     # Prior for mean recapture
    mu[1] ~ dnorm(0, 1)                      # Prior for mean female survival 
    mu[2] ~ dnorm(0, 1)                      # Prior for mean male survival 
    sigma ~ dunif(0, 5)                      # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    
    # Likelihood 
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,f[i]] <- 1
      for (t in (f[i]+1):n.occasions){
        # State process
        z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
        # Observation process
        y[i,t] ~ dbern(p[i,t-1] * z[i,t])
      } #t
    } #i
}
")

# Data
f <- apply(dipper[, 1:7], 1, get.first)
jags.data <- list(y = dipper[,1:7], 
                  f = f, 
                  nind = dim(dipper)[1], 
                  n.occasions = dim(dipper[,1:7])[2], 
                  sex = dipper[, 8])

# Initial values
inits <- function(){list(z = z.inits(dipper[,1:7]), 
                         mu = rnorm(2, 0, 1), 
                         sigma = runif(1, 0, 5), 
                         mean.p = runif(1, 0, 1))}  

# Parameters to be monitored
parameters <- c("mu", "mean.p", "sigma2")

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R
cjs.dipper <- jags(data = jags.data, 
                   inits = inits, 
                   parameters.to.save = parameters, 
                   model.file = "cjs-sex-temp.txt", 
                   n.chains = nc, 
                   n.thin = nt, 
                   n.iter = ni, 
                   n.burnin = nb, 
                   parallel = TRUE)

# Summarize posteriors
print(cjs.dipper, digits = 3)

# Produce histogram for temporal variance, with posterior mean in red
hist(cjs.dipper$sims.list$sigma2, 
     col = "gray", 
     nclass = 35, 
     las = 1, 
     xlab = expression(sigma^2), 
     main = "")
abline(v = mean(cjs.dipper$sims.list$sigma2), 
       col = "red", 
       lwd = 2)

# Estimated survival for females and males, respectively
par(mfrow = c(1, 2))
plot(density(plogis(cjs.dipper$sims.list$mu[,1])), 
     xlab = "", 
     ylab = "", 
     main = "female survival")
plot(density(plogis(cjs.dipper$sims.list$mu[,2])), 
     xlab = "", 
     ylab = "", 
     main = "male survival")


###############
### EXERCISE 3: Write a multistate model w/ two states, 
# in such a way that survival in A is time-dependent (fixed time effect) 
# and survival at site B it is time-dependent with a random temporal effect.
###########################

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
jags.data <- list(y = rch, 
                  f = f, 
                  n.occasions = ncol(rch), 
                  nind = nrow(rch))


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
ms <- jags(data = jags.data, 
           inits = inits, 
           parameters.to.save = parameters, 
           model.file = "state-space-ms.jags", 
           n.chains = nc, 
           n.thin = nt, 
           n.iter = ni, 
           n.burnin = nb,
           parallel = TRUE)

# Inspect results
print(ms, digits = 3)

par(mfrow = c(4,3))
traceplot(ms)


###############
### EXERCISE 4: The purpose of the exercise is to build models incorporating uncertainty in the assignment of states. 
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
# This procedure was implemented in R, see the script make-uncertain.R in the codes/ directory
# Try and fit a model with uncertainty in the state assignment. Consider 2 probabilities to ascertain the breeding status of an individual 
# encountered as Non-Successful-Breeder (deltaNB) or Successful-Breeder (deltaB).
###########################

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

jags.data <- list(nind = nind, 
                  n.occasions = n.occasions, 
                  y = as.matrix(y + 1), f = f)

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
parameters <- c("phiB", 
                "phiNB", 
                "psiNBB", 
                "psiBNB", 
                "piNB", 
                "pB", 
                "pNB", 
                "deltaNB", 
                "deltaB")

# MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 2

# Tadaaaaaaan, call JAGS from R:
ms <- jags(data = jags.data, 
           inits = inits, 
           parameters.to.save = parameters, 
           model.file = "multievent.jags", 
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
