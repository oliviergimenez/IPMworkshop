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


# Specify path
setwd('...')

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
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   }
for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   }

# Define state-transition and observation matrices
for (i in 1:nind){  
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA[t] * (1-psiAB[t])
      ps[1,i,t,2] <- phiA[t] * psiAB[t]
      ps[1,i,t,3] <- 1-phiA[t]
      ps[2,i,t,1] <- phiB[t] * psiBA[t]
      ps[2,i,t,2] <- phiB[t] * (1-psiBA[t])
      ps[2,i,t,3] <- 1-phiB[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-pA[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pB[t]
      po[2,i,t,3] <- 1-pB[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
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
rch <- ch          # Recoded ch
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

inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2, 0, 1), mean.p = runif(2, 0, 1), z = ms.init.z(rch, f))}  

# 1.4. Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

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
ni <- 4000
nt <- 1
nb <- 2000
nc <- 3

# 1.6. Call JAGS from R
ms <- jags(jags.data, inits, parameters, "state-space-ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# 1.7. Inspect results
par(mfrow = c(3,3))
traceplot(ms)

print(ms, digits = 3)


JAGS output for model 'state-space-ms.jags', generated by jagsUI.
Estimates based on 3 chains of 4000 iterations,
burn-in = 2000 iterations and thin rate = 1,
yielding 6000 total samples from the joint posterior. 
MCMC ran for 15.378 minutes at time 2018-10-24 03:39:14.

                mean      sd    2.5%      50%    97.5% overlap0 f  Rhat n.eff
mean.phi[1]    0.777   0.024   0.735    0.774    0.830    FALSE 1 1.030    86
mean.phi[2]    0.678   0.032   0.614    0.679    0.739    FALSE 1 1.098    25
mean.psi[1]    0.309   0.087   0.144    0.301    0.482    FALSE 1 1.216    14
mean.psi[2]    0.572   0.092   0.414    0.563    0.785    FALSE 1 1.270    12
mean.p[1]      0.763   0.094   0.601    0.752    0.955    FALSE 1 1.205    14
mean.p[2]      0.374   0.131   0.212    0.346    0.749    FALSE 1 1.348    12
deviance    1441.412 199.325 950.499 1492.095 1705.852    FALSE 1 1.138    21

**WARNING** Rhat values indicate convergence failure. 
Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
For each parameter, n.eff is a crude measure of effective sample size. 

overlap0 checks if 0 falls in the parameter's 95% credible interval.
f is the proportion of the posterior with the same sign as the mean;
i.e., our confidence that the parameter is positive or negative.

DIC info: (pD = var(deviance)/2) 
pD = 17874.8 and DIC = 19316.25 
DIC is an estimate of expected predictive error (lower is better).



###########################
# Exercise 1:
# Write the model in such a way that survival in A is time-dependent (fixed time effect) and survival at site B it is time-dependent with a random temporal effect.
###########################






###############################################

# 2. Approach: multinomial likelihood

# 2.1. Specify model in BUGS language
cat(file = "multinom-ms.jags", "
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

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   }

for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   }

# Define state-transition and reencounter probabilities
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- phiA[t] * (1-psiAB[t])
      psi[1,t,2] <- phiA[t] * psiAB[t]
      psi[2,t,1] <- phiB[t] * psiBA[t]
      psi[2,t,2] <- phiB[t] * (1-psiBA[t])

      po[1,t] <- pA[t]
      po[2,t] <- pB[t]


      # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities      
      for (s in 1:ns){
         dp[s,t,s] <- po[s,t]
         dq[s,t,s] <- 1-po[s,t]
         } # s
      for (s in 1:(ns-1)){
         for (m in (s+1):ns){
            dp[s,t,m] <- 0
            dq[s,t,m] <- 0
            } # s
         } # m
      for (s in 2:ns){
         for (m in 1:(s-1)){
            dp[s,t,m] <- 0
            dq[s,t,m] <- 0
            } # s
         } # m
      } # t

# Define the multinomial likelihood
for (t in 1:((n.occasions-1)*ns)){
   marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
   }

# Define the cell probabilities of the multistate m-array   
# Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
for (t in 1:(n.occasions-2)){
   U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
   for (j in (t+1):(n.occasions-1)){
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
      }
   }
U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
# Diagonal
for (t in 1:(n.occasions-2)){
   pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
      }
   }
pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]

# Below main diagonal
for (t in 2:(n.occasions-1)){
   for (j in 1:(t-1)){
      pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
      } #j
   } #t

# Last column: probability of non-recapture
for (t in 1:((n.occasions-1)*ns)){
   pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
   } #t
}
")




# 2.2. Data preparation
# Function to create a single or multistate m-array from capture-recapture data
# Input variables
#    ch: matrix with single- or multistate capture histories (0: not captured; 1..X: captured in the 1..X states)
#    unobs: number of unobserved states (default is 0, needs to be given only in specific cases)
#
# Output
#    out: single- or multistate m-array. The last column of each m-array is the number of released individuals that were never recaptured. Thus, the total number of released individuals per occasion and state is the row sum of the m-array.

marray <- function(ch, unobs = 0){
   ns <- length(table(ch)) - 1 + unobs
   no <- ncol(ch)
   out <- matrix(0, ncol = ns*(no-1)+1, nrow = ns*(no-1))
   # Remove capture histories of individuals that are marked at last occasion
   get.first <- function(x) min(which(x!=0))
   first <- apply(ch, 1, get.first)
   last.only <- which(first==no)
   if (length(last.only) > 0) ch <- ch[-last.only,]
   # Compute m-array
   for (i in 1:nrow(ch)){
      cap.occ <- which(ch[i,]!=0)
      state <- ch[i,cap.occ]
      if (length(state) == 1) {
         out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] <- out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] + 1
         }
      if (length(state) > 1) {
         for (t in 2:length(cap.occ)){
            out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] <- out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] + 1
            } # t
         if (max(cap.occ) < no){
            out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] <- out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] + 1
            } # if
         } # if
      } # t
   return(out)
   }    

# Create multistate m-array
ms.arr <- marray(ch)

# Calculate the number of states
ns <- length(unique(as.numeric(ch))) - 1

# Bundle data
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns))


# 2.3. Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1))}  

# 2.4. Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# 2.5. MCMC settings
ni <- 10000; nt <- 1; nb <- 2000; nc <- 3

# 2.6. Call JAGS from R
m.ms <- jags(jags.data, inits, parameters, "multinom-ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


# 2.7. Inspect results
par(mfrow = c(3,3))
traceplot(m.ms)

print(m.ms, digits = 3)



#################################

