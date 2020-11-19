################################################
#
# Capture-recapture models (specifically: Cormack-Jolly-Seber models)
#
# IPM workshop, virtual Montpellier November 2020
#
################################################
# Based on the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" 
# by Marc Kéry & Michael Schaub (2012, Academic Press)
#########################################################################

# Michael Schaub, 24.1.2012, revised 9.7.2012, revised 10.12.2013, revised 12.9.2017, revised 15.10.19
# Olivier Gimenez 19.11.2020

#########################################################################

# Make sure that you have installed JAGS version 3.2.0
library(jagsUI)  # Load package

#########################################################################

setwd("C:/....")     # set working directory (Optional)


####################################################################
#
# 7. Estimation of survival probabilities using capture-recapture data
#
#####################################################################

#####################################################################
#### 7.3. Models with constant parameters
#####################################################################

# Define parameter values
n.occasions <- 6                   # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions-1)
p <- rep(0.4, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break		# If dead, move to next individual 
      # Bernoulli trial: is individual recaptured? 
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
      } #t
    } #i
  return(CH)
  }

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
cat(file = "cjs-c-c.txt", "
model {

    # Constraints
    for (i in 1:nind){
      for (t in 1:(n.occasions-1)){
        phi[i,t] <- mean.phi
        p[i,t] <- mean.p
        } #t
      } #i
    
    # Priors
    mean.phi ~ dunif(0, 1)         # Prior for mean survival
    mean.p ~ dunif(0, 1)           # Prior for mean recapture
    
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH))

# Initial values
# It is vital to specify "good" initial values for the latent state variable z. Here is one possibility that works well.
z.inits <- function(ch){
   state <- ch
   state[state==0] <- 1
   get.first <- function(x) min(which(x!=0))
   f <- apply(ch, 1, get.first)
   for (i in 1:nrow(ch)){
      state[i,1:f[i]] <- NA
      }
   return(state)
   }

inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = z.inits(CH))}


# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs.c.c <- jags(jags.data, inits, parameters, "cjs-c-c.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(cjs.c.c, digits = 3)

# Compare the posterior means and medians of survival and detection 
# to the values we used to simulate data. 

# Note: You may also get posterior distribution for the latent 
# parameters z. 


#####################################################################
#### 7.4. Models with time-variation
#####################################################################

### 7.4.1. Fixed time effects

cat(file = "cjs.txt", "
model {
    # Constraints
    for (i in 1:nind){
      for (t in 1:(n.occasions-1)){
        phi[i,t] <- alpha[t]
        p[i,t] <- beta[t]
        } #t
      } #i

    # Priors
    for (t in 1:(n.occasions-1)){
      alpha[t] ~ dunif(0, 1)        # Priors for time-spec. recapture
      beta[t] ~ dunif(0, 1)         # Priors for time-spec. recapture
      }

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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH))

# Initial values
inits <- function(){list(z = z.inits(CH))}  

# Parameters monitored
parameters <- c("alpha", "beta")

# MCMC settings
ni <- 5000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R
cjs.t.t <- jags(jags.data, inits, parameters, "cjs.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.t.t, digits = 3)

# Again, compare the estimates with the values we used to simulate the data. 


###############
### EXERCISE 1: RECREATE THIS ANALYSIS WITH THE SIMPLEST BUGS MODEL STRUCTURE
###############



### UNDERSTANDING PARAMETER REDUNDANCY

# Look at the parameter estimates: do you see anything conspicuous?
cjs.t.t$summary


# Plot the posteriors
par(mfrow = c(3,5))
for (i in 1:5){
   plot(density(cjs.t.t$sims.list$alpha[,i]), main = paste("Occasion ", i), ylim = c(0,5))
   abline(h = 1, lty = 2)
   }
for (i in 1:5){
   plot(density(cjs.t.t$sims.list$beta[,i]), main = paste("Occasion ", i), ylim = c(0,5))
   abline(h = 1, lty = 2)
   }
for (i in 1:5){
   plot(cjs.t.t$sims.list$alpha[,i], cjs.t.t$sims.list$beta[,i], main = paste("Occasion ", i), pch = 16, cex = 0.5)
   }




### 7.4.2. Random time effects

# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
sigma2.phi <- 1                    # Temporal variance of survival
p <- rep(0.4, n.occasions-1)

# Determine annual survival probabilities (on logistic scale)
logit.phi <- rnorm(n.occasions-1, qlogis(mean.phi), sigma2.phi^0.5)
logit.phi
phi <- plogis(logit.phi)
phi

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
cat(file = "cjs-temp.txt", "
model {
    
    # Constraints
    for (i in 1:nind){
       for (t in 1:(n.occasions-1)){
          logit(phi[i,t]) <- mu + epsilon[t]
          p[i,t] <- mean.p
          } #t
       } #i

    for (t in 1:(n.occasions-1)){
       epsilon[t] ~ dnorm(0, tau)
       
       # calculate annual estimates of survival (optional)
       phi.est[t] <- ilogit(mu + epsilon[t])
       } #t
    
    # Priors
    mean.p ~ dunif(0, 1)                     # Prior for mean recapture
    mean.phi ~ dunif(0, 1)                   # Prior for mean survival
    mu <- logit(mean.phi)                    # Logit transformation
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2])

# Initial values
inits <- function(){list(z = z.inits(CH), mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "phi.est", "mean.p", "sigma2")

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (~ 1 min)
cjs.ran <- jags(jags.data, inits, parameters, "cjs-temp.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.ran, digits = 3)

# You know the drill: take a minute or two to compare your estimates
# with the values we used to simulate the data

# Produce histogram
hist(cjs.ran$sims.list$sigma2, col = "gray", nclass=35, las=1, xlab = expression(sigma^2), main = "")
abline(v = sigma2.phi, col = "red", lwd = 2)



### 7.4.3. Temporal covariates

# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
beta <- -0.3                       # Negative effect of winter severity (regression coeff)
r.sigma2 <- 0.2                    # Residual temporal variance

# Draw annual survival probabilities
winter <- rnorm(n.occasions-1, 0, 1^0.5)
logit.phi <- qlogis(mean.phi) + beta * winter + rnorm(n.occasions-1, 0, r.sigma2^0.5)
phi <- plogis(logit.phi)

# Look at the survival - winter relationship
plot(winter, logit.phi, ylab = 'survival (logit scale)')
abline(a = qlogis(mean.phi), b = beta, col = 'blue', lwd = 2)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
cat(file = "cjs-cov.txt", "
model {
    
    # Constraints
    for (i in 1:nind){
       for (t in 1:(n.occasions-1)){
          logit(phi[i,t]) <- mu + beta * x[t] + epsilon[t]
          p[i,t] <- mean.p
          } #t
       } #i
    for (t in 1:(n.occasions-1)){
       epsilon[t] ~ dnorm(0, tau)
       }

    # Derived parameters (could be calculated offline from the mcmc chains)
    for (t in 1:(n.occasions-1)){
       phi.est[t] <- ilogit(mu + beta * x[t] + epsilon[t]) # Yearly survival
       }

    # Priors
    mu ~ dnorm(0, 1)                     # Prior for logit of mean survival
    mean.phi <- ilogit(mu)                   # Logit transformation
    beta ~ dunif(-10, 10)                   # Prior for slope parameter
    sigma ~ dunif(0, 10)                     # Prior on standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Residual temporal variance
    mean.p ~ dunif(0, 1)                     # Prior for mean recapture
    
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH), x = winter)

# Initial values
inits <- function(){list(z = z.inits(CH), mu = rnorm(1), sigma = runif(1, 0, 5), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "phi.est", "mu", "sigma2", "beta")

# MCMC settings
ni <- 5000
nt <- 1
nb <- 2500
nc <- 3

# Call JAGS from R (~ 3 min)
cjs.cov <- jags(jags.data, inits, parameters, "cjs-cov.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.cov, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.cov$sims.list$beta, nclass = 25, col = "gray", main = "", xlab = expression(beta), ylab = "Frequency")
abline(v = -0.3, col = "red", lwd = 2)
hist(cjs.cov$sims.list$sigma2, nclass = 50, col = "gray", main = "", xlab = expression(sigma^2), ylab = "Frequency", xlim = c(0, 3))
abline(v = 0.2, col = "red", lwd = 2)

# Look at the survival - winter relationship
# True in blue, estimated in green
plot(winter, logit.phi, ylab = 'survival (logit scale)')
abline(a = plogis(mean.phi), b = beta, col = 'blue', lwd = 2)
abline(a = mean(cjs.cov$sims.list$mu), 
       b = mean(cjs.cov$sims.list$beta), col = 'green', lwd = 2)


#####################################################################
# 7.5. Models with individual variation
#####################################################################

### 7.5.1. Fixed group effects

# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
phi.f <- rep(0.65, n.occasions-1)  # Survival of females
p.f <- rep(0.6, n.occasions-1)     # Recapture of females
phi.m <- rep(0.8, n.occasions-1)   # Survival of males
p.m <- rep(0.3, n.occasions-1)     # Reacpture of males

# Define matrices with survival and recapture probabilities
PHI.F <- matrix(phi.f, ncol = n.occasions-1, nrow = sum(marked))
P.F <- matrix(p.f, ncol = n.occasions-1, nrow = sum(marked))
PHI.M <- matrix(phi.m, ncol = n.occasions-1, nrow = sum(marked))
P.M <- matrix(p.m, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

# Merge capture-histories by row
CH <- rbind(CH.F, CH.M)

# Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in JAGS language
cat(file = "cjs-group.txt", "
model {
    
    # Constraints
    for (i in 1:nind){
       for (t in 1:(n.occasions-1)){
          phi[i,t] <- phi.g[group[i]]
          p[i,t] <- beta[t]
          } #t
       } #i

    # Priors
    for (u in 1:g){
       phi.g[u] ~ dunif(0, 1)        # Priors for group-specific survival
       }
    for (t in 1:(n.occasions-1)){
       beta[t] ~ dunif(0, 1)
       }
    
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = z.inits(CH), phi.g = runif(length(unique(group)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g", "beta")

# MCMC settings
ni <- 5000
nt <- 1
nb <- 2500
nc <- 3

# Call JAGS from R (~ 2 min)
cjs.group <- jags(jags.data, inits, parameters, "cjs-group.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.group, digits = 3)

# The estimates are to be compared with the values we used to simulate the data.
# Note that detection is time-dependent, and fully identifiable. 


### 7.5.3. Individual random effects

# Define parameter values
n.occasions <- 12                 # Number of capture occasions
marked <- rep(50, n.occasions-1)  # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
v.ind <- 0.5

# Draw annual survival probabilities
logit.phi <- rnorm(sum(marked), qlogis(mean.phi), v.ind^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = FALSE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
cat(file = "cjs-ind.txt", "
model {
    
    # Constraints
    for (i in 1:nind){
       for (t in 1:(n.occasions-1)){
          logit(phi[i,t]) <- mu + epsilon[i]
          p[i,t] <- mean.p
          } #t
       } #i
    for (i in 1:nind){
       epsilon[i] ~ dnorm(0, tau)   
       }

    # Priors
    mean.phi ~ dunif(0, 1)                   # Prior for mean survival
    mu <- logit(mean.phi)                    # Logit transformation
    sigma ~ dunif(0, 5)                      # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)
    mean.p ~ dunif(0, 1)                     # Prior for mean recapture 
    
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH))

# Initial values 
inits <- function(){list(z = z.inits(CH), mean.phi = runif(1, 0, 1), 
                         mean.p = runif(1, 0, 1), sigma = runif(1, 0, 2))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "sigma2")

# MCMC settings
ni <- 20000
nt <- 2
nb <- 10000
nc <- 3

# since this analysis takes a long time, we use a short run (1 chain) JUST as a demo
# (this is most probably too short for a real analysis)
ni <- 2000
nt <- 1
nb <- 1000
nc <- 1

# Call JAGS from R (NOTE: takes a long time to run!)
cjs.ind <- jags(jags.data, inits, parameters, "cjs-ind.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.ind, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.ind$sims.list$mean.phi, nclass = 25, col = "gray", main = "", xlab = expression(bar(phi)), ylab = "Frequency")
abline(v = mean.phi, col = "red", lwd = 2)
hist(cjs.ind$sims.list$sigma2, nclass = 15, col = "gray", main = "", xlab = expression(sigma^2), ylab = "Frequency", xlim = c(0, 3))
abline(v = v.ind, col = "red", lwd = 2)





#####################################################################
# 7.6. Models with time and group effects
#####################################################################

### 7.6.1. Fixed group and time effects (additive)

# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi.f <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p.f <- rep(0.6, n.occasions-1)
diff <- 0.5     # Difference between male and female survival on logit scale
phi.m <- plogis(qlogis(phi.f) + diff)
p.m <- rep(0.3, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.F <- matrix(rep(phi.f, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P.F <- matrix(rep(p.f, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
PHI.M <- matrix(rep(phi.m, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P.M <- matrix(rep(p.m, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)

# Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

# Merge capture-histories
CH <- rbind(CH.F, CH.M)

# Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
cat(file = "cjs-add.txt", "
model {
    
    # Constraints
    for (i in 1:nind){
       for (t in 1:(n.occasions-1)){
          logit(phi[i,t]) <- beta[group[i]] + gamma[t]
          p[i,t] <- p.g[group[i]]
          } #t
       } #i

    # Priors (and derived parameters)
    # ...for survival parameters
    for (t in 1:(n.occasions-1)){
       gamma[t] ~ dunif(-10, 10)                     # Priors for time effects
       phi.g1[t] <- ilogit(gamma[t])                 # Back-transformed survival of males
       phi.g2[t] <- ilogit(gamma[t] + beta[2])       # Back-transformed survival of females 
       }
    beta[1] <- 0                            # Corner constraint
    beta[2] ~ dunif(-10, 10)                # Prior for difference in male and female survival
    # ...for recapture parameters
    for (u in 1:g){
       p.g[u] ~ dunif(0, 1)                  # Priors for group-spec. recapture
       }
    
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = z.inits(CH), gamma = rnorm(n.occasions-1), beta = c(NA, rnorm(1)), p.g = runif(length(unique(group)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g1", "phi.g2", "p.g", "beta")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (~ 2.5 min)
cjs.add <- jags(jags.data, inits, parameters, "cjs-add.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.add, digits = 3)   


# Figure of male and female survival
lower.f <- upper.f <- lower.m <- upper.m <- numeric()
for (t in 1:(n.occasions-1)){
   lower.f[t] <- quantile(cjs.add$sims.list$phi.g1[,t], 0.025)
   upper.f[t] <- quantile(cjs.add$sims.list$phi.g1[,t], 0.975)
   lower.m[t] <- quantile(cjs.add$sims.list$phi.g2[,t], 0.025)
   upper.m[t] <- quantile(cjs.add$sims.list$phi.g2[,t], 0.975)
   }
plot(x = (1:(n.occasions-1))-0.1, y = cjs.add$mean$phi.g1, type = "b", pch = 16, ylim = c(0.2, 1), ylab = "Survival probability", xlab = "Year", bty = "n", cex = 1.5, axes = FALSE)
axis(1, at = 1:11, labels = rep(NA,11), tcl = -0.25)
axis(1, at = seq(2,10,2), labels = c("2", "4", "6", "8", "10"))
axis(2, at = seq(0.2, 1, 0.1), labels = c("0.2", NA, "0.4", NA, "0.6", NA, "0.8", NA, "1.0"), las = 1)
segments((1:(n.occasions-1))-0.1, lower.f, (1:(n.occasions-1))-0.1, upper.f)
points(x = (1:(n.occasions-1))+0.1, y = cjs.add$mean$phi.g2, type = "b", pch = 1, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))+0.1, lower.m, (1:(n.occasions-1))+0.1, upper.m)



# Code for model with interaction between sex and time

## [change these lines in BUGS code]

# Constraints
for (i in 1:nind){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- eta.phi[group[i],t]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
# Priors 
# ...for survival parameters
for (u in 1:g){
   for (t in 1:(n.occasions-1)){
      eta.phi[u,t] ~ dunif(0, 1)     # Prior for time and group-spec. survival
      } #t
   } #g
# ...for recapture parameters
for (u in 1:g){
   p.g[u] ~ dunif(0, 1)              # Priors for group-spec. recapture
   }


### 7.6.2. Fixed group and random time effects

## [change these lines in BUGS code]

# Constraints
for (i in 1:nind){
   for (t in 1:(n.occasions-1)){
      logit(phi[i,t]) <- eta.phi[group[i],t]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
# survival
for (u in 1:g){
   for (t in 1:(n.occasions-1)){
      eta.phi[u,t] <- mu.phi[u] + epsilon[u,t]
      epsilon[u,t] ~ dnorm(0, tau[u])
      } #t
   } #g 
  
# Priors
# ...for survival parameters
for (u in 1:g){
   for (t in 1:(n.occasions-1)){
      mean.phi[u] ~ dunif(0, 1)      # Priors on mean group-spec. survival     
      mu.phi[u] <- logit(mean.phi[u])
      sigma[u] ~ dunif(0, 10)        # Priors for group-spec. sd 
      tau[u] <- pow(sigma[u], -2)
      sigma2[u] <- pow(sigma[u], 2)
      } #t
   } #g
# ...for recapture parameters
for (u in 1:g){
   p.g[u] ~ dunif(0,1)               # Priors for group-spec. recapture
   }


# Version of the model with correlated random effects between groups

# Specify model in BUGS language
cat(file = "cjs-temp-corr.txt", "
model {
    
    # Priors and constraints
    for (i in 1:nind){
       for (t in 1:(n.occasions-1)){
          logit(phi[i,t]) <- eta.phi[t,group[i]]
          p[i,t] <- p.g[group[i]]
          } #t
       } #i
    # for survival parameters
    for (t in 1:(n.occasions-1)){
       eta.phi[t,1:g] ~ dmnorm(mu.phi[], Omega[,])
       } #t

    for (u in 1:g){      
       mean.phi[u] ~ dunif(0, 1)    # Priors on mean group-spec. survival
       mu.phi[u] <- logit(mean.phi[u])
       } #g
    Omega[1:g, 1:g] ~ dwish(R[,], df)  # Priors for variance-covariance matrix
    Sigma[1:g, 1:g] <- inverse(Omega[,])
    
    # for recapture parameters
    for (u in 1:g){
       p.g[u] ~ dunif(0, 1)            # Priors for group-spec. recapture
       }
    
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH), 
                  g = length(unique(group)), group = group, 
                  R = matrix(c(5, 0, 0, 1), ncol = 2), df = 3)

# Initial values
inits <- function(){list(z = z.inits(CH), p.g = runif(length(unique(group)), 0, 1), 
                         Omega = matrix(c(1, 0, 0, 1), ncol = 2))}  

# Parameters monitored
parameters <- c("eta.phi", "p.g", "Sigma", "mean.phi")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 5 min)
cjs.corr <- jags(jags.data, inits, parameters, "cjs-temp-corr.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = F)

# Summarize posteriors
print(cjs.corr, digits = 3)   



#####################################################################
# 7.7. Models with age effects
#####################################################################

# Define parameter values
n.occasions <- 10                   # Number of capture occasions
marked.j <- rep(200, n.occasions-1) # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)  # Annual number of newly marked adults
phi.juv <- 0.3                      # Juvenile annual survival
phi.ad <- 0.65                      # Adult annual survival
p <- rep(0.5, n.occasions-1)        # Recapture
phi.j <- c(phi.juv, rep(phi.ad, n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
   PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
   }
P.J <- matrix(rep(p, sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)  # individuals marked as juveniles
CH.A <- simul.cjs(PHI.A, P.A, marked.a)  # individuals marked as adults 

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)
f.a <- apply(CH.A, 1, get.first)

# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
x.a <- matrix(NA, ncol = dim(CH.A)[2]-1, nrow = dim(CH.A)[1])
for (i in 1:nrow(CH.J)){
   for (t in f.j[i]:(ncol(CH.J)-1)){
      x.j[i,t] <- 2
      x.j[i,f.j[i]] <- 1   
      } #t
   } #i
for (i in 1:nrow(CH.A)){
   for (t in f.a[i]:(ncol(CH.A)-1)){
      x.a[i,t] <- 2
      } #t
   } #i

CH <- rbind(CH.J, CH.A)
f <- c(f.j, f.a)
x <- rbind(x.j, x.a)

# Specify model in BUGS language
cat(file = "cjs-age.txt", "
model {
    # Constraints
    for (i in 1:nind){
       for (t in f[i]:(n.occasions-1)){
          phi[i,t] <- beta[x[i,t]]
          p[i,t] <- mean.p
          } #t
       } #i

    # Priors
    for (u in 1:2){
       beta[u] ~ dunif(0, 1)      # Priors for age-specific survival
       }
    mean.p ~ dunif(0, 1)         # Prior for mean recapture
    
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

# Bundle data
jags.data <- list(y = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH), x = x)

# Initial values
inits <- function(){list(z = z.inits(CH), beta = runif(2, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 200
nt <- 3
nb <- 100
nc <- 3

# Call JAGS from R (~ 2.5 min)
cjs.age <- jags(jags.data, inits, parameters, "cjs-age.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(cjs.age, digits = 3)



## Version of the code for modelling survival as a linear function of age

# Create matrix X indicating age classes
x <- matrix(NA, ncol = ncol(CH)-1, nrow = nrow(CH))
for (i in 1:nrow(CH)){
   for (t in f[i]:(ncol(CH)-1)){
      x[i,t] <- t-f[i]+1
      } #t 
   } #i

## [change these lines in BUGS code]

# Constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[i,t]
      p[i,t] <- mean.p
      } #t
   } #i

# Priors (and derived parameters)
mu ~ dnorm(0, 0.01)             # Prior for mean of logit survival
beta ~ dnorm(0, 0.01)           # Prior for slope parameter
for (i in 1:(n.occasions-1)){
   phi.age[i] <- ilogit(mu + beta*i)   # Logit back-transformation 
   }
mean.p ~ dunif(0, 1)                # Prior for mean recapture





###############
### EXERCISE 2: Analyse real CR data of European dipper
###############

# data from 293 marked birds with 7 capture occasions 
dipper <- read.csv("dat/dipper.csv", sep = ";")
head(dipper)
dipper$sex <- ifelse(dipper$males == 1, 2, 1)
head(dipper)
dipper <- dipper[, -c(8,9)]
dipper <- as.matrix(dipper)
head(dipper)

# Estimate sex-specific survival, assuming constant capture probability
# Add temporal variance on top of survival in the model






























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

# Bundle data
f <- apply(dipper[, 1:7], 1, get.first)
jags.data <- list(y = dipper[,1:7], f = f, nind = dim(dipper)[1], n.occasions = dim(dipper[,1:7])[2], sex = dipper[, 8])

# Initial values
inits <- function(){list(z = z.inits(dipper[,1:7]), mu = rnorm(2, 0, 1), sigma = runif(1, 0, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mu", "mean.p", "sigma2")

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (~ 1 min)
cjs.dipper <- jags(jags.data, inits, parameters, "cjs-sex-temp.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.dipper, digits = 3)


# Produce histogram for temporal variance, with posterior mean in red
hist(cjs.dipper$sims.list$sigma2, col = "gray", nclass=35, las=1, xlab = expression(sigma^2), main = "")
abline(v = mean(cjs.dipper$sims.list$sigma2), col = "red", lwd = 2)

# Estimated survival for females and males, respectively
par(mfrow = c(1, 2))
plot(density(plogis(cjs.dipper$sims.list$mu[,1])), xlab = "", ylab = "", main = "female survival")
plot(density(plogis(cjs.dipper$sims.list$mu[,2])), xlab = "", ylab = "", main = "male survival")




#####################################################################
# 7.10. Fitting the CJS to data in the m-array format: the multinomial likelihood
#####################################################################

### 7.10.1. Introduction

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- nrow(CH)
  n.occasions <- ncol(CH)
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
     m.array[t,1] <- sum(CH[,t])
     }
  for (i in 1:nind){
     pos <- which(CH[i,]!=0)
     g <- length(pos)
     for (z in 1:(g-1)){
        m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
        } #z
     } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
     m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
     }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}


### 7.10.2. Time-dependent models

# Specify model in BUGS language
cat(file = "cjs-mnl.txt", "
model {

    # Priors and constraints
    for (t in 1:(n.occasions-1)){
       logit(phi[t]) <- logit(mean.phi) + eps[t]         # Priors for survival
       eps[t] ~ dnorm(0, tau)       
       p[t] ~ dunif(0, 1)           # Priors for recapture
       }
    
    mean.phi ~ dunif(0, 1)
    #mu <- logit(mean.phi)
        
    sigma ~ dunif(0, 10)
    tau <- 1/(sigma * sigma)

    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
       marr[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
       }

    # Define the cell probabilities of the m-array
    for (t in 1:(n.occasions-1)){
       q[t] <- 1-p[t]                # Probability of non-recapture

       # Main diagonal
       pr[t,t] <- phi[t]*p[t]
      
       # Above main diagonal
       for (j in (t+1):(n.occasions-1)){
          pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
          } #j

       # Below main diagonal
       for (j in 1:(t-1)){
          pr[t,j] <- 0
          } #j
       } #t

    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
       pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
       } #t
    
    # Assess model fit using Freeman-Tukey statistic [can be done offline]
    
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
       for (j in 1:n.occasions){
          expmarr[t,j] <- rel[t]*pr[t,j]
          E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
          } #j
       } #t

    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
       marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
       for (j in 1:n.occasions){
          E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
          } #j
       } #t
    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
}
")

# Create the m-array from the capture-histories
marr <- marray(CH)

# Bundle data
jags.data <- list(marr = marr, n.occasions = ncol(marr), rel = rowSums(marr))

# Initial values
inits <- function(){list(p = runif(dim(marr)[2]-1, 0, 1))}  

# Parameters monitored
parameters <- c("phi", "p", "mean.phi", "sigma", "fit", "fit.new")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2500
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs <- jags(jags.data, inits, parameters, "cjs-mnl.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(cjs, digits = 3) 

# Evaluation of fit
plot(cjs$sims.list$fit, cjs$sims.list$fit.new, xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1, ylim = c(5, 25), xlim = c(5, 25), bty ="n") 
abline(0, 1, col = "black", lwd = 2)
mean(cjs$sims.list$fit.new > cjs$sims.list$fit)


### 7.10.3. Age-dependent models

# Define parameter values
n.occasions <- 12                    # Number of capture occasions
marked.j <- rep(200, n.occasions-1)  # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)   # Annual number of newly marked adults
phi.juv <- 0.3                       # Juvenile annual survival
phi.ad <- 0.65                       # Adult annual survival
p <- rep(0.5, n.occasions-1)         # Recapture
phi.j <- c(phi.juv, rep(phi.ad,n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:(length(marked.j)-1)){
   PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
   }
P.J <- matrix(rep(p, n.occasions*sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a) 

CH <- rbind(CH.J, CH.A)
age <- c(rep(1,nrow(CH.J)), rep(2,nrow(CH.A)))



#########################################
#
# Function to create age-dependent m-arrays from single-state capture-recapture data
#
# Input variables
#    ch: matrix with capture histories. 
#        Note: the capture history file is a single file that includes the individuals of all age classes
#    age: vector with the age class at first capture for each individual
#    mAge: maximal number of age classes for which m-arrays are constructed. Input is optional and only required if the age matrix has fewer age classes as we want to separate (e.g. CH contains only individuals marked as juveniles, and we want 2 age classes)
#
# Output
#    marr: 3-d array with the m-array. The third dimension is the age class. The last column of each m-array is the number of released individuals that were never recaptured. Thus, the total number of released individuals per occasion is the row sum of the m-array.
#
# Written: 14.3.2016, M.Schaub
#
# Last up-date: 15.08.2018
#
################################################



marray.age <- function(ch, age, mAge = 1){

   # 1. Helper functions
   # 1.1. Function to create a m-array based on capture-histories (ch)
   marray <- function(ch){
      ns <- length(table(ch)) - 1
      no <- ncol(ch)
      out <- matrix(0, ncol = ns*(no-1)+1, nrow = ns*(no-1))
      # Remove capture histories of individuals that are marked at last occasion
      get.first <- function(x) min(which(x!=0))
      first <- apply(ch, 1, get.first)
      last <- which(first==no)
      if (length(last) > 0) ch <- ch[-last,]
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
         } # i
      return(out)
      }  

   # 1.2. Function to remove histories without any capture from a capture-recapture matrix
   clean.ch <- function(ch){
      incl <- which(rowSums(ch)>=1)
      ch <- ch[incl,]
      return(ch)
      }

  # 1.3. Function to remove the first capture in a capture-recapture matrix
   rm.first <- function(ch) {
      get.first <- function(x) min(which(x==1))
      first <- apply(ch, 1, get.first)
      for (i in 1:nrow(ch)){
         ch[i,first[i]] <- 0
         }
      return(ch)
      }

   # 1.4. Function to calculate the occasion of first capture
   get.first <- function(x) min(which(x!=0))


   # 2. Calculations   
   if (is.matrix(ch)==FALSE) ch <- matrix(ch, nrow = 1)   
   maxAge <- max(c(max(age), mAge))
   nind <- nrow(ch)
   n.occasions <- ncol(ch)
   
   first <- apply(ch, 1, get.first)
   age.matrix <- matrix(0, ncol = n.occasions, nrow = nind)
   for (i in 1:nind){
      age.matrix[i,first[i]:n.occasions] <- 1:(n.occasions-first[i]+1)+(age[i]-1)
      }
   age.matrix[age.matrix > maxAge] <- maxAge

   # Recode capture history
   ch.rec <- ch
   for (i in 1:nind){
      h <- which(ch.rec[i,]==1)
      for (j in 1:length(h)){
         ch.rec[i,h[j]] <- j
         } # j
      } # i
   ch.rec[ch.rec > maxAge] <- maxAge
      
   ch.split <- array(0, dim = c(nrow(ch), ncol(ch), maxAge))
   for (a in 1:maxAge){
      for (i in 1:nind){
         j <- which(ch.rec[i,]==a | ch.rec[i,]==(a+1))
         if (length(j)==0) next
         ch.split[i,j[1:2],age.matrix[i,j[1]]] <- 1
         if (length(j)>1){
            ch.split[i,j[2:length(j)],age.matrix[i,j[2]]] <- 1
            }
         } # i
      } # a

   marr <- array(0, dim = c(n.occasions-1, n.occasions, maxAge))
   for (a in 1:(maxAge-1)){
      for (i in 1:nind){
         u <- which(ch.split[i,,a]==1)
         if (length(u)==0) next
         if (u[1]==n.occasions) next
         if (length(u)==1) marr[u,n.occasions,a] <- marr[u,n.occasions,a] + 1
         if (length(u)==2) marr[u[1],u[2]-1,a] <- marr[u[1],u[2]-1,a] + 1
         } # i
      } # a
   a <- maxAge

   if (is.matrix(ch.split[,,a])==FALSE){ 
      ch.split1 <- matrix(ch.split[,,a], nrow = 1)
      marr[,,a] <- marray(ch.split1)
      } # if
   else marr[,,a] <- marray(clean.ch(ch.split[,,a]))      
   return(marr)
   }

marray <- marray.age(CH, age)



# Specify model in BUGS language
cat(file = "cjs-mnl-age.txt", "
model {

    # Constraints
    for (t in 1:(n.occasions-1)){
       phi.juv[t] <- mean.phijuv
       phi.ad[t] <- mean.phiad
       p[t] <- mean.p
       }

    # Priors
    mean.phijuv ~ dunif(0, 1)          # Prior for mean juv. survival
    mean.phiad ~ dunif(0, 1)           # Prior for mean ad. survival
    mean.p ~ dunif(0, 1)               # Prior for mean recapture
    
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
       marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
       marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
       }
    
    # Define the cell probabilities of the m-arrays
    
    for (t in 1:(n.occasions-1)){
       q[t] <- 1-p[t]            # Probability of non-recapture

       # Main diagonal
       pr.j[t,t] <- phi.juv[t]*p[t]
       pr.a[t,t] <- phi.ad[t]*p[t]
      
       # Above main diagonal
       for (j in (t+1):(n.occasions-1)){
          pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q[t:(j-1)])*p[j]
          pr.a[t,j] <- prod(phi.ad[t:j])*prod(q[t:(j-1)])*p[j]
          } #j
      
       # Below main diagonal
       for (j in 1:(t-1)){
          pr.j[t,j] <- 0
          pr.a[t,j] <- 0
          } #j
       } #t

    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
       pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
       pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
       } #t
}
")

# Bundle data
jags.data <- list(marr.j = marray[,,1], marr.a = marray[,,2], n.occasions = dim(marray)[2], rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2]))

# Initial values
inits <- function(){list(mean.phijuv = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phijuv", "mean.phiad", "mean.p")

# MCMC settings
ni <- 3000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT <1 min)
cjs.2 <- jags(jags.data, inits, parameters, "cjs-mnl-age.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)


par(mfrow = c(1, 3), las = 1)
hist(cjs.2$sims.list$mean.phijuv, nclass = 30, col = "gray", main = "", xlab = "Juvenile survival", ylab = "Frequency")
abline(v = phi.juv, col = "red", lwd = 2)
hist(cjs.2$sims.list$mean.phiad, nclass = 30, col = "gray", main = "", xlab = "Adult survival", ylab = "")
abline(v = phi.ad, col = "red", lwd = 2)
hist(cjs.2$sims.list$mean.p, nclass = 30, col = "gray", main = "", xlab = "Recapture", ylab = "")
abline(v = p[1], col = "red", lwd = 2)


###############
### EXERCISE 2bis: Analyse real CR data of female Leisler's bats (7.11)
###############

# m-array obtained from 18 marked female Leisler's bats 
# from Thuringia (Germany) with 19 capture occasions
m.leisleri <- matrix(c(4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                       0,5,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                       0,0,9,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                       0,0,0,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
                       0,0,0,0,10,2,1,0,0,0,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,15,0,0,0,0,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,0,11,2,0,1,0,0,0,0,0,0,0,0,19,
                       0,0,0,0,0,0,0,12,1,1,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,0,0,0,13,2,0,0,0,0,0,0,0,0,4,
                       0,0,0,0,0,0,0,0,0,14,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,0,0,0,0,0,13,1,0,0,0,1,0,0,8,
                       0,0,0,0,0,0,0,0,0,0,0,15,3,1,0,0,0,0,12,
                       0,0,0,0,0,0,0,0,0,0,0,0,12,4,0,1,0,0,7,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,19,2,0,0,0,3,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,28,1,0,0,4,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,7,2,21,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,2,21,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,18), 
                     ncol = 19, nrow = 18, byrow = TRUE)

# Estimate mean annual survival and temporal variance, assuming constant capture probability

# Assess model with Bayesian p-value calculated with a Freeman-Tukey measure of discrepancy

