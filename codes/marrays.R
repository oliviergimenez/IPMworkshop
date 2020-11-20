
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

