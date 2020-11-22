# To artificially generate uncertainty on both the states non-breeder and breeder, 
# we used the R script below to alter the raw capture-recapture data from file titis.inp; the resulting file is titis2.inp.

# read in data
titi <- read.table('dat/titis.txt')
# 1 seen as breeder
# 2 seen as non-breeder
# 0 not seen

# nb of capture occasions
ny <- ncol(titi)
# nb of individuals
nind <- nrow(titi)

titi2 <- titi
for (i in 1:nind)
{
  for (j in 1:ny){
    # 1 seen and ascertained Breeder (with probability .2)
    # 2 seen and ascertained Non-Breeder (with probability .7)
    # 3 seen but not ascertained (Non-Breeders with probability .8 + Breeders with probability .3)
    # 0 not seen
    
    # Breeders are ascertained with probability .2
    if (titi[i,j] == 1)
    {
      temp <- rbinom(1,size=1,prob=.2) 
      if (temp == 1) titi2[i,j] <- 1 # if ascertained B, event = 1
      if (temp == 0) titi2[i,j] <- 3 # if not ascertained, event = 3
    }
    
    # Non-Breeders are ascertained with probability .7 (event = 1), 
    # or not ascertained with probability .3 (event = 2)
    if (titi[i,j] == 2) 
    {
      temp <- rbinom(1,size=1,prob=.7)
      if (temp == 1) titi2[i,j] <- 2 # if ascertained NB, event = 2
      if (temp == 0) titi2[i,j] <- 3 # if not ascertained, event = 3
    }
    
  }
}

# write data in a file
write(t(titi2), file="dat/titis2.txt", ncolumns = ny)
