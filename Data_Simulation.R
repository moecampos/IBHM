source("LIB.R")
source("FUN.R")

combou <- readRDS("common_boundary.rds")
bounds.r <- combou$bounds
coords.r <- combou$coords

sample_sizes <- c(1000,15,60,200) #SDM, pollen, genetic and prediction-only sites, respectively

seed <- 1111
set.seed(seed)
n.m <- sample_sizes[1]
n.p <- sample_sizes[2]
n.g <- sample_sizes[3]
v <- sample_sizes[4]
n <- n.m+n.p+n.g

smp <- sample(1:nrow(coords.r),n)
coo.s <- coords.r[smp,]
sdm.s <- sample(1:n,n.m)
if(length(sdm.s)==0){ #in case we want no SDM
  g.s <- (1:n)[1:n.g]
  p.s <- (1:n)[(n.g+1):(n.g+n.p)]
} else{
  g.s <- (1:n)[-sdm.s][1:n.g]
  p.s <- (1:n)[-sdm.s][(n.g+1):(n.g+n.p)]
  coo.m <- coo.s[sdm.s,]
  
}
coo.g <- coo.s[g.s,]
coo.p <- coo.s[p.s,]

p.smp <- sample((1:nrow(coords.r))[1:nrow(coords.r)%notin%smp],v) #prediction sample
coo.pr <- coords.r[p.smp,]

coo.all <- rbind(coo.s,coo.pr)

rhoX <- 3
sigX <- 0.8 
rhoZ <- 8
sigZ <- 10 
mZ_1 <- 1.2
a.m <- -10
a.p <- -0.8
s2g <- 0.2 
s2p <-0.2
b.m <- 1
b.p <- 1.4
mu <- -0.5 
tau.m <- 0.502

x <- book.rMatern(1, coo.all, range = rhoX, sigma = sigX) #function in FUN.R
enough <- FALSE 
while(!enough){ #this loop makes sure there is at lesat one true point in Sm1
  z <- book.rMatern(1, coo.m, range = rhoZ, sigma = sigZ)
  l.m <- a.m + b.m*(mu + x[sdm.s]) + z
  sm1 <- l.m > log(tau.m/(1-tau.m))
  if(sum(sm1)>=1){
    enough <- TRUE
  }
}
l.m[sm1] <- l.m[sm1]+(mZ_1-1)*z[sm1]
p.m <- exp(l.m)/(1+exp(l.m)) #SDM probabilities are dichotomized based on tau.m
y.g <- mu + x[g.s] + rnorm(n.g, sd = sqrt(s2g)) #genetic data
y.p <- a.p + b.p*(mu + x[p.s]) + rnorm(n.p, sd = sqrt(s2p)) #pollen data

rm(list= ls()[!(ls() %in% c('p.m','y.g', 'y.p'))])


