`%notin%` <- Negate(`%in%`)

book.rMatern <- function(n, coords, sigma=1, range, kappa = sqrt(8*nu)/range, variance = sigma^2, nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}

linefind <- function(s1,s2){
  y1 <- s1[2]
  x1 <- s1[1]
  y2 <- s2[2]
  x2 <- s2[1]
  m <- (y2-y1)/(x2-x1)
  b <- y1-m*x1
  return(c(b,m))
}

plot.b.map <- function(coords,b,pal=hcl.colors(12,"Zissou 1"),pch = 15, ...) {
  cr <- colorRamp(rev(pal))
  legend_image <- as.raster(matrix(pal, ncol=1))
  layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
  plot(coords[,1], coords[,2], col=rgb(cr((b - floor(min(b))) / (ceiling(max(b)) - floor(min(b)))), max=255),
       xlab="Longitude", ylab="Latitude", pch = pch, asp = 1, ...)
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
  text(x=1.5, y = seq(0,1,l=5), labels = round(seq(floor(min(b)),ceiling(max(b)),l=5),3))
  rasterImage(legend_image, 0, 0, 1,1)
  par(mfrow=c(1,1))
}

plot.b.map2 <- function(coords,b,pal=hcl.colors(12,"Zissou 1"),pch = 15, ...) {
  cr <- colorRamp(rev(pal))
  legend_image <- as.raster(matrix(pal, ncol=1))
  layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
  plot(coords[,1], coords[,2], col=rgb(cr((b - min(b)) / (max(b) - min(b))), max=255),
       xlab="Longitude", ylab="Latitude", pch = pch, asp = 1, ...)
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
  text(x=1.5, y = seq(0,1,l=5), labels = round(seq(min(b),max(b),l=5),3))
  rasterImage(legend_image, 0, 0, 1,1)
  par(mfrow=c(1,1))
}


run_INLA <- function(data, bounds, tau.m, PCprior = FALSE, a = 1.1){
  #Retrieve data info
  n.m <- sum(data$source=="SDM")
  n.g <- sum(data$source=="Gen")
  n.p <- sum(data$source=="Pol")
  coords.s <- as.matrix(data[,1:2])
  coord.m <- coords.s[1:n.m,]
  coord.g <- coords.s[n.m+1:n.g,]
  coord.p <- coords.s[n.m+n.g+1:n.p,]
  p.m <- data$Y[1:n.m]
  y.g <- data$Y[n.m+1:n.g]
  y.p <- data$Y[n.m+n.g+1:n.p]
  
  
  max.edge <- 4
  bnd <- inla.nonconvex.hull(bounds, convex = -0.12)
  mesh <- inla.mesh.2d(coords.s, max.edge=c(1,5)*max.edge, cutoff = max.edge/3, boundary = bnd)
  
  #####Define sm1 zone in mesh
  sm1 <- coord.m[p.m >= tau.m,]
  radi <- 1.5
  keep <- rep(FALSE,nrow(mesh$loc))
  if(sum(p.m >= tau.m)==1){
    keep <- (pointDistance(sm1,mesh$loc[,1:2],lonlat = FALSE) <= radi) + keep
  } else{
    for(i in 1:sum(p.m >= tau.m)){
      keep <- (pointDistance(sm1[i,],mesh$loc[,1:2],lonlat = FALSE) <= radi) + keep
    }
  }
  
  zone <- rep(0,length(keep))
  zone[keep>0] <- 1
  
  #####Non-stationary SPDE- eventually include as function arguments
  nu <- 1
  alpha <- 2
  logkappa0 <- log(8 * nu) / 2
  logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 - logkappa0
  
  zspde <- inla.spde2.matern(mesh, 
                             B.tau = cbind(logtau0, -1, -zone, nu),
                             B.kappa = cbind(logkappa0, 0, 0, -1),
                             theta.prior.mean = c(0,1,0), 
                             theta.prior.prec = rep(1, 3))
  spde <- inla.spde2.pcmatern(mesh, prior.range = c(1, 0.01), 
                              prior.sigma = c(3, 0.01)) #alpha= 3/2 for exp, alpha = 2 for matern
  #Priors
  hyper <- list(theta = list(prior = 'normal', param = c(1, 0.001)))
  fam.control <- list(list(),list(),list())
  if(PCprior){
    prior.g <- list(prior = "pc.prec", param = c(a,0.01),
                    fixed = FALSE)
    
    pg <- list(hyper = list(prec = prior.g))
    
    prior.mp <- list(initial = 0, prior = "loggamma", param = c(1, 0.00005),
                     fixed = FALSE)
    
    pmp <- list(hyper = list(prec = prior.mp))
    fam.control <- list(list(),pg,pmp)
  } 
  
  #INLA
  formula <- y ~ 0 + intercept + biasm  + biasp + f(sg, model = spde) + 
    f(sm, copy = "sg", fixed = FALSE, hyper = hyper) +
    f(sp, copy = "sg", fixed = FALSE, hyper = hyper) +
    f(z, model = zspde)
  
  Am <- inla.spde.make.A(mesh, coord.m) 
  Ag <- inla.spde.make.A(mesh, coord.g) 
  Ap <- inla.spde.make.A(mesh, coord.p) 
  
  stackm <- inla.stack(
    data = list(y = cbind(as.vector((p.m>tau.m)*1), NA, NA)),
    A = list(Am), 
    effects = list(list(intercept = 1, biasm = 1, sm = 1:spde$n.spde, z = 1:zspde$n.spde))) 
  
  stackg <- inla.stack(
    data = list(y = cbind(NA, as.vector(y.g), NA)),
    A = list(Ag), 
    effects = list(list(intercept = 1, sg = 1:spde$n.spde)))
  
  stackp <- inla.stack(
    data = list(y = cbind(NA, NA,as.vector(y.p))),
    A = list(Ap), 
    effects = list(list(intercept = 1, biasp = 1, sp = 1:spde$n.spde))) 
  
  stack <- inla.stack(stackm, stackg, stackp) 
  familia <- c("binomial","gaussian","gaussian")
  
  res <- inla(formula,
              family = familia,
              data = inla.stack.data(stack),
              control.predictor = list(link = 1, A = inla.stack.A(stack)),
              control.family = fam.control,
              control.compute = list(config = TRUE),
              control.inla = list(control.correct = list(enable = TRUE, factor = 10)))
  
  return(list(INLA_results = res, mesh = mesh))
}

run_INLA_noSDM <- function(data, bounds, PCprior = FALSE, a = 1.1){
  #Retrieve data info
  n.g <- sum(data$source=="Gen")
  n.p <- sum(data$source=="Pol")
  coords.s <- as.matrix(data[,1:2])
  coord.g <- coords.s[1:n.g,]
  coord.p <- coords.s[n.g+1:n.p,]
  y.g <- data$Y[1:n.g]
  y.p <- data$Y[n.g+1:n.p]
  
  max.edge <- 4
  bnd <- inla.nonconvex.hull(bounds, convex = -0.12)
  mesh <- inla.mesh.2d(coords.s, max.edge=c(1,5)*max.edge, cutoff = max.edge/3, boundary = bnd)
  
  #####Non-stationary SPDE- eventually include as function arguments
  nu <- 1
  alpha <- 2

  spde <- inla.spde2.pcmatern(mesh, prior.range = c(1, 0.01), 
                              prior.sigma = c(3, 0.01)) #alpha= 3/2 for exp, alpha = 2 for matern
  #Priors
  hyper <- list(theta = list(prior = 'normal', param = c(1, 0.001)))
  fam.control <- list(list(),list())
  if(PCprior){
    prior.g <- list(prior = "pc.prec", param = c(a,0.01),
                    fixed = FALSE)
    
    pg <- list(hyper = list(prec = prior.g))
    
    prior.mp <- list(initial = 0, prior = "loggamma", param = c(1, 0.00005),
                     fixed = FALSE)
    
    pmp <- list(hyper = list(prec = prior.mp))
    fam.control <- list(pg,pmp)
  } 
  
  #INLA
  formula <- y ~ 0 + intercept + biasp + f(sg, model = spde) + 
    f(sp, copy = "sg", fixed = FALSE, hyper = hyper) 
  
  Ag <- inla.spde.make.A(mesh, coord.g) 
  Ap <- inla.spde.make.A(mesh, coord.p) 
  
  stackg <- inla.stack(
    data = list(y = cbind(as.vector(y.g), NA)),
    A = list(Ag), 
    effects = list(list(intercept = 1, sg = 1:spde$n.spde)))
  
  stackp <- inla.stack(
    data = list(y = cbind(NA,as.vector(y.p))),
    A = list(Ap), 
    effects = list(list(intercept = 1, biasp = 1, sp = 1:spde$n.spde))) 
  
  stack <- inla.stack(stackg, stackp) 
  familia <- c("gaussian","gaussian")
  
  res <- inla(formula,
              family = familia,
              data = inla.stack.data(stack),
              control.predictor = list(A = inla.stack.A(stack)),
              control.family = fam.control,
              control.compute=list(config = TRUE))
  
  return(list(INLA_results = res, mesh = mesh))
}

p.recover <- function(res,proj,samps=1000){
  x.smp <- inla.posterior.sample(samps, res, selection = list(intercept = 1, sg = 0))
  x.smp <- inla.posterior.sample.eval(function() sg+intercept,x.smp)
  p.smp <- pnorm(x.smp)
  p <- rowMeans(p.smp)
  qs <- apply(p.smp,1,quantile,probs=c(0.025,0.975))
  post.samp <- inla.mesh.project(proj,p.smp)
  p.mean <- inla.mesh.project(proj, p)
  p.q2.5 <- inla.mesh.project(proj, qs[1,])
  p.q97.5 <- inla.mesh.project(proj, qs[2,])
  return(list(means=p.mean,q2.5=p.q2.5,q97.5=p.q97.5,sample=post.samp))
} 

run_INLA_linSDM <- function(data, bounds, tau.m, PCprior = FALSE, a = 1.1){
  #Retrieve data info
  n.m <- sum(data$source=="SDM")
  n.g <- sum(data$source=="Gen")
  n.p <- sum(data$source=="Pol")
  coords.s <- as.matrix(data[,1:2])
  coord.m <- coords.s[1:n.m,]
  coord.g <- coords.s[n.m+1:n.g,]
  coord.p <- coords.s[n.m+n.g+1:n.p,]
  y.m <- data$Y[1:n.m]
  y.g <- data$Y[n.m+1:n.g]
  y.p <- data$Y[n.m+n.g+1:n.p]
  p.m <- pnorm(y.m)
  
  
  max.edge <- 4
  bnd <- inla.nonconvex.hull(bounds, convex = -0.12)
  mesh <- inla.mesh.2d(coords.s, max.edge=c(1,5)*max.edge, cutoff = max.edge/3, boundary = bnd)
  
  #####Define sm1 zone in mesh
  sm1 <- coord.m[p.m >= tau.m,]
  radi <- 1.5
  keep <- rep(FALSE,nrow(mesh$loc))
  if(sum(p.m >= tau.m)==1){
    keep <- (pointDistance(sm1,mesh$loc[,1:2],lonlat = FALSE) <= radi) + keep
  } else{
    for(i in 1:sum(p.m >= tau.m)){
      keep <- (pointDistance(sm1[i,],mesh$loc[,1:2],lonlat = FALSE) <= radi) + keep
    }
  }
  
  zone <- rep(0,length(keep))
  zone[keep>0] <- 1
  
  #####Non-stationary SPDE- eventually include as function arguments
  nu <- 1
  alpha <- 2
  logkappa0 <- log(8 * nu) / 2
  logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 - logkappa0
  
  zspde <- inla.spde2.matern(mesh, 
                             B.tau = cbind(logtau0, -1, -zone, nu),
                             B.kappa = cbind(logkappa0, 0, 0, -1),
                             theta.prior.mean = c(0,1,0), 
                             theta.prior.prec = rep(1, 3))
  spde <- inla.spde2.pcmatern(mesh, prior.range = c(1, 0.01), 
                              prior.sigma = c(3, 0.01)) #alpha= 3/2 for exp, alpha = 2 for matern
  #Priors
  hyper <- list(theta = list(prior = 'normal', param = c(1, 0.001)))
  fam.control <- list(list(),list(),list())
  if(PCprior){
    prior.g <- list(prior = "pc.prec", param = c(a,0.01),
                    fixed = FALSE)
    
    pg <- list(hyper = list(prec = prior.g))
    
    prior.mp <- list(initial = 0, prior = "loggamma", param = c(1, 0.00005),
                     fixed = FALSE)
    
    pmp <- list(hyper = list(prec = prior.mp))
    fam.control <- list(list(),pg,pmp)
  } 
  
  #INLA
  formula <- y ~ 0 + intercept + biasm  + biasp + f(sg, model = spde) + 
    f(sm, copy = "sg", fixed = FALSE, hyper = hyper) +
    f(sp, copy = "sg", fixed = FALSE, hyper = hyper) +
    f(z, model = zspde)
  
  Am <- inla.spde.make.A(mesh, coord.m) 
  Ag <- inla.spde.make.A(mesh, coord.g) 
  Ap <- inla.spde.make.A(mesh, coord.p) 
  
  stackm <- inla.stack(
    data = list(y = cbind(as.vector(y.m), NA, NA)),
    A = list(Am), 
    effects = list(list(intercept = 1, biasm = 1, sm = 1:spde$n.spde, z = 1:zspde$n.spde))) 
  
  stackg <- inla.stack(
    data = list(y = cbind(NA, as.vector(y.g), NA)),
    A = list(Ag), 
    effects = list(list(intercept = 1, sg = 1:spde$n.spde)))
  
  stackp <- inla.stack(
    data = list(y = cbind(NA, NA,as.vector(y.p))),
    A = list(Ap), 
    effects = list(list(intercept = 1, biasp = 1, sp = 1:spde$n.spde))) 
  
  stack <- inla.stack(stackm, stackg, stackp) 
  familia <- c("gaussian","gaussian","gaussian")
  
  res <- inla(formula,
              family = familia,
              data = inla.stack.data(stack),
              control.predictor = list(A = inla.stack.A(stack)),
              control.family = fam.control,
              control.compute = list(config = TRUE))
  
  return(list(INLA_results = res, mesh = mesh))
}