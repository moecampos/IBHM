source("LIB.R")
source("FUN.R")

#db is a n times 4 data.frame where:
##Column 1 and 2 are latitude and longitude, respectively
##Column 3 are the Y values
##Column 4 indicates to what source the Y values belongs to (i.e.SDM, Gen, Pol) - sources should be in that order for now

#bounds.r are the set of coordinates that delimitate the boundary used when creating the mesh

#SDM_treat can be one of three: "None", "Logistic", "Linear"

#when SDM is present and used as logistic, tau.m is the threshold value that separates 0s and 1s

#file_names is a vector of length 4:
##Name 1 is used for saving simulations of the posterior latent field
##Name 2 is used for saving the whole INLA output (very memory intensive)
##Name 3 is used for saving only the posterior summary of the model's parameters
##Name 4 is used for saving the mesh, as it cannot be recovered from the INLA output


IBHM <- function(db, bounds.r, PCa, tau.m = NULL, pcp = TRUE, SDM_treat = "Logistic", 
                 file_names = NULL, save = TRUE){
  ini <- Sys.time()
  
  if(SDM_treat == "None") results <- run_INLA_noSDM(db, bounds.r, PCprior = pcp, a = PCa) 
  if(SDM_treat == "Logistic") results <- run_INLA(db, bounds.r, tau.m, PCprior = pcp, a = PCa)
  if(SDM_treat == "Linear") results <- run_INLA_linSDM(db, bounds.r, tau.m, PCprior = pcp, a = PCa)
  
  res <- results$INLA_results
  mesh <- results$mesh
    
  gproj <- inla.mesh.projector(mesh,coo.all)
  post.p <- p.recover(res,gproj)
    
  parsums <- list(fixed = res$summary.fixed, hyper = res$summary.hyperpar)

  fini <- Sys.time() - ini
  
  if(save){
    saveRDS(post.p$sample, file_names[1])
    saveRDS(res, file_names[2])
    saveRDS(parsums, file_names[3])
    saveRDS(mesh, file_names[4])
  }
}