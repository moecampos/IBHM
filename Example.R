source("Data_Simulation.R")
source("IHBM.R")
combou <- readRDS("common_boundary.rds")
bounds.r <- combou$bounds

#Combine coordinates from all sources
all.coords <- rbind(coo.m,coo.g,coo.p)

#Create data frame for use in IBHM()
db <- data.frame(longitude = all.coords[,1], latitude = all.coords[,2],
                 Y = c(p.m, y.g, y.p), #keep same order from all.coords
                 source = c(rep("SDM",length(p.m)),
                            rep("Gen",length(y.g)),
                            rep("Pol",length(y.p))))

#example of running IBHM() without saving anything
IBHM(db, bounds.r, coo.all, tau.m = 0.502, save = FALSE)
