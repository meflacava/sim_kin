## Run multiple round of forward sim with given parameters
# SET parameters and output file name below
# MUST save before running


#### Load packages and scripts ####
library(kinship2)
Sys.setenv(RGL_USE_NULL = TRUE)
library(optiSel)
source("SimFunctions.R")


#### Set up simulation ####
#NOTE: kept this as separate first step because genIDs for real data are different than
# the sim data, so easier to generate 1st set of spawned pairs here and then send into
# sim loop

## Import pedigree
ped.df <- read.csv("data/pedigree2026.csv", stringsAsFactors=F)
#exclude 1 pair with high kinship that was accidentally crossed
# (to avoid choosing their offspring in the first sim gen)
ped.df <- ped.df[ped.df$ID!=192401 & ped.df$ID!=192402,]

## Create 1st pairs object to simulate progeny
pairs <- data.frame(Sire = ped.df$ID[ped.df$yr==max(ped.df$yr) & (substr(ped.df$ID,6,6)==2)])
pairs$Dam <- pairs$Sire - 1
pairs$kin <- NA
pairs$PC <- as.integer(substr(pairs$Sire,3,5))
pairs$yr <- max(ped.df$yr)


#### Run simulation multiple times ####
n.reps <- 20 #number of repetitons for this set of parameters
for (i in 1:n.reps){
  ped.new <- run_forward_sim(ped.df,pairs,
                             n.generations=50,
                             target_pairs=312, #312=real
                             max.progeny=12,
                             kin_seq=seq(0.005,0.0625,0.001), #c(0.03125,0.0625); seq(0.005,0.0625,0.001)
                             max.fam=4, #4
                             weeks=4, #4 or 8
                             loss=0.25, #0.25
                             early_frac=0.5, #0.5
                             tries=100)
  
  ## Save pedigree & parameters
  file_name <- paste0("sim_50gen_312pc_kinSteps_avail_4weeks_",formatC(i,width=3,flag="0"))
  write.csv(ped.new, paste0("output/iterations/",file_name,".csv"), quote=F, row.names=FALSE)
  
  if (i==1){ #save parameters just once
    at <- attributes(ped.new)[!names(attributes(ped.new)) %in% c("names","row.names","class")]
    params_df <- do.call(
      rbind,
      lapply(names(at), function(nm) {
        data.frame(
          param = nm,
          index = seq_along(at[[nm]]),
          value = at[[nm]],
          stringsAsFactors=FALSE
        )
      })
    )
    write.csv(params_df, paste0("output/iterations/param_",file_name,".csv"), quote=F, row.names=FALSE)
  }
}

print("Done")
