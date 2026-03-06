## Run 1 rep of forward sim
# SET parameters and output file name below
# MUST save before running


#### Load packages and scripts ####
setwd("")
library(kinship2)
library(optiSel)
source("SimFunctions.R")


#### Set up simulation ####
#NOTE: kept this as separate first step because genIDs for real data are different than
# the sim data, so easier to generate 1st set of spawned pairs here and then send into
# sim loop

## Import pedigree
ped.df <- read.csv("data/pedigree2026.csv", stringsAsFactors=F)

## Create 1st pairs object to simulate progeny
pairs <- data.frame(Sire = ped.df$ID[ped.df$yr==max(ped.df$yr) & (substr(ped.df$ID,6,6)==2)])
pairs$Dam <- pairs$Sire - 1
pairs$kin <- NA
pairs$PC <- as.integer(substr(pairs$Sire,3,5))
pairs$yr <- max(ped.df$yr)


#### Run simulation ####
ped.new <- run_forward_sim(ped.df,pairs,
                           n.generations=48,
                           target_pairs=312,
                           max.progeny=12,
                           kin_seq=c(0.03125,0.0625),
                           max.fam=4,
                           weeks=8,
                           loss=0.25,
                           early_frac=0.5,
                           tries=100)
#length(unique(ped.new$yr))-length(unique(ped.df$yr)) #check # gens if stops early

## Save pedigree & parameters
file_name <- "sim_48gen_312pc_varykin_avail_8weeks"
write.csv(ped.new, paste0("output/",file_name,".csv"), quote=F, row.names=FALSE)

#save parameters
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
write.csv(params_df, paste0("output/param_",file_name,".csv"), quote=F, row.names=FALSE)

print("Done")
