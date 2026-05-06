## Optimize fish availability parameters



#### Simulate parameter combos ####
library(kinship2)
Sys.setenv(RGL_USE_NULL = TRUE)
library(optiSel)
source("SimFunctions.R")


##### Set up for sims ####
## Import pedigree
ped.df <- read.csv("data/pedigree2026.csv", stringsAsFactors=F)

## Create 1st pairs object to simulate progeny
pairs <- data.frame(Sire = ped.df$ID[ped.df$yr==max(ped.df$yr) & (substr(ped.df$ID,6,6)==2)])
pairs$Dam <- pairs$Sire - 1
pairs$kin <- NA
pairs$PC <- as.integer(substr(pairs$Sire,3,5))
pairs$yr <- max(ped.df$yr)

ped.sim <- ped.df

# 1. Simulate progeny from previous generation’s pairs
ped.sim <- simulate_progeny(ped.sim, pairs=pairs, max.progeny=12)

# Identify new spawners (the newly created progeny)
spawners <- ped.sim[ped.sim$yr==max(ped.sim$yr),]
# pick one representative progeny per family (PC) for kinship calc
rep_ids <- tapply(spawners$ID, spawners$PC, `[`, 1)


# 2. Kinship for new spawners (needs all their ancestors)
# ancestors of representative progeny only
ancestors <- find_ancestors(ped.sim, rep_ids)
ped.sub <- ped.sim[ped.sim$ID %in% c(rep_ids, ancestors),]

# kinship among representative progeny
K_rep <- kinship(id=ped.sub$ID,
                 dadid=ped.sub$Sire,
                 momid=ped.sub$Dam)
K_rep <- K_rep[rep_ids, rep_ids]

# map each spawner to its representative
rep_for_spawner <- rep_ids[as.character(spawners$PC)]

# expand kinship back to individual spawners
Ks <- K_rep[rep_for_spawner, rep_for_spawner]
rownames(Ks) <- spawners$ID
colnames(Ks) <- spawners$ID


##### Run sims ####

for (i in c(4,8,16)) {        # weeks
  for (j in c(0.25,0.5)) {    # early fraction
    for (k in c(0.25,0.5)) {  # loss
      cat(i, "weeks,", j, "early,", k, "loss\n") #track run status
      p <- choose_pairs_avail(
        spawners,
        Ks,
        target_pairs = 312,
        kin_thresh = 0.01,
        max.fam = 4,
        weeks = i,
        early_frac = j,
        loss = k,
        tries = 20)
      p$DamPC <- as.integer(substr(p$Dam,4,6))
      p$SirePC <- as.integer(substr(p$Sire,4,6))
      fname <- paste0("output/availParams/w",formatC(i,width=2,flag="0"),
                      "_early",j,"_lose",k,".csv")
      write.csv(p, fname, row.names=FALSE)
    }
  }
}



#### Compare sims to real ####

## Import real cross data and transform to 1 row per fish
years <- 2022:2025
real <- vector("list", length(years))
names(real) <- years

for (i in years) {
  rc <- read.csv(paste0("data/source/", i, "_refuge_crosses.csv"))
  rc <- rc[rc$PC.FSG > 0, ]
  rc$Male.PC.FSG <- as.integer(rc$Male.PC.FSG)
  rc$Female.PC.FSG <- as.integer(rc$Female.PC.FSG)
  
  #convert 1 row per PC to 1 row per fish
  real[[as.character(i)]] <- data.frame(
    ID = c(rc$Male, rc$Female),
    sourcePC = c(rc$Male.PC.FSG, rc$Female.PC.FSG),
    newPC = rep(rc$PC.FSG, times = 2),
    year = i
  )
}


## Import sim data
files <- list.files("output/availParams/",pattern="\\.csv$",full.names=TRUE)
simulated <- vector("list", length(files))
names(simulated) <- basename(files)
for (i in seq_along(files)) {
  #simulated[[i]] <- read.csv(files[i])
  rc <- read.csv(files[i])
  simulated[[i]] <- data.frame(
    ID = c(rc$Sire, rc$Dam),
    sourcePC = c(rc$SirePC, rc$DamPC),
    newPC = rep(rc$PC, times = 2),
    year = i
  )
}



##### Quantitative comparison ####

## Exclude param combos that didn't meet target # PCs
#Some parameter combinations were too restrictive and could not identify
# enough pairs to cross, and therefore should be eliminated
sapply(simulated, function(df) nrow(df))
simulated <- simulated[sapply(simulated, nrow) == 624] 
names(simulated)
#only 6 sims remain:
# "w4_early0.25_lose0.25.csv" "w4_early0.25_lose0.5.csv" 
# "w4_early0.5_lose0.25.csv"  "w4_early0.5_lose0.5.csv"  
# "w8_early0.25_lose0.25.csv" "w8_early0.5_lose0.25.csv" 


#set number of bins for comparison
bins <- seq(1, 348, length.out = 17)
#with 624 fish per sim, 16 bins means ~40 fish per bin, which should be enough
# to compare the distribution of PCs between sim bins and real data bins

#score function using full distributions
score_sim <- function(df) {
  df$bin <- cut(df$sourcePC, breaks = bins, include.lowest = TRUE)
  
  real_df <- real_all
  real_df$bin <- cut(real_df$sourcePC, breaks = bins, include.lowest = TRUE)
  
  ks_vals <- sapply(levels(df$bin), function(b) {
    sim_vals  <- df$newPC[df$bin == b]
    real_vals <- real_df$newPC[real_df$bin == b]
    # skip empty bins
    if (length(sim_vals) < 5 || length(real_vals) < 5) return(NA)
    ks.test(sim_vals, real_vals)$statistic
  })
  mean(ks_vals, na.rm = TRUE)
}

#compare
real_all <- do.call(rbind, real)
sim_scores <- sapply(simulated, score_sim)
sort(sim_scores)
# w4_early0.5_lose0.25.csv   w4_early0.5_lose0.5.csv  w8_early0.5_lose0.25.csv 
# 0.3443286                 0.3449455                 0.3493262 
# w8_early0.25_lose0.25.csv w4_early0.25_lose0.25.csv  w4_early0.25_lose0.5.csv 
# 0.3500951                 0.3527990                 0.3589783 

# w4_early0.5_lose0.25 is the best match -> distribution most similar to real data



##### Plot all ####

pdf("output/availParams/ChoosePairs_weeks4-16_0.25-0.5props.pdf",
    width=7, height=7)

par(mfrow = c(4, 4),
    mar = c(3.5, 3.5, 2, 1),
    mgp = c(2, 0.7, 0))

#real data
yrs <- 2022:2025
letters_real <- letters[1:4]

for (i in seq_along(yrs)) {
  df <- real[[as.character(yrs[i])]]
  plot(newPC ~ sourcePC,
       data = df,
       col = rgb(0, 0, 1, 0.3),
       pch = 16,
       xlab = paste0("Source PC (", yrs[i], ")"),
       ylab = paste0("New PC (", yrs[i] + 1, ")"),
       main = paste0("(", letters_real[i], ") ", yrs[i], " real data"),
       cex.main = 0.8)
}

#simulated data
letters_sim <- letters[5:16]
n <- 1

for (nm in names(simulated)) {
  df <- simulated[[nm]]
  
  # extract parameters from object name
  parts <- strsplit(nm, "_")[[1]]
  
  weeks <- as.integer(gsub("w", "", parts[1]))
  early <- gsub("early", "", parts[2])
  loss  <- gsub("lose|\\.csv", "", parts[3])
  
  plot(newPC ~ sourcePC,
       data = df,
       col = rgb(0, 0, 0, 0.2),
       pch = 16,
       xlab = "Source PC",
       ylab = "New PC",
       main = paste0("(", letters_sim[n], ") ",
                     weeks, " blocks, ",
                     early, " early, ",
                     loss, " loss"),
       cex.main = 0.8)
  legend("bottomright",
         legend = paste0(length(unique(df$newPC)), " pairs"),
         cex = 0.7,
         bty = "n")
  n <- n + 1
}

dev.off()
