## Data exploration and wrangling

setwd("")


#### Data prep ####

## Create pedigree input file in excel
# - Used "2026 PmX pedigree" from Mandi
# - Copied first 3 columns (PASTE AS VALUES) into new sheet with headers ID, Sire, Dam
# - Deleted first chunk of rows representing parents of wild fish (ID<100000)
# - Save as pedigree2026_simple.csv

## Import
ped.df <- read.csv("data/source/pedigree2026_simple.csv", stringsAsFactors=F)

## Edit the pedigree as needed
#add sex
ped.df$Sex <- as.integer(substr(ped.df$ID,6,6))
#change wild parent IDs to NA
for (i in 1:nrow(ped.df)){
  if (ped.df$Sire[i]<100000){
    ped.df$Sire[i] <- NA
  }
  if (ped.df$Dam[i]<100000){
    ped.df$Dam[i] <- NA
  }
}


## Validate the pedigree
# NOTE that pedtools wants sire Sex=1 and dam Sex=2, but Mandi uses the opposite in her
#  pedigree, so I need to flip them for this validation
for (i in 1:nrow(ped.df)){
  if (ped.df$Sex[i]==1){
    ped.df$Sex[i] <- 2
  }  else if (ped.df$Sex[i]==2){
    ped.df$Sex[i] <- 1
  }
}
library(pedtools)
validatePed(id=ped.df$ID,fid=ped.df$Sire,mid=ped.df$Dam,sex=ped.df$Sex) #no problems


## Add gen and yr to pedigree
ds.gen <- as.integer(substr(ped.df$ID,1,2))
for (i in 1:nrow(ped.df)){
  ped.df$PC[i]<- as.integer(substr(ped.df$ID[i],3,5))
  if (ds.gen[i] %in% seq(10,90,10)){
    ped.df$yr[i] <- 2007 + as.integer(substr(ds.gen[i],1,1))
  } else {
    ped.df$yr[i] <- 2006 + ds.gen[i]
  }
}

## Save [excluding sex column]
#write.csv(ped.df[, -which(names(ped.df)=="Sex")],"data/pedigree2026.csv",row.names=F,quote=F)




##### FCCL hatchery info ####

## Lindberg et al. 2013
#250 PCs, eggs combined into MFGs with 8 PCs each with 80% or more of the initial families 
# represented at the adult stage. Maintain breeding pop of 500
#  -> 6000 offspring become adults / 250 PCs = ~24 offspring surviving to adulthood per PC 
#     (though this is highly variable in practice)

#DS are reared well in excess of the target annual population of 500 breeding individuals,
# starting with >200,000 eggs to ensure >6,000 adult fish from which to select broodstock
# for the refuge pop and also provide fish for research


## Tsai et al. 2022
#40 larval tanks, 5,954-L larval system, 224,000 larvae capacity
# -> 40 tanks x 8 PC per MFG = ~320 PCs?

#Within these multifamily groups, fish are reared at target densities of 5,600 embryos per 
# incubator, 2,500 late-stage larvae per 400-L tank, 1,500 subjuveniles per 1,100-L tank, 
# 1,000 juveniles per 1,100-L tank, and 600 subadults per 1,100-L tank



###### Number of PCs per year ####
ped <- read.csv("data/pedigree2026.csv")

tapply(X=ped$PC, INDEX=ped$yr, FUN=function(x) length(unique(x)))
#2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 
# 164  247  234  256  281  261  224  243  234  262  256  196  245  297  310  313  311  313 



###### Wild parents per generation ####
ped <- read.csv("data/pedigree2026.csv")

sum(is.na(ped$Sire)) #1070

#tapply(X=ped$Sire,INDEX=ped$yr,FUN=function(x) sum(is.na(x))) #alt method
w <- aggregate(is.na(Sire)~yr,data=ped,FUN=sum)
#      yr is.na(Sire)
# 1  2008         328
# 2  2009          59
# 3  2010          35
# 4  2011          70
# 5  2012          48
# 6  2013          84
# 7  2014          63
# 8  2015          48
# 9  2016          54
# 10 2017          83
# 11 2018          70
# 12 2019          32
# 13 2020          88
# 14 2021           5
# 15 2022           0
# 16 2023           2
# 17 2024           1
# 18 2025           0
mean(w$`is.na(Sire)`[w$yr %in% 2009:2020]) #61
sum(w$`is.na(Sire)`[w$yr %in% 2021:2025]) #8 



###### PC availability ####

#When do PCs from last year become available throughout the course of the spawning
# season? Most accurate would be to look at tagging sheet data, but to keep things
# simple, I am going to use the pedigree to just look at which source PCs from one 
# year contributed to which new PCs in the next year
#  - % of families available week or month 1 [break into diff # sessions]
#  - % of families unavailable as season progresses

#set up output
out <- data.frame(
  yr=integer(0),
  n_spawndays=integer(0),
  n_groups=integer(0),
  group=integer(0),
  n_source=integer(0),
  frac_source=integer(0),
  #n_lost=integer(0),
  frac_lost=numeric(0),
  frac_lost_accum=numeric(0)
)

for(yr in 2022:2025) {
  
  #import refuge crosses
  rc_yr <- read.csv(paste0("data/source/",yr,"_refuge_crosses.csv"))
  rc <- rc[rc$PC.FSG>0,]
  rc$Male.PC.FSG <- as.integer(rc$Male.PC.FSG)
  rc$Female.PC.FSG <- as.integer(rc$Female.PC.FSG)
  
  all_fams <- unique(c(rc_yr$Male.PC.FSG, rc_yr$Female.PC.FSG))
  
  #test diff # groups
  group_sizes <- c(4,8,16)
  for(n_groups in group_sizes) {
    
    grp <- cut(seq_len(nrow(rc_yr)),
               breaks=n_groups,
               labels=seq_len(n_groups))
    
    src_by_group <- lapply(levels(grp), function(g) {
      rows <- which(grp == g)
      unique(c(rc_yr$Male.PC.FSG[rows],
               rc_yr$Female.PC.FSG[rows]))
    })
    
    for(g in 1:(n_groups - 1)) {
      
      early <- unique(unlist(src_by_group[1:g]))
      later <- unique(unlist(src_by_group[(g + 1):n_groups]))
      prev <- unique(unlist(src_by_group[g]))
      upcoming <- unique(unlist(src_by_group[g + 1]))
      
      out <- rbind(out,
                   data.frame(
                     yr=yr,
                     n_groups=n_groups,
                     group=g,
                     n_source=length(early),
                     frac_source=round(length(early)/length(all_fams),2),
                     #n_lost=length(setdiff(early, later)),
                     frac_lost=round(length(setdiff(prev, upcoming)) / length(prev),2),
                     frac_lost_accum=round(length(setdiff(early, later)) / length(early),2)
                   )
      )
    }
  }
}

## early_frac AKA % of families available in first session
out$frac_source[out$n_groups==16 & out$group==1]
#4 sessions =  0.39 0.45 0.46 0.45 -> rep with early_frac=0.5
#8 sessions =  0.20 0.25 0.24 0.27 -> rep with early_frac=0.25
#16 sessions = 0.12 0.14 0.13 0.14 -> rep with early_frac=0.125

## loss AKA % of families lost each session
aggregate(frac_lost~n_groups,data=out,FUN=median)
#   n_groups frac_lost
# 1        4     0.635 -> rep with loss=0.5
# 2        8     0.800 -> rep with loss=0.75
# 3       16     0.890 -> rep with loss=0.875
aggregate(frac_lost~yr+n_groups,data=out,FUN=median)
#      yr n_groups frac_lost
# 1  2022        4      0.56
# 2  2023        4      0.63
# 3  2024        4      0.66
# 4  2025        4      0.57
# 5  2022        8      0.79
# 6  2023        8      0.80
# 7  2024        8      0.85
# 8  2025        8      0.77
# 9  2022       16      0.89
# 10 2023       16      0.87
# 11 2024       16      0.92
# 12 2025       16      0.88
