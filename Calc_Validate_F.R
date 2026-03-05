## Validate R-based F values match Anne's software output

setwd("")

#### Validate F values ####
# Check to make sure that F values calculated in R match those calculated by Anne
# (both use the table method, so they should match) - need to make sure they match
# so that I can forward simulate the pedigree and calc F values to track them


##### Calc F in R ####

## Load required packages
#install.packages("optiSel")
Sys.setenv(RGL_USE_NULL = "TRUE")
library(optiSel)

## Import pedigree
ped.df <- read.csv("data/pedigree2026.csv")

## Calculate F using tabular method
# F(fish) = kinship(sire,dam) AKA the F value is the F value for the fish
#  in the pedigree that resulted from its parents' shared ancestry
ped.df$F_R <- pedInbreeding(ped.df)[,2]



##### Match up F values between R and Anne's software ####

## Generate genID to add to Anne's df
# Exported crosses from Anne do not have genIDs, so I need to re-create them
anne <- read.csv("data/source/refuge_crosses_2008-2025.csv",stringsAsFactors=F)
for (i in 1:nrow(anne)){
  if (anne$yr[i] %in% 2008:2016){
    gen <- anne$yr[i] - 2007
    anne$Sire[i] <- paste0(gen,formatC(anne$PC.FSG[i], width = 4, format = "d", flag = "0"),2)
    anne$Dam[i] <- paste0(gen,formatC(anne$PC.FSG[i], width = 4, format = "d", flag = "0"),1)
  } else {
    gen <- anne$yr[i] - 2006
    anne$Sire[i] <- paste0(gen,formatC(anne$PC.FSG[i], width = 3, format = "d", flag = "0"),2)
    anne$Dam[i] <- paste0(gen,formatC(anne$PC.FSG[i], width = 3, format = "d", flag = "0"),1)
  }
}
anne$Sire <- as.integer(anne$Sire)
anne$Dam <- as.integer(anne$Dam)

#save
#write.csv(anne,"data/refuge_crosses_2008-2025_withGenID.csv",quote=F, row.names=F)



## Add Anne's F values to pedigree
# Anne's F values are actually the relatedness values between pair crosses, so they are all
#  2x the pair cross kinship coefficient and therefore 2x the offspring's inbreeding
#  coefficient. I will 1/2 Anne's values to make them comparable to the pedigree-based ones
#  I calculated using optiSel (which are also how the forward simulations will go)
anne <- read.csv("data/refuge_crosses_2008-2025_withGenID.csv",stringsAsFactors=F)

#Create new object
pedF <- ped.df

#Reduce Anne's values by 1/2 and add to pedigree
for (i in 1:nrow(pedF)){
  if (!is.na(pedF$Sire[i])){
    pedF$F_anne[i] <- anne$F[anne$Sire==pedF$Sire[i]]/2
  } else {
    pedF$F_anne[i] <- 0
  }
}
#write.csv(pedF, "output/pedigree2026_Fvalues.csv", row.names=FALSE)


##### Compare F values ####

pedF <- read.csv("output/pedigree2026_Fvalues.csv")
range(pedF$F_R) #0.00000000 0.06367883
range(pedF$F_anne) #0.000000 0.063635

## Look at weird values
#pedF[pedF$F_R==0.25,]
#         ID   Sire    Dam index gen   yr F_R F_anne
#3096 701222 600082 600081  3096  70 2014 0.25      0
# In the pedigree, these parents are siblings (both from parents 501051/2),
#  but in refuge_crosses, 600081 is a wild fish. So the pedigree is wrong, but Anne's
#  calculated F value is correct because she used the refuge cross lists, not the
#  pedigree. I updated the pedigree to list NA for 600081 parents and re-generated files


## Compare mismatched values
c <- pedF[pedF$F_anne>0,]
c$diff <- round(c$F_R/c$F_anne,1)
hist(c$diff)
range(c$diff)
table(c$yr,c$diff)
nrow(c[c$diff<1 & c$diff>0.8,]) #79 rows with diff 0.8-0.9
nrow(c[c$diff<0.8,]) #28 rows with diff 0.4-0.7
c[c$diff<0.8,]
