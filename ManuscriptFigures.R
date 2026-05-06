#Code for figures


#### Color scheme ####
#library(RColorBrewer)
#display.brewer.all(type="qual")
#barplot(1:8,col=brewer.pal(8,"Dark2"))
#cols <- brewer.pal(8,"Dark2")
#cols <- cols[c(1,2,3,6)]
#barplot(1:4,col=cols)

#F thresholds
cols.f <- c("#4D4D4D","#7F7F7F","#B3B3B3")
lty.f <- c(5,2,4)
# "#4D4D4D" darkest gray, lty=5 = 1st cousin (0.0625)
# "#7F7F7F" medium gray, lty=2 = 2nd cousin (0.03125)
# "#B3B3B3" lightest gray, lty=3 = 3rd cousin (0.015625)
plot(NA, xlim=c(0,3), ylim=c(0,3))
for(i in 1:3){abline(h=i, col=cols.f[i], lty=lty.f[i], lwd=3)}

#Fig 1b - random vs. genetically informed pairing
cols.gr <- c("#1F78B4","#E6AB02")
# "#1F78B4" blue = main sim
# "#E6AB02" gold = random mating
barplot(1:2,col=cols.gr)

#Fig ? - varying # PCs
cols.pc <- c("#1B9E77","#D95F02","#7570B3")
#cols[1] "#1B9E77" green = 312 PCs
#cols[2] "#D95F02" orange = 400 PCs 
#cols[3] "#7570B3" purple = 500 PCs
barplot(1:3,col=cols.pc)

barplot(1:5,col=c(cols.gr,cols.pc))


#### Fig 1 - real F values ####

## Import pedigree
pedF <- read.csv("output/pedigree2026_Fvalues.csv")

## Fig 1a
pdf("output/Fig1.pdf",height=4,width=8)
par(mfrow=c(1,2),
    mar=c(3,3.5,1,0.4)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)
plot(F_R~yr, data=pedF,ylim=c(0,0.0625),pch=16,col=rgb(0,0,0,alpha=0.5),
     main="",xlab="",ylab="",xaxt="n",yaxt="n")
axis(1,at=seq(2010,2025,5),labels=F,tick=T)
axis(1,at=seq(2010,2025,5),line=-0.5,tick=F)
mtext("year", side = 1, line = 1.6) 
axis(2,at=seq(0,0.0625,0.02),las=1,hadj=0.75)
mtext("Inbreeding coefficient (F)", side = 2, line = 2.5)
abline(h=0.0625,lty=lty.f[1],col=cols.f[1])
abline(h=0.03125,lty=lty.f[2],lwd=1,col=cols.f[2])
#add median lines for each year
meds <- aggregate(F_R~yr,data=pedF,median,na.rm=T)
segments(meds$yr-0.3,meds$F,meds$yr+0.3,meds$F,col="red",lwd=2) #0.3=length of tick
legend(2008,0.055,c("Median F","Cousins","2nd cousins"),cex=0.8,bty="n",
       lty=c(1,lty.f[1:2]),lwd=c(2,1,1),col=c("red",cols.f[1:2]))
mtext("(a)",font=2,side=2,line=2.5,adj=1,las=1,padj=-14) #panel label

## Fig 1b
plot(F_R~yr, data=pedF[pedF$F_R<0.06,],pch=16,col=rgb(0,0,0,alpha=0.5),
     main="",xlab="",ylab="",xaxt="n",yaxt="n")
axis(1,at=seq(2010,2025,5),labels=F,tick=T)
axis(1,at=seq(2010,2025,5),line=-0.5,tick=F)
mtext("year", side = 1, line = 1.6) 
axis(2,at=seq(0,0.02,0.005),las=1,hadj=0.75)
#mtext("Inbreeding coefficient (F)", side = 2, line = 2.7)
abline(h=0.03125,lty=2,lwd=1,col="green3")
#add median lines for each year
meds <- aggregate(F_R~yr,data=pedF,median,na.rm=T)
segments(meds$yr-0.3,meds$F,meds$yr+0.3,meds$F,col="red",lwd=2) #0.3=length of tick
#legend(2008,0.013,c("2nd cousins","Median F"),cex=0.8,
#       lty=c(2,1),lwd=c(1,2),col=c("green3","red"))
mtext("(b)",font=2,side=2,line=2.6,adj=1,las=1,padj=-14) #panel label
dev.off()




#### Fig 3a - realistic sim ####
ped.plot <- read.csv("output/sim_50gen_312pc_kinSteps_avail_4weeks.csv")
table(ped.plot$yr)
length(unique(ped.plot$yr))-length(unique(ped.plot$yr[ped.plot$yr<2026])) #check # successful sim gens

pdf("output/Fig3.pdf",height=4,width=9)
par(mar=c(3,3.5,1,0.5)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2), nrow = 1),widths=c(1.5,1))
#layout.show(2)
#par(def.par)
plot(F_R~yr, data=ped.plot[ped.plot$yr<2026,],
     xlim=c(min(ped.plot$yr),max(ped.plot$yr)),ylim=c(0,0.0625),
     pch=16,col="gray",xaxt="n",yaxt="n",
     xlab="",ylab="")
points(F_R~yr, data=ped.plot[ped.plot$yr>2025,],pch=16,col="black")
axis(1,at=seq(2010,2070,10),labels=F,tick=T)
axis(1,at=seq(2010,2070,10),line=-0.5,tick=F)
mtext("year", side = 1, line = 1.6) 
axis(2,at=seq(0,0.0625,0.02),las=1,hadj=0.75)
mtext("Inbreeding coefficient (F)", side = 2, line = 2.5) 
abline(h=0.0625,lty=lty.f[1],lwd=1,col=cols.f[1])
abline(h=0.03125,lty=lty.f[2],lwd=1,col=cols.f[2])
#add median lines for each year
meds <- aggregate(F_R~yr,data=ped.plot,median,na.rm=T)
segments(meds$yr-0.3,meds$F,meds$yr+0.3,meds$F,col="red",lwd=2) #0.3=length of tick
legend("bottomright",c("Real","Simulated","Median F","Cousins","2nd cousins"),cex=0.8,
       pch=c(16,16,NA,NA,NA),lty=c(NA,NA,1,lty.f[1:2]),lwd=c(NA,NA,2,1,1),
       col=c("gray","black","red",cols.f[1:2]),bty="n")
mtext("(a)",font=2,side=2,line=2,adj=1,las=1,padj=-13) #panel label


#### Fig 3b - gen man vs. random ####
summ <- function(ped){
  yrs <- sort(unique(ped$yr))
  data.frame(
    yr = yrs,
    med = sapply(yrs,function(y) median(ped$F_R[ped$yr==y])),
    lo  = sapply(yrs,function(y) quantile(ped$F_R[ped$yr==y],0.25)),
    hi  = sapply(yrs,function(y) quantile(ped$F_R[ped$yr==y],0.75))
  )
}

## List of sim files you want to analyze
sims <- list(gen=read.csv("output/sim_50gen_312pc_kinSteps_avail_4weeks.csv"),
             random=read.csv("output/sim_50gen_312pc_random_avail_4weeks.csv"))
sims <- lapply(sims, function(df) df[df$yr>2025,])

#calc medians and IQRs
S <- lapply(sims, summ)

## Gen when median F exceeded a threshold
lapply(S, function(df) df$yr[df$med > 0.03125][1])
#genetically managed = 2064
#random = 2054

## Plot
par(mar=c(3,1,1,0.1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)
#ylim <- range(unlist(lapply(S,function(x) c(x$lo,x$hi))))
ylim <- c(0,0.0625)
plot(S[[1]]$yr, S[[1]]$med, type="n",xaxt="n",yaxt="n",
     ylim=ylim, xlab="year", ylab="")
for(i in seq_along(S)){
  x <- S[[i]]
  # polygon(c(x$yr, rev(x$yr)),
  #         c(x$lo, rev(x$hi)),
  #         col=adjustcolor(cols.gr[i], alpha.f=0.2),
  #         border=NA)
  lines(x$yr, x$med, col=cols.gr[i], lwd=2)
}
axis(1,at=seq(2030,2070,10),labels=F,tick=T)
axis(1,at=seq(2030,2070,10),line=-0.5,tick=F)
mtext("year", side = 1, line = 1.6)
#axis(2,at=seq(0,0.0625,0.02),las=1,hadj=0.75)
abline(h=0.0625,lty=lty.f[1],col=cols.f[1])
abline(h=0.03125,lty=lty.f[2],col=cols.f[2])
legend("bottomright", bty="n",
       legend=c("Genetically informed pairing","Random pairing"),
       col=cols.gr, lty=1, lwd=2,cex=0.8)
       #legend=c("Genetically informed pairing","Random pairing","Cousins","2nd cousins"),
       #col=c(cols.gr,cols.f[1:2]), lty=c(1,1,lty.f[1:2]), lwd=c(2,2,1,1),cex=0.8)
mtext("(b)",font=2,side=2,line=0.15,adj=1,las=1,padj=-13) #panel label
dev.off()


##### Slope of median F ####
slopes <- data.frame(sim=names(S),slope=NA,compto1=NA)
for(i in seq_along(S)){
  slopes$slope[i] <- coef(lm(med ~ yr, data = S[[i]]))["yr"]
  slopes$compto1[i] <- 1 - slopes$slope[i]/slopes$slope[1]
}
slopes
#     sim        slope    compto1
#1    gen 0.0007807341  0.0000000
#2 random 0.0010346194 -0.3251879



#### Fig 4 - vary # PCs ####
summ <- function(ped){
  yrs <- sort(unique(ped$yr))
  data.frame(
    yr = yrs,
    med = sapply(yrs,function(y) median(ped$F_R[ped$yr==y])),
    lo  = sapply(yrs,function(y) quantile(ped$F_R[ped$yr==y],0.25)),
    hi  = sapply(yrs,function(y) quantile(ped$F_R[ped$yr==y],0.75))
  )
}

## List of sim files you want to analyze
sims <- list(pc312=read.csv("output/sim_35gen_312pc_kinSteps_avail_4weeks.csv"),
             pc400=read.csv("output/sim_35gen_400pc_kinSteps_avail_4weeks.csv"),
             pc500=read.csv("output/sim_35gen_500pc_kinSteps_avail_4weeks.csv"))
sims <- lapply(sims, function(df) df[df$yr>2025,])

#calc medians and IQRs
S <- lapply(sims, summ)

## Gen when median F exceeded a threshold
lapply(S, function(df) df$yr[df$med > 0.03125][1])
#312 = 2054, 400 = NA, 500 = NA
lapply(S, function(df) df$yr[df$med > 0.01625][1])
#312 = 2042, 400 = 2045, 500 = 2047

## Plot
pdf("output/Fig4.pdf",height=4,width=5)
par(mar=c(3,3.5,1,0.5)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)
ylim <- range(unlist(lapply(S,function(x) c(x$lo,x$hi))))
#ylim <- c(0,0.0625)
plot(S[[1]]$yr, S[[1]]$med, type="n",
     ylim=ylim, xlab="", ylab="",xaxt="n",yaxt="n")
for(i in seq_along(S)){
  x <- S[[i]]
  # polygon(c(x$yr, rev(x$yr)),
  #         c(x$lo, rev(x$hi)),
  #         col=adjustcolor(cols.pc[i], alpha.f=0.2),
  #         border=NA)
  lines(x$yr, x$med, col=cols.pc[i], lwd=2)
}
axis(1,at=seq(2025,2060,5),labels=F,tick=T)
axis(1,at=seq(2025,2060,5),line=-0.5,tick=F)
mtext("year", side = 1, line = 1.6) 
axis(2,at=seq(0,0.03125,0.01),las=1,hadj=0.75)
mtext("Inbreeding coefficient (F)", side = 2, line = 2.5) 
#abline(h=0.0625,lty=lty.f[1],col=cols.f[1])
abline(h=0.03125,lty=lty.f[2],col=cols.f[2])
abline(h=0.015625,lty=lty.f[3],col=cols.f[3])
legend("bottomright", legend=c("312 crosses","400 crosses","500 crosses","2nd cousins","3rd cousins"),cex=0.8,
       col=c(cols.pc,cols.f[2:3]), lty=c(1,1,1,lty.f[2:3]), lwd=c(2,2,2,1,1),bty="n")
dev.off()


##### Slope of median F ####
slopes <- data.frame(sim=names(S),slope=NA,compto1=NA)
for(i in seq_along(S)){
  slopes$slope[i] <- coef(lm(med ~ yr, data = S[[i]]))["yr"]
  slopes$compto1[i] <- 1 - slopes$slope[i]/slopes$slope[1]
}
slopes
#    sim        slope   compto1
#1 pc312 0.0007725320 0.0000000
#2 pc400 0.0006065004 0.2149188
#3 pc500 0.0004985754 0.3546217




##### Identical parameters comparison ####
#Verify that running multiple iterations of identical parameters produces 
# similar F trends
its <- list.files(path="./output/iterations",
                  pattern="^sim_50gen_312pc_kinSteps_avail_4weeks")
meds <- data.frame(file=character(0),
                   rep=integer(),
                   yr=integer(),
                   median=numeric())
for (i in its){
  x <- read.csv(paste0("output/iterations/",i))
  x <- x[x$yr>2025,]
  m <- aggregate(F_R~yr,data=x,FUN=median)
  meds <- rbind(meds, 
                data.frame(file=i,
                           rep=as.integer(sub(".*_(\\d{3})\\.csv$", "\\1", i)),
                           yr=m$yr,
                           median=m$F_R))
}

#plot trend in median values across years for each sim rep
plot(median~yr,data=meds,type="l",col=rgb(0,0,0,0.2),
     main="",ylab="median F among reps")

## Differences among reps for each year
res <- do.call(rbind, lapply(split(meds, meds$yr), function(d) {
  x <- d$median
  
  diffs <- abs(outer(x, x, "-"))
  diffs <- diffs[upper.tri(diffs)]
  
  rel <- diffs / median(x)
  
  c(
    median_abs_diff = median(diffs),
    mean_abs_diff = mean(diffs),
    min_abs_diff  = min(diffs),
    max_abs_diff  = max(diffs),
    median_rel_diff = median(rel),
    mean_rel_diff = mean(rel),
    min_rel_diff  = min(rel),
    max_rel_diff  = max(rel)
  )
}))
res <- data.frame(yr = as.numeric(rownames(res)), res, row.names = NULL)
#could report:
options(scipen = 999)
c(
  median_of_medians = median(res$median_abs_diff),
  mean_of_means = mean(res$mean_abs_diff),
  max_of_max    = max(res$max_abs_diff),
  median_rel = median(res$median_rel_diff),
  mean_rel      = mean(res$mean_rel_diff),
  max_rel       = max(res$max_rel_diff)
)
options(scipen = 0) #restore
#median_of_medians     mean_of_means        max_of_max        median_rel          mean_rel           max_rel 
#    0.00007329675     0.00019024396     0.00103180481     0.00353065116     0.00985725843     0.12260857280 

#plot
boxplot(lapply(split(meds, meds$yr), function(d) {
  x <- d$median
  diffs <- abs(outer(x, x, "-"))
  diffs <- diffs[upper.tri(diffs)]
  diffs / median(x)
}),
names = sort(unique(meds$yr)),
xlab = "Year", ylab = "Relative difference")
