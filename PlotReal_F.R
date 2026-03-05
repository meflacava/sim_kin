## Plot real F values from hatchery

setwd("")


#### Plot real F values ####
## Import pedigree
pedF <- read.csv("output/pedigree2026_Fvalues.csv")

plot(F_R~yr, data=pedF,ylim=c(0,0.0625),pch=16,col=rgb(0,0,0,alpha=0.5),
     main="",xlab="year",ylab="Inbreeding coefficient (F)")
abline(h=0.0625,lty=2,col="blue")
abline(h=0.03125,lty=2,lwd=1,col="green3")
#add median lines for each year
meds <- aggregate(F_R~yr,data=pedF,median,na.rm=T)
segments(meds$yr-0.3,meds$F,meds$yr+0.3,meds$F,col="red",lwd=2) #0.3=length of tick
legend(2008,0.055,c("Cousins","2nd cousins","Median F"),cex=0.75,
       lty=c(2,2,1),lwd=c(1,1,2),col=c("blue","green3","red"))


## Without outliers
plot(F_R~yr, data=pedF[pedF$F_R<0.06,],pch=16,col=rgb(0,0,0,alpha=0.5),
     main="",xlab="year",ylab="Inbreeding coefficient (F)")
abline(h=0.03125,lty=2,lwd=1,col="green3")
#add median lines for each year
meds <- aggregate(F_R~yr,data=pedF,median,na.rm=T)
segments(meds$yr-0.3,meds$F,meds$yr+0.3,meds$F,col="red",lwd=2) #0.3=length of tick
#legend(2008,0.013,c("2nd cousins","Median F"),cex=0.8,
#       lty=c(2,1),lwd=c(1,2),col=c("green3","red"))


##### High F fish ####

cuz <- pedF[pedF$F_R>0.06,]
nrow(cuz) #22
length(unique(cuz$Sire)) #10 pair crosses between cousins
