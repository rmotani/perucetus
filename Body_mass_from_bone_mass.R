# Set your working directly to the unzipped folder containing this code

# BEFORE RUNNING THE SCRIPT
#
# We could not include three of the data files in this distribution because of 
# copyright issues. You need to download them. Make sure to place in the 
# working directory.
#
# -From https://github.com/n8upham/MamPhy_v1/tree/master/_DATA,
# downlaod MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre
#
# -From https://www.nature.com/articles/s41586-023-06381-1#Sec17, 
# download 
# -AquaticAndEstimations.csv
# -TerrestrialExtant.csv
# which are contained in Supplementary Data 2. Place these two files directly 
# in your working directory, without the folder structure in the zip file.

library(ape)
library(phytools)
library(scales)
library(nlme)

scl <- NULL
for(i in -5:11) scl <- c(scl,seq(10^i,9*10^i,10^i))
log10scl <- log10(scl)
logscl <- log(scl)
SCL <- 10^seq(-5,11)
log10SCL <- log10(SCL)
logSCL <- log(SCL)
bmscl <- seq(10000,400000,10000)
log10bmscl <- log10(bmscl)


###
### DATA
###

### Terrestrial data from Bianucci et al.(2023)
LAND <- read.csv("TerrestrialExtant.csv")

### Aquatic data from Bianucci et al. (2023)
WATER <- read.csv("AquaticAndEstimations.csv")

### Perucetus skeletal mass range given by Bianucci et al.
peru <- data.frame(SkeletalMass=log10(c(5302.278152,6465.529349,7628.780546)))
peru.dom <- data.frame(boneM=log10(c(5302.278152,6465.529349,7628.780546)))

### Manatee data from Domning and Buffrenil (1991)
Dom <- DOM <- read.csv("Domning_Buffrenil_1991.csv")
Dom[,sapply(DOM,is.numeric)] <- log10(DOM[,sapply(DOM,is.numeric)])
dom <- na.omit(Dom[,c(1,2,4,10)])

### Selecting data columns in old fashion
Land <- LAND[,colnames(LAND)=="Taxon" | colnames(LAND)=="BodyMass" | 
              colnames(LAND)=="SkeletalMass"]
Water <- WATER[WATER$Cet4Reg=="1",colnames(WATER)=="Taxon" | 
              colnames(WATER)=="BodyMass" | colnames(WATER)=="SkeletalMass"]
Odonto <- WATER[WATER$Clade=="1" & WATER$Cet4Reg=="1",colnames(WATER)=="Taxon" | 
              colnames(WATER)=="BodyMass" | colnames(WATER)=="SkeletalMass"]
Mysti <- WATER[WATER$Clade=="2" & WATER$Cet4Reg=="1",colnames(WATER)=="Taxon" | 
              colnames(WATER)=="BodyMass" | colnames(WATER)=="SkeletalMass"]
Sirenia <- SIRENIA <- WATER[WATER$Clade=="3",colnames(WATER)=="Taxon" | 
              colnames(WATER)=="BodyMass" | colnames(WATER)=="SkeletalMass"]

### Species means
land <- Land <- aggregate(.~Taxon,data=Land,FUN=mean)
water <- Water <- aggregate(.~Taxon,data=Water,FUN=mean)
odonto <- Odonto <- aggregate(.~Taxon,data=Odonto,FUN=mean)
mysti <- Mysti <- aggregate(.~Taxon,data=Mysti,FUN=mean)
sirenia <- Sirenia <- aggregate(.~Taxon,data=Sirenia,FUN=mean)
manatee <- Manatee <- SIRENIA[-1,]

### Log10 transformation of numeric variables
land[,sapply(Land,is.numeric)] <- log10(Land[,sapply(land,is.numeric)])
water[,sapply(Water,is.numeric)] <- log10(Water[,sapply(water,is.numeric)])
odonto[,sapply(Odonto,is.numeric)] <- log10(Odonto[,sapply(Odonto,is.numeric)])
mysti[,sapply(Mysti,is.numeric)] <- log10(Mysti[,sapply(Mysti,is.numeric)])
sirenia[,sapply(Sirenia,is.numeric)] <- log10(Sirenia[,sapply(Sirenia,is.numeric)])
manatee[,sapply(Manatee,is.numeric)] <- log10(Manatee[,sapply(Manatee,is.numeric)])


###
### Non-phylogenetic Ordinary Least Square Regression
###

### Linear models and their summary
lm.land <- lm(BodyMass~SkeletalMass,land)
lm.water <- lm(BodyMass~SkeletalMass,water)
lm.odonto <- lm(BodyMass~SkeletalMass,odonto)
lm.mysti <- lm(BodyMass~SkeletalMass,mysti)
lm.sirenia <- lm(BodyMass~SkeletalMass,sirenia)
lm.manatee <- lm(BodyMass~SkeletalMass,manatee)
lm.dom <- lm(bodyM~boneM,dom)
summary(lm.land)
summary(lm.water)
summary(lm.odonto)
summary(lm.mysti)
summary(lm.sirenia)
summary(lm.manatee)
summary(lm.dom)

### Mean absolute errors
mean(abs(10^predict(lm.land) - Land$BodyMass)/Land$BodyMass)
mean(abs(10^predict(lm.water) - Water$BodyMass)/Water$BodyMass)
mean(abs(10^predict(lm.sirenia) - Sirenia$BodyMass)/Sirenia$BodyMass)
mean(abs(10^predict(lm.odonto) - Odonto$BodyMass)/Odonto$BodyMass)
mean(abs(10^predict(lm.manatee) - Manatee$BodyMass)/Manatee$BodyMass)
mean(abs(10^predict(lm.dom) - DOM$bodyM)/DOM$bodyM)

### Perucetus mass estimates
10^predict(lm.water,peru)
10^predict(lm.land,peru)
10^predict(lm.sirenia,peru)
10^predict(lm.odonto,peru)
10^predict(lm.manatee,peru,interval="confidence")
10^predict(lm.manatee,peru,interval="prediction")
10^predict(lm.dom,peru.dom,interval="confidence")
10^predict(lm.dom,peru.dom,interval="prediction")

    
###
### Phylogenetic Generalized Least Square Regression
###

# The tree below is from https://github.com/n8upham/MamPhy_v1/tree/master/_DATA
#<---- Copying the scripts from PGLS.R of Bianucci et al.(2023)

#### Tree loading and tip adjustments ####

treeAll <- read.nexus("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")# from Upham et al. (2019): http://vertlife.org/data/mammals/
allLabels <- treeAll$tip.label; shortLabels <- c()
for (i in 1:length(allLabels)){
  shortLabels[i] <- paste(unlist(strsplit(allLabels[i], split = "_"))[1], unlist(strsplit(allLabels[i], split = "_"))[2],sep="_")
}

#Species match corrections
shortLabels[grep("Sorex_ari",shortLabels)] <- "Soricid_sp"
shortLabels[grep("Peromyscus_fraterculus",shortLabels)] <- "Peromyscus_sp"
shortLabels[grep("Tamias_striatus",shortLabels)] <- "Eutamias_sp"
shortLabels[grep("Spermophilus_pygmaeus",shortLabels)] <- "Spermophilus_sp"
shortLabels[grep("Thomomys_talpoides",shortLabels)] <- "Thomomys_sp"
shortLabels[grep("Castor_fiber",shortLabels)] <- "Castor_sp"
shortLabels[grep("Canis_lupus",shortLabels)] <- "Canis_sp"
shortLabels[grep("Mustela_strigidorsa",shortLabels)] <- "Mustela_sp"
shortLabels[grep("Spilogale_putorius",shortLabels)] <- "Spilogale_sp"
shortLabels[grep("Felis_nigripes",shortLabels)] <- "Felis_sp"

#Update tip names  
treeAll$tip.label <- shortLabels
#-----> End copying the scripts from PGLS.R of Bianucci et al. (2023)


#Trees cropped with sample--modified from Bianucci et al. (2023)
tre.land <-  keep.tip(treeAll,land$Taxon)
tre.water <-  keep.tip(treeAll,water$Taxon)
tre.water <- read.newick("cetacea.nwk")
tre.water$tip.label[tre.water$tip.label=="Physeter_catodon"] <- "Physeter_macrocephalus"
tre.odonto <- keep.tip(tre.water,odonto$Taxon)
tre.mysti <- keep.tip(tre.water,mysti$Taxon)

### Terrestrial phylogenetic GLS, Brownian
gls1 <- gls(BodyMass~SkeletalMass,data = land,correlation=corPagel(value=1,phy=tre.land, fixed=T, form=~Taxon),method="ML")
summary(gls1)
mean(abs(10^predict(gls1)-Land$BodyMass)/Land$BodyMass)
range(abs(10^predict(gls1)-Land$BodyMass)/Land$BodyMass)
10^predict(gls1,peru)

### Cetacea phylogenetic GLS, Brownian
gls2 <- gls(BodyMass~SkeletalMass,data = water,correlation=corPagel(value=1,phy=tre.water, fixed=T, form=~Taxon),method="ML")
summary(gls2)
mean(abs(10^predict(gls2)-10^water$BodyMass)/(10^water$BodyMass))
range(abs(10^predict(gls2)-10^water$BodyMass)/(10^water$BodyMass))
10^predict(gls2,peru)

### Odontoceti phylogenetic GLS, Brownian
gls3 <- gls(BodyMass~SkeletalMass,data = odonto,correlation=corPagel(value=1,phy=tre.odonto, fixed=T, form=~Taxon),method="ML")
summary(gls3)
mean(abs(10^predict(gls3)-10^odonto$BodyMass)/(10^odonto$BodyMass))
range(abs(10^predict(gls3)-10^odonto$BodyMass)/(10^odonto$BodyMass))
10^predict(gls3,peru)

### Mysticeti phylogenetic GLS, Brownian
gls4 <- gls(BodyMass~SkeletalMass,data = mysti,correlation=corPagel(value=1,phy=tre.mysti, fixed=T, form=~Taxon),method="ML")
summary(gls4)
mean(abs(10^predict(gls4)-10^mysti$BodyMass)/(10^mysti$BodyMass))
range(abs(10^predict(gls4)-10^mysti$BodyMass)/(10^mysti$BodyMass))
10^predict(gls4,peru)

### Cetacea isometry assuming the mean skeleton/body mass ratio
cetmean <- colMeans(water[,sapply(water,is.numeric)])
cet.iso <- function(x) x+cetmean[2]-cetmean[1]
mean(abs(10^cet.iso(water$SkeletalMass)-10^water$BodyMass)/(10^water$BodyMass))
range(abs(10^cet.iso(water$SkeletalMass)-10^water$BodyMass)/(10^water$BodyMass))
10^cet.iso(peru)

### Cetacea isometry assuming the minimum skeleton/body mass ratio, 
### resulting in maximum body mass estimates
max.iso <- function(x) x+water[8,3]-water[8,2]
mean(abs(10^max.iso(water$SkeletalMass)-10^water$BodyMass)/(10^water$BodyMass))
range(abs(10^max.iso(water$SkeletalMass)-10^water$BodyMass)/(10^water$BodyMass))
10^max.iso(peru)

### Cetacea isometry assuming the maximum skeleton/body mass ratio, 
### resulting in minimum body mass estimates
min.iso <- function(x) x+water[11,3]-water[11,2]
mean(abs(10^min.iso(water$SkeletalMass)-10^water$BodyMass)/(10^water$BodyMass))
range(abs(10^min.iso(water$SkeletalMass)-10^water$BodyMass)/(10^water$BodyMass))
10^min.iso(peru)

### PGLS Results Summary
rslt <- rbind(c(coef(gls1),10^predict(gls1,peru)),
      c(coef(gls2),10^predict(gls2,peru)),
      c(coef(gls3),10^predict(gls3,peru)),
      c(coef(gls4),10^predict(gls4,peru))
      )
rownames(rslt) <- c("Terrestrial","Cetacea","Odontoceti","Mysticeti")
colnames(rslt) <- c("Intercept","Slope","BM form Min SM","BM from mean SM","BM from Max SM")
write.csv(rslt,file="PGLS_summary.csv")

### Plotting
x11(h=9,w=8)
#pdf(h=9,w=8,file="PGLS_9.pdf")
layout(matrix(c(1,2,3,3),2,2))
  plot(BodyMass~SkeletalMass,land,xlim=c(-3,5),ylim=c(-2,6),typ="n",axes=F,
         xlab="Skeletal Mass (kg)", ylab="Body Mass (kg)")
    axis(1,at=log10SCL,label=SCL,las=3)
    axis(2,at=log10SCL,label=SCL,las=1)
    box()
    rect(log10(5302.278152),-3,log10(7628.780546),7,col=scales::alpha("blue",0.2),border=F)
    rect(-4,log10(98000),5,log10(114000),col=scales::alpha("red",0.2),border=F)
    rect(-4,log10(71500),5,log10(82900),col=scales::alpha("red",0.2),border=F)
    rect(-4,log10(41300),5,log10(48000),col=scales::alpha("red",0.2),border=F)
    rect(log10(4000),log10(30000),log10(10000),log10(390000),lwd=2)
    points(BodyMass~SkeletalMass,land,pch=18)
    #    points(BodyMass~SkeletalMass,water,col=4,pch=16)
    points(BodyMass~SkeletalMass,mysti,col=2,pch=17)
    points(BodyMass~SkeletalMass,odonto,col=3,pch=17)
    # points(BodyMass~SkeletalMass,manatee,col=5,pch=17)
    # points(bodyM~boneM,dom,col=6,pch=17)
    abline(gls1,lwd=2,col="Maroon")
    abline(gls2,lwd=3,col=1)
    abline(gls3,lwd=2,col=3)
    abline(gls4,lwd=2,col=2)
    abline(lm.manatee,lwd=2,col=5)
    abline(lm.dom,lwd=2,col=6)
    abline(b=1,a=cetmean[2]-cetmean[1],lty=5,col=4,lwd=2)
    abline(b=1,a=water[8,3]-water[8,2],lty=5,col=1,lwd=2)
    abline(v=log10(6465.529349),col=4)
    abline(h=log10(c(30000,390000)))
  plot(bodyM~boneM,dom,xlim=c(-3,5),ylim=c(-2,6),typ="n",axes=F,
         xlab="Skeletal Mass (kg)", ylab="Body Mass (kg)")
    axis(1,at=log10SCL,label=SCL,las=3)
    axis(2,at=log10SCL,label=SCL,las=1)
    box()
    points(bodyM~boneM,dom,col=6,pch=17)
    points(BodyMass~SkeletalMass,manatee,col=5,pch=17,cex=0.6)
    abline(lm.manatee,lwd=2,col=5)
    abline(lm.dom,lwd=2,col=6)
  plot(BodyMass~SkeletalMass,land,xlim=c(3.6,4),ylim=c(4.5,5.55),axes=F,
         xlab="Skeletal Mass (kg)", ylab="Body Mass (kg)")
    axis(1,at=log10scl,label=scl,las=3)
    axis(2,at=log10bmscl,label=bmscl,las=1)
    box()
    abline(h=log10(bmscl),col="gray")
    rect(log10(5302.278152),-3,log10(7628.780546),7,col=scales::alpha("blue",0.2),border=F)
    rect(-4,log10(98000),5,log10(114000),col=scales::alpha("red",0.2),border=F)
    rect(-4,log10(71500),5,log10(82900),col=scales::alpha("red",0.2),border=F)
    rect(-4,log10(41300),5,log10(48000),col=scales::alpha("red",0.2),border=F)
    points(BodyMass~SkeletalMass,land,pch=18)
    #    points(BodyMass~SkeletalMass,water,col=4,pch=16)
    # points(BodyMass~SkeletalMass,mysti,col=2,pch=17)
    # points(BodyMass~SkeletalMass,odonto,col=3,pch=17)
    points(bodyM~boneM,dom,col=6,pch=17)
    abline(gls1,lwd=2,col="Maroon")
    abline(gls2,lwd=3,col=1)
    abline(gls3,lwd=2,col=3)
    abline(gls4,lwd=2,col=2)
    abline(lm.manatee,lwd=2,col=5)
    abline(lm.dom,lwd=2,col=6)
    cetmean <- colMeans(water[,sapply(water,is.numeric)])
    abline(b=1,a=cetmean[2]-cetmean[1],lty=5,col=4,lwd=2)
    abline(b=1,a=water[8,3]-water[8,2],lty=5,col=1,lwd=2)
    abline(v=log10(6465.529349),col=4)
#dev.off()
    