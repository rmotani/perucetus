# Set your working directly to the unzipped folder containing this code

library(ape)
library(phytools)
library(scales)
library(nlme)

### DATA

Dat <- DAT <- read.csv("Whale_mass_length.csv")
Dat[,sapply(DAT,is.numeric)] <- log10(DAT[,sapply(DAT,is.numeric)])
dat <- Dat[substr(Dat$Source,1,4) == "Ohno" | 
             substr(Dat$Source,1,4) == "Omur" |
             substr(Dat$Source,1,4) == "Nish" |
             substr(Dat$Source,1,4) == "Robi" |
             substr(Dat$Source,1,4) == "Buff" |
             substr(Dat$Source,1,4) == "Bian" |
             substr(Dat$Source,1,4) == "Fuji",]


### Adult Delphinus delphis means 

Buf <- BUF <- read.csv("Buffrenil_EA_1985.csv")
Buf[,sapply(BUF,is.numeric)] <- log10(BUF[,sapply(BUF,is.numeric)])
buf <- Buf[(BUF$Age >= 12 & BUF$Sex == "M")|(BUF$Age >= 7 & BUF$Sex == "F"),5:6]
buf.means <- colMeans(na.omit(buf))


### PGLS

tre.dry <- read.newick(file="Cetacea_dry.nwk")
tre.wet <- read.newick(file="Cetacea_wet.nwk")
tre.dry$tip.label[tre.dry$tip.label=="Physeter_catodon"] <- "Physeter_macrocephalus"
tre.wet$tip.label[tre.wet$tip.label=="Physeter_catodon"] <- "Physeter_macrocephalus"

Dry <- na.omit(dat[dat$Bone == "dry",c(2,10,11)])
Wet <- na.omit(dat[dat$Bone == "wet",c(2,10,11)])
Dry.spmeans <- aggregate(Dry, by=list(Dry$Species), FUN=mean)
Wet.spmeans <- aggregate(Wet, by=list(Wet$Species), FUN=mean)
rownames(Dry.spmeans) <- Dry.spmeans$Group.1
rownames(Wet.spmeans) <- Wet.spmeans$Group.1
dry <- Dry.spmeans[tre.dry$tip.label,-c(1,2)]
dry["Delphinus_delphis",] <- buf.means
wet <- Wet.spmeans[tre.wet$tip.label,-c(1,2)]

gls.dry <- gls(boneM~bodyM,data = dry,correlation=corPagel(value=1,phy=tre.dry, 
          fixed=T, form = ~1),method="ML")
gls.wet <- gls(boneM~bodyM,data = wet,correlation=corPagel(value=1,phy=tre.wet, 
          fixed=T, form = ~1),method="ML")
summary(gls.dry)
summary(gls.wet)


### Linear Model

dat$grp <- as.factor(paste0(dat$Species,"_",dat$Bone,"_",dat$Geography))

lm1 <- lm(boneM~bodyM+grp,dat)
summary(lm1)


### Plot

pscale <- NULL; for(i in -3:2) pscale <- c(pscale,seq(10^i,9*10^i,10^i))
log10pscale <- log10(pscale)
cs.base <- palette.colors(9)
ss <- c(3,16,17,13,15,18,5,1,19,4,2,17,6,7,8,9,0,15,16,10,11,12)
cs <- cs.base[c(1,2,3,4,4,5,9,6,6,1,7,7,1,1,1,1,8,8,8,1,1,1)]
x11(w=8,h=8)
#pdf(file="Bone-Body.pdf",h=6.5,w=6.5)
  plot(boneM~bodyM,dat, type="n",axes=F)
    axis(1,at=log10pscale,labels=pscale*1000,las=3)
    axis(2,at=log10pscale,labels=pscale*1000,las=1)
    box()
    for(i in 1:nlevels(dat$grp)){
      dd <- na.omit(dat[dat$grp==levels(dat$grp)[i],10:11])
      if(!is.vector(dd)){
        if(nrow(dd)>1) polygon(dd[chull(dd),],border=cs[i],lwd=2,
                               col=scales::alpha(cs[i],0.2))
      } 
    }
    for(i in 1:nlevels(dat$grp)){
      dd <- na.omit(dat[dat$grp==levels(dat$grp)[i],10:11])
      points(boneM~bodyM,dd,col=cs[i],pch=ss[i])
    }
  legend("topleft",levels(dat$grp),col=cs,pch=ss)
#dev.off()


