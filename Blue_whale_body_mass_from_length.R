# set your working directly to the unzipped folder containing this code

dat <- DAT <- read.csv("Whale_mass_length.csv")
dat[,sapply(DAT,is.numeric)] <- log10(DAT[,sapply(DAT,is.numeric)])

## Selecting blue whales
Blue <- dat[dat$grp2=="Blue" & dat$Pregnant=="N",]

## Selecting data with a uniform collection method
blue <- Blue[substr(Blue$Source,1,4)=="Ohno" | substr(Blue$Source,1,4)=="Nish",]
#blue <- Blue

## OLS
lm.blue <- lm(bodyM~bodyL,blue)
summary(lm.blue)
10^predict(lm.blue,interval="confidence")
10^predict(lm.blue,interval="prediction")
10^predict(lm.blue,newdata=data.frame(bodyL=log10(c(22.73,23.37))),interval="prediction")

## Body mass estimation for an individual with a body length of 30 m
largest <- 10^predict(lm.blue,data.frame(bodyL=log10(30)),interval="confidence")
largest <- c(largest,10^predict(lm.blue,data.frame(bodyL=log10(30)),interval="prediction")[,-1])
names(largest) <- c("mean","ci-low","ci-high","pi-low","pi-high")

## Assuming blood loss of 7%
largest/0.93

## Assuming blood loss of 14%
largest/0.86
