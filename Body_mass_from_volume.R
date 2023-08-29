# Set your working directly to the unzipped folder containing this code
# Make sure that subfolders "Perucetus" and "Balaenoptera_musculus" are present,
# containing body silhouette PNG images.

library(paleomass)
library(rgl)

### Blue whale

# Based on the lateral and ventral views
paleomass(Folder="Balaenoptera_musculus", body.axis.l=30,n1=2,n2=2.3,
          ffin.onset=0.29,ffin.adj.up=60,ffin.adj.lat=-20,
          cfin.roll=pi/2,cfin.onset=0.92,cfin.adj.up=-105,
          dfin.onset=0.75,dfin.adj.up=5)

# Based on the lateral view only
paleomass(Folder="Balaenoptera_musculus",BodyV.File="BodyL",body.axis.l=30,n1=2,n2=2.3,
          ffin.onset=0.29,ffin.adj.up=60,ffin.adj.lat=-20,
          cfin.roll=pi/2,cfin.onset=0.92,cfin.adj.up=-105,
          dfin.onset=0.75,dfin.adj.up=5)

### Perucetus

paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=20,n1=2,n2=2.3,
          cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
          ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
          hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=19,n1=2,n2=2.3,
          cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
          ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
          hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55,Save.Total.Mesh=F)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=18,n1=2,n2=2.3,
          cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
          ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
          hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55,Save.Total.Mesh=F)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=17,n1=2,n2=2.3,
          cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
          ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
          hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55,Save.Total.Mesh=F)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=16,n1=2,n2=2.3,
          cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
          ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
          hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55,Save.Total.Mesh=F)
paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=15,n1=2,n2=2.3,
          cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
          ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
          hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55,Save.Total.Mesh=F)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
rgl::close3d(dev = rgl::cur3d(), silent = TRUE)


mean.Manatee <- 74181.55/0.8
mean.Manatee <- 45292.4/0.8
ff <- function(x){
  pm.out <- paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=x,n1=2,n2=2.3,
                      cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
                      ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
                      hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55,Save.Total.Mesh=F)
  rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
  rgl::close3d(dev = rgl::cur3d(), silent = TRUE)
  abs(pm.out[14,6]-mean.Manatee)
}

# optimize(ff,c(15,20))
# The optimum is at 17.99741

paleomass(Folder="Perucetus",BodyV.File="BodyL",body.axis.l=18.0,n1=2,n2=2.3,
          cfin.onset=0.88,cfin.roll=pi/2,cfin.adj.up=55,
          ffin.onset=0.117,ffin.adj.lat=-10,ffin.adj.up=30,ffin.roll=pi/6,
          hfin.onset=0.7,hfin.adj.lat=-200,hfin.adj.up=55,Save.Total.Mesh=F)
