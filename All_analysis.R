
## @knitr Loading packages, results='hide',message=FALSE

require(raster)
require(rgdal)
require(maptools)
require(vegan)

source("anetwork.R")
source("MidD.R")
source("simunetwork.R")

Custom.Palette<- colorRampPalette(colors=c('lightblue','yellow','red'),bias=1,space='rgb')
color.select<-function(x){Custom.Palette(101)[1+round(((x-min(x))/diff(range(x)))*100)]}

op<-par(no.readonly=TRUE)


## @knitr Importing data

mammal.data=read.csv('NC5.csv')# Data from each locality




## @knitr Importing maps, results="hide"

folder<-"Shapefiles/"

brazil<-readShapeSpatial(paste(folder,"BR_Contorno.shp",sep=""))
biomes<-readShapeSpatial(paste(folder,"bioma.shp",sep=""))



## @knitr Importing environmental layers from worldclim, results="hide"

wclim <- getData("worldclim", var="bio", res=2.5, path="")
wclim.clipped<-crop(wclim,extent(-70,-30,-35,5))
#rm(wclim)



## @knitr Transforming data, results="hide"

### Calculate the distance from the coast

poly<-(brazil@polygons)[[1]]@Polygons

LAT.coast<-unlist(lapply(poly,function(x){x@coords[,2]}))
LONG.coast<-unlist(lapply(poly,function(x){x@coords[,1]}))
long.chui<--(53+22/60+25/3600)

LONG.coast.east<-LONG.coast[LONG.coast>long.chui]
LAT.coast.east<-LAT.coast[LONG.coast>long.chui]

LONG.matrix.east<-t(matrix(LONG.coast.east,nrow=length(LONG.coast.east),ncol=length(mammal.data$Lat)))
LAT.matrix.east<-t(matrix(LAT.coast.east,nrow=length(LAT.coast.east),ncol=length(mammal.data$Lat)))

Longitude.matrix<-matrix(mammal.data$Long,nrow=length(mammal.data$Lat),ncol=length(LONG.coast.east))
Latitude.matrix<-matrix(mammal.data$Lat,nrow=length(mammal.data$Lat),ncol=length(LONG.coast.east))

LONG.diff<-LONG.matrix.east-Longitude.matrix
LAT.diff<-LAT.matrix.east-Latitude.matrix

Distance.coast.matrix<-sqrt(LONG.diff^2+LAT.diff^2)
mammal.data$Distance.coast<-apply(Distance.coast.matrix,1,min)

mammal.data$Long.closest<-apply(Distance.coast.matrix,1,function(x){LONG.coast.east[x==min(x)]})
mammal.data$Lat.closest<-apply(Distance.coast.matrix,1,function(x){LAT.coast.east[x==min(x)]})

rm(Distance.coast.matrix,LAT.diff,LONG.diff,Longitude.matrix,Latitude.matrix,
   LONG.matrix.east,LAT.matrix.east,LONG.coast.east,LAT.coast.east,LAT.coast,LONG.coast,
   long.chui,poly)

#Grouping the latitude and longitude in 2x2 degrees (create new variables in the dataframe)
mammal.data$Long2<-floor(mammal.data$Long/2)*2+1
mammal.data$Lat2<-floor(mammal.data$Lat/2)*2+1

attach(mammal.data)

mammal.short.PA<-tapply(rep(1,length(LOCALIDADE)),list(LOCALIDADE,ESPECIE),sum)
mammal.short.PA[is.na(mammal.short.PA)]<-0

#head(mammal.data)



## @knitr Creating_2x2_degrees_for_the_whole_AF

AF.Long<-unlist(lapply((biomes[6,]@polygons)[[1]]@Polygons,function(x)x@coords[,1]))
AF.Lat<-unlist(lapply((biomes[6,]@polygons)[[1]]@Polygons,function(x)x@coords[,2]))

AF.Lat2<-ceiling(c(AF.Lat,Lat)/2)*2-1
AF.Long2<-ceiling(c(AF.Long,Long)/2)*2-1

#Create a dataframe with the unique locations representing the AF (the point 53 is an isolated island))
environment.AF.2d<-unique(data.frame(AF.Long2=AF.Long2,AF.Lat2=AF.Lat2))



## @knitr Summarizing environmental data for each quadrant

species.2d<-tapply(rep(1,length(Lat2)),list(paste(Long2,Lat2),ESPECIE),sum)
species.2d[is.na(species.2d)]<-0
species.2d[species.2d>0]<-1

environment.2d<-data.frame(unique(cbind(Long2,Lat2)))
environment.2d<-environment.2d[order(factor(paste(environment.2d[,1],environment.2d[,2]))),]

environment.2d$conservation<-tapply(as.numeric(CONSERVACAO),paste(Long2,Lat2),mean)
environment.2d$vegetation<-tapply(as.numeric(FISIONOMIA),paste(Long2,Lat2),mean)

environment.2d$matchin.AF.2d<-match(paste(environment.2d[,1],environment.2d[,2]),paste(environment.AF.2d[,1],environment.AF.2d[,2]))

# r extracting environmental data from worldclim layers, results="hide"}

environment.2d[paste("bio",1:19,sep="")]<-NA

for(i in 1:nrow(environment.2d)){
  temp.clim<-extract(wclim.clipped,extent(environment.2d[i,1]-1,environment.2d[i,1]+1,environment.2d[i,2]-1,environment.2d[i,2]+1))
  environment.2d[i,paste("bio",1:19,sep="")]<-colMeans(temp.clim,na.rm=TRUE)}

environment.2d[,paste("PCA.wclim.",1:3,sep="")]<-prcomp(decostand(environment.2d[,paste("bio",1:19,sep="")],"standardize"))$x[,1:3]


# Creating observed statistics

environment.2d$rich<-rowSums(species.2d)

environment.2d[,paste("pcoa.jac.",1:3,sep="")]<-cmdscale(vegdist(species.2d,"jaccard"),k=3)

#plot(pcoa.jac.1~Lat2,data=environment.2d)



## @knitr plotting maps, fig.width=10,fig.height=10,echo=FALSE


par(xpd=NA)

#par (new=T)

plot(biomes[6,],col='lightblue',border="0",xlim=c(-45,-38))

#plot(subset(wclim.clipped,"bio1"),add=T)
plot(brazil,add=T,lwd=2)
#plot(biomes[6,],add=T,col=adjustcolor(4,alpha=.2),border=0)

rect(environment.AF.2d$AF.Long2-1,environment.AF.2d$AF.Lat2-1,environment.AF.2d$AF.Long2+1,environment.AF.2d$AF.Lat2+1,col=adjustcolor("darkgrey",alpha=.1),border=adjustcolor("darkgrey",alpha=1))

rect(environment.2d$Long2-1,environment.2d$Lat2-1,environment.2d$Long2+1,environment.2d$Lat2+1,col=adjustcolor("green",alpha=.2),border=adjustcolor("darkgrey",alpha=1))
points(mammal.data$Long,mammal.data$Lat,cex=log(ifelse(is.na(AREA..ha.),100,AREA..ha.)+1)/3,pch=21,bg=adjustcolor(2,alpha=.1),col=0)


range.size<-range(log(ifelse(is.na(AREA..ha.),100,AREA..ha.)+1))
divs=5
sizes<-seq(range.size[1]+(diff(range.size)/divs/2),range.size[2]-(diff(range.size)/divs/2),length=divs)

legend(-30,-5, legend=round(exp(sizes)) , pch = 21, col=0,pt.bg=adjustcolor(2,alpha=.5),pt.cex=sizes/3,y.intersp=exp(seq(log(1.2),log(1.5),length=divs)),x.intersp=3,bty="n")

text(-28,-4,"Area (ha)",cex=2)

#points(Long2,Lat2)
#rect(Long2-1,Lat2-1,Long2+1,Lat2+1,col=adjustcolor(3,alpha=.01))


#text(tapply(Long2+1,paste(Lat2,Long2),min),tapply(Lat2+1,paste(Lat2,Long2),min),
#tapply(Lat2,paste(Lat2,Long2),length),font=2)

#text(tapply(Long2+1,paste(Lat2,Long2),min),tapply(Lat2+1,paste(Lat2,Long2),min),
#tapply(AUTOR.E.ANO,paste(Lat2,Long2),function(x)nlevels(factor(as.character(x)))))

#text(environment.2d$Long2,environment.2d$Lat2,rowSums(species.2d))



## @knitr Connectivity, echo=FALSE

#Create a matrix of euclidean distances between all pairs of quadrants
dist.AF.2d<-as.matrix(dist(environment.AF.2d[,1:2]))

#Set as 1 the migration rate when the distance is below 4 and above 2 degrees (all adjacent quadrants, including diagonal but not self loops)

connect.AF.2d<-ifelse(dist.AF.2d<4&dist.AF.2d>0,1,0)
connect.2d<-connect.AF.2d[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d]

#rownames(connect.2d)<-1:26




## @knitr plotconnectivity, echo=FALSE, fig.width=16,fig.height=8

par(mfrow=c(1,2))

plot(brazil,border="darkgrey")
rect(environment.2d$Long2-1,environment.2d$Lat2-1,environment.2d$Long2+1,environment.2d$Lat2+1,border=adjustcolor("darkgrey",alpha=1))

points(environment.2d$Long2,environment.2d$Lat2,pch=21,col=0,bg=2)

nnodes.2d<-nrow(environment.2d)

segments(matrix(environment.2d$Long2,nnodes.2d,nnodes.2d),
         matrix(environment.2d$Lat2,nnodes.2d,nnodes.2d),
         t(matrix(environment.2d$Long2,nnodes.2d,nnodes.2d)),
         t(matrix(environment.2d$Lat2,nnodes.2d,nnodes.2d)),col=connect.2d*2)


par(op)


## @knitr Mid-Domain, cache=FALSE

reps=999

similarity.MidD<-array(NA,dim=c(nrow(environment.2d),nrow(environment.2d),reps))
rich.MidD<-array(NA,dim=c(nrow(environment.2d),reps))

for(i in 1:reps){
  simu1<-MidD(connect.2d,colSums(species.2d),simplify=T)
  similarity.MidD[,,i]<-as.matrix(vegdist(simu1,"jaccard"))
  rich.MidD[,i]<-rowSums(simu1)
  }


#sum(colSums(simu1)-colSums(species.2d))

environment.2d$rich.MidD.mean<-rowMeans(rich.MidD)
environment.2d$rich.MidD.sd<-apply(rich.MidD,1,sd)

similarity.MidD.mean<-apply(similarity.MidD,2,function(x)rowMeans(x,na.rm=TRUE))
similarity.MidD.sd<-apply(similarity.MidD,2,function(x)apply(x,1,sd,na.rm=TRUE))

environment.2d[,paste("pcoa.MidD.jac.",1:3,sep="")]<-(
  prcomp(similarity.MidD.mean)$x)[,1:3]

#plot(pcoa.jac.1~pcoa.MidD.jac.1,data=environment.2d)



## @knitr ,fig.width=18,fig.height=9,echo=FALSE

par(mfrow=c(1,2))

plot(brazil)
rect(environment.2d$Long2-1,environment.2d$Lat2-1,environment.2d$Long2+1,environment.2d$Lat2+1,
     col=adjustcolor(color.select(environment.2d$rich.MidD.mean),alpha=.8))
title("A)",cex=3,adj=0)


plot(brazil)
rect(environment.2d$Long2-1,environment.2d$Lat2-1,environment.2d$Long2+1,environment.2d$Lat2+1,
     col=adjustcolor(color.select(environment.2d$pcoa.MidD.jac.1),alpha=.8))
title("B)",cex=3,adj=0)

par(op)


## @knitr plot.connectivity, echo=FALSE, fig.width=16,fig.height=8

par(mfrow=c(1,2))

plot(brazil,border="darkgrey")
rect(environment.AF.2d$AF.Long2-1,environment.AF.2d$AF.Lat2-1,environment.AF.2d$AF.Long2+1,environment.AF.2d$AF.Lat2+1,border=adjustcolor("darkgrey",alpha=1))

points(environment.AF.2d$AF.Long2,environment.AF.2d$AF.Lat2,pch=21,col=0,bg=2)

nnodes.AF.2d<-nrow(environment.AF.2d)

segments(matrix(environment.AF.2d$AF.Long2,nnodes.AF.2d,nnodes.AF.2d),
         matrix(environment.AF.2d$AF.Lat2,nnodes.AF.2d,nnodes.AF.2d),
         t(matrix(environment.AF.2d$AF.Long2,nnodes.AF.2d,nnodes.AF.2d)),
         t(matrix(environment.AF.2d$AF.Lat2,nnodes.AF.2d,nnodes.AF.2d)),col=connect.AF.2d*2)

par(op)


## @knitr Optimization, cache=TRUE, results='hide',fig.show='hide'

#Set initial parameters

m=0.01790858  #migration rate
v=0.0017789240  #speciation rate
N=100 #Number of individuals in each node

#Create migration matrix (the sum of all migration rates in a node equals 1)
#M<-connect.AF.2d*(m/rowSums(connect.AF.2d))
#M[is.na(M)]<-0

#diag(M)<-1-(rowSums(M)-diag(M))

#rowSums(M)

#neutral.AF<-anetwork(N,M,v)
#F1=neutral.AF$finalF

#environment.AF.2d$Hill<-1/(diag(neutral.AF$finalF))
#morisita.AF.2d<-extractMH(neutral.AF$finalF)


optimized<-optim(c(m,v),lower=0.001,upper=0.999,method="L-BFGS-B",function(x){
  
  M<-connect.AF.2d*(x[1]/rowSums(connect.AF.2d))
  M[is.na(M)]<-0
  
  diag(M)<-1-(rowSums(M)-diag(M))
  
  #rowSums(M)
  
  neutral.AF<-anetwork(N,M,x[2],nruns=1000)
  F1=neutral.AF$finalF
  
  morisita.AF.2d<-extractMH(neutral.AF$finalF)
  
  sum((decostand(prcomp(morisita.AF.2d[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d])$x[,1],'range')-
         decostand(environment.2d$pcoa.jac.1,'range'))^2)
  
  })




## @knitr 
m<-optimized$par[1]
v<-optimized$par[2]

M<-connect.AF.2d*(m/rowSums(connect.AF.2d))
M[is.na(M)]<-0

diag(M)<-1-(rowSums(M)-diag(M))

#rowSums(M)

neutral.AF<-anetwork(N,M,v)
F1=neutral.AF$finalF

environment.AF.2d$Hill<-1/(diag(neutral.AF$finalF))
morisita.AF.2d<-extractMH(neutral.AF$finalF)

environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]<-prcomp((morisita.AF.2d))$x[,1:3]
environment.AF.2d[,paste("pcoa.mor.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA

environment.AF.2d[environment.2d$matchin.AF.2d,paste("pcoa.mor.",1:3,sep="")]<-#PCA.AF[,1:3]
  prcomp(morisita.AF.2d[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d])$x[,1:3]

#plot(environment.2d$Lat2,decostand(environment.2d$pcoa.jac.1,'range'),bg="gold",pch=22)
#points((environment.AF.2d$AF.Lat2),decostand(environment.AF.2d$pcoa.mor.1,'range',na.rm=T),bg=1,pch=21)
#points(environment.2d$Lat2,decostand(environment.2d$pcoa.MidD.jac.1,'range'),bg=2,pch=23)

#plot(decostand(prcomp(morisita.AF.2d[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d])$x[,1],'range'),decostand(environment.2d$pcoa.jac.1,'range'))



## @knitr neutral_theory_real_simulation

N<-matrix(N,nrow(M))

if(file.exists("resu.simu10k.txt")){
  
  resu.simu<- read.table("resu.simu10k.txt",header=TRUE)
  
  }else{
    
    #m=0.01790858  #migration rate
    #v=0.0017789240  #speciation rate
    #N=100 #Number of individuals in each node
    
    reps2<-10000
    initbuffer<-10000
    
    sp.abund2<-simunetwork(N,M,v,nruns=initbuffer,simplify=FALSE)
    
    resu.simu<-matrix(NA,sum(N),reps2)
    
    for(rep in 1:reps2){
      
      sp.abund2<-simunetwork(N,M,v,nruns=300,simplify=FALSE,sp.abund=sp.abund2)
      
      resu.simu[,rep]<-sp.abund2
      
      #print(rep)
      }
    
    write.table(resu.simu,file="resu.simu10k.txt")
    
    #tapply(rep(1,sum(N)),list(rep(1:length(N),N),resu.simu[,1]),sum)
    
    }


## @knitr Calculating_statistics_from_neutral_real_simu

similarity.simu.total<-array(NA,dim=c(nrow(environment.AF.2d),nrow(environment.AF.2d),ncol(resu.simu)))
similarity.simu.jac.total<-array(NA,dim=c(nrow(environment.AF.2d),nrow(environment.AF.2d),ncol(resu.simu)))
similarity.simu<-array(NA,dim=c(nrow(environment.2d),nrow(environment.2d),ncol(resu.simu)))
similarity.simu.jac<-array(NA,dim=c(nrow(environment.2d),nrow(environment.2d),ncol(resu.simu)))
hill.simu<-matrix(NA,nrow(environment.AF.2d),ncol(resu.simu))
rich.simu<-matrix(NA,nrow(environment.AF.2d),ncol(resu.simu))

for(i in 1:ncol(resu.simu)){
  
  ver<-tapply(rep(1,sum(N)),list(rep(1:length(N),N),resu.simu[,i]),sum)
  ver[is.na(ver)]<-0
  
  similarity.simu.total[,,i]<-as.matrix(vegdist(ver,"morisita"))
  similarity.simu.jac.total[,,i]<-as.matrix(vegdist(ver,"jaccard"))
  
  similarity.simu[,,i]<-as.matrix(vegdist(ver,"morisita"))[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d]
  similarity.simu.jac[,,i]<-as.matrix(vegdist(ver,"jaccard"))[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d]
  
  hill.simu[,i]<-(1/(rowSums((ver/rowSums(ver))^2)*(rowSums(ver)-1)/(rowSums(ver))))
  rich.simu[,i]<-rowSums(ver>0)
  }


environment.AF.2d$rich.simu.mean<-rowMeans(rich.simu)
environment.AF.2d$rich.simu.sd<-apply(rich.simu,1,sd)

environment.AF.2d$hill.simu.mean<-rowMeans(hill.simu)
environment.AF.2d$hill.simu.sd<-apply(hill.simu,1,sd)

environment.AF.2d[,paste("pcoa.mor.sim.total.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA

environment.AF.2d[,paste("pcoa.jac.sim.total.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA

environment.AF.2d[,paste("pcoa.mor.sim.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA

environment.AF.2d[,paste("pcoa.jac.sim.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA


environment.AF.2d[,paste("pcoa.mor.sim.total.",1:3,sep="")]<-#PCA.AF[,1:3]
  prcomp(apply(similarity.simu.total[,,round(.6*ncol(resu.simu)):ncol(resu.simu)],2,function(x)rowMeans(x)))$x[,1:3]

environment.AF.2d[,paste("pcoa.jac.sim.total.",1:3,sep="")]<-#PCA.AF[,1:3]
  prcomp(apply(similarity.simu.jac.total[,,round(.6*ncol(resu.simu)):ncol(resu.simu)],2,function(x)rowMeans(x)))$x[,1:3]


environment.AF.2d[environment.2d$matchin.AF.2d,paste("pcoa.mor.sim.",1:3,sep="")]<-#PCA.AF[,1:3]
  prcomp(apply(similarity.simu[,,round(.6*ncol(resu.simu)):ncol(resu.simu)],2,function(x)rowMeans(x)))$x[,1:3]

environment.AF.2d[environment.2d$matchin.AF.2d,paste("pcoa.jac.sim.",1:3,sep="")]<-#PCA.AF[,1:3]
  prcomp(apply(similarity.simu.jac[,,round(.6*ncol(resu.simu)):ncol(resu.simu)],2,function(x)rowMeans(x)))$x[,1:3]

#cor(environment.AF.2d$pcoa.mor.1,environment.AF.2d$pcoa.mor.sim.1,use="c")




## @knitr ,fig.width=18,fig.height=9,echo=FALSE

par(mfrow=c(1,2))

plot(brazil)
rect(environment.AF.2d$AF.Long2-1,environment.AF.2d$AF.Lat2-1,environment.AF.2d$AF.Long2+1,environment.AF.2d$AF.Lat2+1,
     col=adjustcolor(color.select(environment.AF.2d$rich.simu.mean),alpha=.8))
title("A)",cex=3,adj=0)


plot(brazil)
rect(environment.AF.2d$AF.Long2-1,environment.AF.2d$AF.Lat2-1,environment.AF.2d$AF.Long2+1,environment.AF.2d$AF.Lat2+1,
     col=adjustcolor(color.select(environment.AF.2d$pcoa.jac.sim.total.1),alpha=.8))
title("B)",cex=3,adj=0)

par(op)


## @knitr , fig.width=16,fig.height=8

par(mfrow=c(1,2))

plot(environment.2d$Lat2,decostand(environment.2d$pcoa.jac.1,'range'),bg="gold",pch=22)

points(environment.2d$Lat2,decostand(environment.2d$pcoa.MidD.jac.1,'range'),bg=2,pch=23)

points((environment.AF.2d$AF.Lat2),decostand(environment.AF.2d$pcoa.mor.1,'range',na.rm=T),bg=4,pch=24)

points((environment.AF.2d$AF.Lat2),decostand(environment.AF.2d$pcoa.mor.sim.1,'range',na.rm=T),bg=1,pch=3,cex=2)

plot.new()

legend("topleft",c("Observed","MidD","Neutral-ana","Neutral-simu"),pch=c(22,23,24,3),pt.bg=c("gold",2,4,3),bty="n",y.intersp=1.2,cex=2)

#plot(decostand(environment.AF.2d$pcoa.mor.1,'range',na.rm=T),decostand(environment.AF.2d$pcoa.mor.sim.1,'range',na.rm=T),bg=1,pch=3,cex=2)
#plot(environment.AF.2d$pcoa.mor.1,environment.AF.2d$pcoa.mor.sim.1,bg=1,pch=3,cex=2)



## @knitr using_sdm-logistic_regression_of_each_species, warning=FALSE

#Calculates the glm coefficients (slope and intercept) for each species
#glm.results gives the glm results (with p values) for all species
#glm.coef gives just the coefficients for each species and variable
#the other variables are meaningless and are necessary just to create the others

glm.results.2d<-list()
bysp.glm.results.2d<-list()

glm.coef.2d<-list()
bysp.glm.coef.2d<-matrix(NA,ncol(species.2d),2,dimnames=list(colnames(species.2d),c('a','b')))

glm.variables.2d<-c("Long2","Lat2","conservation","vegetation",paste("bio",1:19,sep=""),"PCA.wclim.1","PCA.wclim.2")

for(k in glm.variables.2d){
  for(i in colnames(species.2d)){
    
    bysp.glm.results.2d[[i]]<-glm(as.formula(paste('species.2d[,i]~',paste(k))),family=binomial,data=environment.2d)
    bysp.glm.coef.2d[i,]<-coefficients(bysp.glm.results.2d[[i]])
    
    }
  glm.coef.2d[[k]]<-list(bysp.glm.coef.2d,environment.2d[k])
  }


################################
#Estimate the probability of occurrence for each species in each site given the environmental variables

probs.2d<-lapply(glm.coef.2d,function(x){(exp(x[[1]][,1]+x[[1]][,2]%*%t(x[[2]])))/(1+exp(x[[1]][,1]+x[[1]][,2]%*%t(x[[2]])))}) 

################################

Hill.logis<-data.frame(Hill.logis=lapply(probs.2d,function(x)1/(1-(colSums((x/colSums(x))^2)))))
rich.logis<-data.frame(rich.logis=lapply(probs.2d,function(x)colSums(x)))

environment.AF.2d[,as.character(colnames(Hill.logis))]<-NA
environment.AF.2d[,as.character(colnames(rich.logis))]<-NA

environment.AF.2d[environment.2d$matchin.AF.2d,colnames(Hill.logis)]<-Hill.logis
environment.AF.2d[environment.2d$matchin.AF.2d,colnames(rich.logis)]<-rich.logis



## @knitr species_similarity_in_sdms

species.2d.logis<-species.2d*0

reps=100

similarity.logis.jac<-array(NA,dim=c(nrow(environment.2d),nrow(environment.2d),reps))

dim(similarity.logis.jac)

for(j in glm.variables.2d){
  
  j=glm.variables.2d[[25]]
  
  for(k in 1:reps){
    for(i in 1:nrow(species.2d.logis)){
      species.2d.logis[i,sample(ncol(species.2d.logis),sum(species.2d[i,]),prob=probs.2d[[j]][,i])] <- 1
      }
    similarity.logis.jac[,,k]<-as.matrix(vegdist(species.2d.logis,"jaccard"))
    }
  
  environment.AF.2d[,paste("pcoa.logis.jac.",j,".",1:3,sep="")]<-NA
  
  environment.AF.2d[environment.2d$matchin.AF.2d,paste("pcoa.logis.jac.",j,".",1:3,sep="")]<-
    prcomp(apply(similarity.logis.jac,2,rowMeans))$x[,1:3]
  
  }




## @knitr plot_sdms,fig.align='center',fig.width=16,fig.height=16

par(op)

par(mgp=c(.2,0,0))
par(mar=c(2,2,0,0))

par(mfrow=c(ceiling(length(Hill.logis)/ceiling(sqrt(length(Hill.logis)))),ceiling(sqrt(length(Hill.logis)))))
for (i in glm.variables.2d){plot(Hill.logis[[paste("Hill.logis.",i,sep="")]]~environment.2d[[i]],pch=21,bg=color.select(Hill.logis[[paste("Hill.logis.",i,sep="")]]),xlab=i,ylab="Diversity",axes=FALSE);box()}

par(mfrow=c(ceiling(length(rich.logis)/ceiling(sqrt(length(rich.logis)))),ceiling(sqrt(length(rich.logis)))))
for (i in glm.variables.2d){plot(rich.logis[[paste("rich.logis.",i,sep="")]]~environment.2d[[i]],pch=21,bg=color.select(rich.logis[[paste("rich.logis.",i,sep="")]]),xlab=i,ylab="Diversity",axes=FALSE);box()}




## @knitr 
################################
# Plotting the diversity predictions on the map

attach(environment.2d)

par(mfrow=c(ceiling(length(Hill.logis)/ceiling(sqrt(length(Hill.logis)))),ceiling(sqrt(length(Hill.logis)))),mar=c(0,0,0,0))
for (i in glm.variables.2d){
  plot(brazil)
  title(main=i,line=-3)
  points(Long2,Lat2,pch=22,bg='darkgrey',col=0,cex=.5)
  rect(Long2-1,Lat2-1,Long2+1,Lat2+1,col=color.select(Hill.logis[[paste("Hill.logis.",i,sep="")]]))
  plot(brazil,add=T)
  }


par(mfrow=c(ceiling(length(rich.logis)/ceiling(sqrt(length(rich.logis)))),ceiling(sqrt(length(rich.logis)))),mar=c(0,0,0,0))
for (i in glm.variables.2d){
  plot(brazil)
  title(main=i,line=-3)
  points(Long2,Lat2,pch=22,bg='darkgrey',col=0,cex=.5)
  rect(Long2-1,Lat2-1,Long2+1,Lat2+1,col=color.select(rich.logis[[paste("rich.logis.",i,sep="")]]))
  plot(brazil,add=T)
  }

detach(environment.2d)



## @knitr binding_all_data_frames

environment.AF.2d[,colnames(environment.2d)]<-NA
environment.AF.2d[environment.2d$match,colnames(environment.2d)]<-environment.2d


## @knitr comparisons of the models

# Species richness

# Simple linear regressions among the predictions and the real data

#################
# Raw environmental variables (spline) + logistic regression for each species individually

vars<-c("Lat2","Long2","conservation","vegetation",paste("bio",1:19,sep=""),"PCA.wclim.1","PCA.wclim.2")

par(mgp=c(.2,0,0))
par(mar=c(2,2,0,0))
par(mfrow=c(ceiling(length(vars)/ceiling(sqrt(length(vars)))),ceiling(sqrt(length(vars)))))

for (i in vars){
  plot(as.formula(paste("rich~",i)),data=environment.AF.2d[!is.na(environment.AF.2d$bio19),],axes=FALSE)
  lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),i],environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),paste("rich.logis.",i,sep="")],nknots=5),lwd=2,col=4)
  lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),i],environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),paste("rich",sep="")],nknots=5),lwd=2,col=2)
  #points(environment.AF.2d[,i],environment.AF.2d[,paste("rich.logis.",i,sep="")],pch=21,cex=.2,bg=4,col=0)
  box()
  #abline(lm(as.formula(paste("rich~",i)),data=environment.AF.2d))
  }

#Long was excluded for being strongly correlated to Lat

summary(step(lm(rich~.,data=environment.AF.2d[,c("rich","Lat2","conservation","vegetation","PCA.wclim.1","PCA.wclim.2")])))
summary(step(lm(rich~.,data=environment.AF.2d[,c("rich",paste("rich.logis.",c("Lat2","conservation","vegetation","PCA.wclim.1","PCA.wclim.2"),sep=""))])))

environment.AF.2d[,c("rich",paste("rich.logis.",c("Lat2","conservation","vegetation","PCA.wclim.1","PCA.wclim.2"),sep=""))]

(environment.AF.2d)[,"rich.logis.Lat2"]


###############
# Neutral model

summary(lm(rich~rich.simu.mean,data=environment.AF.2d))

par(op)
plot(rich~rich.simu.mean,data=environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),])
lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),"rich.simu.mean"],environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),"rich"],nknots=5),lwd=2,col=2)


###############
# Mid-Domain

summary(step(lm(rich~rich.MidD.mean,data=environment.2d)))
par(op)
plot(rich~rich.MidD.mean,data=environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),])
lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),"rich.MidD.mean"],environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),"rich"],nknots=5),lwd=2,col=2)


#############################
# Comparing species composition

# Raw environmental variables (spline) + logistic regression for each species individually

vars<-c("Lat2","Long2","conservation","vegetation",paste("bio",1:19,sep=""),"PCA.wclim.1","PCA.wclim.2")

par(mgp=c(.2,0,0))
par(mar=c(2,2,0,0))
par(mfrow=c(ceiling(length(vars)/ceiling(sqrt(length(vars)))),ceiling(sqrt(length(vars)))))

for (i in vars){
  plot(as.formula(paste("pcoa.jac.1~",i)),data=environment.AF.2d[!is.na(environment.AF.2d$bio19),],axes=FALSE,ylim=range(environment.AF.2d[,c(paste("pcoa.logis.jac.",i,".1",sep=""),"pcoa.jac.1")],na.rm=T))
  lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),i],environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),paste("pcoa.logis.jac.",i,".1",sep="")],nknots=5),lwd=2,col=4)
  lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,"pcoa.jac.1"]),i],environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),paste("pcoa.jac.1",sep="")],nknots=5),lwd=2,col=2)
  #points(environment.AF.2d[,i],environment.AF.2d[,paste("rich.logis.",i,sep="")],pch=21,cex=.2,bg=4,col=0)
  box()
  #abline(lm(as.formula(paste("rich~",i)),data=environment.AF.2d))
  }

summary(step(lm(pcoa.jac.1~.,data=environment.AF.2d[,c("pcoa.jac.1","Lat2","conservation","vegetation","PCA.wclim.1","PCA.wclim.2")])))

#summary(step(lm(pcoa.jac.1~.,data=environment.AF.2d[,c("pcoa.jac.1",paste("pcoa.logis.jac.",paste(c("PCA.wclim.1","PCA.wclim.2"),".1",sep=""),sep=""))])))

#############
# Neutral model

summary(lm(pcoa.jac.1~pcoa.jac.sim.1,data=environment.AF.2d))

colnames(environment.AF.2d)

par(op)
plot(rich~rich.simu.mean,data=environment.AF.2d[!is.na(environment.AF.2d[,"rich"]),])


#################
# Mid-Domain







## @knitr detach_data

detach(mammal.data)



## @knitr extracting_R_chunks


purl("All_analysis.Rmd")



## @knitr knit2html, echo=FALSE, eval=FALSE
## 
## library(knitr)
## 
## knit2html("All_analysis.Rmd")
## 


