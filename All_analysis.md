Comparing dispersal, environmental and mid-domain effects on species distribution
============================




## By CSDambros


>html updated in 2014-02-12 09:14:32

### Load required packages, import auxiliary functions and set default graphical parameters



```r

require(raster)
require(rgdal)
require(maptools)
require(vegan)

source("anetwork.R")
source("MidD.R")
source("simunetwork.R")

Custom.Palette <- colorRampPalette(colors = c("lightblue", "yellow", "red"), 
    bias = 1, space = "rgb")
color.select <- function(x) {
    Custom.Palette(101)[1 + round(((x - min(x))/diff(range(x))) * 100)]
}

op <- par(no.readonly = TRUE)
```


### Importing the main data file


```r

mammal.data = read.csv("NC5.csv")  # Data from each locality

```



### Importing maps from shapefiles
>The maps must be in a folder called "Shapefiles/" before running the code


```r

folder <- "Shapefiles/"

brazil <- readShapeSpatial(paste(folder, "BR_Contorno.shp", sep = ""))
biomes <- readShapeSpatial(paste(folder, "bioma.shp", sep = ""))
```


### Importing data from the WorldClim (automatically download data if necessary)


```r

wclim <- getData("worldclim", var = "bio", res = 2.5, path = "")
wclim.clipped <- crop(wclim, extent(-70, -30, -35, 5))
# rm(wclim)
```


### Transform data and extract variables from data sets


```r

### Calculate the distance from the coast

poly <- (brazil@polygons)[[1]]@Polygons

LAT.coast <- unlist(lapply(poly, function(x) {
    x@coords[, 2]
}))
LONG.coast <- unlist(lapply(poly, function(x) {
    x@coords[, 1]
}))
long.chui <- -(53 + 22/60 + 25/3600)

LONG.coast.east <- LONG.coast[LONG.coast > long.chui]
LAT.coast.east <- LAT.coast[LONG.coast > long.chui]

LONG.matrix.east <- t(matrix(LONG.coast.east, nrow = length(LONG.coast.east), 
    ncol = length(mammal.data$Lat)))
LAT.matrix.east <- t(matrix(LAT.coast.east, nrow = length(LAT.coast.east), ncol = length(mammal.data$Lat)))

Longitude.matrix <- matrix(mammal.data$Long, nrow = length(mammal.data$Lat), 
    ncol = length(LONG.coast.east))
Latitude.matrix <- matrix(mammal.data$Lat, nrow = length(mammal.data$Lat), ncol = length(LONG.coast.east))

LONG.diff <- LONG.matrix.east - Longitude.matrix
LAT.diff <- LAT.matrix.east - Latitude.matrix

Distance.coast.matrix <- sqrt(LONG.diff^2 + LAT.diff^2)
mammal.data$Distance.coast <- apply(Distance.coast.matrix, 1, min)

mammal.data$Long.closest <- apply(Distance.coast.matrix, 1, function(x) {
    LONG.coast.east[x == min(x)]
})
mammal.data$Lat.closest <- apply(Distance.coast.matrix, 1, function(x) {
    LAT.coast.east[x == min(x)]
})

rm(Distance.coast.matrix, LAT.diff, LONG.diff, Longitude.matrix, Latitude.matrix, 
    LONG.matrix.east, LAT.matrix.east, LONG.coast.east, LAT.coast.east, LAT.coast, 
    LONG.coast, long.chui, poly)

# Grouping the latitude and longitude in 2x2 degrees (create new variables
# in the dataframe)
mammal.data$Long2 <- floor(mammal.data$Long/2) * 2 + 1
mammal.data$Lat2 <- floor(mammal.data$Lat/2) * 2 + 1

attach(mammal.data)

mammal.short.PA <- tapply(rep(1, length(LOCALIDADE)), list(LOCALIDADE, ESPECIE), 
    sum)
mammal.short.PA[is.na(mammal.short.PA)] <- 0

# head(mammal.data)
```


### Create a data with the points grouped in 2x2 quadrants


```r

AF.Long <- unlist(lapply((biomes[6, ]@polygons)[[1]]@Polygons, function(x) x@coords[, 
    1]))
AF.Lat <- unlist(lapply((biomes[6, ]@polygons)[[1]]@Polygons, function(x) x@coords[, 
    2]))

AF.Lat2 <- ceiling(c(AF.Lat, Lat)/2) * 2 - 1
AF.Long2 <- ceiling(c(AF.Long, Long)/2) * 2 - 1

# Create a dataframe with the unique locations representing the AF (the
# point 53 is an isolated island))
environment.AF.2d <- unique(data.frame(AF.Long2 = AF.Long2, AF.Lat2 = AF.Lat2))
```


### Summarizing environmental data and groupping species in each quadrant


```r

species.2d <- tapply(rep(1, length(Lat2)), list(paste(Long2, Lat2), ESPECIE), 
    sum)
species.2d[is.na(species.2d)] <- 0
species.2d[species.2d > 0] <- 1

environment.2d <- data.frame(unique(cbind(Long2, Lat2)))
environment.2d <- environment.2d[order(factor(paste(environment.2d[, 1], environment.2d[, 
    2]))), ]

environment.2d$conservation <- tapply(as.numeric(CONSERVACAO), paste(Long2, 
    Lat2), mean)
environment.2d$vegetation <- tapply(as.numeric(FISIONOMIA), paste(Long2, Lat2), 
    mean)

environment.2d$matchin.AF.2d <- match(paste(environment.2d[, 1], environment.2d[, 
    2]), paste(environment.AF.2d[, 1], environment.AF.2d[, 2]))

# r extracting environmental data from worldclim layers, results='hide'}

environment.2d[paste("bio", 1:19, sep = "")] <- NA

for (i in 1:nrow(environment.2d)) {
    temp.clim <- extract(wclim.clipped, extent(environment.2d[i, 1] - 1, environment.2d[i, 
        1] + 1, environment.2d[i, 2] - 1, environment.2d[i, 2] + 1))
    environment.2d[i, paste("bio", 1:19, sep = "")] <- colMeans(temp.clim, na.rm = TRUE)
}

environment.2d[, paste("PCA.wclim.", 1:3, sep = "")] <- prcomp(decostand(environment.2d[, 
    paste("bio", 1:19, sep = "")], "standardize"))$x[, 1:3]


# Creating observed statistics

environment.2d$rich <- rowSums(species.2d)

environment.2d[, paste("pcoa.jac.", 1:3, sep = "")] <- cmdscale(vegdist(species.2d, 
    "jaccard"), k = 3)

# plot(pcoa.jac.1~Lat2,data=environment.2d)
```


### Plot a map of the original data and data grouped in quadrants

![plot of chunk plotting maps](figure/plotting maps.png) 


Mid-Domain model
===============

Simulate the range of distribution for each species
-----------------------------

>The model picks one species at a time. Then measures how in how many quadrants the species was present
>Select one quadrant at random and spread the species from there until the number of occupied quadrants equals the original number
>Repeat for all species several times. See the complementary function MidD (in "MidD.R") for details.


* Creating a network (matrix) describing the connectivity between quadrants

>In this model, migration can just occur between quadrants in contact to each other




![plot of chunk plotconnectivity](figure/plotconnectivity.png) 


### Run the Mid-Domain model


```r

reps = 999

similarity.MidD <- array(NA, dim = c(nrow(environment.2d), nrow(environment.2d), 
    reps))
rich.MidD <- array(NA, dim = c(nrow(environment.2d), reps))

for (i in 1:reps) {
    simu1 <- MidD(connect.2d, colSums(species.2d), simplify = T)
    similarity.MidD[, , i] <- as.matrix(vegdist(simu1, "jaccard"))
    rich.MidD[, i] <- rowSums(simu1)
}


# sum(colSums(simu1)-colSums(species.2d))

environment.2d$rich.MidD.mean <- rowMeans(rich.MidD)
environment.2d$rich.MidD.sd <- apply(rich.MidD, 1, sd)

similarity.MidD.mean <- apply(similarity.MidD, 2, function(x) rowMeans(x, na.rm = TRUE))
similarity.MidD.sd <- apply(similarity.MidD, 2, function(x) apply(x, 1, sd, 
    na.rm = TRUE))

environment.2d[, paste("pcoa.MidD.jac.", 1:3, sep = "")] <- (prcomp(similarity.MidD.mean)$x)[, 
    1:3]

# plot(pcoa.jac.1~pcoa.MidD.jac.1,data=environment.2d)
```


### Plot species diversity (Richness; A) and species composition (B) predicted by the Mid-Domain model 

>Areas with similar colors in B represent areas with similar composition of species


![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 



Neutral model
===============
>based in Economo and Keitt (2008)

>This model spreads a given number of individuals (*N*) in each quadrant. Initially all individuals in all quadrants are from a single species. Then the model is run for several generations. In each generation all individuals reproduce and die. Based on the parameters *v* and *m*, the indviduals in the offspring can speciate (point speciation) and migrate. In a single quadrant, two individuals can be from the same species if none of them have speciated with rate *v* in the previous generation, AND if they had the same parent or different parents from the same species, either from the same quadrant or not (none, one or both have migrated). Eventually the diversity within and between quadrants reach a steady-state. At this point the model stops. See the complementary function anetwork (in "anetwork.R") for details. This function actually does not simulate the whole process described here, but instead uses a method borrowed from the coalescent theory of population genetics. See Neutral model 2 for a real simulation which gives the same results.


![plot of chunk plot.connectivity](figure/plot.connectivity.png) 


This code optimizes the migration (*m*) and speciation (*v*) (single) parameters in order to better fit the species composition change in the observed data.

>In other words: the model runs hundreds of times for thousands of generations until the best combination of *m* and *v* is reached.
>This procedure is analogue to fit a regression line in a linear regression (when using Maximum Likelihood).


```r

# Set initial parameters

m = 0.01790858  #migration rate
v = 0.001778924  #speciation rate
N = 100  #Number of individuals in each node

# Create migration matrix (the sum of all migration rates in a node equals
# 1) M<-connect.AF.2d*(m/rowSums(connect.AF.2d)) M[is.na(M)]<-0

# diag(M)<-1-(rowSums(M)-diag(M))

# rowSums(M)

# neutral.AF<-anetwork(N,M,v) F1=neutral.AF$finalF

# environment.AF.2d$Hill<-1/(diag(neutral.AF$finalF))
# morisita.AF.2d<-extractMH(neutral.AF$finalF)


optimized <- optim(c(m, v), lower = 0.001, upper = 0.999, method = "L-BFGS-B", 
    function(x) {
        
        M <- connect.AF.2d * (x[1]/rowSums(connect.AF.2d))
        M[is.na(M)] <- 0
        
        diag(M) <- 1 - (rowSums(M) - diag(M))
        
        # rowSums(M)
        
        neutral.AF <- anetwork(N, M, x[2], nruns = 1000)
        F1 = neutral.AF$finalF
        
        morisita.AF.2d <- extractMH(neutral.AF$finalF)
        
        sum((decostand(prcomp(morisita.AF.2d[environment.2d$matchin.AF.2d, environment.2d$matchin.AF.2d])$x[, 
            1], "range") - decostand(environment.2d$pcoa.jac.1, "range"))^2)
        
    })

```


Uses the optimized parameters to recreate the neutral simulation (run one more time)


```r
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
```



Neutral model 2
===============
>based simulation of individuals (takes more than 3 hours to run in a regular computer)


```r

N <- matrix(N, nrow(M))

if (file.exists("resu.simu10k.txt")) {
    
    resu.simu <- read.table("resu.simu10k.txt", header = TRUE)
    
} else {
    
    # m=0.01790858 #migration rate v=0.0017789240 #speciation rate N=100
    # #Number of individuals in each node
    
    reps2 <- 10000
    initbuffer <- 10000
    
    sp.abund2 <- simunetwork(N, M, v, nruns = initbuffer, simplify = FALSE)
    
    resu.simu <- matrix(NA, sum(N), reps2)
    
    for (rep in 1:reps2) {
        
        sp.abund2 <- simunetwork(N, M, v, nruns = 300, simplify = FALSE, sp.abund = sp.abund2)
        
        resu.simu[, rep] <- sp.abund2
        
        # print(rep)
    }
    
    write.table(resu.simu, file = "resu.simu10k.txt")
    
    # tapply(rep(1,sum(N)),list(rep(1:length(N),N),resu.simu[,1]),sum)
    
}
```




```r

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

```


### Plot species diversity (Richness; A) and species composition (B) predicted by the Neutral Model (2nd) 

>Areas with similar colors in B represent areas with similar composition of species


![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


### Comparison of the species composition observed and predicted by the Mid Domain and Neutral (x2) models along the Latitudinal gradient


```r

par(mfrow = c(1, 2))

plot(environment.2d$Lat2, decostand(environment.2d$pcoa.jac.1, "range"), bg = "gold", 
    pch = 22)

points(environment.2d$Lat2, decostand(environment.2d$pcoa.MidD.jac.1, "range"), 
    bg = 2, pch = 23)

points((environment.AF.2d$AF.Lat2), decostand(environment.AF.2d$pcoa.mor.1, 
    "range", na.rm = T), bg = 4, pch = 24)

points((environment.AF.2d$AF.Lat2), decostand(environment.AF.2d$pcoa.mor.sim.1, 
    "range", na.rm = T), bg = 1, pch = 3, cex = 2)

plot.new()

legend("topleft", c("Observed", "MidD", "Neutral-ana", "Neutral-simu"), pch = c(22, 
    23, 24, 3), pt.bg = c("gold", 2, 4, 3), bty = "n", y.intersp = 1.2, cex = 2)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```r

#
# plot(decostand(environment.AF.2d$pcoa.mor.1,'range',na.rm=T),decostand(environment.AF.2d$pcoa.mor.sim.1,'range',na.rm=T),bg=1,pch=3,cex=2)
# plot(environment.AF.2d$pcoa.mor.1,environment.AF.2d$pcoa.mor.sim.1,bg=1,pch=3,cex=2)
```


### Using GLMs to estimate the probability of each species occurring along each environmental gradient

>Based on these coefficients, calculates 
>> The expected number of species in each quadrant
>> The probability of two individuals being from different species (PIE)
>> The morisita-horn index of species similarity


```r

# Calculates the glm coefficients (slope and intercept) for each species
# glm.results gives the glm results (with p values) for all species
# glm.coef gives just the coefficients for each species and variable the
# other variables are meaningless and are necessary just to create the
# others

glm.results.2d <- list()
bysp.glm.results.2d <- list()

glm.coef.2d <- list()
bysp.glm.coef.2d <- matrix(NA, ncol(species.2d), 2, dimnames = list(colnames(species.2d), 
    c("a", "b")))

glm.variables.2d <- c("Long2", "Lat2", "conservation", "vegetation", paste("bio", 
    1:19, sep = ""), "PCA.wclim.1", "PCA.wclim.2")

for (k in glm.variables.2d) {
    for (i in colnames(species.2d)) {
        
        bysp.glm.results.2d[[i]] <- glm(as.formula(paste("species.2d[,i]~", 
            paste(k))), family = binomial, data = environment.2d)
        bysp.glm.coef.2d[i, ] <- coefficients(bysp.glm.results.2d[[i]])
        
    }
    glm.coef.2d[[k]] <- list(bysp.glm.coef.2d, environment.2d[k])
}


################################ Estimate the probability of occurrence
################################ for each species in each site given the
################################ environmental variables

probs.2d <- lapply(glm.coef.2d, function(x) {
    (exp(x[[1]][, 1] + x[[1]][, 2] %*% t(x[[2]])))/(1 + exp(x[[1]][, 1] + x[[1]][, 
        2] %*% t(x[[2]])))
})

################################

Hill.logis <- data.frame(Hill.logis = lapply(probs.2d, function(x) 1/(1 - (colSums((x/colSums(x))^2)))))
rich.logis <- data.frame(rich.logis = lapply(probs.2d, function(x) colSums(x)))

environment.AF.2d[, as.character(colnames(Hill.logis))] <- NA
environment.AF.2d[, as.character(colnames(rich.logis))] <- NA

environment.AF.2d[environment.2d$matchin.AF.2d, colnames(Hill.logis)] <- Hill.logis
environment.AF.2d[environment.2d$matchin.AF.2d, colnames(rich.logis)] <- rich.logis
```


#### Simulate the distribution of individuals using the probabilities and then extract the jaccard index


```r

species.2d.logis <- species.2d * 0

reps = 100

similarity.logis.jac <- array(NA, dim = c(nrow(environment.2d), nrow(environment.2d), 
    reps))

dim(similarity.logis.jac)
```

```
## [1]  26  26 100
```

```r

for (j in glm.variables.2d) {
    
    j = glm.variables.2d[[25]]
    
    for (k in 1:reps) {
        for (i in 1:nrow(species.2d.logis)) {
            species.2d.logis[i, sample(ncol(species.2d.logis), sum(species.2d[i, 
                ]), prob = probs.2d[[j]][, i])] <- 1
        }
        similarity.logis.jac[, , k] <- as.matrix(vegdist(species.2d.logis, "jaccard"))
    }
    
    environment.AF.2d[, paste("pcoa.logis.jac.", j, ".", 1:3, sep = "")] <- NA
    
    environment.AF.2d[environment.2d$matchin.AF.2d, paste("pcoa.logis.jac.", 
        j, ".", 1:3, sep = "")] <- prcomp(apply(similarity.logis.jac, 2, rowMeans))$x[, 
        1:3]
    
}

```


### Species diversity (Hill numbers and Richness) along environmental gradients predicted by the SDMs


```r

par(op)

par(mgp = c(0.2, 0, 0))
par(mar = c(2, 2, 0, 0))

par(mfrow = c(ceiling(length(Hill.logis)/ceiling(sqrt(length(Hill.logis)))), 
    ceiling(sqrt(length(Hill.logis)))))
for (i in glm.variables.2d) {
    plot(Hill.logis[[paste("Hill.logis.", i, sep = "")]] ~ environment.2d[[i]], 
        pch = 21, bg = color.select(Hill.logis[[paste("Hill.logis.", i, sep = "")]]), 
        xlab = i, ylab = "Diversity", axes = FALSE)
    box()
}
```

![plot of chunk plot_sdms](figure/plot_sdms1.png) 

```r

par(mfrow = c(ceiling(length(rich.logis)/ceiling(sqrt(length(rich.logis)))), 
    ceiling(sqrt(length(rich.logis)))))
for (i in glm.variables.2d) {
    plot(rich.logis[[paste("rich.logis.", i, sep = "")]] ~ environment.2d[[i]], 
        pch = 21, bg = color.select(rich.logis[[paste("rich.logis.", i, sep = "")]]), 
        xlab = i, ylab = "Diversity", axes = FALSE)
    box()
}
```

![plot of chunk plot_sdms](figure/plot_sdms2.png) 

```r

```


### Species diversity (Hill numbers and Richness) plotted in the map as predicted by the SDMs of each environmental gradient


```r
################################ Plotting the diversity predictions on the
################################ map

attach(environment.2d)
```

```
## The following objects are masked from mammal.data:
## 
##     Lat2, Long2
```

```r

par(mfrow = c(ceiling(length(Hill.logis)/ceiling(sqrt(length(Hill.logis)))), 
    ceiling(sqrt(length(Hill.logis)))), mar = c(0, 0, 0, 0))
for (i in glm.variables.2d) {
    plot(brazil)
    title(main = i, line = -3)
    points(Long2, Lat2, pch = 22, bg = "darkgrey", col = 0, cex = 0.5)
    rect(Long2 - 1, Lat2 - 1, Long2 + 1, Lat2 + 1, col = color.select(Hill.logis[[paste("Hill.logis.", 
        i, sep = "")]]))
    plot(brazil, add = T)
}
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) 

```r


par(mfrow = c(ceiling(length(rich.logis)/ceiling(sqrt(length(rich.logis)))), 
    ceiling(sqrt(length(rich.logis)))), mar = c(0, 0, 0, 0))
for (i in glm.variables.2d) {
    plot(brazil)
    title(main = i, line = -3)
    points(Long2, Lat2, pch = 22, bg = "darkgrey", col = 0, cex = 0.5)
    rect(Long2 - 1, Lat2 - 1, Long2 + 1, Lat2 + 1, col = color.select(rich.logis[[paste("rich.logis.", 
        i, sep = "")]]))
    plot(brazil, add = T)
}
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) 

```r

detach(environment.2d)
```



```r

environment.AF.2d[, colnames(environment.2d)] <- NA
environment.AF.2d[environment.2d$match, colnames(environment.2d)] <- environment.2d
```



```r

# Species richness

# Simple linear regressions among the predictions and the real data

################# Raw environmental variables (spline) + logistic
################# regression for each species individually

vars <- c("Lat2", "Long2", "conservation", "vegetation", paste("bio", 1:19, 
    sep = ""), "PCA.wclim.1", "PCA.wclim.2")

par(mgp = c(0.2, 0, 0))
par(mar = c(2, 2, 0, 0))
par(mfrow = c(ceiling(length(vars)/ceiling(sqrt(length(vars)))), ceiling(sqrt(length(vars)))))

for (i in vars) {
    plot(as.formula(paste("rich~", i)), data = environment.AF.2d[!is.na(environment.AF.2d$bio19), 
        ], axes = FALSE)
    lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), 
        i], environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), paste("rich.logis.", 
        i, sep = "")], nknots = 5), lwd = 2, col = 4)
    lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), 
        i], environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), paste("rich", 
        sep = "")], nknots = 5), lwd = 2, col = 2)
    #
    # points(environment.AF.2d[,i],environment.AF.2d[,paste('rich.logis.',i,sep='')],pch=21,cex=.2,bg=4,col=0)
    box()
    # abline(lm(as.formula(paste('rich~',i)),data=environment.AF.2d))
}
```

![plot of chunk comparisons of the models](figure/comparisons of the models1.png) 

```r

# Long was excluded for being strongly correlated to Lat

summary(step(lm(rich ~ ., data = environment.AF.2d[, c("rich", "Lat2", "conservation", 
    "vegetation", "PCA.wclim.1", "PCA.wclim.2")])))
```

```
## Start:  AIC=98.78
## rich ~ Lat2 + conservation + vegetation + PCA.wclim.1 + PCA.wclim.2
## 
##                Df Sum of Sq RSS   AIC
## - Lat2          1       8.2 740  97.1
## - PCA.wclim.1   1      22.6 755  97.6
## <none>                      732  98.8
## - vegetation    1      83.9 816  99.6
## - PCA.wclim.2   1     131.2 863 101.1
## - conservation  1     155.2 887 101.8
## 
## Step:  AIC=97.07
## rich ~ conservation + vegetation + PCA.wclim.1 + PCA.wclim.2
## 
##                Df Sum of Sq RSS  AIC
## - PCA.wclim.1   1      30.7 771 96.1
## <none>                      740 97.1
## - vegetation    1      80.9 821 97.8
## - PCA.wclim.2   1     125.0 865 99.1
## - conservation  1     150.7 891 99.9
## 
## Step:  AIC=96.12
## rich ~ conservation + vegetation + PCA.wclim.2
## 
##                Df Sum of Sq RSS   AIC
## <none>                      771  96.1
## - vegetation    1      63.1 834  96.2
## - PCA.wclim.2   1     136.8 908  98.4
## - conservation  1     202.6 973 100.2
```

```
## 
## Call:
## lm(formula = rich ~ conservation + vegetation + PCA.wclim.2, 
##     data = environment.AF.2d[, c("rich", "Lat2", "conservation", 
##         "vegetation", "PCA.wclim.1", "PCA.wclim.2")])
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -8.586 -3.841 -0.921  3.182 11.931 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)   
## (Intercept)    19.957      6.029    3.31   0.0032 **
## conservation   -2.541      1.057   -2.40   0.0251 * 
## vegetation      2.486      1.852    1.34   0.1933   
## PCA.wclim.2    -1.024      0.518   -1.98   0.0608 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.92 on 22 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.413,	Adjusted R-squared:  0.333 
## F-statistic: 5.17 on 3 and 22 DF,  p-value: 0.00744
```

```r
summary(step(lm(rich ~ ., data = environment.AF.2d[, c("rich", paste("rich.logis.", 
    c("Lat2", "conservation", "vegetation", "PCA.wclim.1", "PCA.wclim.2"), sep = ""))])))
```

```
## Start:  AIC=104.3
## rich ~ rich.logis.Lat2 + rich.logis.conservation + rich.logis.vegetation + 
##     rich.logis.PCA.wclim.1 + rich.logis.PCA.wclim.2
## 
##                           Df Sum of Sq  RSS AIC
## - rich.logis.PCA.wclim.1   1       9.6  917 103
## - rich.logis.vegetation    1      22.9  930 103
## - rich.logis.PCA.wclim.2   1      26.5  933 103
## - rich.logis.Lat2          1      29.4  936 103
## <none>                                  907 104
## - rich.logis.conservation  1     156.0 1063 106
## 
## Step:  AIC=102.6
## rich ~ rich.logis.Lat2 + rich.logis.conservation + rich.logis.vegetation + 
##     rich.logis.PCA.wclim.2
## 
##                           Df Sum of Sq  RSS AIC
## - rich.logis.PCA.wclim.2   1      18.0  935 101
## - rich.logis.Lat2          1      25.4  942 101
## - rich.logis.vegetation    1      52.9  969 102
## <none>                                  917 103
## - rich.logis.conservation  1     150.9 1067 105
## 
## Step:  AIC=101.1
## rich ~ rich.logis.Lat2 + rich.logis.conservation + rich.logis.vegetation
## 
##                           Df Sum of Sq  RSS   AIC
## - rich.logis.Lat2          1      26.8  961  99.9
## - rich.logis.vegetation    1      57.0  992 100.7
## <none>                                  935 101.1
## - rich.logis.conservation  1     133.3 1068 102.6
## 
## Step:  AIC=99.87
## rich ~ rich.logis.conservation + rich.logis.vegetation
## 
##                           Df Sum of Sq  RSS   AIC
## - rich.logis.vegetation    1      76.2 1038  99.9
## <none>                                  961  99.9
## - rich.logis.conservation  1     117.0 1078 100.9
## 
## Step:  AIC=99.85
## rich ~ rich.logis.conservation
## 
##                           Df Sum of Sq  RSS   AIC
## <none>                                 1038  99.9
## - rich.logis.conservation  1       276 1314 104.0
```

```
## 
## Call:
## lm(formula = rich ~ rich.logis.conservation, data = environment.AF.2d[, 
##     c("rich", paste("rich.logis.", c("Lat2", "conservation", 
##         "vegetation", "PCA.wclim.1", "PCA.wclim.2"), sep = ""))])
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -9.34  -5.08  -1.52   3.88  11.49 
## 
## Coefficients:
##                         Estimate Std. Error t value Pr(>|t|)  
## (Intercept)                1.278      5.225    0.24    0.809  
## rich.logis.conservation    0.909      0.360    2.53    0.018 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 6.58 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.21,	Adjusted R-squared:  0.177 
## F-statistic: 6.39 on 1 and 24 DF,  p-value: 0.0185
```

```r

environment.AF.2d[, c("rich", paste("rich.logis.", c("Lat2", "conservation", 
    "vegetation", "PCA.wclim.1", "PCA.wclim.2"), sep = ""))]
```

```
##       rich rich.logis.Lat2 rich.logis.conservation rich.logis.vegetation
## 1        5           12.95                  10.821                 11.99
## 11      20           12.95                  18.491                 28.01
## 88       6           12.95                  13.963                 11.99
## 136     NA              NA                      NA                    NA
## 342     20           12.39                  10.821                 11.99
## 894      8           12.39                  12.745                 11.99
## 2286     6           12.45                  10.821                 11.99
## 2852    12           12.45                  22.064                 11.99
## 3639    NA              NA                      NA                    NA
## 4620    NA              NA                      NA                    NA
## 4654    NA              NA                      NA                    NA
## 4657    NA              NA                      NA                    NA
## 4703    NA              NA                      NA                    NA
## 4791    NA              NA                      NA                    NA
## 4902    NA              NA                      NA                    NA
## 4997    NA              NA                      NA                    NA
## 5109    NA              NA                      NA                    NA
## 5500    NA              NA                      NA                    NA
## 5706    NA              NA                      NA                    NA
## 5806    NA              NA                      NA                    NA
## 5929    21           12.45                  18.942                 18.01
## 5954    NA              NA                      NA                    NA
## 6076    22           12.88                  11.208                 16.10
## 6239     8           13.43                  10.821                 11.99
## 6461    15           13.43                  10.821                 11.99
## 6594    NA              NA                      NA                    NA
## 6647    20           14.05                  13.475                 11.99
## 6714    NA              NA                      NA                    NA
## 6815    NA              NA                      NA                    NA
## 6857    NA              NA                      NA                    NA
## 6859    NA              NA                      NA                    NA
## 7379    NA              NA                      NA                    NA
## 7485    12           17.68                  19.361                 15.75
## 7600    NA              NA                      NA                    NA
## 7793     8           18.60                  10.821                 11.99
## 7949    NA              NA                      NA                    NA
## 8199    15           19.46                  15.383                 11.99
## 8276    NA              NA                      NA                    NA
## 8466    NA              NA                      NA                    NA
## 11791   NA              NA                      NA                    NA
## 13076   19           16.72                  15.249                 15.75
## 13801   16           15.82                  19.361                 15.75
## 14805   NA              NA                      NA                    NA
## 15180   19           14.05                  15.480                 13.51
## 15410   28           13.43                  18.880                 14.22
## 17270   NA              NA                      NA                    NA
## 18664   20           12.88                  10.500                 13.86
## 21898   11           12.88                  11.449                 15.75
## 26283    9           12.45                  13.963                 15.75
## 27490   29           12.39                  17.851                 15.75
## 30837   NA              NA                      NA                    NA
## 33569   NA              NA                      NA                    NA
## 34294   NA              NA                      NA                    NA
## 37269    6           17.05                  12.745                 11.99
## 37287    3           14.39                   9.981                 11.99
## 37533    8           13.43                   9.981                 11.99
##       rich.logis.PCA.wclim.1 rich.logis.PCA.wclim.2
## 1                      12.74                  12.58
## 11                     17.98                  12.60
## 88                     12.60                  12.69
## 136                       NA                     NA
## 342                    12.21                  12.63
## 894                    15.62                  12.55
## 2286                   12.27                  12.87
## 2852                   11.70                  12.55
## 3639                      NA                     NA
## 4620                      NA                     NA
## 4654                      NA                     NA
## 4657                      NA                     NA
## 4703                      NA                     NA
## 4791                      NA                     NA
## 4902                      NA                     NA
## 4997                      NA                     NA
## 5109                      NA                     NA
## 5500                      NA                     NA
## 5706                      NA                     NA
## 5806                      NA                     NA
## 5929                   11.96                  13.45
## 5954                      NA                     NA
## 6076                   11.75                  15.74
## 6239                   12.59                  17.00
## 6461                   12.35                  18.85
## 6594                      NA                     NA
## 6647                   13.44                  16.39
## 6714                      NA                     NA
## 6815                      NA                     NA
## 6857                      NA                     NA
## 6859                      NA                     NA
## 7379                      NA                     NA
## 7485                   16.58                  13.61
## 7600                      NA                     NA
## 7793                   19.52                  13.24
## 7949                      NA                     NA
## 8199                   18.55                  12.55
## 8276                      NA                     NA
## 8466                      NA                     NA
## 11791                     NA                     NA
## 13076                  16.57                  14.00
## 13801                  17.71                  14.03
## 14805                     NA                     NA
## 15180                  17.02                  13.07
## 15410                  14.15                  13.19
## 17270                     NA                     NA
## 18664                  13.17                  13.75
## 21898                  11.69                  17.19
## 26283                  11.99                  14.01
## 27490                  12.48                  12.85
## 30837                     NA                     NA
## 33569                     NA                     NA
## 34294                     NA                     NA
## 37269                  14.27                  14.07
## 37287                  12.29                  14.67
## 37533                  12.82                  15.87
```

```r

(environment.AF.2d)[, "rich.logis.Lat2"]
```

```
##  [1] 12.95 12.95 12.95    NA 12.39 12.39 12.45 12.45    NA    NA    NA
## [12]    NA    NA    NA    NA    NA    NA    NA    NA    NA 12.45    NA
## [23] 12.88 13.43 13.43    NA 14.05    NA    NA    NA    NA    NA 17.68
## [34]    NA 18.60    NA 19.46    NA    NA    NA 16.72 15.82    NA 14.05
## [45] 13.43    NA 12.88 12.88 12.45 12.39    NA    NA    NA 17.05 14.39
## [56] 13.43
```

```r


############### Neutral model

summary(lm(rich ~ rich.simu.mean, data = environment.AF.2d))
```

```
## 
## Call:
## lm(formula = rich ~ rich.simu.mean, data = environment.AF.2d)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -9.84  -6.79  -1.04   5.10  13.80 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)
## (Intercept)      -13.30      18.15   -0.73     0.47
## rich.simu.mean     3.99       2.64    1.51     0.14
## 
## Residual standard error: 7.07 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.087,	Adjusted R-squared:  0.049 
## F-statistic: 2.29 on 1 and 24 DF,  p-value: 0.143
```

```r

par(op)
plot(rich ~ rich.simu.mean, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ])
lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), "rich.simu.mean"], 
    environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), "rich"], nknots = 5), 
    lwd = 2, col = 2)
```

![plot of chunk comparisons of the models](figure/comparisons of the models2.png) 

```r


############### Mid-Domain

summary(step(lm(rich ~ rich.MidD.mean, data = environment.2d)))
```

```
## Start:  AIC=104
## rich ~ rich.MidD.mean
## 
##                  Df Sum of Sq  RSS AIC
## - rich.MidD.mean  1      95.3 1314 104
## <none>                        1219 104
## 
## Step:  AIC=104
## rich ~ 1
```

```
## 
## Call:
## lm(formula = rich ~ 1, data = environment.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -11.077  -6.077  -0.577   5.923  14.923 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    14.08       1.42     9.9  3.9e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.25 on 25 degrees of freedom
```

```r
par(op)
plot(rich ~ rich.MidD.mean, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ])
lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), "rich.MidD.mean"], 
    environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), "rich"], nknots = 5), 
    lwd = 2, col = 2)
```

![plot of chunk comparisons of the models](figure/comparisons of the models3.png) 

```r


############################# Comparing species composition

# Raw environmental variables (spline) + logistic regression for each
# species individually

vars <- c("Lat2", "Long2", "conservation", "vegetation", paste("bio", 1:19, 
    sep = ""), "PCA.wclim.1", "PCA.wclim.2")

par(mgp = c(0.2, 0, 0))
par(mar = c(2, 2, 0, 0))
par(mfrow = c(ceiling(length(vars)/ceiling(sqrt(length(vars)))), ceiling(sqrt(length(vars)))))

for (i in vars) {
    plot(as.formula(paste("pcoa.jac.1~", i)), data = environment.AF.2d[!is.na(environment.AF.2d$bio19), 
        ], axes = FALSE, ylim = range(environment.AF.2d[, c(paste("pcoa.logis.jac.", 
        i, ".1", sep = ""), "pcoa.jac.1")], na.rm = T))
    lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), 
        i], environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), paste("pcoa.logis.jac.", 
        i, ".1", sep = "")], nknots = 5), lwd = 2, col = 4)
    lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[, "pcoa.jac.1"]), 
        i], environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), paste("pcoa.jac.1", 
        sep = "")], nknots = 5), lwd = 2, col = 2)
    #
    # points(environment.AF.2d[,i],environment.AF.2d[,paste('rich.logis.',i,sep='')],pch=21,cex=.2,bg=4,col=0)
    box()
    # abline(lm(as.formula(paste('rich~',i)),data=environment.AF.2d))
}
```

```
## Error: undefined columns selected
```

```r

summary(step(lm(pcoa.jac.1 ~ ., data = environment.AF.2d[, c("pcoa.jac.1", "Lat2", 
    "conservation", "vegetation", "PCA.wclim.1", "PCA.wclim.2")])))
```

```
## Start:  AIC=-100.2
## pcoa.jac.1 ~ Lat2 + conservation + vegetation + PCA.wclim.1 + 
##     PCA.wclim.2
## 
##                Df Sum of Sq   RSS    AIC
## - Lat2          1    0.0010 0.349 -102.1
## - conservation  1    0.0077 0.355 -101.6
## <none>                      0.348 -100.2
## - vegetation    1    0.0969 0.445  -95.8
## - PCA.wclim.2   1    0.1304 0.478  -93.9
## - PCA.wclim.1   1    0.1316 0.479  -93.8
## 
## Step:  AIC=-102.1
## pcoa.jac.1 ~ conservation + vegetation + PCA.wclim.1 + PCA.wclim.2
## 
##                Df Sum of Sq   RSS    AIC
## - conservation  1     0.008 0.357 -103.5
## <none>                      0.349 -102.1
## - vegetation    1     0.099 0.447  -97.6
## - PCA.wclim.2   1     0.136 0.484  -95.6
## - PCA.wclim.1   1     1.047 1.396  -68.0
## 
## Step:  AIC=-103.5
## pcoa.jac.1 ~ vegetation + PCA.wclim.1 + PCA.wclim.2
## 
##               Df Sum of Sq   RSS    AIC
## <none>                     0.357 -103.5
## - PCA.wclim.2  1     0.129 0.486  -97.5
## - vegetation   1     0.181 0.538  -94.8
## - PCA.wclim.1  1     1.178 1.534  -67.6
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ vegetation + PCA.wclim.1 + PCA.wclim.2, 
##     data = environment.AF.2d[, c("pcoa.jac.1", "Lat2", "conservation", 
##         "vegetation", "PCA.wclim.1", "PCA.wclim.2")])
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2707 -0.0796 -0.0015  0.0652  0.3421 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.17535    0.05812    3.02   0.0063 ** 
## vegetation  -0.11620    0.03477   -3.34   0.0030 ** 
## PCA.wclim.1  0.07368    0.00865    8.52    2e-08 ***
## PCA.wclim.2  0.02962    0.01050    2.82   0.0100 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.127 on 22 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.795,	Adjusted R-squared:  0.767 
## F-statistic: 28.5 on 3 and 22 DF,  p-value: 9.22e-08
```

```r

#
# summary(step(lm(pcoa.jac.1~.,data=environment.AF.2d[,c('pcoa.jac.1',paste('pcoa.logis.jac.',paste(c('PCA.wclim.1','PCA.wclim.2'),'.1',sep=''),sep=''))])))

############# Neutral model

summary(lm(pcoa.jac.1 ~ pcoa.jac.sim.1, data = environment.AF.2d))
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ pcoa.jac.sim.1, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.1826 -0.1044  0.0059  0.0728  0.3556 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    -2.28e-17   2.63e-02    0.00        1    
## pcoa.jac.sim.1 -5.85e-01   6.85e-02   -8.55  9.6e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.134 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.753,	Adjusted R-squared:  0.742 
## F-statistic:   73 on 1 and 24 DF,  p-value: 9.62e-09
```

```r

colnames(environment.AF.2d)
```

```
##   [1] "AF.Long2"                     "AF.Lat2"                     
##   [3] "Hill"                         "pcoa.mor.AF1"                
##   [5] "pcoa.mor.AF2"                 "pcoa.mor.AF3"                
##   [7] "pcoa.mor.1"                   "pcoa.mor.2"                  
##   [9] "pcoa.mor.3"                   "rich.simu.mean"              
##  [11] "rich.simu.sd"                 "hill.simu.mean"              
##  [13] "hill.simu.sd"                 "pcoa.mor.sim.total.1"        
##  [15] "pcoa.mor.sim.total.2"         "pcoa.mor.sim.total.3"        
##  [17] "pcoa.jac.sim.total.1"         "pcoa.jac.sim.total.2"        
##  [19] "pcoa.jac.sim.total.3"         "pcoa.mor.sim.1"              
##  [21] "pcoa.mor.sim.2"               "pcoa.mor.sim.3"              
##  [23] "pcoa.jac.sim.1"               "pcoa.jac.sim.2"              
##  [25] "pcoa.jac.sim.3"               "Hill.logis.Long2"            
##  [27] "Hill.logis.Lat2"              "Hill.logis.conservation"     
##  [29] "Hill.logis.vegetation"        "Hill.logis.bio1"             
##  [31] "Hill.logis.bio2"              "Hill.logis.bio3"             
##  [33] "Hill.logis.bio4"              "Hill.logis.bio5"             
##  [35] "Hill.logis.bio6"              "Hill.logis.bio7"             
##  [37] "Hill.logis.bio8"              "Hill.logis.bio9"             
##  [39] "Hill.logis.bio10"             "Hill.logis.bio11"            
##  [41] "Hill.logis.bio12"             "Hill.logis.bio13"            
##  [43] "Hill.logis.bio14"             "Hill.logis.bio15"            
##  [45] "Hill.logis.bio16"             "Hill.logis.bio17"            
##  [47] "Hill.logis.bio18"             "Hill.logis.bio19"            
##  [49] "Hill.logis.PCA.wclim.1"       "Hill.logis.PCA.wclim.2"      
##  [51] "rich.logis.Long2"             "rich.logis.Lat2"             
##  [53] "rich.logis.conservation"      "rich.logis.vegetation"       
##  [55] "rich.logis.bio1"              "rich.logis.bio2"             
##  [57] "rich.logis.bio3"              "rich.logis.bio4"             
##  [59] "rich.logis.bio5"              "rich.logis.bio6"             
##  [61] "rich.logis.bio7"              "rich.logis.bio8"             
##  [63] "rich.logis.bio9"              "rich.logis.bio10"            
##  [65] "rich.logis.bio11"             "rich.logis.bio12"            
##  [67] "rich.logis.bio13"             "rich.logis.bio14"            
##  [69] "rich.logis.bio15"             "rich.logis.bio16"            
##  [71] "rich.logis.bio17"             "rich.logis.bio18"            
##  [73] "rich.logis.bio19"             "rich.logis.PCA.wclim.1"      
##  [75] "rich.logis.PCA.wclim.2"       "pcoa.logis.jac.PCA.wclim.2.1"
##  [77] "pcoa.logis.jac.PCA.wclim.2.2" "pcoa.logis.jac.PCA.wclim.2.3"
##  [79] "Long2"                        "Lat2"                        
##  [81] "conservation"                 "vegetation"                  
##  [83] "matchin.AF.2d"                "bio1"                        
##  [85] "bio2"                         "bio3"                        
##  [87] "bio4"                         "bio5"                        
##  [89] "bio6"                         "bio7"                        
##  [91] "bio8"                         "bio9"                        
##  [93] "bio10"                        "bio11"                       
##  [95] "bio12"                        "bio13"                       
##  [97] "bio14"                        "bio15"                       
##  [99] "bio16"                        "bio17"                       
## [101] "bio18"                        "bio19"                       
## [103] "PCA.wclim.1"                  "PCA.wclim.2"                 
## [105] "PCA.wclim.3"                  "rich"                        
## [107] "pcoa.jac.1"                   "pcoa.jac.2"                  
## [109] "pcoa.jac.3"                   "rich.MidD.mean"              
## [111] "rich.MidD.sd"                 "pcoa.MidD.jac.1"             
## [113] "pcoa.MidD.jac.2"              "pcoa.MidD.jac.3"
```

```r

par(op)
plot(rich ~ rich.simu.mean, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ])
```

![plot of chunk comparisons of the models](figure/comparisons of the models4.png) 

```r


################# Mid-Domain

```




```r

detach(mammal.data)
```



```r


purl("All_analysis.Rmd")
```

```
## 
## 
## processing file: All_analysis.Rmd
```

```
##   |                                                                         |                                                                 |   0%  |                                                                         |.                                                                |   2%  |                                                                         |..                                                               |   4%  |                                                                         |...                                                              |   5%  |                                                                         |.....                                                            |   7%  |                                                                         |......                                                           |   9%  |                                                                         |.......                                                          |  11%  |                                                                         |........                                                         |  12%  |                                                                         |.........                                                        |  14%  |                                                                         |..........                                                       |  16%  |                                                                         |...........                                                      |  18%  |                                                                         |.............                                                    |  19%  |                                                                         |..............                                                   |  21%  |                                                                         |...............                                                  |  23%  |                                                                         |................                                                 |  25%  |                                                                         |.................                                                |  26%  |                                                                         |..................                                               |  28%  |                                                                         |...................                                              |  30%  |                                                                         |.....................                                            |  32%  |                                                                         |......................                                           |  33%  |                                                                         |.......................                                          |  35%  |                                                                         |........................                                         |  37%  |                                                                         |.........................                                        |  39%  |                                                                         |..........................                                       |  40%  |                                                                         |...........................                                      |  42%  |                                                                         |.............................                                    |  44%  |                                                                         |..............................                                   |  46%  |                                                                         |...............................                                  |  47%  |                                                                         |................................                                 |  49%  |                                                                         |.................................                                |  51%  |                                                                         |..................................                               |  53%  |                                                                         |...................................                              |  54%  |                                                                         |....................................                             |  56%  |                                                                         |......................................                           |  58%  |                                                                         |.......................................                          |  60%  |                                                                         |........................................                         |  61%  |                                                                         |.........................................                        |  63%  |                                                                         |..........................................                       |  65%  |                                                                         |...........................................                      |  67%  |                                                                         |............................................                     |  68%  |                                                                         |..............................................                   |  70%  |                                                                         |...............................................                  |  72%  |                                                                         |................................................                 |  74%  |                                                                         |.................................................                |  75%  |                                                                         |..................................................               |  77%  |                                                                         |...................................................              |  79%  |                                                                         |....................................................             |  81%  |                                                                         |......................................................           |  82%  |                                                                         |.......................................................          |  84%  |                                                                         |........................................................         |  86%  |                                                                         |.........................................................        |  88%  |                                                                         |..........................................................       |  89%  |                                                                         |...........................................................      |  91%  |                                                                         |............................................................     |  93%  |                                                                         |..............................................................   |  95%  |                                                                         |...............................................................  |  96%  |                                                                         |................................................................ |  98%  |                                                                         |.................................................................| 100%
```

```
## output file: All_analysis.R
```

```
## [1] "All_analysis.R"
```





