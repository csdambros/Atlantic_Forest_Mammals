Comparing dispersal, environmental and mid-domain effects on species distribution
============================



## By CSDambros


>html updated in 2014-04-07 21:50:22

### Load required packages, import auxiliary functions and set default graphical parameters



```r

require(raster)
require(rgdal)
require(maptools)
require(vegan)
require(mgcv)

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

wclim <- raster::getData("worldclim", var = "bio", res = 2.5, path = "")
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

# Grouping the latitude and longitude in 2x2 degrees (create new variables
# in the dataframe)
mammal.data$Long2 <- floor(mammal.data$Long/2) * 2 + 1
mammal.data$Lat2 <- floor(mammal.data$Lat/2) * 2 + 1

Longitude.matrix2 <- matrix(mammal.data$Long2, nrow = length(mammal.data$Lat2), 
    ncol = length(LONG.coast.east))
Latitude.matrix2 <- matrix(mammal.data$Lat2, nrow = length(mammal.data$Lat2), 
    ncol = length(LONG.coast.east))

LONG.diff2 <- LONG.matrix.east - Longitude.matrix2
LAT.diff2 <- LAT.matrix.east - Latitude.matrix2

Distance.coast.matrix2 <- sqrt(LONG.diff2^2 + LAT.diff2^2)
mammal.data$Distance.coast2 <- apply(Distance.coast.matrix2, 1, min)

mammal.data$Long.closest2 <- apply(Distance.coast.matrix2, 1, function(x) {
    LONG.coast.east[x == min(x)]
})
mammal.data$Lat.closest2 <- apply(Distance.coast.matrix2, 1, function(x) {
    LAT.coast.east[x == min(x)]
})

rm(Distance.coast.matrix, LAT.diff, LONG.diff, Longitude.matrix, Latitude.matrix, 
    LONG.matrix.east, LAT.matrix.east, LONG.coast.east, LAT.coast.east, LAT.coast, 
    LONG.coast, long.chui, poly, Distance.coast.matrix2, LAT.diff2, LONG.diff2, 
    Longitude.matrix2, Latitude.matrix2)

attach(mammal.data)

#
# mammal.short.PA<-tapply(rep(1,length(LOCALIDADE)),list(LOCALIDADE,ESPECIE),sum)
# mammal.short.PA[is.na(mammal.short.PA)]<-0

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



```r
#
# species.2d<-tapply(rep(1,length(Lat2)),list(paste(Long2,Lat2),ESPECIE),sum)
# 
# species.ab<-tapply(ABUNDANCIA,list(LOCALIDADE,ESPECIE),sum)# CAUTION:
# SOME SPECIES WERE NOT DETECTED IN ANY STUDY RECORDING ABUNDANCE
# locality.species.AB[is.na(locality.species.AB)]<-0
species.ab.loc <- tapply(rep(1, length(Lat)), list(paste(Long, Lat), ESPECIE), 
    sum)
species.ab.loc[is.na(species.ab.loc)] <- 0
# dim(species.ab)
# 


species.loc <- species.ab.loc
species.loc[species.loc > 0] <- 1

#
```


### Summarizing environmental data


```r

environment <- data.frame(unique(cbind(Long, Lat)))
environment <- environment[order(factor(paste(environment[, 1], environment[, 
    2]))), ]

environment$dist.coast <- tapply(as.numeric(Distance.coast), paste(Long, Lat), 
    mean)
environment$conservation <- tapply(as.numeric(CONSERVACAO), paste(Long, Lat), 
    mean)
environment$vegetation <- tapply(as.numeric(FISIONOMIA), paste(Long, Lat), mean)

###

environment[paste("bio", 1:19, sep = "")] <- NA

# Take the mean climatic values of four points around the observed area

for (i in 1:nrow(environment)) {
    temp.clim <- extract(wclim.clipped, extent(environment[i, 1] - 0.05, environment[i, 
        1] + 0.05, environment[i, 2] - 0.05, environment[i, 2] + 0.05))
    environment[i, paste("bio", 1:19, sep = "")] <- colMeans(temp.clim, na.rm = TRUE)
}

environment[, paste("PCA.wclim.", 1:3, sep = "")] <- prcomp(decostand(environment[, 
    paste("bio", 1:19, sep = "")], "standardize"))$x[, 1:3]

# Creating observed statistics

environment$rich <- rowSums(species.loc)
similarity.jac.loc <- vegdist(species.loc, "jaccard")
environment[, paste("pcoa.jac.", 1:3, sep = "")] <- cmdscale(similarity.jac.loc, 
    k = 3)

environment <- data.frame(environment, std = decostand(environment, "standardize"))
```


Regressions of the richness and species composition against the climatic variables


```r

summary(step(lm(rich ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment)))
```

```
## Start:  AIC=226.3
## rich ~ std.dist.coast + std.conservation + std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                    Df Sum of Sq  RSS AIC
## - std.conservation  1       1.2 1308 224
## - std.vegetation    1       6.6 1313 225
## - std.dist.coast    1       9.5 1316 225
## <none>                          1307 226
## - std.PCA.wclim.1   1     126.6 1433 231
## - std.PCA.wclim.2   1     147.3 1454 232
## 
## Step:  AIC=224.4
## rich ~ std.dist.coast + std.vegetation + std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                   Df Sum of Sq  RSS AIC
## - std.vegetation   1       9.5 1317 223
## - std.dist.coast   1      10.7 1318 223
## <none>                         1308 224
## - std.PCA.wclim.1  1     139.4 1447 230
## - std.PCA.wclim.2  1     147.0 1455 230
## 
## Step:  AIC=222.9
## rich ~ std.dist.coast + std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                   Df Sum of Sq  RSS AIC
## - std.dist.coast   1       4.6 1322 221
## <none>                         1317 223
## - std.PCA.wclim.2  1     139.9 1457 228
## - std.PCA.wclim.1  1     168.9 1486 230
## 
## Step:  AIC=221.2
## rich ~ std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                   Df Sum of Sq  RSS AIC
## <none>                         1322 221
## - std.PCA.wclim.2  1       145 1467 227
## - std.PCA.wclim.1  1       207 1529 230
```

```
## 
## Call:
## lm(formula = rich ~ std.PCA.wclim.1 + std.PCA.wclim.2, data = environment)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -7.995 -2.528 -0.613  1.394 17.772 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        8.160      0.495   16.49   <2e-16 ***
## std.PCA.wclim.1    1.674      0.498    3.36   0.0012 ** 
## std.PCA.wclim.2   -1.400      0.498   -2.81   0.0064 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.28 on 72 degrees of freedom
## Multiple R-squared:  0.21,	Adjusted R-squared:  0.189 
## F-statistic:  9.6 on 2 and 72 DF,  p-value: 0.000202
```

```r

summary(step(lm(pcoa.jac.1 ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment)))
```

```
## Start:  AIC=-280.8
## pcoa.jac.1 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                    Df Sum of Sq  RSS  AIC
## - std.PCA.wclim.2   1     0.015 1.53 -282
## - std.vegetation    1     0.017 1.53 -282
## <none>                          1.51 -281
## - std.conservation  1     0.054 1.57 -280
## - std.dist.coast    1     0.167 1.68 -275
## - std.PCA.wclim.1   1     1.415 2.93 -233
## 
## Step:  AIC=-282.1
## pcoa.jac.1 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1
## 
##                    Df Sum of Sq  RSS  AIC
## - std.vegetation    1     0.031 1.56 -283
## <none>                          1.53 -282
## - std.conservation  1     0.047 1.57 -282
## - std.dist.coast    1     0.174 1.70 -276
## - std.PCA.wclim.1   1     1.416 2.94 -235
## 
## Step:  AIC=-282.6
## pcoa.jac.1 ~ std.dist.coast + std.conservation + std.PCA.wclim.1
## 
##                    Df Sum of Sq  RSS  AIC
## - std.conservation  1     0.026 1.58 -283
## <none>                          1.56 -283
## - std.dist.coast    1     0.144 1.70 -278
## - std.PCA.wclim.1   1     1.811 3.37 -227
## 
## Step:  AIC=-283.4
## pcoa.jac.1 ~ std.dist.coast + std.PCA.wclim.1
## 
##                   Df Sum of Sq  RSS  AIC
## <none>                         1.58 -283
## - std.dist.coast   1     0.147 1.73 -279
## - std.PCA.wclim.1  1     2.144 3.73 -221
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ std.dist.coast + std.PCA.wclim.1, data = environment)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.3072 -0.0912 -0.0107  0.1178  0.2803 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      1.89e-16   1.71e-02    0.00    1.000    
## std.dist.coast  -4.70e-02   1.82e-02   -2.59    0.012 *  
## std.PCA.wclim.1 -1.79e-01   1.82e-02   -9.88    5e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.148 on 72 degrees of freedom
## Multiple R-squared:  0.649,	Adjusted R-squared:  0.64 
## F-statistic: 66.7 on 2 and 72 DF,  p-value: <2e-16
```

```r

summary(step(lm(pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment)))
```

```
## Start:  AIC=-271
## pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                    Df Sum of Sq  RSS  AIC
## - std.vegetation    1     0.011 1.74 -272
## - std.dist.coast    1     0.015 1.74 -272
## - std.conservation  1     0.016 1.74 -272
## <none>                          1.72 -271
## - std.PCA.wclim.1   1     0.049 1.77 -271
## - std.PCA.wclim.2   1     0.876 2.60 -242
## 
## Step:  AIC=-272.5
## pcoa.jac.2 ~ std.dist.coast + std.conservation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                    Df Sum of Sq  RSS  AIC
## - std.conservation  1     0.028 1.76 -273
## - std.dist.coast    1     0.037 1.77 -273
## <none>                          1.74 -272
## - std.PCA.wclim.1   1     0.075 1.81 -271
## - std.PCA.wclim.2   1     0.896 2.63 -243
## 
## Step:  AIC=-273.3
## pcoa.jac.2 ~ std.dist.coast + std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                   Df Sum of Sq  RSS  AIC
## - std.dist.coast   1     0.037 1.80 -274
## <none>                         1.76 -273
## - std.PCA.wclim.1  1     0.054 1.82 -273
## - std.PCA.wclim.2  1     0.880 2.64 -245
## 
## Step:  AIC=-273.8
## pcoa.jac.2 ~ std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                   Df Sum of Sq  RSS  AIC
## - std.PCA.wclim.1  1     0.033 1.83 -274
## <none>                         1.80 -274
## - std.PCA.wclim.2  1     0.856 2.66 -246
## 
## Step:  AIC=-274.4
## pcoa.jac.2 ~ std.PCA.wclim.2
## 
##                   Df Sum of Sq  RSS  AIC
## <none>                         1.83 -274
## - std.PCA.wclim.2  1     0.856 2.69 -248
```

```
## 
## Call:
## lm(formula = pcoa.jac.2 ~ std.PCA.wclim.2, data = environment)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.3786 -0.1060  0.0023  0.1086  0.3449 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      1.42e-17   1.83e-02    0.00        1    
## std.PCA.wclim.2 -1.08e-01   1.84e-02   -5.84  1.3e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.158 on 73 degrees of freedom
## Multiple R-squared:  0.319,	Adjusted R-squared:  0.309 
## F-statistic: 34.1 on 1 and 73 DF,  p-value: 1.34e-07
```


Mantel tests


```r

geodist.loc <- dist(environment[, c("Lat", "Long")])
envdist.loc <- dist(environment[grep("std", colnames(environment))][3:24])

plot(1 - similarity.jac.loc ~ geodist.loc, col = 2)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r
plot(1 - similarity.jac.loc ~ envdist.loc, col = 2)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 

```r

mantel(similarity.jac.loc, geodist.loc)
```

```
## 
## Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel(xdis = similarity.jac.loc, ydis = geodist.loc) 
## 
## Mantel statistic r: 0.479 
##       Significance: 0.001 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0530 0.0654 0.0773 0.0895 
## 
## Based on 999 permutations
```

```r
mantel(similarity.jac.loc, envdist.loc)
```

```
## 
## Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel(xdis = similarity.jac.loc, ydis = envdist.loc) 
## 
## Mantel statistic r: 0.321 
##       Significance: 0.001 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0522 0.0701 0.0869 0.1059 
## 
## Based on 999 permutations
```

```r

mantel.partial(similarity.jac.loc, geodist.loc, envdist.loc)
```

```
## 
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = similarity.jac.loc, ydis = geodist.loc,      zdis = envdist.loc) 
## 
## Mantel statistic r: 0.378 
##       Significance: 0.001 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0599 0.0748 0.0930 0.1076 
## 
## Based on 999 permutations
```

```r
mantel.partial(similarity.jac.loc, envdist.loc, geodist.loc)
```

```
## 
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = similarity.jac.loc, ydis = envdist.loc,      zdis = geodist.loc) 
## 
## Mantel statistic r: 0.0457 
##       Significance: 0.17 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0624 0.0801 0.0936 0.1205 
## 
## Based on 999 permutations
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

environment.2d$dist.coast <- tapply(as.numeric(Distance.coast2), paste(Long2, 
    Lat2), mean)
environment.2d$conservation <- tapply(as.numeric(CONSERVACAO), paste(Long2, 
    Lat2), mean)
environment.2d$vegetation <- tapply(as.numeric(FISIONOMIA), paste(Long2, Lat2), 
    mean)

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
similarity.jac <- vegdist(species.2d, "jaccard")
environment.2d[, paste("pcoa.jac.", 1:3, sep = "")] <- cmdscale(similarity.jac, 
    k = 3)

# plot(pcoa.jac.1~Lat2,data=environment.2d)

environment.2d <- data.frame(environment.2d, std = decostand(environment.2d, 
    "standardize"))

environment.2d$matchin.AF.2d <- match(paste(environment.2d[, 1], environment.2d[, 
    2]), paste(environment.AF.2d[, 1], environment.AF.2d[, 2]))
```


Regressions of the richness and species composition against the climatic variables (2d)


```r

summary(step(lm(rich ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment.2d)))
```

```
## Start:  AIC=99.07
## rich ~ std.dist.coast + std.conservation + std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                    Df Sum of Sq RSS   AIC
## - std.dist.coast    1       0.0 740  97.1
## - std.PCA.wclim.1   1      23.2 763  97.9
## - std.vegetation    1      59.0 799  99.1
## <none>                          740  99.1
## - std.PCA.wclim.2   1     116.8 857 100.9
## - std.conservation  1     133.1 873 101.4
## 
## Step:  AIC=97.07
## rich ~ std.conservation + std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                    Df Sum of Sq RSS  AIC
## - std.PCA.wclim.1   1      30.7 771 96.1
## <none>                          740 97.1
## - std.vegetation    1      80.9 821 97.8
## - std.PCA.wclim.2   1     125.0 865 99.1
## - std.conservation  1     150.7 891 99.9
## 
## Step:  AIC=96.12
## rich ~ std.conservation + std.vegetation + std.PCA.wclim.2
## 
##                    Df Sum of Sq RSS   AIC
## <none>                          771  96.1
## - std.vegetation    1      63.1 834  96.2
## - std.PCA.wclim.2   1     136.8 908  98.4
## - std.conservation  1     202.6 973 100.2
```

```
## 
## Call:
## lm(formula = rich ~ std.conservation + std.vegetation + std.PCA.wclim.2, 
##     data = environment.2d)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -8.586 -3.841 -0.921  3.182 11.931 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         14.08       1.16   12.13  3.2e-11 ***
## std.conservation    -3.46       1.44   -2.40    0.025 *  
## std.vegetation       1.83       1.37    1.34    0.193    
## std.PCA.wclim.2     -2.49       1.26   -1.98    0.061 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.92 on 22 degrees of freedom
## Multiple R-squared:  0.413,	Adjusted R-squared:  0.333 
## F-statistic: 5.17 on 3 and 22 DF,  p-value: 0.00744
```

```r
summary(step(lm(pcoa.jac.1 ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment.2d)))
```

```
## Start:  AIC=-102
## pcoa.jac.1 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                    Df Sum of Sq   RSS    AIC
## - std.conservation  1     0.019 0.344 -102.5
## - std.dist.coast    1     0.024 0.349 -102.1
## <none>                          0.324 -102.0
## - std.vegetation    1     0.036 0.360 -101.3
## - std.PCA.wclim.2   1     0.157 0.481  -93.7
## - std.PCA.wclim.1   1     0.679 1.003  -74.6
## 
## Step:  AIC=-102.5
## pcoa.jac.1 ~ std.dist.coast + std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                   Df Sum of Sq   RSS    AIC
## - std.dist.coast   1     0.013 0.357 -103.5
## <none>                         0.344 -102.5
## - std.vegetation   1     0.118 0.461  -96.8
## - std.PCA.wclim.2  1     0.138 0.482  -95.7
## - std.PCA.wclim.1  1     0.902 1.246  -71.0
## 
## Step:  AIC=-103.5
## pcoa.jac.1 ~ std.vegetation + std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                   Df Sum of Sq   RSS    AIC
## <none>                         0.357 -103.5
## - std.PCA.wclim.2  1     0.129 0.486  -97.5
## - std.vegetation   1     0.181 0.538  -94.8
## - std.PCA.wclim.1  1     1.178 1.534  -67.6
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2, data = environment.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2707 -0.0796 -0.0015  0.0652  0.3421 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     -4.77e-17   2.50e-02    0.00    1.000    
## std.vegetation  -8.58e-02   2.57e-02   -3.34    0.003 ** 
## std.PCA.wclim.1  2.18e-01   2.56e-02    8.52    2e-08 ***
## std.PCA.wclim.2  7.19e-02   2.55e-02    2.82    0.010 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.127 on 22 degrees of freedom
## Multiple R-squared:  0.795,	Adjusted R-squared:  0.767 
## F-statistic: 28.5 on 3 and 22 DF,  p-value: 9.22e-08
```

```r
summary(step(lm(pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment.2d)))
```

```
## Start:  AIC=-94.54
## pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                    Df Sum of Sq   RSS   AIC
## - std.PCA.wclim.2   1    0.0320 0.464 -94.7
## <none>                          0.432 -94.5
## - std.dist.coast    1    0.0394 0.471 -94.3
## - std.conservation  1    0.0684 0.500 -92.7
## - std.PCA.wclim.1   1    0.0772 0.509 -92.3
## - std.vegetation    1    0.0916 0.524 -91.5
## 
## Step:  AIC=-94.68
## pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1
## 
##                    Df Sum of Sq   RSS   AIC
## <none>                          0.464 -94.7
## - std.conservation  1    0.0428 0.507 -94.4
## - std.PCA.wclim.1   1    0.0611 0.525 -93.5
## - std.dist.coast    1    0.0650 0.529 -93.3
## - std.vegetation    1    0.1283 0.592 -90.3
```

```
## 
## Call:
## lm(formula = pcoa.jac.2 ~ std.dist.coast + std.conservation + 
##     std.vegetation + std.PCA.wclim.1, data = environment.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.3014 -0.0785  0.0149  0.1250  0.1673 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      -1.51e-18   2.91e-02    0.00    1.000  
## std.dist.coast   -6.15e-02   3.59e-02   -1.71    0.101  
## std.conservation  5.09e-02   3.66e-02    1.39    0.178  
## std.vegetation   -9.58e-02   3.97e-02   -2.41    0.025 *
## std.PCA.wclim.1  -5.75e-02   3.46e-02   -1.66    0.111  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.149 on 21 degrees of freedom
## Multiple R-squared:  0.54,	Adjusted R-squared:  0.453 
## F-statistic: 6.17 on 4 and 21 DF,  p-value: 0.0019
```


Mantel tests


```r

geodist.2d <- dist(environment.2d[, c("Lat2", "Long2")])
envdist.2d <- dist(environment.2d[grep("std", colnames(environment.2d))][3:24])

plot(1 - similarity.jac ~ geodist.2d, col = 2)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-61.png) 

```r
plot(1 - similarity.jac ~ envdist.2d, col = 2)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-62.png) 

```r

mantel(similarity.jac, geodist.2d)
```

```
## 
## Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel(xdis = similarity.jac, ydis = geodist.2d) 
## 
## Mantel statistic r: 0.493 
##       Significance: 0.001 
## 
## Upper quantiles of permutations (null model):
##   90%   95% 97.5%   99% 
## 0.113 0.149 0.180 0.216 
## 
## Based on 999 permutations
```

```r
mantel(similarity.jac, envdist.2d)
```

```
## 
## Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel(xdis = similarity.jac, ydis = envdist.2d) 
## 
## Mantel statistic r: 0.482 
##       Significance: 0.001 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0955 0.1200 0.1387 0.1572 
## 
## Based on 999 permutations
```

```r

mantel.partial(similarity.jac, geodist.2d, envdist.2d)
```

```
## 
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = similarity.jac, ydis = geodist.2d, zdis = envdist.2d) 
## 
## Mantel statistic r: 0.214 
##       Significance: 0.01 
## 
## Upper quantiles of permutations (null model):
##   90%   95% 97.5%   99% 
## 0.111 0.139 0.171 0.212 
## 
## Based on 999 permutations
```

```r
mantel.partial(similarity.jac, envdist.2d, geodist.2d)
```

```
## 
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = similarity.jac, ydis = envdist.2d, zdis = geodist.2d) 
## 
## Mantel statistic r: 0.177 
##       Significance: 0.016 
## 
## Upper quantiles of permutations (null model):
##   90%   95% 97.5%   99% 
## 0.102 0.128 0.154 0.194 
## 
## Based on 999 permutations
```



### Plot a map of the original data and data grouped in quadrants

![plot of chunk plotting maps](figure/plotting maps.png) 




# Running models of species distribution

### There are four different models in this project: Mid Domain, Neutral model (1 and 2), direct effect of environment on species richness and GLMs for each species along environmental gradients

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

MidD.sim <- array(NA, dim = c(nrow(connect.2d), ncol(species.2d), reps))

for (i in 1:reps) {
    MidD.sim[, , i] <- MidD(connect.2d, colSums(species.2d), simplify = T)
}

#### Extracting summary statistics

rich.MidD.vec <- apply(MidD.sim, 3, rowSums)

similarity.MidD.vec <- apply(MidD.sim, 3, function(mat) as.vector(as.matrix(vegdist(mat, 
    "jaccard"))))


similarity.MidD.mean <- matrix(rowMeans(similarity.MidD.vec), nrow(species.2d), 
    nrow(species.2d))

environment.2d$rich.MidD.mean <- rowMeans(rich.MidD.vec)
environment.2d$rich.MidD.sd <- apply(rich.MidD.vec, 1, sd)

similarity.MidD.mean <- matrix(rowMeans(similarity.MidD.vec), nrow(species.2d), 
    nrow(species.2d))
similarity.MidD.sd <- matrix(apply(similarity.MidD.vec, 1, sd), nrow(species.2d), 
    nrow(species.2d))

environment.2d[, paste("pcoa.MidD.jac.", 1:3, sep = "")] <- (prcomp(similarity.MidD.mean)$x)[, 
    1:3]

# plot(pcoa.jac.1~pcoa.MidD.jac.1,data=environment.2d)
```


### Plot species diversity (Richness; A) and species composition (B) predicted by the Mid-Domain model 

>Areas with similar colors in B represent areas with similar composition of species


![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 



Neutral model
===============
>based in Economo and Keitt (2008)

>This model spreads a given number of individuals (*N*) in each quadrant. Initially all individuals in all quadrants are from a single species. Then the model is run for several generations. In each generation all individuals reproduce and die. Based on the parameters *v* and *m*, the indviduals in the offspring can speciate (point speciation) and migrate. In a single quadrant, two individuals can be from the same species if none of them have speciated with rate *v* in the previous generation, AND if they had the same parent or different parents from the same species, either from the same quadrant or not (none, one or both have migrated). Eventually the diversity within and between quadrants reach a steady-state. At this point the model stops. See the complementary function anetwork (in "anetwork.R") for details. This function actually does not simulate the whole process described here, but instead uses a method borrowed from the coalescent theory of population genetics. See Neutral model 2 for a real simulation which gives the same results.


![plot of chunk plot.connectivity](figure/plot.connectivity.png) 


### This code optimizes the migration (*m*) and speciation (*v*) (single) parameters in order to better fit the species composition change in the observed data.

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


### Uses the optimized parameters to recreate the neutral simulation (run one more time)


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
    
    neutral.sim.ind <- read.table("resu.simu10k.txt", header = TRUE)
    
} else {
    
    # m=0.01790858 #migration rate v=0.0017789240 #speciation rate N=100
    # #Number of individuals in each node
    
    reps2 <- 10000
    initbuffer <- 10000
    
    sp.abund2 <- simunetwork(N, M, v, nruns = initbuffer, simplify = FALSE)
    
    neutral.sim.ind <- matrix(NA, sum(N), reps2)
    
    for (rep in 1:reps2) {
        
        sp.abund2 <- simunetwork(N, M, v, nruns = 300, simplify = FALSE, sp.abund = sp.abund2)
        
        neutral.sim.ind[, rep] <- sp.abund2
        
        # print(rep)
    }
    
    write.table(neutral.sim.ind, file = "neutral.sim.ind10k.txt")
    
    # tapply(rep(1,sum(N)),list(rep(1:length(N),N),neutral.sim.ind[,1]),sum)
    
}

```




```r


similarity.neutral.total<-array(NA,dim=c(nrow(environment.AF.2d),nrow(environment.AF.2d),ncol(neutral.sim.ind)))
similarity.neutral.jac.total<-array(NA,dim=c(nrow(environment.AF.2d),nrow(environment.AF.2d),ncol(neutral.sim.ind)))
similarity.neutral<-array(NA,dim=c(nrow(environment.2d),nrow(environment.2d),ncol(neutral.sim.ind)))
similarity.neutral.jac<-array(NA,dim=c(nrow(environment.2d),nrow(environment.2d),ncol(neutral.sim.ind)))
hill.neutral<-matrix(NA,nrow(environment.AF.2d),ncol(neutral.sim.ind))
rich.neutral<-matrix(NA,nrow(environment.AF.2d),ncol(neutral.sim.ind))

for(i in 1:ncol(neutral.sim.ind)){
  
  ver<-tapply(rep(1,sum(N)),list(rep(1:length(N),N),neutral.sim.ind[,i]),sum)
  ver[is.na(ver)]<-0
  
  similarity.neutral.total[,,i]<-as.matrix(vegdist(ver,"morisita"))
  similarity.neutral.jac.total[,,i]<-as.matrix(vegdist(ver,"jaccard"))
  
  similarity.neutral[,,i]<-as.matrix(vegdist(ver,"morisita"))[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d]
  similarity.neutral.jac[,,i]<-as.matrix(vegdist(ver,"jaccard"))[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d]
  
  hill.neutral[,i]<-(1/(rowSums((ver/rowSums(ver))^2)*(rowSums(ver)-1)/(rowSums(ver))))
  rich.neutral[,i]<-rowSums(ver>0)
  }


environment.AF.2d$rich.neutral.mean<-rowMeans(rich.neutral)
environment.AF.2d$rich.neutral.sd<-apply(rich.neutral,1,sd)

environment.AF.2d$hill.neutral.mean<-rowMeans(hill.neutral)
environment.AF.2d$hill.neutral.sd<-apply(hill.neutral,1,sd)

environment.AF.2d[,paste("pcoa.mor.neutral.total.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA

environment.AF.2d[,paste("pcoa.jac.neutral.total.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA

environment.AF.2d[,paste("pcoa.mor.neutral.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA

environment.AF.2d[,paste("pcoa.jac.neutral.",1:3,sep="")]<-environment.AF.2d[,paste("pcoa.mor.AF",1:3,sep="")]*NA


environment.AF.2d[,paste("pcoa.mor.neutral.total.",1:3,sep="")]<-#PCA.AF[,1:3]
  prcomp(apply(similarity.neutral.total[,,round(.6*ncol(neutral.sim.ind)):ncol(neutral.sim.ind)],2,function(x)rowMeans(x)))$x[,1:3]

environment.AF.2d[,paste("pcoa.jac.neutral.total.",1:3,sep="")]<-#PCA.AF[,1:3]
  -prcomp(apply(similarity.neutral.jac.total[,,round(.6*ncol(neutral.sim.ind)):ncol(neutral.sim.ind)],2,function(x)rowMeans(x)))$x[,1:3]

similarity.neutral.mean<-apply(similarity.neutral[,,round(.6*ncol(neutral.sim.ind)):ncol(neutral.sim.ind)],2,function(x)rowMeans(x))

environment.AF.2d[environment.2d$matchin.AF.2d,paste("pcoa.mor.neutral.",1:3,sep="")]<-#PCA.AF[,1:3]
  prcomp(similarity.neutral.mean)$x[,1:3]

similarity.neutral.jac.mean<-apply(similarity.neutral.jac[,,round(.6*ncol(neutral.sim.ind)):ncol(neutral.sim.ind)],2,function(x)rowMeans(x))

environment.AF.2d[environment.2d$matchin.AF.2d,paste("pcoa.jac.neutral.",1:3,sep="")]<-#PCA.AF[,1:3]
  -prcomp(similarity.neutral.jac.mean)$x[,1:3]

#cor(environment.AF.2d$pcoa.mor.1,environment.AF.2d$pcoa.mor.neutral.1,use="c")

```


### Plot species diversity (Richness; A) and species composition (B) predicted by the Neutral Model (2nd) 

>Areas with similar colors in B represent areas with similar composition of species


![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


### Comparison of the species composition observed and predicted by the Mid Domain and Neutral (x2) models along the Latitudinal gradient

>The Latitudinal gradient is one of the best known gradients to affect the species diversity in the area.


```r

par(mfrow = c(1, 2))

plot(environment.2d$Lat2, decostand(environment.2d$pcoa.jac.1, "range"), bg = "gold", 
    pch = 22)

points(environment.2d$Lat2, decostand(environment.2d$pcoa.MidD.jac.1, "range"), 
    bg = 2, pch = 23)

points((environment.AF.2d$AF.Lat2), decostand(environment.AF.2d$pcoa.mor.1, 
    "range", na.rm = T), bg = 4, pch = 24)

points((environment.AF.2d$AF.Lat2), decostand(environment.AF.2d$pcoa.mor.neutral.1, 
    "range", na.rm = T), bg = 1, pch = 3, cex = 2)

plot.new()

legend("topleft", c("Observed", "MidD", "Neutral-ana", "Neutral-simu"), pch = c(22, 
    23, 24, 3), pt.bg = c("gold", 2, 4, 3), bty = "n", y.intersp = 1.2, cex = 2)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

```r

#
# plot(decostand(environment.AF.2d$pcoa.mor.1,'range',na.rm=T),decostand(environment.AF.2d$pcoa.mor.neutral.1,'range',na.rm=T),bg=1,pch=3,cex=2)
# plot(environment.AF.2d$pcoa.mor.1,environment.AF.2d$pcoa.mor.neutral.1,bg=1,pch=3,cex=2)
```


Using GLMs to estimate the probability of each species occurring along each environmental gradient
==================================

>Based on these coefficients, calculates 
>> The expected number of species in each quadrant  
>> The probability of two individuals being from different species (PIE)  
>> The morisita-horn index of species similarity  

> Calculates the glm coefficients (slope and intercept) for each species
> glm.results gives the glm results (with p values) for all species
> glm.coef gives just the coefficients for each species and variable the other variables are meaningless and are necessary just to create the others


```r

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

# plot(species.2d[,22]~unlist(glm.coef.2d[[2]][[2]]),type='n')

# for(i in 1:64){
# points(probs.2d[[2]][i,]~unlist(glm.coef.2d[[2]][[2]]),pch=21,bg=2,type='l')
# }

# plot(colSums(probs.2d[[2]])~unlist(glm.coef.2d[[2]][[2]]))


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

reps = 999

similarity.logis.jac <- array(NA, dim = c(nrow(environment.2d), nrow(environment.2d), 
    reps))

# dim(similarity.logis.jac)
# 
# for(j in glm.variables.2d){
# 
# # j=glm.variables.2d[[25]]
# 
# for(k in 1:reps){
# 
# species.2d.logis<-species.2d*0
# 
# for(i in 1:nrow(species.2d.logis)){
# species.2d.logis[i,sample(ncol(species.2d.logis),sum(species.2d[i,]),prob=probs.2d[[j]][,i])]
# <- 1 }
# similarity.logis.jac[,,k]<-as.matrix(vegdist(species.2d.logis,'jaccard'))
# }
# 
# environment.AF.2d[,paste('pcoa.logis.jac.',j,'.',1:3,sep='')]<-NA
# 
# environment.AF.2d[environment.2d$matchin.AF.2d,paste('pcoa.logis.jac.',j,'.',1:3,sep='')]<-
# prcomp(apply(similarity.logis.jac,2,rowMeans))$x[,1:3]
# 
# }

species.2d.logis.ls <- list()
similarity.logis.jac.ls <- list()

for (j in glm.variables.2d) {
    
    # j=glm.variables.2d[[25]]
    
    species.2d.logis.reps <- list()
    
    for (k in 1:reps) {
        
        species.2d.logis <- species.2d * 0
        
        for (i in 1:ncol(species.2d.logis)) {
            
            species.2d.logis[sample(nrow(species.2d.logis), sum(species.2d[, 
                i]), prob = probs.2d[[j]][i, ]), i] <- 1
            
        }
        
        similarity.logis.jac[, , k] <- as.matrix(vegdist(species.2d.logis, "jaccard"))
        species.2d.logis.reps[[k]] <- species.2d.logis
        
    }
    
    similarity.logis.jac.ls[[j]] <- similarity.logis.jac
    species.2d.logis.ls[[j]] <- species.2d.logis.reps
    
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
    
    ver <- do.call(cbind, lapply(species.2d.logis.ls[[i]], rowSums))
    
    plot(rowSums(species.2d) ~ environment.2d[[i]], xlab = i, ylab = "Diversity", 
        axes = FALSE, ylim = c(0, max(rowSums(species.2d))), pch = 21, col = 0, 
        bg = "dark grey")
    
    points(rowMeans(ver) ~ environment.2d[[i]], xlab = i, ylab = "Diversity", 
        axes = FALSE, ylim = c(0, max(rowSums(species.2d))), pch = 21, col = 0, 
        bg = 1)
    
    points(rich.logis[[paste("rich.logis.", i, sep = "")]] ~ environment.2d[[i]], 
        pch = 21, bg = color.select(rich.logis[[paste("rich.logis.", i, sep = "")]]))
    
    box()
    
}
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

```
## Warning: "axes" is not a graphical parameter
```

![plot of chunk plot_sdms](figure/plot_sdms2.png) 


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

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-111.png) 

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

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-112.png) 

```r

detach(environment.2d)
```



```r

environment.AF.2d[, colnames(environment.2d)] <- NA
environment.AF.2d[environment.2d$match, colnames(environment.2d)] <- environment.2d
```




```r

MR1.1 <- step(lm(rich ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment.AF.2d))
```

```
## Start:  AIC=99.07
## rich ~ std.dist.coast + std.conservation + std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                    Df Sum of Sq RSS   AIC
## - std.dist.coast    1       0.0 740  97.1
## - std.PCA.wclim.1   1      23.2 763  97.9
## - std.vegetation    1      59.0 799  99.1
## <none>                          740  99.1
## - std.PCA.wclim.2   1     116.8 857 100.9
## - std.conservation  1     133.1 873 101.4
## 
## Step:  AIC=97.07
## rich ~ std.conservation + std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                    Df Sum of Sq RSS  AIC
## - std.PCA.wclim.1   1      30.7 771 96.1
## <none>                          740 97.1
## - std.vegetation    1      80.9 821 97.8
## - std.PCA.wclim.2   1     125.0 865 99.1
## - std.conservation  1     150.7 891 99.9
## 
## Step:  AIC=96.12
## rich ~ std.conservation + std.vegetation + std.PCA.wclim.2
## 
##                    Df Sum of Sq RSS   AIC
## <none>                          771  96.1
## - std.vegetation    1      63.1 834  96.2
## - std.PCA.wclim.2   1     136.8 908  98.4
## - std.conservation  1     202.6 973 100.2
```

```r

MR1.1 <- (lm(rich ~ 1, data = environment.AF.2d))

summary(MR1.2 <- update(MR1.1, ~. + rich.neutral.mean))
```

```
## 
## Call:
## lm(formula = rich ~ rich.neutral.mean, data = environment.AF.2d)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -9.84  -6.79  -1.04   5.10  13.80 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)
## (Intercept)         -13.30      18.15   -0.73     0.47
## rich.neutral.mean     3.99       2.64    1.51     0.14
## 
## Residual standard error: 7.07 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.087,	Adjusted R-squared:  0.049 
## F-statistic: 2.29 on 1 and 24 DF,  p-value: 0.143
```

```r
summary(MR1.3 <- update(MR1.1, ~. + rich.MidD.mean))
```

```
## 
## Call:
## lm(formula = rich ~ rich.MidD.mean, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -11.349  -6.197   0.209   5.255  13.155 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)
## (Intercept)       6.218      5.745    1.08     0.29
## rich.MidD.mean    0.558      0.396    1.41     0.17
## 
## Residual standard error: 7.11 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.0765,	Adjusted R-squared:  0.038 
## F-statistic: 1.99 on 1 and 24 DF,  p-value: 0.171
```

```r

AIC(MR1.1, MR1.2, MR1.3)
```

```
##       df   AIC
## MR1.1  2 179.8
## MR1.2  3 179.4
## MR1.3  3 179.7
```

```r

summary(MC1.1 <- step(lm(pcoa.jac.1 ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment.AF.2d)))
```

```
## Start:  AIC=-102
## pcoa.jac.1 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                    Df Sum of Sq   RSS    AIC
## - std.conservation  1     0.019 0.344 -102.5
## - std.dist.coast    1     0.024 0.349 -102.1
## <none>                          0.324 -102.0
## - std.vegetation    1     0.036 0.360 -101.3
## - std.PCA.wclim.2   1     0.157 0.481  -93.7
## - std.PCA.wclim.1   1     0.679 1.003  -74.6
## 
## Step:  AIC=-102.5
## pcoa.jac.1 ~ std.dist.coast + std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2
## 
##                   Df Sum of Sq   RSS    AIC
## - std.dist.coast   1     0.013 0.357 -103.5
## <none>                         0.344 -102.5
## - std.vegetation   1     0.118 0.461  -96.8
## - std.PCA.wclim.2  1     0.138 0.482  -95.7
## - std.PCA.wclim.1  1     0.902 1.246  -71.0
## 
## Step:  AIC=-103.5
## pcoa.jac.1 ~ std.vegetation + std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                   Df Sum of Sq   RSS    AIC
## <none>                         0.357 -103.5
## - std.PCA.wclim.2  1     0.129 0.486  -97.5
## - std.vegetation   1     0.181 0.538  -94.8
## - std.PCA.wclim.1  1     1.178 1.534  -67.6
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2707 -0.0796 -0.0015  0.0652  0.3421 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      1.34e-17   2.50e-02    0.00    1.000    
## std.vegetation  -8.58e-02   2.57e-02   -3.34    0.003 ** 
## std.PCA.wclim.1  2.18e-01   2.56e-02    8.52    2e-08 ***
## std.PCA.wclim.2  7.19e-02   2.55e-02    2.82    0.010 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.127 on 22 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.795,	Adjusted R-squared:  0.767 
## F-statistic: 28.5 on 3 and 22 DF,  p-value: 9.22e-08
```

```r

summary(MC1.2 <- update(MC1.1, ~. + pcoa.jac.neutral.1))
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2 + pcoa.jac.neutral.1, data = environment.AF.2d)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.20944 -0.07443  0.00627  0.04848  0.26666 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)  
## (Intercept)        -1.89e-17   2.26e-02    0.00    1.000  
## std.vegetation     -6.69e-02   2.45e-02   -2.73    0.013 *
## std.PCA.wclim.1     4.17e-02   7.66e-02    0.54    0.592  
## std.PCA.wclim.2    -2.65e-02   4.68e-02   -0.57    0.578  
## pcoa.jac.neutral.1  5.32e-01   2.20e-01    2.42    0.025 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.115 on 21 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.84,	Adjusted R-squared:  0.809 
## F-statistic: 27.5 on 4 and 21 DF,  p-value: 4.35e-08
```

```r
summary(MC1.3 <- update(MC1.1, ~. + pcoa.MidD.jac.1))
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ std.vegetation + std.PCA.wclim.1 + 
##     std.PCA.wclim.2 + pcoa.MidD.jac.1, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2334 -0.0549 -0.0220  0.0330  0.3144 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      2.35e-17   2.39e-02    0.00    1.000  
## std.vegetation  -6.17e-02   2.81e-02   -2.19    0.040 *
## std.PCA.wclim.1  8.29e-02   8.14e-02    1.02    0.320  
## std.PCA.wclim.2  1.16e-02   4.23e-02    0.27    0.787  
## pcoa.MidD.jac.1  1.88e-01   1.08e-01    1.74    0.096 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.122 on 21 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.821,	Adjusted R-squared:  0.787 
## F-statistic: 24.1 on 4 and 21 DF,  p-value: 1.36e-07
```

```r

summary(MC1.4 <- update(MC1.1, ~pcoa.jac.neutral.1))
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ pcoa.jac.neutral.1, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.1826 -0.1044  0.0059  0.0728  0.3556 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        -2.28e-17   2.63e-02    0.00        1    
## pcoa.jac.neutral.1  5.85e-01   6.85e-02    8.55  9.6e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.134 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.753,	Adjusted R-squared:  0.742 
## F-statistic:   73 on 1 and 24 DF,  p-value: 9.62e-09
```

```r
summary(MC1.5 <- update(MC1.1, ~pcoa.MidD.jac.1))
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ pcoa.MidD.jac.1, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2096 -0.0752 -0.0287  0.0377  0.3618 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     2.80e-17   2.52e-02    0.00        1    
## pcoa.MidD.jac.1 2.85e-01   3.16e-02    9.01  3.6e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.129 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.772,	Adjusted R-squared:  0.762 
## F-statistic: 81.2 on 1 and 24 DF,  p-value: 3.6e-09
```

```r

AIC(MC1.1, MC1.2, MC1.3, MC1.4, MC1.5)
```

```
##       df    AIC
## MC1.1  5 -27.72
## MC1.2  6 -32.10
## MC1.3  6 -29.24
## MC1.4  3 -26.80
## MC1.5  3 -28.90
```

```r

summary(MC2.1 <- step(lm(pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
    std.PCA.wclim.1 + std.PCA.wclim.2, data = environment.AF.2d)))
```

```
## Start:  AIC=-94.54
## pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1 + std.PCA.wclim.2
## 
##                    Df Sum of Sq   RSS   AIC
## - std.PCA.wclim.2   1    0.0320 0.464 -94.7
## <none>                          0.432 -94.5
## - std.dist.coast    1    0.0394 0.471 -94.3
## - std.conservation  1    0.0684 0.500 -92.7
## - std.PCA.wclim.1   1    0.0772 0.509 -92.3
## - std.vegetation    1    0.0916 0.524 -91.5
## 
## Step:  AIC=-94.68
## pcoa.jac.2 ~ std.dist.coast + std.conservation + std.vegetation + 
##     std.PCA.wclim.1
## 
##                    Df Sum of Sq   RSS   AIC
## <none>                          0.464 -94.7
## - std.conservation  1    0.0428 0.507 -94.4
## - std.PCA.wclim.1   1    0.0611 0.525 -93.5
## - std.dist.coast    1    0.0650 0.529 -93.3
## - std.vegetation    1    0.1283 0.592 -90.3
```

```
## 
## Call:
## lm(formula = pcoa.jac.2 ~ std.dist.coast + std.conservation + 
##     std.vegetation + std.PCA.wclim.1, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.3014 -0.0785  0.0149  0.1250  0.1673 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      -3.46e-17   2.91e-02    0.00    1.000  
## std.dist.coast   -6.15e-02   3.59e-02   -1.71    0.101  
## std.conservation  5.09e-02   3.66e-02    1.39    0.178  
## std.vegetation   -9.58e-02   3.97e-02   -2.41    0.025 *
## std.PCA.wclim.1  -5.75e-02   3.46e-02   -1.66    0.111  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.149 on 21 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.54,	Adjusted R-squared:  0.453 
## F-statistic: 6.17 on 4 and 21 DF,  p-value: 0.0019
```

```r

summary(MC2.2 <- update(MC2.1, ~. + pcoa.jac.neutral.2))
```

```
## 
## Call:
## lm(formula = pcoa.jac.2 ~ std.dist.coast + std.conservation + 
##     std.vegetation + std.PCA.wclim.1 + pcoa.jac.neutral.2, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.3564 -0.0526  0.0206  0.0965  0.1717 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)  
## (Intercept)        -3.92e-17   2.76e-02    0.00    1.000  
## std.dist.coast     -5.44e-02   3.41e-02   -1.59    0.126  
## std.conservation    7.46e-02   3.68e-02    2.02    0.056 .
## std.vegetation     -7.97e-02   3.85e-02   -2.07    0.052 .
## std.PCA.wclim.1    -3.38e-02   3.51e-02   -0.96    0.346  
## pcoa.jac.neutral.2  2.10e-01   1.12e-01    1.87    0.076 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.14 on 20 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.609,	Adjusted R-squared:  0.511 
## F-statistic: 6.23 on 5 and 20 DF,  p-value: 0.00123
```

```r
summary(MC2.3 <- update(MC2.1, ~. + pcoa.MidD.jac.2))
```

```
## 
## Call:
## lm(formula = pcoa.jac.2 ~ std.dist.coast + std.conservation + 
##     std.vegetation + std.PCA.wclim.1 + pcoa.MidD.jac.2, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.3488 -0.0501  0.0221  0.1007  0.1719 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      -3.78e-17   2.79e-02    0.00    1.000  
## std.dist.coast   -6.04e-02   3.43e-02   -1.76    0.093 .
## std.conservation  6.89e-02   3.65e-02    1.89    0.073 .
## std.vegetation   -7.98e-02   3.91e-02   -2.04    0.054 .
## std.PCA.wclim.1  -3.82e-02   3.49e-02   -1.10    0.286  
## pcoa.MidD.jac.2   9.57e-02   5.52e-02    1.73    0.099 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.142 on 20 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:   0.6,	Adjusted R-squared:   0.5 
## F-statistic: 6.01 on 5 and 20 DF,  p-value: 0.0015
```

```r

AIC(MC2.1, MC2.2, MC2.3)
```

```
##       df    AIC
## MC2.1  6 -18.90
## MC2.2  7 -21.09
## MC2.3  7 -20.53
```

```r

```




Comparing the models fit
========================

> For this step I calculated the Mean Square Error of each model as described in Gotelli et al. 2009
> The mean square error combines the measurement of precision and bias of the model


### Species Richness

#### Using the raw environmental data as predictors

>Long was excluded for being strongly correlated to Lat


```r
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


#### Using the environmental data as predictors of individual species (in a logistic regression)

>The individual environmental variables were regressed against each species. This generated a probability of a given species to occur in a given site taking into account just one variable. The sum of the probability for all species is the expected species richness in a site. The following code shows how the expected richness given the environmental variables in a logistic regression match the observed richness.

>This code needs to be modified (include logistic prediction as the fitted)


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


### Neutral model


```r

summary(lm(rich ~ rich.neutral.mean, data = environment.AF.2d))
```

```
## 
## Call:
## lm(formula = rich ~ rich.neutral.mean, data = environment.AF.2d)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -9.84  -6.79  -1.04   5.10  13.80 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)
## (Intercept)         -13.30      18.15   -0.73     0.47
## rich.neutral.mean     3.99       2.64    1.51     0.14
## 
## Residual standard error: 7.07 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.087,	Adjusted R-squared:  0.049 
## F-statistic: 2.29 on 1 and 24 DF,  p-value: 0.143
```


### Mid-Domain



```r

summary(step(lm(rich ~ rich.MidD.mean, data = environment.2d)))
```

```
## Start:  AIC=103.9
## rich ~ rich.MidD.mean
## 
##                  Df Sum of Sq  RSS AIC
## <none>                        1213 104
## - rich.MidD.mean  1       100 1314 104
```

```
## 
## Call:
## lm(formula = rich ~ rich.MidD.mean, data = environment.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -11.349  -6.197   0.209   5.255  13.155 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)
## (Intercept)       6.218      5.745    1.08     0.29
## rich.MidD.mean    0.558      0.396    1.41     0.17
## 
## Residual standard error: 7.11 on 24 degrees of freedom
## Multiple R-squared:  0.0765,	Adjusted R-squared:  0.038 
## F-statistic: 1.99 on 1 and 24 DF,  p-value: 0.171
```

```r

```


### Figures comparing the models


```r

# Environment

vars <- c("Lat2", "Long2", "conservation", "vegetation", paste("bio", 1:19, 
    sep = ""), "PCA.wclim.1", "PCA.wclim.2")

par(mgp = c(0.2, 0, 0))
par(mar = c(2, 2, 0, 0))
par(mfrow = c(ceiling(length(vars)/ceiling(sqrt(length(vars)))), ceiling(sqrt(length(vars)))))

# for (i in vars){
# plot(as.formula(paste('rich~',i)),data=environment.AF.2d[!is.na(environment.AF.2d$bio19),],axes=FALSE)
# lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,'rich']),i],environment.AF.2d[!is.na(environment.AF.2d[,'rich']),paste('rich.logis.',i,sep='')],nknots=5),lwd=2,col=4)
# lines(smooth.spline(environment.AF.2d[!is.na(environment.AF.2d[,'rich']),i],environment.AF.2d[!is.na(environment.AF.2d[,'rich']),paste('rich',sep='')],nknots=5),lwd=2,col=2)
# #points(environment.AF.2d[,i],environment.AF.2d[,paste('rich.logis.',i,sep='')],pch=21,cex=.2,bg=4,col=0)
# box() #abline(lm(as.formula(paste('rich~',i)),data=environment.AF.2d)) }


for (i in vars) {
    plot(as.formula(paste("rich~", i)), data = environment.AF.2d[!is.na(environment.AF.2d$bio19), 
        ], axes = FALSE)
    abline(lm(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), paste("rich", 
        sep = "")] ~ environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), 
        i]), lwd = 2, col = 2)
    abline(lm(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), paste("rich.logis.", 
        i, sep = "")] ~ environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), 
        i]), lwd = 2, col = 4)
    points(environment.AF.2d[, i], environment.AF.2d[, paste("rich.logis.", 
        i, sep = "")], pch = 21, bg = 2, col = 0)
    box()
    # abline(lm(as.formula(paste('rich~',i)),data=environment.AF.2d))
}
```

![plot of chunk comparisons of the models](figure/comparisons of the models1.png) 

```r

# Neutral

par(op)
plot(rich ~ rich.neutral.mean, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ], pch = 21, col = NULL, bg = 4)
abline(lm(rich ~ rich.neutral.mean, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ]))
```

![plot of chunk comparisons of the models](figure/comparisons of the models2.png) 

```r

# Mid Domain

par(op)
plot(rich ~ rich.MidD.mean, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ], col = NULL, bg = 2, pch = 21)
abline(lm(rich ~ rich.MidD.mean, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ]))
```

![plot of chunk comparisons of the models](figure/comparisons of the models3.png) 


### Species Composition

#### Using the raw environmental data as predictors

>Long was excluded for being strongly correlated to Lat


```r

vars <- c("Lat2", "Long2", "conservation", "vegetation", paste("bio", 1:19, 
    sep = ""), "PCA.wclim.1", "PCA.wclim.2")

par(mgp = c(0.2, 0, 0))
par(mar = c(2, 2, 0, 0))
par(mfrow = c(ceiling(length(vars)/ceiling(sqrt(length(vars)))), ceiling(sqrt(length(vars)))))

for (i in vars) {
    plot(as.formula(paste("pcoa.jac.1~", i)), data = decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19), 
        ], "range"), axes = FALSE)
    abline(lm(as.formula(paste("pcoa.jac.1~", i)), data = decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19), 
        ], "range")), lwd = 2, col = 4)
    abline(lm(as.formula(paste("pcoa.logis.jac.", i, ".1~", i, sep = "")), data = decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19), 
        ], "range")), lwd = 2, col = 2)
    points(decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19), i], 
        "range"), decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19), 
        paste("pcoa.logis.jac.", i, ".1", sep = "")], "range"), pch = 21, bg = 2, 
        col = 0)
    box()
    # abline(lm(as.formula(paste('rich~',i)),data=environment.AF.2d))
}
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-171.png) 

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

summary(lm(pcoa.jac.1 ~ pcoa.jac.neutral.1, data = environment.AF.2d))
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ pcoa.jac.neutral.1, data = environment.AF.2d)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.1826 -0.1044  0.0059  0.0728  0.3556 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        -2.28e-17   2.63e-02    0.00        1    
## pcoa.jac.neutral.1  5.85e-01   6.85e-02    8.55  9.6e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.134 on 24 degrees of freedom
##   (30 observations deleted due to missingness)
## Multiple R-squared:  0.753,	Adjusted R-squared:  0.742 
## F-statistic:   73 on 1 and 24 DF,  p-value: 9.62e-09
```

```r

par(op)
plot(pcoa.jac.1 ~ pcoa.jac.neutral.1, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ])

abline(lm(environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), "pcoa.jac.1"] ~ 
    environment.AF.2d[!is.na(environment.AF.2d[, "rich"]), "pcoa.jac.neutral.1"]), 
    lwd = 2, col = 4)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-172.png) 

```r

# colnames(environment.AF.2d)

par(op)

################# Mid-Domain

par(op)

summary(lm(pcoa.jac.1 ~ pcoa.MidD.jac.1, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ]))
```

```
## 
## Call:
## lm(formula = pcoa.jac.1 ~ pcoa.MidD.jac.1, data = environment.AF.2d[!is.na(environment.AF.2d[, 
##     "rich"]), ])
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.2096 -0.0752 -0.0287  0.0377  0.3618 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     2.80e-17   2.52e-02    0.00        1    
## pcoa.MidD.jac.1 2.85e-01   3.16e-02    9.01  3.6e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.129 on 24 degrees of freedom
## Multiple R-squared:  0.772,	Adjusted R-squared:  0.762 
## F-statistic: 81.2 on 1 and 24 DF,  p-value: 3.6e-09
```

```r

plot(pcoa.jac.1 ~ pcoa.MidD.jac.1, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ])

abline(lm(pcoa.jac.1 ~ pcoa.MidD.jac.1, data = environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ]), lwd = 2, col = 2)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-173.png) 

```r

```




# Comparison of the observed and predicted by the models


```r

par(mfrow = c(1, 3))

plot(0:1, 0:1, type = "n", xlab = "Predicted", ylab = "Observed")

abline(0, 1, lwd = 2, lty = 2)

title("Richness")

abline(lm(rich ~ rich.neutral.mean, data = decostand(environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ], "range")), lwd = 2, col = 4)

abline(lm(rich ~ rich.MidD.mean, data = decostand(environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ], "range")), lwd = 2, col = 2)

vars <- c("conservation", "PCA.wclim.1", "PCA.wclim.2")

# for (i in vars){

#
# var.model<-lm(as.formula(paste('rich~',i)),data=decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19),],'range'))

#
# abline(lm(rich~predict(var.model),data=decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19),],'range')))

#
# abline(lm(as.formula(paste('rich~',paste('rich.logis.',i,sep=''))),data=decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19),],'range')),lwd=2,col=6)
# }


# Species composition (PCA)

plot(0:1, 0:1, type = "n", xlab = "Predicted", ylab = "Observed")
abline(0, 1, lwd = 2, lty = 2)

title("Species composition (PCA)")

abline(lm(pcoa.jac.1 ~ pcoa.jac.neutral.1, data = decostand(environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ], "range")), col = 4)

#
# points(rich~rich.MidD.mean,data=decostand(environment.AF.2d[!is.na(environment.AF.2d[,'rich']),],'range'),pch=21,bg=2)
abline(lm(pcoa.jac.1 ~ pcoa.MidD.jac.1, data = decostand(environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ], "range")), col = 2)

#
# vars<-c('Lat2','Long2','conservation','vegetation',paste('bio',1:19,sep=''),'PCA.wclim.1','PCA.wclim.2')

# vars<-c('vegetation','PCA.wclim.1','PCA.wclim.2')

# for (i in vars){
# abline(lm(as.formula(paste('pcoa.jac.1~',i)),data=decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19),],'range')),lwd=2,col=3)
# abline(lm(as.formula(paste('pcoa.jac.1~',paste('pcoa.logis.jac.',i,'.1',sep=''))),data=decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19),],'range')),lwd=2,col=6)
# }


# Species composition (Jaccard index )

plot(0:1, 0:1, type = "n", xlab = "Predicted", ylab = "Observed")

abline(0, 1, lwd = 2, lty = 2)

title("Species composition (Jaccard)")

#
# points(1-as.dist(morisita.AF.2d[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d]),similarity.jac,col=4)

# library(vegan)

points(1 - decostand(as.vector(as.dist(morisita.AF.2d[environment.2d$matchin.AF.2d, 
    environment.2d$matchin.AF.2d])), "range"), decostand(as.vector(similarity.jac), 
    "range"), col = 4)

# points(similarity.jac~as.dist(similarity.MidD.mean),col=2)

points(decostand(as.vector(as.dist(similarity.MidD.mean)), "range"), decostand(as.vector(similarity.jac), 
    "range"), col = 2)

#
# summary(lm(similarity.jac~1-as.dist(morisita.AF.2d[environment.2d$matchin.AF.2d,environment.2d$matchin.AF.2d])))

# summary(lm(similarity.jac~as.dist(similarity.MidD.mean)))

# points(as.dist(similarity.neutral.jac.mean),similarity.jac,col=4)

# points(as.dist(similarity.MidD.mean),similarity.jac,col=2)

# abline(lm(similarity.jac~as.dist(1-similarity.neutral.jac.mean)),col=4)

abline(lm(decostand(as.vector(similarity.jac), "range") ~ decostand(-as.vector(as.dist(morisita.AF.2d[environment.2d$matchin.AF.2d, 
    environment.2d$matchin.AF.2d])), "range")), col = 4)


#
# points(rich~rich.MidD.mean,data=decostand(environment.AF.2d[!is.na(environment.AF.2d[,'rich']),],'range'),pch=21,bg=2)
abline(lm(pcoa.jac.1 ~ pcoa.MidD.jac.1, data = decostand(environment.AF.2d[!is.na(environment.AF.2d[, 
    "rich"]), ], "range")), col = 2)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 

```r

# points(similarity.jac~as.dist(similarity.MidD.mean),col=2)

#
# points(decostand(as.vector(as.dist(similarity.MidD.mean)),'range'),decostand(as.vector(similarity.jac),'range'),col=2)


#
# vars<-c('Lat2','Long2','conservation','vegetation',paste('bio',1:19,sep=''),'PCA.wclim.1','PCA.wclim.2')

# vars<-c('vegetation','PCA.wclim.1','PCA.wclim.2')

# for (i in vars){
# abline(lm(as.formula(paste('pcoa.jac.1~',i)),data=decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19),],'range')),lwd=2,col=3)
# abline(lm(as.formula(paste('pcoa.jac.1~',paste('pcoa.logis.jac.',i,'.1',sep=''))),data=decostand(environment.AF.2d[!is.na(environment.AF.2d$bio19),],'range')),lwd=2,col=6)
# }
```



Compare the models using the statistics in Gotelli et al. (2009)

Mean Square Error (MSE = Bias + Variance)

Richness


```r

O <- decostand(environment.2d[!is.na(environment.2d[, "rich"]), "rich"], "standardize")

S <- decostand(rich.MidD.vec, "standardize")
E <- rowMeans(S)

s.bias.sq.MidD <- sum((O - E)^2)
s.var.MidD <- (sum(colSums((S - E)^2)))/(ncol(S) - 1)
s.mse.MidD = s.bias.sq.MidD + s.var.MidD

S <- decostand(rich.neutral[!is.na(environment.AF.2d[, "rich"]), ], "standardize")
E <- rowMeans(S)

s.bias.sq.neutral <- sum((O - E)^2)
s.var.neutral <- (sum(colSums((S - E)^2)))/(ncol(S) - 1)
s.mse.neutral = s.bias.sq.neutral + s.var.neutral

s.bias.sq.logis <- {
}
s.var.logis <- {
}
s.mse.logis = {
}

for (i in glm.variables.2d) {
    
    S <- decostand(do.call(cbind, lapply(species.2d.logis.ls[[i]], rowSums)), 
        "standardize")
    
    E <- rowMeans(S)
    
    s.bias.sq.logis[i] <- sum((O - E)^2)
    s.var.logis[i] <- (sum(colSums((S - E)^2)))/(ncol(S) - 1)
    s.mse.logis[i] = s.bias.sq.logis[i] + s.var.logis[i]
    
}

data.frame(BIASsq = c(MidD = s.bias.sq.MidD, Neutral = s.bias.sq.neutral, s.bias.sq.logis), 
    VAR = c(MidD = s.var.MidD, Neutral = s.var.neutral, s.var.logis), sMSE = c(MidD = s.mse.MidD, 
        Neutral = s.mse.neutral, s.mse.logis))
```

```
##              BIASsq    VAR  sMSE
## MidD          29.91  9.103 39.01
## Neutral       24.82 23.164 47.98
## Long2         31.17 16.666 47.84
## Lat2          34.04 19.578 53.62
## conservation  21.16 11.691 32.85
## vegetation    22.17 13.978 36.15
## bio1          31.30 17.171 48.47
## bio2          28.40 22.385 50.79
## bio3          35.17 17.210 52.38
## bio4          32.73 15.561 48.29
## bio5          24.22 23.022 47.24
## bio6          34.61 20.225 54.84
## bio7          29.31 21.380 50.69
## bio8          29.77 13.779 43.55
## bio9          24.69 15.609 40.30
## bio10         25.17 21.752 46.92
## bio11         34.08 16.550 50.63
## bio12         29.39 20.258 49.65
## bio13         32.15 21.324 53.47
## bio14         28.88 17.508 46.39
## bio15         31.98 17.758 49.74
## bio16         32.20 22.070 54.27
## bio17         29.06 17.940 47.00
## bio18         30.44 22.901 53.34
## bio19         25.95 17.535 43.49
## PCA.wclim.1   31.30 16.990 48.29
## PCA.wclim.2   31.30 20.625 51.93
```


Similarity


```r

O <- decostand(as.vector(as.matrix(similarity.jac)), "standardize")

S <- decostand(similarity.MidD.vec, "standardize")
E <- rowMeans(S)

s.bias.sq.MidD.jac <- sum((O - E)^2)
s.var.MidD.jac <- (sum(colSums((S - E)^2)))/(ncol(S) - 1)
s.mse.MidD.jac = s.bias.sq.MidD.jac + s.var.MidD.jac

S <- decostand(apply(similarity.neutral.jac, 3, as.vector), "standardize")
E <- rowMeans(S)


s.bias.sq.neu.jac <- sum((O - E)^2)
s.var.neu.jac <- (sum(colSums((S - E)^2)))/(ncol(S) - 1)
s.mse.neu.jac <- s.bias.sq.neu.jac + s.var.neu.jac

s.bias.sq.logis.jac <- {
}
s.var.logis.jac <- {
}
s.mse.logis.jac <- {
}

for (i in glm.variables.2d) {
    
    S <- decostand(apply(similarity.logis.jac.ls[[i]], 3, as.vector), "standardize")
    E <- rowMeans(S)
    
    s.bias.sq.logis.jac[i] <- sum((O - E)^2)
    s.var.logis.jac[i] <- (sum(colSums((S - E)^2)))/(ncol(S) - 1)
    s.mse.logis.jac[i] = s.bias.sq.logis.jac[i] + s.var.logis.jac[i]
    
}

data.frame(BIASsq = c(MidD = s.bias.sq.MidD.jac, Neutral = s.bias.sq.neu.jac, 
    s.bias.sq.logis.jac), VAR = c(MidD = s.var.MidD.jac, Neutral = s.var.neu.jac, 
    s.var.logis.jac), sMSE = c(MidD = s.mse.MidD.jac, Neutral = s.mse.neu.jac, 
    s.mse.logis.jac))
```

```
##              BIASsq    VAR  sMSE
## MidD          353.1  65.33 418.5
## Neutral       274.4  78.02 352.4
## Long2         222.6 142.56 365.1
## Lat2          212.2 148.17 360.4
## conservation  316.4 160.76 477.2
## vegetation    286.0 161.73 447.7
## bio1          252.9 152.43 405.3
## bio2          312.6 161.24 473.8
## bio3          245.9 151.51 397.4
## bio4          221.4 143.41 364.8
## bio5          334.2 163.56 497.7
## bio6          261.1 152.18 413.3
## bio7          270.4 151.89 422.3
## bio8          246.7 142.56 389.3
## bio9          312.1 157.34 469.4
## bio10         316.8 160.20 477.0
## bio11         225.6 148.63 374.2
## bio12         307.1 159.97 467.1
## bio13         292.7 161.21 453.9
## bio14         276.6 153.63 430.2
## bio15         263.2 155.68 418.9
## bio16         293.3 161.19 454.5
## bio17         282.0 154.09 436.1
## bio18         310.5 161.61 472.2
## bio19         303.4 159.37 462.8
## PCA.wclim.1   240.4 147.31 387.7
## PCA.wclim.2   306.9 161.27 468.1
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
##   |                                                                         |                                                                 |   0%  |                                                                         |.                                                                |   1%  |                                                                         |.                                                                |   2%  |                                                                         |..                                                               |   3%  |                                                                         |...                                                              |   5%  |                                                                         |....                                                             |   6%  |                                                                         |....                                                             |   7%  |                                                                         |.....                                                            |   8%  |                                                                         |......                                                           |   9%  |                                                                         |.......                                                          |  10%  |                                                                         |.......                                                          |  11%  |                                                                         |........                                                         |  13%  |                                                                         |.........                                                        |  14%  |                                                                         |..........                                                       |  15%  |                                                                         |..........                                                       |  16%  |                                                                         |...........                                                      |  17%  |                                                                         |............                                                     |  18%  |                                                                         |.............                                                    |  20%  |                                                                         |.............                                                    |  21%  |                                                                         |..............                                                   |  22%  |                                                                         |...............                                                  |  23%  |                                                                         |................                                                 |  24%  |                                                                         |................                                                 |  25%  |                                                                         |.................                                                |  26%  |                                                                         |..................                                               |  28%  |                                                                         |...................                                              |  29%  |                                                                         |...................                                              |  30%  |                                                                         |....................                                             |  31%  |                                                                         |.....................                                            |  32%  |                                                                         |......................                                           |  33%  |                                                                         |......................                                           |  34%  |                                                                         |.......................                                          |  36%  |                                                                         |........................                                         |  37%  |                                                                         |.........................                                        |  38%  |                                                                         |.........................                                        |  39%  |                                                                         |..........................                                       |  40%  |                                                                         |...........................                                      |  41%  |                                                                         |............................                                     |  43%  |                                                                         |............................                                     |  44%  |                                                                         |.............................                                    |  45%  |                                                                         |..............................                                   |  46%  |                                                                         |...............................                                  |  47%  |                                                                         |...............................                                  |  48%  |                                                                         |................................                                 |  49%  |                                                                         |.................................                                |  51%  |                                                                         |..................................                               |  52%  |                                                                         |..................................                               |  53%  |                                                                         |...................................                              |  54%  |                                                                         |....................................                             |  55%  |                                                                         |.....................................                            |  56%  |                                                                         |.....................................                            |  57%  |                                                                         |......................................                           |  59%  |                                                                         |.......................................                          |  60%  |                                                                         |........................................                         |  61%  |                                                                         |........................................                         |  62%  |                                                                         |.........................................                        |  63%  |                                                                         |..........................................                       |  64%  |                                                                         |...........................................                      |  66%  |                                                                         |...........................................                      |  67%  |                                                                         |............................................                     |  68%  |                                                                         |.............................................                    |  69%  |                                                                         |..............................................                   |  70%  |                                                                         |..............................................                   |  71%  |                                                                         |...............................................                  |  72%  |                                                                         |................................................                 |  74%  |                                                                         |.................................................                |  75%  |                                                                         |.................................................                |  76%  |                                                                         |..................................................               |  77%  |                                                                         |...................................................              |  78%  |                                                                         |....................................................             |  79%  |                                                                         |....................................................             |  80%  |                                                                         |.....................................................            |  82%  |                                                                         |......................................................           |  83%  |                                                                         |.......................................................          |  84%  |                                                                         |.......................................................          |  85%  |                                                                         |........................................................         |  86%  |                                                                         |.........................................................        |  87%  |                                                                         |..........................................................       |  89%  |                                                                         |..........................................................       |  90%  |                                                                         |...........................................................      |  91%  |                                                                         |............................................................     |  92%  |                                                                         |.............................................................    |  93%  |                                                                         |.............................................................    |  94%  |                                                                         |..............................................................   |  95%  |                                                                         |...............................................................  |  97%  |                                                                         |................................................................ |  98%  |                                                                         |................................................................ |  99%  |                                                                         |.................................................................| 100%
```

```
## output file: All_analysis.R
```

```
## [1] "All_analysis.R"
```





