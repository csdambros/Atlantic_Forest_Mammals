Comparing dispersal, environmental and mid-domain effects on species distribution
============================



##### By CSDambros

>html updated at 2014-04-16 10:38:38

### Methods

#### Study site
 
  The study was conducted in the whole area recognized as the Atlantic forest biome (Fig. 1). This biome encompasses an extent of 102,012 km², but today representing only 7.9% of its original forest cover. It has several vegetation types such as rainforests, mixed (Araucarian) moist forests, semideciduous forests, dry forests and upland grasslands. Rainforests tend to occur near the coast and semideciduous and dry forests far from coast; mixed forests are common in the south of AF (SOS Mata Atlântica 2007). AF shows moist tropical and subtropical climates, without well-defined dry season, and annual mean temperatures above 15ºC (Leite, 2002). 








#### Data collection

  We compiled a database based on the literature in which diverse authors have sampled small mammals along the AF (Appendix S1 in Supporting Information). We have used Google @ search tool for sampling articles, and the main keywords used in combination were “small mammal”, “marsupial”, “rodent”, “community”, “composition”, “richness”, “diversity”, and “Atlantic Forest”. Unpublished data sampled by NCC was also included (Appendix S1). There was some variation in the area covered by each surveyed study (see Appendix S1), but we have established a minimum of sampling effort for a given article to be accepted in our database: at least 1000 trap-nights, 6 months of field work, and use of wire type and/or Sherman live-traps installed on the ground or understory level of the forest. Besides that, and for independence assumptions, we have chosen surveyed areas distant at least 10 km from each other; most of the field surveys available in the literature sampled only one geographical location under this condition. From each selected survey, we obtained local species composition and abundance data when available. 

  Data on conservation were taken directly from information available on the articles consulted, and classified as categories of forest conservation (degrees 1 to 5), where primary conserved forest means 5 and disturbed secondary forest, including with clearings, means 1. Physiognomy data were taken either from information available on the area descriptions and from general sources such as (Veloso et al., 1991). We compiled the 19 environmental variables available in Bioclim (www.worldclim.org/bioclim) in a scale of 2.5 arc minutes: annual mean temperature (1), mean diurnal range (2), isothermality (3), temperature seasonality (4), maximum and minimum temperature of the warmest and coldest months (5 and 6), temperature annual range (7), mean temperature of the wettest, driest, warmest and coldest quarters (8-11), annual precipitation (12), precipitation of the wettest and driest months (13 and 14), precipitation seasonality (15), and precipitation of the wettest, direst, warmest and coldest quarters (16-19).  








#### Analysis

  The environmental predictor variables were used as predictors of the species diversity and included in simulation models (see below). Most of the 19 WorldClim variables are correlated to each other, so we grouped these variables in two Principal Component Axes (PCA1 and PCA2) (But see Supplementary material S4 for individual comparisons). We standardized all the variables previously to the analyses. We preferred to use the described environmental variables as predictors and excluded the highly correlated variables of latitude and longitude. We used them to calculate the pairwise distance among sampling units and the linear distance of each sampling site to the coast.









  We analyzed the data using the original local sampling units and grouping the original data points in 26 2x2 degrees quadrants (Fig. 1). Most localities were surveyed by a single study, then most of the sampling points represent a single study (see data collection). We used the number of species encountered in each locality (S), and the similarity in species composition as response variables.
  We assessed the similarity in species composition by calculating the Jaccard similarity index among all pairs of samples. The matrix of similarities was summarized using a Principal Coordinates Analysis (PCoA). We used only the first two ordination axes because all axes in higher dimensions captured alone less than 10% of the variation from the similarity matrix.








![plot of chunk Fig. 1 Study site](figure/Fig. 1 Study site.png) 

#### Fig. 1. Map of the Atlantic Forest (AF; blue shade) showing the orinal sampling points (red circles) and the quadrants encompassing the entire AF. The size of the circles represent the size o the original sampling area in the log scale. Green shade represent the quadrants where small mammal data was available. 


  For estimating the influence of dispersal limitation on species local diversity (S) and turnover (Jaccard index; PCoA1 and PCoA2), we simulated the AF as a network of interconnected quadrants where the species or individuals could move between adjacent quadrants (Fig. 2). We used two models to recreate the species distribution under dispersal alone: the Mid-Domain spreading-dye model (Colwell and Hurtt 1994; Colwell and Lees 2000), and the analytical neutral approach borrowed from the population genetics (Nagylaki, 1980) and presented by Economo & Keitt (2008) for community ecology. These models were created just using the grouped data in quadrants because of the completeness of data, and due to computational limitation.
  In the Mid-Domain spreading dye model, we recorded the number of quadrants occupied by each small mammal species. For each species, one of the 26 quadrants was randomly selected and the species occurrence was spread from the selected quadrant to neighboring quadrants until the original number of quadrants was occupied. This procedure was repeated 10,000 times for the 64 species. Each quadrant had up to eight neighbors (Moore neighborhood), and the model was bounded by the domain where actual small mammal data was recorded (26 quadrants).








  In the neutral model, the whole area comprising the AF was divided in 56 2x2 degrees quadrants from which 26 had small mammal data available (Fig. 1). The model started with a single species occupying all the 56 quadrants. In each generation, new species were added in each quadrant by point speciation with rate v, set the same for all quadrants (see Economo & Keitt (2010) and Muneepeerakul et al. (2008) for more details and other uses). v represents the probability of an individual to become a new species but could also represent the addition of new species by immigration from a larger species pool outside the AF (eg. the Cerrado or Amazonian forests). To recreate the dispersal of individuals, we determined that a quadrant could just be colonized by a neighbor (Moore neighborhood), and that all quadrants had the same migration rate (parameter m). The local community size (number of individuals) was set the same for all quadrants (N = 100). The model was run for multiple generations, until the diversity within (PIE and HillPIE) and among quadrants (Morisita-Horn similarity) reached a steady-state (usually more than 30,000 generations).


![plot of chunk Fig. 2 - Connectivity map 26 and 56](figure/Fig. 2 - Connectivity map 26 and 56.png) 

#### Fig. 2. Map representing the connectivity of quadrants used to simulate the individual dispersal in the Mid Domain (A) and neutral (B) models.


  Differently from the Mid-domain model, the values of m, v, and N were initially set independently of the observed data. Because different combinations of these parameters can create the same patterns of species distribution (Etienne 2006?; Economo and Keitt 2010), and because estimates for these parameters based on empirical data are not available, we did not attempt to estimate realistic parameters, but focused on the general patterns that the model could produce (see discussion for implications of this procedure). N was arbitrarily defined as 100 individuals for all quadrants, and m and v were used as knobs to best fit the model to the observed similarity in species composition. These parameters were simultaneously tunned using the L-BFSG-S optimization algorithm. In summary, m and v were adjusted to minimize the differences between the predictions of the neutral model to the observed data, giving the best possible explanation of dispersal to the observed data (this procedure is conceptually analog to a regression analysis using maximum likelihood optimization). Because the number of individuals per quadrant was kept constant during the optimization, the differences in diversity were due to the location of a quadrant in the network and the number of neighbor quadrants, but not affected by the differences in abundance.











  The model proposed by Economo and Keitt (2008) is probabilistic and does not require the simulation of each individual in the metacommunity, being less computationally intensive. This model allowed us to investigate thoroughly the parameter space of m and v to find the best explanatory model for the observed data. However, because the model is based on probabilities, it does not allow one to calculate statistics based on the raw data (eg. the Jaccard similarity index). Because not all original studies recorded species abundances, we compared the Morisita-Horn from the neutral model and the Jaccard indexis from the data for optimization. These indexes are usulally highly correlated and we believe this procedure did not have profound effects on the analysis (CITATION - CHAO?). 
  To calculate the Jaccard similarity index from the neutral model, We  simulated each individual in the metacommunity using the previously optimized parameters. We run an initial buffer of 30,000 generations. Then the model was run for additional 1000 generations 5,000 times, representing 5,000 simulations. The mean of all the 5,000 simulations was used to calculate the remain summary statistics (see below). The correlation between the probabilistic and the simulated model was 99.98% when comparing the Morisita-Horn index of pairwise similarity, so we believe the simulations run long enough to capture the final predictions of the neutral model.







To test the association of the species diversity with the environment, we fit individual logistic regressions for all species against the predictor variables, allowing each species to vary independently of one another.




We then summed the probability of occurrence generated by the logistic model for each species to estimate the expected number of species in each quadrant:
  
$$
\begin{equation}
\begin{split}
\hat S_i = \sum_{j=1}^{S}{\frac{e^{X_i\beta_j}}{1+e^{X_i\beta_j}}}
\end{split}
\end{equation}
$$
  
Where $X_i$ denotes a vector of the environmental variables in the quadrant $i$, and $\beta_j$ is a vector with the coefficients from the logistic regression for the species $j$.




To calculate the raw statistics based on these probabilities, we simulated the distribution of each species based on the probabilities of occurrence generated by the logistic curve, similarly to what was done for the neutral model. For each species, we counted the number of quadrants where the species was observed. Thus we spread the species into the map based on its probability of occurrence until the number of quadrants occupied in the simulation matched the number observed. This allowed each species to vary independently, and to the relationship between the environmental variables and the response variables (eg. number of species) to be non-linear. Differently from the non-linearities that could be used in other statistical tests, such as Generalized additive models, the deterministic part of this model has an underlying mechanistic function (logistic).








We compared the three simulation models (Mid Domain, Neutral, and based on logistic regressions) by their Mean Square Error (MSE; Gotelli et al. 2009). The MSE was calculated as the sum of the squared bias and model variance, so the MSE is a balance between precision and accuracy (see Gotelli et al. 2009 and Appendix S3 for detailed description).







Additionally, we tested for the direct association of the response variables against the environmental predictors. For this test, we used a linear regression model (OLS) starting with all environmental predictors in the model (conservation, vegetation, PCAbclim1, PCAbclim2). For the grouped data, we created three sets of models: (1) Using the neutral and Mid-domain predictions as a predictor variables along with the environmental variables in the model, as suggested by Letten et al. (2013); (2) using the neutral and Mid-domain predictions as a null hypothesis and testing for association of residuals with the environmental variables; and (3) a regression with only the environmental variables as predictors. We used both the raw environmental data and the prediction from the logistic simulation as predictor variables. This gave a total of 10 models (Table S1;S3). The models were compared by their AIC values.
  



  Finally, we tested for patterns of distance-decay in the species similarity using Mantel and Partial Mantel tests correlating the matrix of Jaccard similarity (M) to the geographical distance (D) and environmental dissimilarity (E) matrices (Euclidean distance) (see Thompson and Townsend 2006). Legendre et al. (2005) and Tuomisto and Ruokolainen (2006; 2008) suggest using multiple regression on similarity matrices to separate the effects of niche and neutral processes. However, there is a strong debate about the validity of these models (Legendre et al. 2008; Tuomisto and Ruokolainen 2008). Some authors also suggest using partial Redundancy Analysis (Borcard et al. 1992) for this purpose (Gilbert and Lechowicz 2004), but this method also has serious limitations (Smith and Lundholm 2010). We preferred to base our discussion on the mantel tests and multiple regressions using the summarized data (each quadrant as a unit), as described above. 









  All the analyses were conducted in the R program (R Development Core Team, 2013). The models and most of the summary statistics calculations were implemented by the authors, and are available as a supplementary material (see Appendix S3). We used the packages vegan (Oksanen et al., 2008) for the remaining analyses.

  	
### Results

  All the models had a very poor fit to species richness (Fig. S1-S2). The model with the lowest mean square error and bias was the logistic using only conservation status as a predictor followed by the logistic model based on vegetation type (Table 1). The model with the lowest variance was the neutral model, but the difference among the models was much smaller than for the bias (Table 1).

#### Table 1. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for species richness. The Mean Square Error is the sum of the Bias and Variance.


<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Apr 16 10:48:11 2014 -->
<TABLE border=0, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 1278.09 </TD> <TD align="right"> 196.17 </TD> <TD align="right"> 1474.26 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 1518.38 </TD> <TD align="right"> 145.41 </TD> <TD align="right"> 1663.79 </TD> </TR>
  <TR> <TD align="right"> conservation </TD> <TD align="right"> 1038.67 </TD> <TD align="right"> 202.79 </TD> <TD align="right"> 1241.46 </TD> </TR>
  <TR> <TD align="right"> vegetation </TD> <TD align="right"> 1100.87 </TD> <TD align="right"> 199.00 </TD> <TD align="right"> 1299.87 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 1362.20 </TD> <TD align="right"> 181.26 </TD> <TD align="right"> 1543.46 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 1413.74 </TD> <TD align="right"> 199.67 </TD> <TD align="right"> 1613.41 </TD> </TR>
   </TABLE>


  The regression models gave very similar results. The model with the lowest AIC was the linear regression of the species richness against the raw environmental predictors (Distance to the coast, vegetation, conservation and PCA axes of edaphic conditions) (Table 2). However, the inclusion of neutral and Mid-Domain predictions along the environmental variables gave similar results (Table 2). The environmental variables were able to explain 33% of the variation on species richness, and conservation status had the strongest effect (Fig. SXX).
  Besides the poor fit for the data, the Mid Domain and Neutral models created the well known richness peak in the central areas of the Atlantic forest (Fig. S1).
  
#### Table 2. AIC values comparing the 10 regression models tested. MR1.1: Environmental predictors; MR1.2: Environmental + Neutral predictions; MR1.3: ; MR1.4: ; MR1.5: ; MR1.6: ; MR1.7: ; MR1.8: ; MR1.9: ; MR1.10: .

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Apr 16 10:48:11 2014 -->
<TABLE border=0, bgcolor=#989898>
<TR> <TH>  </TH> <TH> df </TH> <TH> AIC </TH>  </TR>
  <TR> <TD align="right"> MR1.1 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 69.84 </TD> </TR>
  <TR> <TD align="right"> MR1.2 </TD> <TD align="right"> 7.00 </TD> <TD align="right"> 71.78 </TD> </TR>
  <TR> <TD align="right"> MR1.3 </TD> <TD align="right"> 7.00 </TD> <TD align="right"> 70.64 </TD> </TR>
  <TR> <TD align="right"> MR1.4 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 76.32 </TD> </TR>
  <TR> <TD align="right"> MR1.5 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 76.64 </TD> </TR>
  <TR> <TD align="right"> MR1.6 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 75.45 </TD> </TR>
  <TR> <TD align="right"> MR1.7 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 75.03 </TD> </TR>
  <TR> <TD align="right"> MR1.8 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 75.96 </TD> </TR>
  <TR> <TD align="right"> MR1.9 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 76.32 </TD> </TR>
  <TR> <TD align="right"> MR1.10 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 76.64 </TD> </TR>
   </TABLE>


  The similarity in species composition was well fit for either the neutral, the Mid Domain, and the environmental models (Fig. S3-S4). The mean square error in the models and the difference between the models was much smaller for the species composition than for the richness, besides the much larger number of comparisons for the species composition (325 for the pairwise similarity matrix against 26 for the vector of species richness). The lowest mean square error and bias were obtained using the edaphic conditions as predicctor, while the Mid Domain model had the lowest variance. The neutral model had the highest mean square error due to the largest variance compared to the remain models (Table 3).
  
  
#### Table 3. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for the jaccard pairwise similarity. The Mean Square Error is the sum of the Bias and Variance.

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Apr 16 10:48:11 2014 -->
<TABLE border=0, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 8.80 </TD> <TD align="right"> 1.88 </TD> <TD align="right"> 10.68 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 7.07 </TD> <TD align="right"> 9.43 </TD> <TD align="right"> 16.51 </TD> </TR>
  <TR> <TD align="right"> conservation </TD> <TD align="right"> 5.99 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 8.25 </TD> </TR>
  <TR> <TD align="right"> vegetation </TD> <TD align="right"> 5.52 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 7.77 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 4.65 </TD> <TD align="right"> 2.18 </TD> <TD align="right"> 6.84 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 5.76 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 8.00 </TD> </TR>
   </TABLE>

 
 Converselly, three of the four the best regression models for the first ordination axis were obtained including the neutral results as a predictor variable (Table 4). The fourth best model was obtained using just the logistic model as a predictor. For the second axis, the best models included the environmental predictors alone, or in combination with the Neutral or Mid Domain results. The neutral and Mid Domain models alone were able to explain 78% and 77% of the variation in the first axis of species composition. The environmental variables combined explained 77% of the variation is the first axis. The combination of environmental and neutral predictors gave an explanation of 81% for the first ordination axis. The environmental variables combined with the neutral model explained 49% of the variation in the second ordination axis. The use of the Mid Domain instead of the neutral predictor gave similar results (R2 = 0.48; Table 4).
  The Mid Domain and Neutral models were able to explain a non-linear change in species composition along the latitudinal gradient (Fig. S3;S5).  
  
#### Table 4. AIC values comparing the 10 regression models tested. MR1.1: Environmental predictors; MR1.2: Environmental + Neutral predictions; MR1.3: ; MR1.4: ; MR1.5: ; MR1.6: ; MR1.7: ; MR1.8: ; MR1.9: ; MR1.10: .

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Apr 16 10:48:11 2014 -->
<TABLE border=0, bgcolor=#989898>
<TR> <TH>  </TH> <TH> df </TH> <TH> AIC </TH>  </TR>
  <TR> <TD align="right"> MC1.1 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 42.92 </TD> </TR>
  <TR> <TD align="right"> MC1.2 </TD> <TD align="right"> 7.00 </TD> <TD align="right"> 37.40 </TD> </TR>
  <TR> <TD align="right"> MC1.3 </TD> <TD align="right"> 7.00 </TD> <TD align="right"> 40.31 </TD> </TR>
  <TR> <TD align="right"> MC1.4 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 39.28 </TD> </TR>
  <TR> <TD align="right"> MC1.5 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 40.48 </TD> </TR>
  <TR> <TD align="right"> MC1.6 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 43.45 </TD> </TR>
  <TR> <TD align="right"> MC1.7 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 45.35 </TD> </TR>
  <TR> <TD align="right"> MC1.8 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 39.12 </TD> </TR>
  <TR> <TD align="right"> MC1.9 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 39.28 </TD> </TR>
  <TR> <TD align="right"> MC1.10 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 40.48 </TD> </TR>
  <TR> <TD align="right"> MC2.1 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 64.97 </TD> </TR>
  <TR> <TD align="right"> MC2.2 </TD> <TD align="right"> 7.00 </TD> <TD align="right"> 65.43 </TD> </TR>
  <TR> <TD align="right"> MC2.3 </TD> <TD align="right"> 7.00 </TD> <TD align="right"> 66.46 </TD> </TR>
  <TR> <TD align="right"> MC2.4 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 75.94 </TD> </TR>
  <TR> <TD align="right"> MC2.5 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 74.98 </TD> </TR>
  <TR> <TD align="right"> MC2.6 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 79.45 </TD> </TR>
  <TR> <TD align="right"> MC2.7 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 78.39 </TD> </TR>
  <TR> <TD align="right"> MC2.8 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 76.15 </TD> </TR>
  <TR> <TD align="right"> MC2.9 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 75.94 </TD> </TR>
  <TR> <TD align="right"> MC2.10 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 74.98 </TD> </TR>
   </TABLE>


The mantel test results show a stronger environmental effect when grouping the points in quadrants than when using the original data points (Table 5). The geographical distance was highly correlated with the species composition in both cases. Geographical distance was correlated with the species composition when the environment was partialed out, but the correlation was much weaker than in the unconstrained model.

#### Table 5. Mantel and partial Mantel test results comparing the correlation of species similarity against geographical distance and environmental dissimilarity. The minus sign indicate partial results (Geo - Env: Geographical distance without the effect of environmenal dissimilarity).

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Apr 16 10:48:11 2014 -->
<TABLE border=0, bgcolor=#989898>
<TR> <TH>  </TH> <TH> statistic </TH> <TH> signif </TH>  </TR>
  <TR> <TD align="right"> Local Geo </TD> <TD align="right"> 0.48 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> Local Env </TD> <TD align="right"> 0.32 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 2d Geo </TD> <TD align="right"> 0.49 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 2d Env </TD> <TD align="right"> 0.48 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> Local Geo-Env </TD> <TD align="right"> 0.38 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> Local Env-Geo </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.17 </TD> </TR>
  <TR> <TD align="right"> 2d Geo-Env </TD> <TD align="right"> 0.21 </TD> <TD align="right"> 0.01 </TD> </TR>
  <TR> <TD align="right"> 2d Env-Geo </TD> <TD align="right"> 0.18 </TD> <TD align="right"> 0.02 </TD> </TR>
   </TABLE>


Both neutal and Mid Domain models had very similar predictions for the distance-decay pattern in species similarity (Fig. 5)


![plot of chunk Fig 3. distance decay real.midDomain and neutral mor](figure/Fig 3. distance decay real.midDomain and neutral mor.png) 

#### Fig. 3. Distance-decay in species similarity using the Jaccard similarity index. The Mid Domain (red) and Neutral model (blue) have similar fit for the data, but the variance was much smaller than observed (grey).

### Discussion



----------------------------

![plot of chunk Fig S1. Map Richness Mid Domain vs Neutral](figure/Fig S1. Map Richness Mid Domain vs Neutral.png) 

#### Fig S1. Projection of the species richness predicted by dispersal models on the map and comparison with the observed species distribution. Warmer colors represent areas with hiher species richness. (A) Observed species richness; (B) Predicted by the Mid Domain model; (C) Predicted by the neutral model. Shadded quadrants represent areas included in the neutral model but where actual data is not available.



```
## The following objects are masked from mammal.data:
## 
##     Lat2, Long2
```

![plot of chunk Fig S2. rich.logis in the map](figure/Fig S2. rich.logis in the map.png) 

#### Fig S2. Projection of the species richness predicted by logistic models on the map and comparison with the observed species distribution. Warmer colors represent areas with hiher species richness. (A) Observed species richness; (B-Z) Predicted by logistic regressions of individual species against environmental gradients. 


![plot of chunk Fig S3. Map Composition Mid Domain vs Neutral](figure/Fig S3. Map Composition Mid Domain vs Neutral.png) 

#### Fig S3. Projection of the species composition predicted by dispersal models on the map and comparison with the observed species distribution. The colors represent the similarity in species composition measured by the pairwise jaccard similarity index between sites summirized in one axis of a principal coordinates analysis (pcoa). Those quadrants with similar colors have a similar composition of species. (A) Observed species composition; (B) Species composition predicted by the Mid Domain model; (\C) Species composition predicted by the neutral model. Shadded quadrants represent areas included in the neutral model but where actual data is not available.


![plot of chunk Fig S4. jac.logis in the map](figure/Fig S4. jac.logis in the map.png) 

#### Fig S4. Projection of the species composition predicted by dispersal models on the map and comparison with the observed species distribution. The colors represent the similarity in species composition measured by the pairwise jaccard similarity index between sites summirized in one axis of a principal coordinates analysis (pcoa). Those quadrants with similar colors have a similar composition of species. (A) Observed species composition; (B-Z) Species compostition predicted by logistic regressions of individual species against the environmental gradients.

![plot of chunk Fig S5. plot all predictions against lat](figure/Fig S5. plot all predictions against lat.png) 

#### Fig S5. Comparison of the change in species composition along the latitudinal gradient predicted by the Neutral and Mid Domain models. The composition was measured by the Jaccard similarity index between all pairs of sites and summarized by the first axis of a Principal Coordinates Analysis (PCoA).


#### Table S1. List of the sites, authors, etc. used in this manuscript

#### Table S2. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for species richness. The Mean Square Error is the sum of the Bias and Variance.

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Apr 16 10:48:16 2014 -->
<TABLE border=0, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 1278.09 </TD> <TD align="right"> 196.17 </TD> <TD align="right"> 1474.26 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 1518.38 </TD> <TD align="right"> 145.41 </TD> <TD align="right"> 1663.79 </TD> </TR>
  <TR> <TD align="right"> Long2 </TD> <TD align="right"> 1357.45 </TD> <TD align="right"> 179.67 </TD> <TD align="right"> 1537.12 </TD> </TR>
  <TR> <TD align="right"> Lat2 </TD> <TD align="right"> 1459.36 </TD> <TD align="right"> 188.67 </TD> <TD align="right"> 1648.03 </TD> </TR>
  <TR> <TD align="right"> conservation </TD> <TD align="right"> 1038.67 </TD> <TD align="right"> 202.79 </TD> <TD align="right"> 1241.46 </TD> </TR>
  <TR> <TD align="right"> vegetation </TD> <TD align="right"> 1100.87 </TD> <TD align="right"> 199.00 </TD> <TD align="right"> 1299.87 </TD> </TR>
  <TR> <TD align="right"> bio1 </TD> <TD align="right"> 1366.41 </TD> <TD align="right"> 180.91 </TD> <TD align="right"> 1547.31 </TD> </TR>
  <TR> <TD align="right"> bio2 </TD> <TD align="right"> 1362.95 </TD> <TD align="right"> 196.31 </TD> <TD align="right"> 1559.26 </TD> </TR>
  <TR> <TD align="right"> bio3 </TD> <TD align="right"> 1441.70 </TD> <TD align="right"> 189.88 </TD> <TD align="right"> 1631.58 </TD> </TR>
  <TR> <TD align="right"> bio4 </TD> <TD align="right"> 1377.84 </TD> <TD align="right"> 175.20 </TD> <TD align="right"> 1553.04 </TD> </TR>
  <TR> <TD align="right"> bio5 </TD> <TD align="right"> 1261.78 </TD> <TD align="right"> 203.83 </TD> <TD align="right"> 1465.61 </TD> </TR>
  <TR> <TD align="right"> bio6 </TD> <TD align="right"> 1454.76 </TD> <TD align="right"> 187.13 </TD> <TD align="right"> 1641.89 </TD> </TR>
  <TR> <TD align="right"> bio7 </TD> <TD align="right"> 1355.76 </TD> <TD align="right"> 186.80 </TD> <TD align="right"> 1542.57 </TD> </TR>
  <TR> <TD align="right"> bio8 </TD> <TD align="right"> 1301.47 </TD> <TD align="right"> 191.31 </TD> <TD align="right"> 1492.78 </TD> </TR>
  <TR> <TD align="right"> bio9 </TD> <TD align="right"> 1177.35 </TD> <TD align="right"> 184.58 </TD> <TD align="right"> 1361.93 </TD> </TR>
  <TR> <TD align="right"> bio10 </TD> <TD align="right"> 1270.14 </TD> <TD align="right"> 197.19 </TD> <TD align="right"> 1467.32 </TD> </TR>
  <TR> <TD align="right"> bio11 </TD> <TD align="right"> 1431.17 </TD> <TD align="right"> 177.59 </TD> <TD align="right"> 1608.76 </TD> </TR>
  <TR> <TD align="right"> bio12 </TD> <TD align="right"> 1355.68 </TD> <TD align="right"> 193.86 </TD> <TD align="right"> 1549.54 </TD> </TR>
  <TR> <TD align="right"> bio13 </TD> <TD align="right"> 1436.17 </TD> <TD align="right"> 198.81 </TD> <TD align="right"> 1634.98 </TD> </TR>
  <TR> <TD align="right"> bio14 </TD> <TD align="right"> 1324.93 </TD> <TD align="right"> 189.42 </TD> <TD align="right"> 1514.34 </TD> </TR>
  <TR> <TD align="right"> bio15 </TD> <TD align="right"> 1379.40 </TD> <TD align="right"> 187.84 </TD> <TD align="right"> 1567.25 </TD> </TR>
  <TR> <TD align="right"> bio16 </TD> <TD align="right"> 1437.55 </TD> <TD align="right"> 198.11 </TD> <TD align="right"> 1635.66 </TD> </TR>
  <TR> <TD align="right"> bio17 </TD> <TD align="right"> 1315.26 </TD> <TD align="right"> 189.81 </TD> <TD align="right"> 1505.06 </TD> </TR>
  <TR> <TD align="right"> bio18 </TD> <TD align="right"> 1410.36 </TD> <TD align="right"> 202.59 </TD> <TD align="right"> 1612.95 </TD> </TR>
  <TR> <TD align="right"> bio19 </TD> <TD align="right"> 1247.06 </TD> <TD align="right"> 196.66 </TD> <TD align="right"> 1443.72 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 1362.20 </TD> <TD align="right"> 181.26 </TD> <TD align="right"> 1543.46 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 1413.74 </TD> <TD align="right"> 199.67 </TD> <TD align="right"> 1613.41 </TD> </TR>
   </TABLE>


#### Table S3. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for the jaccard pairwise similarity. The Mean Square Error is the sum of the Bias and Variance.

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Apr 16 10:48:16 2014 -->
<TABLE border=0, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 8.80 </TD> <TD align="right"> 1.88 </TD> <TD align="right"> 10.68 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 7.07 </TD> <TD align="right"> 9.43 </TD> <TD align="right"> 16.51 </TD> </TR>
  <TR> <TD align="right"> Long2 </TD> <TD align="right"> 4.38 </TD> <TD align="right"> 2.15 </TD> <TD align="right"> 6.53 </TD> </TR>
  <TR> <TD align="right"> Lat2 </TD> <TD align="right"> 4.22 </TD> <TD align="right"> 2.19 </TD> <TD align="right"> 6.41 </TD> </TR>
  <TR> <TD align="right"> conservation </TD> <TD align="right"> 5.99 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 8.25 </TD> </TR>
  <TR> <TD align="right"> vegetation </TD> <TD align="right"> 5.52 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 7.77 </TD> </TR>
  <TR> <TD align="right"> bio1 </TD> <TD align="right"> 4.94 </TD> <TD align="right"> 2.23 </TD> <TD align="right"> 7.17 </TD> </TR>
  <TR> <TD align="right"> bio2 </TD> <TD align="right"> 5.91 </TD> <TD align="right"> 2.22 </TD> <TD align="right"> 8.13 </TD> </TR>
  <TR> <TD align="right"> bio3 </TD> <TD align="right"> 4.75 </TD> <TD align="right"> 2.21 </TD> <TD align="right"> 6.96 </TD> </TR>
  <TR> <TD align="right"> bio4 </TD> <TD align="right"> 4.29 </TD> <TD align="right"> 2.15 </TD> <TD align="right"> 6.44 </TD> </TR>
  <TR> <TD align="right"> bio5 </TD> <TD align="right"> 6.27 </TD> <TD align="right"> 2.27 </TD> <TD align="right"> 8.54 </TD> </TR>
  <TR> <TD align="right"> bio6 </TD> <TD align="right"> 5.01 </TD> <TD align="right"> 2.21 </TD> <TD align="right"> 7.22 </TD> </TR>
  <TR> <TD align="right"> bio7 </TD> <TD align="right"> 5.13 </TD> <TD align="right"> 2.19 </TD> <TD align="right"> 7.32 </TD> </TR>
  <TR> <TD align="right"> bio8 </TD> <TD align="right"> 4.72 </TD> <TD align="right"> 2.16 </TD> <TD align="right"> 6.88 </TD> </TR>
  <TR> <TD align="right"> bio9 </TD> <TD align="right"> 5.90 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 8.15 </TD> </TR>
  <TR> <TD align="right"> bio10 </TD> <TD align="right"> 5.96 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 8.22 </TD> </TR>
  <TR> <TD align="right"> bio11 </TD> <TD align="right"> 4.46 </TD> <TD align="right"> 2.21 </TD> <TD align="right"> 6.67 </TD> </TR>
  <TR> <TD align="right"> bio12 </TD> <TD align="right"> 5.78 </TD> <TD align="right"> 2.23 </TD> <TD align="right"> 8.01 </TD> </TR>
  <TR> <TD align="right"> bio13 </TD> <TD align="right"> 5.53 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 7.76 </TD> </TR>
  <TR> <TD align="right"> bio14 </TD> <TD align="right"> 5.32 </TD> <TD align="right"> 2.22 </TD> <TD align="right"> 7.54 </TD> </TR>
  <TR> <TD align="right"> bio15 </TD> <TD align="right"> 5.10 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 7.35 </TD> </TR>
  <TR> <TD align="right"> bio16 </TD> <TD align="right"> 5.61 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 7.86 </TD> </TR>
  <TR> <TD align="right"> bio17 </TD> <TD align="right"> 5.37 </TD> <TD align="right"> 2.23 </TD> <TD align="right"> 7.60 </TD> </TR>
  <TR> <TD align="right"> bio18 </TD> <TD align="right"> 5.83 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 8.09 </TD> </TR>
  <TR> <TD align="right"> bio19 </TD> <TD align="right"> 5.75 </TD> <TD align="right"> 2.23 </TD> <TD align="right"> 7.98 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 4.65 </TD> <TD align="right"> 2.18 </TD> <TD align="right"> 6.84 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 5.76 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 8.00 </TD> </TR>
   </TABLE>


----------------------------






![plot of chunk plot logis hill](figure/plot logis hill.png) 

#### Fig S6.

![plot of chunk plot logis rich and srich](figure/plot logis rich and srich.png) 

#### Fig S7.



### Figures comparing the models




### Fig sp comp




# Comparison of the observed and predicted by the models















