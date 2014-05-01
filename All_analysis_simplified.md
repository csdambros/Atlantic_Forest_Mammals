---
output: html_document
---


Comparing dispersal, environmental and mid-domain effects on species distribution
============================







##### By CSDambros

>html updated at 2014-04-30 10:07:31

### Methods

#### Study site
 
  The study was conducted in the whole area recognized as the Atlantic forest biome (Fig. 1). This biome encompasses an extent of 102,012 km², but today representing only 7.9% of its original forest cover. It has several vegetation types such as rainforests, mixed (Araucarian) moist forests, semideciduous forests, dry forests and upland grasslands. Rainforests tend to occur near the coast and semideciduous and dry forests far from coast; mixed forests are common in the south of AF (SOS Mata Atlântica 2007). AF shows moist tropical and subtropical climates, without well-defined dry season, and annual mean temperatures above 15ºC (Leite, 2002). 








#### Data collection

  We compiled a database based on the literature in which diverse authors have sampled small mammals along the AF (Appendix S1 in Supporting Information). We have used Google @ search tool for sampling articles, and the main keywords used in combination were “small mammal”, “marsupial”, “rodent”, “community”, “composition”, “richness”, “diversity”, and “Atlantic Forest”. Unpublished data sampled by NCC was also included (Appendix S1). There was some variation in the area covered by each surveyed study (see Appendix S1), but we have established a minimum of sampling effort for a given article to be accepted in our database: at least 1000 trap-nights, 6 months of field work, and use of wire type and/or Sherman live-traps installed on the ground or understory level of the forest. Besides that, and for independence assumptions, we have chosen surveyed areas distant at least 10 km from each other; most of the field surveys available in the literature sampled only one geographical location under this condition. From each selected survey, we obtained local species composition and abundance data when available. 
  We compiled the 19 environmental variables available in Bioclim (www.worldclim.org/bioclim) in a scale of 2.5 arc minutes: annual mean temperature (1), mean diurnal range (2), isothermality (3), temperature seasonality (4), maximum and minimum temperature of the warmest and coldest months (5 and 6), temperature annual range (7), mean temperature of the wettest, driest, warmest and coldest quarters (8-11), annual precipitation (12), precipitation of the wettest and driest months (13 and 14), precipitation seasonality (15), and precipitation of the wettest, direst, warmest and coldest quarters (16-19). Because most of the variables are correlated to each other, the climatic variables were summarized in two Principal Component axes. These axes were then used as predictor variables in all models. We present the results using individual variables as supplements (see below).  
  Data on habitat quality were taken directly from information available on the articles consulted, and classified as categories of forest conservation (degrees 1 to 5), where a score of 5 means primary pristine forest, and disturbed secondary forest, including with clearings, has a score of 1.








#### Analysis

  The standardized environmental variables were used as predictors of the species diversity and included in simulation models (see below). Latitude and Longitude were highly correlated with the climatic varialbes, so they were not included in the models as predictor variables, but were used to calculate the pairwise distance among sampling units.









  We analyzed the data using the original local sampling units and grouping the original data points in 26 2x2 degrees grid cells (Fig. 1). Most localities were surveyed by a single study, then most of the sampling points represent a single study (see data collection). The number of species encountered in each locality (S; alpha diversity), and the similarity in species composition (beta diversity) were used as response variables.
  The similarity in species composition was acessed by calculating the Jaccard similarity index among all pairs of samples. The matrix of similarities was summarized using a Principal Coordinates Analysis (PCoA). We used only the first two ordination axes because all axes in higher dimensions captured alone less than 10% of the variation from the similarity matrix.








![plot of chunk Fig_1_Study_site](figure/Fig_1_Study_site.png) 

#### Fig. 1. Map of the Atlantic Forest (AF; blue shade) showing the orinal sampling points (red circles) and the grid cells encompassing the entire AF. The size of the circles represent the size o the original sampling area in the log scale. Green shade represent the grid cells where small mammal data was available. 

####

  For estimating the influence of dispersal limitation on species local diversity (S) and turnover (Jaccard index; PCoA1 and PCoA2), we created a network of interconnected grid cells representing the AF. This network was used to determine the flux of species or individuals between grid cells (Fig. S2). Two models were used to recreate the species distribution under dispersal alone: the mid domain spreading-dye model (Colwell and Hurtt 1994; Colwell and Lees 2000), and the analytical neutral approach borrowed from the population genetics (Nagylaki, 1980) and presented by Economo & Keitt (2008) for community ecology. These models were created just using the grouped data in grid cells because of the completeness of data, and due to computational limitation.
  In the mid domain spreading dye model, the number of grid cells occupied by each small mammal species was recorded. For each species, one of the 26 grid cells was randomly selected and the species occurrence was spread from the selected cell to neighboring cells until the original number of grid cells was occupied. This procedure was repeated 10,000 times for the 64 species. Each cell had up to eight neighbors (Moore neighborhood; Fig. S2), and the model was bounded by the domain where actual small mammal data was recorded (26 grid cells).








  In the neutral model, the whole area comprising the AF was divided in 56 2x2 degrees grid cells from which 26 had small mammal data available (Fig. 1). The model started with a single species occupying all the 56 cells. In each generation, new species were added in each cell by point speciation with rate $\nu$, set the same for al cells (see Economo & Keitt (2010) and Muneepeerakul et al. (2008) for more details and other uses). $\nu$ represents the probability of an individual to become a new species but could also represent the addition of new species by immigration from a larger species pool outside the AF (eg. the Cerrado or Amazonian forests). To recreate the dispersal of individuals, we determined that a cell could just be colonized by a neighbor (Moore neighborhood; Fig S2), and that all grid cells had the same migration rate (parameter $m$). The local community size (number of individuals) was set the same for all grid cells ($N$ = 100). The model was run for multiple generations, until the diversity within (PIE and HillPIE) and between all cells (Morisita-Horn similarity) reached a steady-state (usually more than 30,000 generations).


  Differently from the Mid-domain model, the values of $m$, $\nu$, and $N$ were initially set independently of the observed data. Because different combinations of these parameters can create the same patterns of species distribution (Etienne 2006?; Economo and Keitt 2010), and because estimates for these parameters based on empirical data are not available, we did not attempt to estimate realistic parameters, but focused on the general patterns that the model could produce (see discussion for implications of this procedure). $N$ was arbitrarily defined as 100 individuals for all grid cells, and $m$ and $\nu$ were used as knobs to best fit the model to the observed similarity in species composition. These parameters were simultaneously tunned using the L-BFSG-S optimization algorithm. In summary, $m$ and $\nu$ were adjusted to minimize the differences between the predictions of the neutral model to the observed data, giving the best possible explanation of dispersal to the observed data (this procedure is conceptually analog to a regression analysis using maximum likelihood optimization). Because the number of individuals per grid cell was kept constant during the optimization, the differences in diversity were due to the location of a cell in the network and the number of neighbor grid cells, but not affected by the differences in abundance.











  The model proposed by Economo and Keitt (2008) is probabilistic and does not require the simulation of each individual in the metacommunity, being less computationally intensive. This model allowed us to investigate thoroughly the parameter space of $m$ and $\nu$ to find the best explanatory model to the observed data. However, because the model is based on probabilities, it does not allow one to calculate statistics based on the raw data (eg. the Jaccard similarity index). Because not all original studies recorded species abundances, we compared the Morisita-Horn similarity matrix from the neutral model with the observed Jaccard similarity for optimization. These indexes are usually highly correlated (Chao et al., 2006; Krasnov et al., 2005) and we believe this procedure did not have profound effects on the neutral fit.
  To calculate the Jaccard similarity index from the neutral model, each individual in the metacommunity was simulated using the previously optimized parameters. We run an initial buffer of 30,000 generations. Then the model was run for additional 1000 generations 5,000 times, representing 5,000 simulations. The mean of all the 5,000 simulations was used to calculate the remaining summary statistics (see below). The correlation between the probabilistic and the simulated model was 99.98% when comparing the Morisita-Horn index of pairwise similarity, so we believe the simulations run long enough to capture the final predictions of the neutral model.







To test the association of the species diversity with the climatic and habitat quality variables, individual logistic regressions were fit for all species against these predictor variables, allowing each species to vary independently of one another.




 The probability of occurrence generated by the logistic model for each species was summed to estimate the expected number of species in each quadrant:
  
$$
\hat S_i = \sum_{j=1}^{S}{\frac{e^{X_i\beta_j}}{1+e^{X_i\beta_j}}}
$$
  
Where $X_i$ denotes a vector of the environmental variables in the quadrant $i$, and $\beta_j$ is a vector with the coefficients from the logistic regression for the species $j$.




To calculate the raw statistics based on these probabilities, the distribution of each species was simulated based on the probabilities of occurrence generated by the logistic curve, similarly to what was done for the neutral model. For each species, the number of grid cells where the species was observed was counted. The model spread each species into the map based on its probability of occurrence until the number of grid cells occupied in the simulation matched the number observed. This allowed each species to vary independently, and to the relationship between the environmental variables and the response variables (eg. number of species) to be non-linear. Differently from the non-linearities that could be used in other statistical tests, such as Generalized Additive Models, the deterministic part of this model has an underlying mechanistic function (logistic). Note that this model does not require the species to have continous ranges as in the mid domain model.








We compared the three simulation models (Mid Domain, Neutral, and based on logistic regressions) by their Mean Square Error (MSE; Gotelli et al. 2009). The MSE was calculated as the sum of the squared bias and model variance, so the MSE is a balance between precision and accuracy (see Gotelli et al. 2009 and Appendix S3 for detailed description).







Additionally, we tested for the direct association of the response variables against the environmental predictors. For this test, linear regression models (OLS) were fit for species richness and the two ordination axes of species composition. For the original samplin units, the OLS models had the climatic and habitat quality variables as predictors. For the grouped data, we created three sets of models: (1) Using just the climatic and habitat quality variables as predictors; (2) using just the neutral and mid domain models as predictors; (3) using the neutral and mid domain predictions as predictor variables along with the climatic and habitat quality variables, as suggested by Letten et al. (2013); and (4) using the neutral and mid domain predictions as null hypotheses, and testing for the association of residuals with the climatic and habitat quality variables. Both the raw environmental variables and the prediction from the logistic simulations were used as predictor variables. This gave a total of 20 models (REnv; Neutral; MidD; REnv+Neutral; REnv+MidD; LEnv; REnv+Neutral; REnv+MidD; ResNeutral+Env;  ResMidD+Env; Table 2; Table 4). We assumed normality in the residuals of all models. The models were compared by their AIC values.
  The regression coefficients were set to 1 for the neutal, mid domain, and logistic simulation predictors, so this parameter was not allowed to vary. We did not include an intercept in the models using only these variables as predictors as well, then there were no extra coefficients to be estimated besides the standard deviation in the residuals (one extra parameter). Note however that an intercept and slope were calculated for each species in the logistic models. These extra parameters are added as penalties for the AIC calculations.




  



  Finally, we tested for patterns of distance-decay in the species similarity using simple and partial Mantel tests correlating the matrix of Jaccard similarity (M) to the geographical distance (D), and environmental dissimilarity (E) matrices (Euclidean distance) (see Thompson and Townsend 2006). Legendre et al. (2005) and Tuomisto and Ruokolainen (2006; 2008) suggest using multiple regression on similarity matrices to separate the effects of niche and neutral processes. However, there is a strong debate about the validity of these models (Legendre et al. 2008; Tuomisto and Ruokolainen 2008). Some authors also suggest using partial Redundancy Analysis (Borcard et al. 1992) for this purpose (Gilbert and Lechowicz 2004), but this method also has serious limitations (Smith and Lundholm 2010). We preferred to base our discussion on the mantel tests and multiple regressions using the summarized data (each grid cell as a unit), as described above. We used (simple) generalized linear models (GLMs) with log links to estimate the relationship between the similarity in species composition and geographical or environmental distances (Millar et al. 2011). Because the Jaccard similarity is a proportion (proportion of shared species), the error of this model was fit with a binomial distribution (Millar et al. 2011). The R-squares of the GLMs was calculated by using the McFadden's apporach (McFadden 1974).













  All the analyses were conducted in the R program (R Development Core Team, 2013). The models and most of the summary statistics calculations were implemented by the authors, and are available as a supplementary material (see Appendix SXX). We used the packages vegan (Oksanen et al., 2008) for the remaining analyses.

  	
### Results

  All the models had a very poor fit to species richness (Fig. 2; Fig. S2). The model with the lowest mean square error and bias was the logistic simulation using habitat quality as a predictor of species presence (Table 1). The model with the lowest variance was the neutral model, but the difference among the models was much smaller than for the bias (Table 1).

#### Table 1. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for species richness. The Mean Square Error is the sum of the Bias and Variance.

####

<TABLE border=1, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 1277.60 </TD> <TD align="right"> 199.85 </TD> <TD align="right"> 1477.44 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 1518.38 </TD> <TD align="right"> 145.41 </TD> <TD align="right"> 1663.79 </TD> </TR>
  <TR> <TD align="right"> habitat </TD> <TD align="right"> 1033.74 </TD> <TD align="right"> 201.66 </TD> <TD align="right"> 1235.41 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 1363.23 </TD> <TD align="right"> 180.92 </TD> <TD align="right"> 1544.15 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 1403.98 </TD> <TD align="right"> 203.76 </TD> <TD align="right"> 1607.73 </TD> </TR>
   </TABLE>


###

Similarly, the regression model with the lowest AIC was the linear regression of the species richness against habitat quality (Table 2). However, the inclusion of the neutral model prediction in the same model produced a similar AIC value and increased the explanatory power in about 7% (Table 2). Habitat quality alone was able to explain 24% of the variation in species richness. The climatic variables were able to explain 7% of the variation on species richness. The logistic models had similar explanatory power to the models using the raw environmental predictors.
  Besides the poor fit for the data, the Mid Domain and Neutral models created the well known richness peak in the central areas of the Atlantic forest (Fig. 2). The neutral and Mid Domain models were able to explain alone 9% and 8% of the variation on species richness, respectively.
  In the small grain (original sampling units), the habitat quality and climatic variables were able to explain 21% and 2% of the variation in species richness.

#### Table 2. AIC values comparing the 10 regression models tested for species Richness. MR1.1: Environmental variables as predictors; MR1.2: Environmental variables + Neutral as predictors; MR1.3: Environmental variables + Mid Domain as predictors; MR1.4: Neutral; MR1.5: Mid Domain; MR1.6: Environmental variables as predictors of residuals from Neutral; MR1.7:  Environmental variables as predictors of residuals from Mid Domain; MR1.8: Logistic models as predictors; MR1.9: Logistic model + Neutral as predictors; MR1.10: Logistic model + Mid Domain as predictors.


####

<TABLE border=1, bgcolor=#989898>
<TR> <TH>  </TH> <TH> model </TH> <TH> Var </TH> <TH> logLik </TH> <TH> npar </TH> <TH> AIC </TH> <TH> rsquares </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> MR1.1 </TD> <TD> Clim </TD> <TD align="right"> -35.44 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 78.88 </TD> <TD align="right"> 0.07 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> MR1.2 </TD> <TD> Clim+Neu </TD> <TD align="right"> -34.51 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 81.02 </TD> <TD align="right"> 0.15 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> MR1.3 </TD> <TD> Clim+Mid </TD> <TD align="right"> -35.40 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 78.80 </TD> <TD align="right"> 0.16 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> MR1.4 </TD> <TD> Neu </TD> <TD align="right"> -37.19 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 80.39 </TD> <TD align="right"> 0.09 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> MR1.5 </TD> <TD> Mid </TD> <TD align="right"> -41.14 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 84.28 </TD> <TD align="right"> 0.08 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> MR1.6 </TD> <TD> ResN~Clim </TD> <TD align="right"> -14.43 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 40.86 </TD> <TD align="right"> 0.14 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> MR1.7 </TD> <TD> ResN~Clim </TD> <TD align="right"> -15.42 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 38.85 </TD> <TD align="right"> 0.02 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> MR1.8 </TD> <TD> LogisC </TD> <TD align="right"> -36.29 </TD> <TD align="right"> 129.00 </TD> <TD align="right"> 330.57 </TD> <TD align="right"> 0.01 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> MR1.9 </TD> <TD> LogisC+Neu </TD> <TD align="right"> -34.31 </TD> <TD align="right"> 131.00 </TD> <TD align="right"> 330.62 </TD> <TD align="right"> 0.19 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> MR1.10 </TD> <TD> LogisC+Mid </TD> <TD align="right"> -35.88 </TD> <TD align="right"> 129.00 </TD> <TD align="right"> 329.76 </TD> <TD align="right"> 0.18 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> MR1.11 </TD> <TD> Hab </TD> <TD align="right"> -32.79 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 71.57 </TD> <TD align="right"> 0.24 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> MR1.12 </TD> <TD> Hab+Neu </TD> <TD align="right"> -32.88 </TD> <TD align="right"> 5.00 </TD> <TD align="right"> 75.77 </TD> <TD align="right"> 0.30 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> MR1.13 </TD> <TD> Hab+Mid </TD> <TD align="right"> -36.21 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 78.42 </TD> <TD align="right"> 0.31 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> MR1.16 </TD> <TD> ResN~Hab </TD> <TD align="right"> -15.23 </TD> <TD align="right"> 5.00 </TD> <TD align="right"> 40.45 </TD> <TD align="right"> 0.01 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> MR1.17 </TD> <TD> ResN~Hab </TD> <TD align="right"> -15.51 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 37.01 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> MR1.18 </TD> <TD> LogisH </TD> <TD align="right"> -33.31 </TD> <TD align="right"> 65.00 </TD> <TD align="right"> 196.63 </TD> <TD align="right"> 0.21 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> MR1.19 </TD> <TD> LogisH+Neu </TD> <TD align="right"> -34.05 </TD> <TD align="right"> 67.00 </TD> <TD align="right"> 202.11 </TD> <TD align="right"> 0.27 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> MR1.20 </TD> <TD> LogisH+Mid </TD> <TD align="right"> -36.95 </TD> <TD align="right"> 65.00 </TD> <TD align="right"> 203.90 </TD> <TD align="right"> 0.27 </TD> </TR>
   </TABLE>


###


![plot of chunk Fig_2_Map_Richness_Mid_Domain_vs_Neutral](figure/Fig_2_Map_Richness_Mid_Domain_vs_Neutral1.png) ![plot of chunk Fig_2_Map_Richness_Mid_Domain_vs_Neutral](figure/Fig_2_Map_Richness_Mid_Domain_vs_Neutral2.png) ![plot of chunk Fig_2_Map_Richness_Mid_Domain_vs_Neutral](figure/Fig_2_Map_Richness_Mid_Domain_vs_Neutral3.png) 

#### Fig 2. Projection of the species richness predicted by dispersal models on the map and comparison with the observed species distribution. Warmer colors represent areas with hiher species richness. (A) Observed species richness; (B) Predicted by the climatic variables in the logistic models; (B) Predicted by habitat quality in the logistic models; (D) Predicted by the Mid Domain model; (E) Predicted by the neutral model. Blank cells represent areas included in the neutral model but where actual data is not available.


####


  The similarity in species composition was well fit for both the neutral, the Mid Domain, and the climatic models. The MSE in the models, and the difference in MSE between the models were much smaller for the species composition than for richness (Table 3), besides the much larger number of comparisons for the species composition (325 for the pairwise similarity matrix against 26 for the vector of species richness). The lowest mean square error and bias were obtained using the climatic variables as predictors of the species occurrence in the logistic models, while the mid domain model had the lowest variance (Table 3). The neutral model had the highest mean square error due to its largest variance (Table 3).

#### Table 3. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for the jaccard pairwise similarity. The Mean Square Error is the sum of the Bias and Variance.

####

<TABLE border=1, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 8.83 </TD> <TD align="right"> 1.88 </TD> <TD align="right"> 10.71 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 7.07 </TD> <TD align="right"> 9.43 </TD> <TD align="right"> 16.51 </TD> </TR>
  <TR> <TD align="right"> habitat </TD> <TD align="right"> 5.95 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 8.20 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 4.64 </TD> <TD align="right"> 2.20 </TD> <TD align="right"> 6.84 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 5.77 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 8.02 </TD> </TR>
   </TABLE>


####

Conversely, the best regression models for the first ordination axis were obtained including the mid domain or neutral model results as predictor variables (Table 4). For the second ordination axis, the best model included habitat quality, and habitat quality along the mid domain model. The neutral and mid domain models alone were able to explain 78% and 77% of the variation in the first axis of species composition. The climatic variables combined explained 69% of the variation is the first axis. Habitat quality alone was a poor predictor of the species composition in the first ordination axis, but the best explanatory model was obtained by combining habitat quality and one of the dispersal models (R²=82%). 
For the second axis, the mid domain model along habitat quality explained about 20% of the variation in species composition (Table 4). Besides the largest AICs due to the high number of parameters, the combination of the dispersal models with the climatic variables in the logistic regressions gave the best explanatory power for the second ordination axis (R²=68%; Table 4). 

#### Table 4. AIC values comparing the 10 regression models tested for species composition. Specias composition was measured as the first (MC1) and second (MC2) axes of the principal component analysis summarizing the jaccard similarity matrix. MC1.1 and MC2.1: Environmental variables as predictors; MC1.2 and MC2.2: Environmental variables + Neutral as predictors; MC1.3 and MC2.3: Environmental variables + Mid Domain as predictors; MC1.4 and MC2.4: Neutral; MC1.5 and MC2.5: Mid Domain; MC1.6 and MC2.6: Environmental variables as predictors of residuals from Neutral; MC1.7 and MC2.7:  Environmental variables as predictors of residuals from Mid Domain; MR1.8 and MC2.8: Logistic models as predictors; MC1.9 and MC2.9: Logistic model + Neutral as predictors; MC1.10 and MC2.10: Logistic model + Mid Domain as predictors.

####

<TABLE border=1, bgcolor=#989898>
<TR> <TH>  </TH> <TH> model </TH> <TH> Var </TH> <TH> logLik </TH> <TH> npar </TH> <TH> AIC </TH> <TH> rsquares </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> MC1.1 </TD> <TD> Clim </TD> <TD align="right"> -21.10 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 50.20 </TD> <TD align="right"> 0.69 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> MC1.2 </TD> <TD> Clim+Neu </TD> <TD align="right"> -16.36 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 44.72 </TD> <TD align="right"> 0.79 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> MC1.3 </TD> <TD> Clim+Mid </TD> <TD align="right"> -16.72 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 41.44 </TD> <TD align="right"> 0.78 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> MC1.4 </TD> <TD> Neu </TD> <TD align="right"> -17.42 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 40.84 </TD> <TD align="right"> 0.78 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> MC1.5 </TD> <TD> Mid </TD> <TD align="right"> -18.07 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 38.13 </TD> <TD align="right"> 0.77 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> MC1.6 </TD> <TD> ResN~Clim </TD> <TD align="right"> -6.18 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 24.35 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> MC1.7 </TD> <TD> ResN~Clim </TD> <TD align="right"> -4.91 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 17.82 </TD> <TD align="right"> 0.05 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> MC1.8 </TD> <TD> LogisC </TD> <TD align="right"> -21.59 </TD> <TD align="right"> 129.00 </TD> <TD align="right"> 301.19 </TD> <TD align="right"> 0.68 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> MC1.9 </TD> <TD> LogisC+Neu </TD> <TD align="right"> -16.39 </TD> <TD align="right"> 131.00 </TD> <TD align="right"> 294.79 </TD> <TD align="right"> 0.79 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> MC1.10 </TD> <TD> LogisC+Mid </TD> <TD align="right"> -17.16 </TD> <TD align="right"> 129.00 </TD> <TD align="right"> 292.32 </TD> <TD align="right"> 0.77 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> MC1.11 </TD> <TD> Hab </TD> <TD align="right"> -35.46 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 76.93 </TD> <TD align="right"> 0.07 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> MC1.12 </TD> <TD> Hab+Neu </TD> <TD align="right"> -14.99 </TD> <TD align="right"> 5.00 </TD> <TD align="right"> 39.98 </TD> <TD align="right"> 0.82 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> MC1.13 </TD> <TD> Hab+Mid </TD> <TD align="right"> -15.33 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 36.66 </TD> <TD align="right"> 0.82 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> MC1.16 </TD> <TD> ResN~Hab </TD> <TD align="right"> -3.49 </TD> <TD align="right"> 5.00 </TD> <TD align="right"> 16.97 </TD> <TD align="right"> 0.39 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> MC1.17 </TD> <TD> ResN~Hab </TD> <TD align="right"> -3.80 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 13.60 </TD> <TD align="right"> 0.22 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> MC1.18 </TD> <TD> LogisH </TD> <TD align="right"> -35.28 </TD> <TD align="right"> 65.00 </TD> <TD align="right"> 200.55 </TD> <TD align="right"> 0.08 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> MC1.19 </TD> <TD> LogisH+Neu </TD> <TD align="right"> -14.34 </TD> <TD align="right"> 67.00 </TD> <TD align="right"> 162.68 </TD> <TD align="right"> 0.83 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> MC1.20 </TD> <TD> LogisH+Mid </TD> <TD align="right"> -14.77 </TD> <TD align="right"> 65.00 </TD> <TD align="right"> 159.54 </TD> <TD align="right"> 0.83 </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD> MC2.1 </TD> <TD> Clim </TD> <TD align="right"> -34.00 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 76.01 </TD> <TD align="right"> 0.17 </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD> MC2.2 </TD> <TD> Clim+Neu </TD> <TD align="right"> -33.69 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 79.37 </TD> <TD align="right"> 0.21 </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD> MC2.3 </TD> <TD> Clim+Mid </TD> <TD align="right"> -33.95 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 75.90 </TD> <TD align="right"> 0.20 </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD> MC2.4 </TD> <TD> Neu </TD> <TD align="right"> -40.36 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 86.72 </TD> <TD align="right"> 0.10 </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD> MC2.5 </TD> <TD> Mid </TD> <TD align="right"> -39.51 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 81.01 </TD> <TD align="right"> 0.13 </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD> MC2.6 </TD> <TD> ResN~Clim </TD> <TD align="right"> -15.28 </TD> <TD align="right"> 6.00 </TD> <TD align="right"> 42.56 </TD> <TD align="right"> 0.02 </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD> MC2.7 </TD> <TD> ResN~Clim </TD> <TD align="right"> -14.81 </TD> <TD align="right"> 4.00 </TD> <TD align="right"> 37.61 </TD> <TD align="right"> 0.05 </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD> MC2.8 </TD> <TD> LogisC </TD> <TD align="right"> -33.18 </TD> <TD align="right"> 129.00 </TD> <TD align="right"> 324.35 </TD> <TD align="right"> 0.22 </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD> MC2.9 </TD> <TD> LogisC+Neu </TD> <TD align="right"> -31.78 </TD> <TD align="right"> 131.00 </TD> <TD align="right"> 325.55 </TD> <TD align="right"> 0.35 </TD> </TR>
  <TR> <TD align="right"> 28 </TD> <TD> MC2.10 </TD> <TD> LogisC+Mid </TD> <TD align="right"> -32.09 </TD> <TD align="right"> 129.00 </TD> <TD align="right"> 322.19 </TD> <TD align="right"> 0.36 </TD> </TR>
  <TR> <TD align="right"> 29 </TD> <TD> MC2.11 </TD> <TD> Hab </TD> <TD align="right"> -33.80 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 73.60 </TD> <TD align="right"> 0.18 </TD> </TR>
  <TR> <TD align="right"> 30 </TD> <TD> MC2.12 </TD> <TD> Hab+Neu </TD> <TD align="right"> -33.69 </TD> <TD align="right"> 5.00 </TD> <TD align="right"> 77.37 </TD> <TD align="right"> 0.21 </TD> </TR>
  <TR> <TD align="right"> 31 </TD> <TD> MC2.13 </TD> <TD> Hab+Mid </TD> <TD align="right"> -33.95 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 73.90 </TD> <TD align="right"> 0.20 </TD> </TR>
  <TR> <TD align="right"> 32 </TD> <TD> MC2.16 </TD> <TD> ResN~Hab </TD> <TD align="right"> -15.29 </TD> <TD align="right"> 5.00 </TD> <TD align="right"> 40.58 </TD> <TD align="right"> 0.02 </TD> </TR>
  <TR> <TD align="right"> 33 </TD> <TD> MC2.17 </TD> <TD> ResN~Hab </TD> <TD align="right"> -15.00 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 36.00 </TD> <TD align="right"> 0.02 </TD> </TR>
  <TR> <TD align="right"> 34 </TD> <TD> MC2.18 </TD> <TD> LogisH </TD> <TD align="right"> -36.32 </TD> <TD align="right"> 65.00 </TD> <TD align="right"> 202.64 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 35 </TD> <TD> MC2.19 </TD> <TD> LogisH+Neu </TD> <TD align="right"> -40.36 </TD> <TD align="right"> 67.00 </TD> <TD align="right"> 214.72 </TD> <TD align="right"> 0.10 </TD> </TR>
  <TR> <TD align="right"> 36 </TD> <TD> MC2.20 </TD> <TD> LogisH+Mid </TD> <TD align="right"> -39.49 </TD> <TD align="right"> 65.00 </TD> <TD align="right"> 208.98 </TD> <TD align="right"> 0.13 </TD> </TR>
   </TABLE>


###

![plot of chunk Fig_3_Map_Composition_Mid_Domain_vs_Neutral](figure/Fig_3_Map_Composition_Mid_Domain_vs_Neutral1.png) ![plot of chunk Fig_3_Map_Composition_Mid_Domain_vs_Neutral](figure/Fig_3_Map_Composition_Mid_Domain_vs_Neutral2.png) ![plot of chunk Fig_3_Map_Composition_Mid_Domain_vs_Neutral](figure/Fig_3_Map_Composition_Mid_Domain_vs_Neutral3.png) 

#### Fig 3. Projection of the species composition predicted by dispersal models on the map and comparison with the observed species distribution. The colors represent the similarity in species composition measured by the pairwise jaccard similarity index between sites summirized in one axis of a principal coordinates analysis (pcoa). Those grid cells with similar colors have a similar composition of species. (A) Observed species composition; (B) Species composition predicted by the climatic variables in the logistic models; (C) Species composition predicted by habitat quality in the logistic models; (D) Species composition predicted by the Mid Domain model; (E) Species composition predicted by the neutral model. Shadded grid cells represent areas included in the neutral model but where actual data is not available.

####


  The mantel test results showed a stronger climatic effect on species composition when grouping the data in grid cells (Table 5). The geographical distance was highly correlated with the species composition in both cases. Geographical distance was correlated with the species composition when the climatic effect was partialed out, but the correlation was stronger in the unconstrained model (Table 5). The similarity in species composition was also weakly associated with the climatic variables when controlling for geographical distance in the grouped data. Habitat quality was a poor predictor of species composition in all models (Table 5).

#### Table 5. Mantel and partial Mantel test results comparing the correlation of species similarity against geographical distance and environmental dissimilarity. The minus sign indicate partial results (Geo - Env: Geographical distance without the effect of environmenal dissimilarity).

####


<TABLE border=1, bgcolor=#989898>
<TR> <TH>  </TH> <TH> statistic </TH> <TH> signif </TH>  </TR>
  <TR> <TD align="right"> Local Geo </TD> <TD align="right"> 0.48 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> Local Clim </TD> <TD align="right"> 0.35 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> Local Hab </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.41 </TD> </TR>
  <TR> <TD align="right"> 2d Geo </TD> <TD align="right"> 0.49 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 2d Env </TD> <TD align="right"> 0.42 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 2d Hab </TD> <TD align="right"> 0.02 </TD> <TD align="right"> 0.36 </TD> </TR>
  <TR> <TD align="right"> Local Geo-Env </TD> <TD align="right"> 0.38 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> Local Env-Geo </TD> <TD align="right"> 0.14 </TD> <TD align="right"> 0.01 </TD> </TR>
  <TR> <TD align="right"> Local Hab-Geo </TD> <TD align="right"> -0.02 </TD> <TD align="right"> 0.76 </TD> </TR>
  <TR> <TD align="right"> 2d Geo-Env </TD> <TD align="right"> 0.31 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 2d Env-Geo </TD> <TD align="right"> 0.11 </TD> <TD align="right"> 0.09 </TD> </TR>
  <TR> <TD align="right"> 2d Hab-Geo </TD> <TD align="right"> 0.11 </TD> <TD align="right"> 0.09 </TD> </TR>
   </TABLE>


###

Both neutral and Mid Domain models had very similar predictions for the distance-decay pattern in species similarity (Fig. 4). However, the decay predicted the dispersal models was much steeper in shorter distances than estimated by the glm model of species similarity against geographical distance (Fig. 4). The climatic variables predicted an almost linear decay in species similarity against geographical distance (Fig. 4).


![plot of chunk Fig_4_distance_decay_real_midDomain_and_neutral_mor](figure/Fig_4_distance_decay_real_midDomain_and_neutral_mor.png) 

#### Fig. 4. Distance-decay in species similarity using the Jaccard similarity index. The Mid Domain (red) and Neutral model (blue) have similar fit for the data, but the variance was much smaller than observed (grey). Habitat quality and the climatic variables did not predict the exponential decay with geographical distance.

####

----------------------------


![plot of chunk Fig_S1-Connectivity_map_26_and_56](figure/Fig_S1-Connectivity_map_26_and_56.png) 

#### Fig. S1. Map representing the connectivity of grid cells used to simulate the individual dispersal in the Mid Domain (A) and neutral (B) models.

####

![plot of chunk Fig_S2_rich_logis_in_the_map](figure/Fig_S2_rich_logis_in_the_map.png) 

#### Fig S2. Projection of the species richness predicted by logistic models on the map and comparison with the observed species distribution. Warmer colors represent areas with hiher species richness. (A) Observed species richness; (B-Z) Predicted by logistic regressions of individual species against environmental gradients. 


####

![plot of chunk Fig_S3_jac_logis_in_the_map](figure/Fig_S3_jac_logis_in_the_map.png) 

#### Fig S3. Projection of the species composition predicted by dispersal models on the map and comparison with the observed species distribution. The colors represent the similarity in species composition measured by the pairwise jaccard similarity index between sites summarized in one axis of a principal coordinates analysis (pcoa). Those grid cells with similar colors have a similar composition of species. (A) Observed species composition; (B-Z) Species composition predicted by logistic regressions of individual species against the environmental gradients.


####

![plot of chunk Fig_S4_plot_all_predictions_against_lat](figure/Fig_S4_plot_all_predictions_against_lat.png) 

#### Fig S4. Change in species composition along the latitudinal gradient as predicted by the Neutral and Mid Domain models. The composition was measured by the Jaccard similarity index between all pairs of sites and summarized by the first axis of a Principal Coordinates Analysis (PCoA).


![plot of chunk Fig_S5_Residuals_similarity](figure/Fig_S5_Residuals_similarity.png) 

#### Fig S5. Residual plots of the the fitted models using the pairwise Jaccard similarity as a response variable. The plots represent the models using (A) geographical distance; (B) environmental distance; (C) Neutral; and (D) Mid-Domain models as predictors. For A and B, the fitted model was a Generalized Linear Model with a log link and binomial distribution. For C and D, the graphs represent the deviations predicted by simulations (see data analysis for details).

###

![plot of chunk FigS6_residuals_regressions_richness](figure/FigS6_residuals_regressions_richness.png) 

#### Fig S6. Residual plots of the the fitted models using the pairwise Jaccard similarity as a response variable. The plots represent the models using (A) geographical distance; (B) environmental distance; (C) Neutral; and (D) Mid-Domain models as predictors. For A and B, the fitted model was a Generalized Linear Model with a log link and binomial distribution. For C and D, the graphs represent the deviations predicted by simulations (see data analysis for details).

###


![plot of chunk FigS7_residuals_regressions_composition](figure/FigS7_residuals_regressions_composition.png) 

#### Fig S7. Residual plots of the the fitted models using the pairwise Jaccard similarity as a response variable. The plots represent the models using (A) geographical distance; (B) environmental distance; (C) Neutral; and (D) Mid-Domain models as predictors. For A and B, the fitted model was a Generalized Linear Model with a log link and binomial distribution. For C and D, the graphs represent the deviations predicted by simulations (see data analysis for details).


####



#### Table S1. List of the sites, authors, etc. used in this manuscript

###

#### Table S2. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for species richness. The Mean Square Error is the sum of the Bias and Variance.

####

<TABLE border=1, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 1277.60 </TD> <TD align="right"> 199.85 </TD> <TD align="right"> 1477.44 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 1518.38 </TD> <TD align="right"> 145.41 </TD> <TD align="right"> 1663.79 </TD> </TR>
  <TR> <TD align="right"> Long2 </TD> <TD align="right"> 1353.43 </TD> <TD align="right"> 183.81 </TD> <TD align="right"> 1537.24 </TD> </TR>
  <TR> <TD align="right"> Lat2 </TD> <TD align="right"> 1452.79 </TD> <TD align="right"> 182.48 </TD> <TD align="right"> 1635.28 </TD> </TR>
  <TR> <TD align="right"> habitat </TD> <TD align="right"> 1033.74 </TD> <TD align="right"> 201.66 </TD> <TD align="right"> 1235.41 </TD> </TR>
  <TR> <TD align="right"> vegetation </TD> <TD align="right"> 1089.84 </TD> <TD align="right"> 195.53 </TD> <TD align="right"> 1285.38 </TD> </TR>
  <TR> <TD align="right"> bio1 </TD> <TD align="right"> 1366.20 </TD> <TD align="right"> 180.76 </TD> <TD align="right"> 1546.95 </TD> </TR>
  <TR> <TD align="right"> bio2 </TD> <TD align="right"> 1355.39 </TD> <TD align="right"> 196.31 </TD> <TD align="right"> 1551.70 </TD> </TR>
  <TR> <TD align="right"> bio3 </TD> <TD align="right"> 1455.14 </TD> <TD align="right"> 188.38 </TD> <TD align="right"> 1643.52 </TD> </TR>
  <TR> <TD align="right"> bio4 </TD> <TD align="right"> 1384.24 </TD> <TD align="right"> 178.17 </TD> <TD align="right"> 1562.41 </TD> </TR>
  <TR> <TD align="right"> bio5 </TD> <TD align="right"> 1274.38 </TD> <TD align="right"> 202.78 </TD> <TD align="right"> 1477.16 </TD> </TR>
  <TR> <TD align="right"> bio6 </TD> <TD align="right"> 1463.32 </TD> <TD align="right"> 188.82 </TD> <TD align="right"> 1652.14 </TD> </TR>
  <TR> <TD align="right"> bio7 </TD> <TD align="right"> 1369.02 </TD> <TD align="right"> 187.30 </TD> <TD align="right"> 1556.32 </TD> </TR>
  <TR> <TD align="right"> bio8 </TD> <TD align="right"> 1302.98 </TD> <TD align="right"> 189.30 </TD> <TD align="right"> 1492.28 </TD> </TR>
  <TR> <TD align="right"> bio9 </TD> <TD align="right"> 1161.33 </TD> <TD align="right"> 185.87 </TD> <TD align="right"> 1347.20 </TD> </TR>
  <TR> <TD align="right"> bio10 </TD> <TD align="right"> 1283.26 </TD> <TD align="right"> 196.89 </TD> <TD align="right"> 1480.15 </TD> </TR>
  <TR> <TD align="right"> bio11 </TD> <TD align="right"> 1431.97 </TD> <TD align="right"> 180.08 </TD> <TD align="right"> 1612.05 </TD> </TR>
  <TR> <TD align="right"> bio12 </TD> <TD align="right"> 1340.50 </TD> <TD align="right"> 197.49 </TD> <TD align="right"> 1537.99 </TD> </TR>
  <TR> <TD align="right"> bio13 </TD> <TD align="right"> 1419.80 </TD> <TD align="right"> 199.53 </TD> <TD align="right"> 1619.33 </TD> </TR>
  <TR> <TD align="right"> bio14 </TD> <TD align="right"> 1306.61 </TD> <TD align="right"> 189.95 </TD> <TD align="right"> 1496.55 </TD> </TR>
  <TR> <TD align="right"> bio15 </TD> <TD align="right"> 1377.16 </TD> <TD align="right"> 187.30 </TD> <TD align="right"> 1564.46 </TD> </TR>
  <TR> <TD align="right"> bio16 </TD> <TD align="right"> 1423.43 </TD> <TD align="right"> 200.08 </TD> <TD align="right"> 1623.51 </TD> </TR>
  <TR> <TD align="right"> bio17 </TD> <TD align="right"> 1321.05 </TD> <TD align="right"> 191.03 </TD> <TD align="right"> 1512.08 </TD> </TR>
  <TR> <TD align="right"> bio18 </TD> <TD align="right"> 1417.65 </TD> <TD align="right"> 201.31 </TD> <TD align="right"> 1618.96 </TD> </TR>
  <TR> <TD align="right"> bio19 </TD> <TD align="right"> 1233.77 </TD> <TD align="right"> 196.68 </TD> <TD align="right"> 1430.45 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 1363.23 </TD> <TD align="right"> 180.92 </TD> <TD align="right"> 1544.15 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 1403.98 </TD> <TD align="right"> 203.76 </TD> <TD align="right"> 1607.73 </TD> </TR>
   </TABLE>


###

#### Table S3. Bias, Variance and Mean Square Error of the Mid Domain, Neutral and Logistic simulation models for the jaccard pairwise similarity. The Mean Square Error is the sum of the Bias and Variance.

<TABLE border=1, bgcolor=#989898>
<TR> <TH>  </TH> <TH> BIASsq </TH> <TH> VAR </TH> <TH> sMSE </TH>  </TR>
  <TR> <TD align="right"> MidD </TD> <TD align="right"> 8.83 </TD> <TD align="right"> 1.88 </TD> <TD align="right"> 10.71 </TD> </TR>
  <TR> <TD align="right"> Neutral </TD> <TD align="right"> 7.07 </TD> <TD align="right"> 9.43 </TD> <TD align="right"> 16.51 </TD> </TR>
  <TR> <TD align="right"> Long2 </TD> <TD align="right"> 4.36 </TD> <TD align="right"> 2.16 </TD> <TD align="right"> 6.52 </TD> </TR>
  <TR> <TD align="right"> Lat2 </TD> <TD align="right"> 4.20 </TD> <TD align="right"> 2.15 </TD> <TD align="right"> 6.35 </TD> </TR>
  <TR> <TD align="right"> habitat </TD> <TD align="right"> 5.95 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 8.20 </TD> </TR>
  <TR> <TD align="right"> vegetation </TD> <TD align="right"> 5.47 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 7.73 </TD> </TR>
  <TR> <TD align="right"> bio1 </TD> <TD align="right"> 4.92 </TD> <TD align="right"> 2.23 </TD> <TD align="right"> 7.15 </TD> </TR>
  <TR> <TD align="right"> bio2 </TD> <TD align="right"> 5.90 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 8.14 </TD> </TR>
  <TR> <TD align="right"> bio3 </TD> <TD align="right"> 4.75 </TD> <TD align="right"> 2.21 </TD> <TD align="right"> 6.96 </TD> </TR>
  <TR> <TD align="right"> bio4 </TD> <TD align="right"> 4.32 </TD> <TD align="right"> 2.16 </TD> <TD align="right"> 6.49 </TD> </TR>
  <TR> <TD align="right"> bio5 </TD> <TD align="right"> 6.28 </TD> <TD align="right"> 2.27 </TD> <TD align="right"> 8.55 </TD> </TR>
  <TR> <TD align="right"> bio6 </TD> <TD align="right"> 5.06 </TD> <TD align="right"> 2.22 </TD> <TD align="right"> 7.28 </TD> </TR>
  <TR> <TD align="right"> bio7 </TD> <TD align="right"> 5.13 </TD> <TD align="right"> 2.18 </TD> <TD align="right"> 7.30 </TD> </TR>
  <TR> <TD align="right"> bio8 </TD> <TD align="right"> 4.71 </TD> <TD align="right"> 2.15 </TD> <TD align="right"> 6.87 </TD> </TR>
  <TR> <TD align="right"> bio9 </TD> <TD align="right"> 5.85 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 8.11 </TD> </TR>
  <TR> <TD align="right"> bio10 </TD> <TD align="right"> 5.99 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 8.24 </TD> </TR>
  <TR> <TD align="right"> bio11 </TD> <TD align="right"> 4.46 </TD> <TD align="right"> 2.21 </TD> <TD align="right"> 6.67 </TD> </TR>
  <TR> <TD align="right"> bio12 </TD> <TD align="right"> 5.78 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 8.02 </TD> </TR>
  <TR> <TD align="right"> bio13 </TD> <TD align="right"> 5.54 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 7.78 </TD> </TR>
  <TR> <TD align="right"> bio14 </TD> <TD align="right"> 5.29 </TD> <TD align="right"> 2.21 </TD> <TD align="right"> 7.51 </TD> </TR>
  <TR> <TD align="right"> bio15 </TD> <TD align="right"> 5.08 </TD> <TD align="right"> 2.24 </TD> <TD align="right"> 7.33 </TD> </TR>
  <TR> <TD align="right"> bio16 </TD> <TD align="right"> 5.59 </TD> <TD align="right"> 2.25 </TD> <TD align="right"> 7.84 </TD> </TR>
  <TR> <TD align="right"> bio17 </TD> <TD align="right"> 5.35 </TD> <TD align="right"> 2.22 </TD> <TD align="right"> 7.57 </TD> </TR>
  <TR> <TD align="right"> bio18 </TD> <TD align="right"> 5.84 </TD> <TD align="right"> 2.27 </TD> <TD align="right"> 8.11 </TD> </TR>
  <TR> <TD align="right"> bio19 </TD> <TD align="right"> 5.66 </TD> <TD align="right"> 2.23 </TD> <TD align="right"> 7.89 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.1 </TD> <TD align="right"> 4.64 </TD> <TD align="right"> 2.20 </TD> <TD align="right"> 6.84 </TD> </TR>
  <TR> <TD align="right"> PCA.wclim.2 </TD> <TD align="right"> 5.77 </TD> <TD align="right"> 2.26 </TD> <TD align="right"> 8.02 </TD> </TR>
   </TABLE>


###


----------------------------











