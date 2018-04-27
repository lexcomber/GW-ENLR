# GW-ENLR: Geographically weighted elastic net logistic regression - code and data
The code and data are all described in the file `Comber_Harris_GW-ENLR_code_simulation.R`. To run the simulation data you will need to download the zip file from this repository and unzip the files to a directory of your choice. This should then be set as the working directory in R using the `setwd()` function. 

Please contact Lex Comber [a.comber@leeds.ac.uk](a.comber@leeds.ac.uk) if you have any questions.

## Paper title: Geographically weighted elastic net logistic regression
Alexis Comber<sup>1</sup> and Paul Harris<sup>2</sup>

<sup>1</sup>School of Geography, University of Leeds, Leeds, UK LS2 9JT\
<sup>2</sup>Rothamsted Research, North Wyke, Okehampton, Devon, UK EX20 2SB

## Abstract
This paper develops and applies a localized approach to elastic net logistic regression. It uses a geographically weighted framework to subset and weight data for case studies to predict a) whether counties in the USA voted for Trump in the 2016 presidential elections using socio-economic data, and b) species presence-absence from environmental data. The rationale for the approach relates to: the efficiency of techniques such as penalized regression at model selection and information reduction and overcoming predictor variable collinearity; the lack of sensitivity to and obfuscation of local trends in the data by global, whole map models; and the potential for high local collinearity amongst predictor variables in local data subsets in geographically weighted models even when none is observed globally. Collinearity can result in a loss of precision and power in the coefficient estimates leading to poor inferences. The GW-ELNR was compared with other regression models and improved prediction for the election case study only. This exhibited much greater spatial heterogeneity in the binary response than the species case study. However, model comparisons for both case studies showed that standard geographically weighted logistic regression over-estimated relationship non-stationarity because it failed to adequately deal with collinearity and model selection. Mapping the GW-ELNR coefficient estimates and indicating the locations of local collinearity amongst predictor variables provides a deeper understanding of the data structure and how best to develop a potential predictive model by. It also confirms or not perceptions of relationship non-stationarity. GW-ELNR was evaluated through a simulation experiment which endorsed this value of the new model. On-going work is investigating locally-derived local elastic net parameters.


April 2018.
