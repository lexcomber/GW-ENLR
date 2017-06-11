# GW-ENLR: Geographically weighted elastic net logistic regression - code and data
This analysis is described in a paper soon to be submitted. The code and data are all described in the file `Comber_Harris_GW-ENLR_code.R`. Some of the file paths (`setwd()`) will have to be changed to reflect your local directory.

Please contact Lex Comber [a.comber@leeds.ac.uk](a.comber@leeds.ac.uk) if you have any questions.

## Paper title: Geographically weighted elastic net logistic regression
Alexis Comber<sup>1</sup> and Paul Harris<sup>2</sup>

<sup>1</sup>School of Geography, University of Leeds, Leeds, UK LS2 9JT\
<sup>2</sup>Rothamsted Research, North Wyke, Okehampton, Devon, UK EX20 2SB

## Abstract
This paper proposes and applies a Geographically Weighted Elastic Net Logistic Regression (GW-ENLR) model. It extends previous research that has suggested the notion of geographically weighted elastic net, by including the case of logistic or generalized linear regression. It applies GW-ENLR to 2 case studies: patterns of voting at the county level in the 2016 USA presidential election, specifically examining the spatial structure of socio-economic factors associated with voting for Trump; and a species presence-absence dataset, linked to explanatory environmental and climatic factors at locations on regular grid covering mainland USA. GW-ENLR is used to construct predictive models for the 2 case studies and the results are compared with the results from standard generalized linear regression (GLM), generalized geographically weighted regression (GGWR) and a standard a-spatial elastic net logistic regression model. In all case studies prediction rates were improved under the GW-ENLR. The results are discussed in the context of covariate collinearity, model selection and commonly observed assumptions in statistical modelling when spatial autocorrelation and heterogeneity are present in data.

June 2017.
