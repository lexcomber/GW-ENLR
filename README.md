# GW-ENLR: Geographically weighted elastic net logistic regression - code and data
This analysis is described in a paper soon to be submitted. The code and data are all described in the file `Comber_Harris_GW-ENLR_code.R`. Some of the file paths (`setwd()`) will have to be changed to reflect your local directory.

Please contact Lex Comber [a.comber@leeds.ac.uk](a.comber@leeds.ac.uk) if you have any questions.

## Paper title: Geographically weighted elastic net logistic regression
Alexis Comber<sup>1</sup> and Paul Harris<sup>2</sup>

<sup>1</sup>School of Geography, University of Leeds, Leeds, UK LS2 9JT\
<sup>2</sup>Rothamsted Research, North Wyke, Okehampton, Devon, UK EX20 2SB

## Abstract
This paper proposes and applies a Geographically Weighted Elastic Net Logistic Regression (GW-ENLR) model. It extends previous research that has suggested the notion of a localized elastic net, as an extension to a localized ridge regression or a localized lasso.  All such models have the objective to capture data relationships that vary across space, where this study’s elastic net extension allows for a robust local model selection component that also alleviates local collinearity. The paper applies GW-ENLR to two case studies: patterns of voting at the county level in the 2016 USA presidential election, specifically examining the spatial structure of socio-economic factors associated with voting for Trump; and a species presence-absence dataset linked to explanatory environmental and climatic factors at locations on regular grid covering mainland USA. GW-ENLR is used to construct predictive classification models for the two case studies and the results are compared with those from an a-spatial logistic regression, an a-spatial elastic net logistic regression and a geographically weighted logistic regression. In all cases, classification rates were improved under the new GW-ENLR model. Results are discussed in the context of predictor variable collinearity, model (predictor variable) selection and the relationship heterogeneities that are present in data.

July 2017.
