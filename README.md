## Replication_MBSTS_TL

This is the replication code for Paper Statistical Machine Learning Meets High-Dimensional Spatiotemporal Challenges - A Case Study of COVID-19 Modeling.

## File Description

**`Real_data_comparison.R`** includes the source code for implementing mbsts-tl approach and bsts approach.

**`Real_data_comparison.R`** includes the source code for providing additional function in R package 'mbsts' to analyze the mbsts model in covid data.

**`mbsts_tl.R`** provides an R function for implementing mbsts-tl in our paper with some examples.

**`Alldata_weeklyAvg_state_2020.xlsx`** is the covid data.

## Details for 'mbsts_tl' function

Here we provide more details of mbsts-tl method in the paper.

**Description**

**Input**

-   'X': A (n\*K)-dimensional matrix containing all candidate predictor series for each target series. K=∑ k_i is the number of all candidate predictors for all target series. The first k_1 variables are the set of candidate predictors for the first target series, and the next k_2 variables are the set of candidate predictors for the second target series, etc.

-   'Y': A (n\*m)-dimensional matrix containing multiple target series, where n is the number of observations and m is the number of target series. Note that usually m should not be too large.

-   'array.number': A vector of integer values denoting the acumulated number of predictors for target series. For example, if there are three target series where the first has 8 predictors, the second has 6 predictors, and the third has 10 predictors, then the vector is c(8,14,24).

-   'rho.space/season.space/vrho.space/lambda.space': candidate for the hyperparameter in mbsts model, more details for rho, season, vrho and lambda can be found in mbsts_function in R package mbsts.

-   'time_begin/time_end': the starting time and ending time for each time period.

**output**

-   'error': error is a vector for calculating normalized absolute error for each time period.

-   'optimal.parameter': optimal.parameter is a m\*4-dimension matrix that the optimal hyperparameter we use in mbsts model, where m is the number of period we explore.

-   'beta': beta is a m\*dimension matrix where dimension is the dimension of covariate for all time series.

## Reference
