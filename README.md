# Zero-Inflated Poisson INGARCH Model for Earthquake Prediction

## Introduction

This repository contains the R scripts and data for the implementation and evaluation of an innovative Zero-Inflated Poisson Integer-Valued Generalized Autoregressive Conditional Heteroskedasticity (ZIP-INGARCH) model designed for analyzing and predicting the monthly frequency of earthquakes with a magnitude above 5.

## Project Abstract

Earthquakes are natural disasters with the potential to impact millions globally. Analyzing earthquake data is paramount to understanding seismic mechanisms and enhancing public safety. This project adopts the ZIP-INGARCH model to address the discrete, sporadic, and zero-inflated nature of major earthquake data. The Bayesian approach, leveraging a Markov Chain Monte Carlo (MCMC) algorithm, is employed for parameter estimation.

Various orders of ZIP-INGARCH models are fitted to the data, with subsequent one-step-ahead point and interval predictions made for earthquake occurrences. Predictive accuracy is assessed using the Average Absolute Error and the Average Relative Error metrics. The ZIP-INGARCH(2,2) model demonstrates superior performance, with the median of its predictive distribution achieving an Average Absolute Error of 0.125 and an Average Relative Error of 0.0625. Additionally, the prediction intervals for all models reliably encompass the actual observed values, affirming the precision of the proposed forecasting method.

## Keywords

- Earthquake Prediction
- INGARCH Model
- Zero-Inflated Poisson Distribution
- Bayesian Method
- MCMC

## Repository Contents

- `MCMC_InGarch_independent_sampling_order(1,1).R`: This R script meticulously applies the Zero-Inflated Poisson Integer-Valued GARCH (ZIP-INGARCH) model of order (1,1) to earthquake data. It includes MCMC algorithm implementation for parameter estimation, diagnostic checks for convergence, forecasting future earthquake occurrences with point predictions, and evaluating the forecast accuracy against actual data.

- `Forcasting_Model_Evaluation.R`: R script for evaluating earthquake frequency forecasts. It includes data preprocessing, convergence diagnostics, parameter estimation for the ZIP INGARCH(1,1) model, and generating point and interval forecasts using various statistical measures..
- `/Source Data`: Directory containing the datasets used in the models.
- `/Forecasting Results`: Directory with the output results and predictions from the models.

## Models and Predictions

The repository includes scripts for:

- Fitting ZIP-INGARCH models of different orders to earthquake frequency data.
- Generating one-step-ahead point predictions and prediction intervals.
- Computing error metrics to evaluate predictive precision.




## Contact

For any queries regarding this project, please open an issue in the repository or contact Vincent Li at liwenhao.vincent@gmail.com.

---

*This README is a brief overview of the comprehensive earthquake prediction research documented in the accompanying thesis.*
