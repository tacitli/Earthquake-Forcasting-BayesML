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

- `ZIP_INGARCH_Modeling.R`: R script implementing the ZIP-INGARCH model fitting and predictions.
- `Earthquake_Prediction_Analysis.R`: R script for data analysis and visualization of earthquake frequency predictions.
- `/data`: Directory containing the datasets used in the models.
- `/results`: Directory with the output results and predictions from the models.

## Models and Predictions

The repository includes scripts for:

- Fitting ZIP-INGARCH models of different orders to earthquake frequency data.
- Generating one-step-ahead point predictions and prediction intervals.
- Computing error metrics to evaluate predictive precision.

## Usage

To run the analysis and predictions:

1. Ensure R and the necessary packages (`coda`, `rjags`, etc.) are installed.
2. Set the working directory to your local path where the scripts are stored.
3. Run the `ZIP_INGARCH_Modeling.R` and `Earthquake_Prediction_Analysis.R` scripts sequentially.

## Contributing

Contributions to this project are welcome. To contribute:

1. Fork the repository.
2. Create a new branch for your feature.
3. Commit your changes.
4. Push to the branch.
5. Create a new Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any queries regarding this project, please open an issue in the repository or contact [Your Name/Contact Information].

---

*This README is a brief overview of the comprehensive earthquake prediction research documented in the accompanying thesis.*
