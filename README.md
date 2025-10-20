# NonLinearRC

This repository contains files for our R package "NonLinearRC", which applies regression calibration to estimate parameters in a Cox regression model when the exposure variable is measured with error, and has a nonlinear effect on the hazard function

> **Estimating Non-Linear Exposure and Event Time Assocication in the Presence of Exposure Measurement Error**  
> Authors: Trevor J Thomson, and Ying Huang
> Submitted to "The New England Journal of Statistics in Data Science"

## 📌 Overview
This repository provides:
- Functions used to simulate data and estimate model parameters
- Code to illustrate our estimation procedure with a simulated dataset

## 🛠️ Requirements
The code was written in **R** and uses the following packages:
- `MASS`
- `pbapply`
- `pbapply`
- `Rcpp`
- `RcppArmadillo`
- `stringr`
- `numDeriv`
- `cubature`
- `survival`
- `splines2`

The **R** script calls **C++** files to execute our proposed estimation strategy.

## 📧 Contact
For questions, feedback, or collaboration inquiries, feel free to reach out to:

👤 Trevor Thomson

📍 Fred Hutchinson Cancer Center

📬 Email: tthomson@fredhutch.org
