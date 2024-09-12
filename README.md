# Ridge and Lasso Regression Analysis

This repository contains an implementation of Ridge and Lasso regression models for analyzing the relationship between OTU data and indoxyl sulfate levels, along with a binomial classification model for GCB and ABC subtypes.

## Repository Structure

- input data files.
  - `OTU_data_spike.rds`: OTU data used for analysis.
  - `indoxylSulfate.csv`: Indoxyl sulfate data.
  - `mmml_vsn.rda`: Molecular data for classification.
- `scripts/`: Contains R scripts for the analysis.
  - `analysis.R`: The main analysis script implementing Ridge, Lasso, and cross-validation techniques.

## Instructions

1. Clone the repository:

2. Install the necessary R packages:
```R
install.packages("glmnet")
install.packages("ROCR")
