# Ecotoxicity prediction of chemical compounds using machine learning and different molecular structure representations

This repository contains scripts for prediction HC50 values (hazardous concentration) using a QSPR model. This model is built in R, with RDKit used for computing molecular descriptors via reticulate.  
The best model was trained using 100 moleculars descriptors scaled between 0 and 1, MACCS keys as molecular fingerprints and the `xgbTree` algorithm.

# 1. Setup
## 1.1. Install software
Ensure you have the following installed on your system:
- Anaconda 23.7.4
- R version 4.3.2 (2023-10-31 ucrt)
- RStudio 2024.04.2+764

## 1.2. Create Conda Environment with RDKit
1. Open Anaconda Prompt and run:
  > conda create -n rdkit-env python=3.10 rdkit -c conda-forge
2. Activate the environment:
  > conda activate rdkit-env

The installed version of RDKit is 2022.03.05.

# 2. Required R packages
## 2.1 Run the following in RStudio:
> install.packages("webchem", version = "1.3.0")  
> install.packages("readxl", version = "1.4.3")  
> install.packages("reticulate", version = "1.34.0")  
> install.packages("caret", version = "6.0.94")  
> install.packages("dplyr", version = "1.1.4")

# 3. Set up the Python path in R scripts
## 3.1. Update Python Path
Before running the scripts, update the Python path in `2_compute_descriptors.R` and `4_predict_HC50.R`:
>  use_python("path\\to\\python.exe")

# 4. Running the scripts
## 4.1. Overview
The project consists of two pipelines:  
  1. Model training pipeline:  
    - `1_cas_to_smiles.R`: converts CAS numbers to SMILES notation  
    - `2_compute_descriptors.R`: computers molecular descriptors and MACCS keys  
    - `3_train_model.R`: trains the QSPR model  
  
  2. Prediction Pipeline:  
    - `4_precict_HC50.R`: uses the best trained model to predict HC50 value of a user provided SMILES string

## 4.2. Input SMILES
In `4_predict_HC50.R` find the following section
> #Enter desired SMILES notation (without space)  
smiles_codes <- "CCCCCC"

Replace "CCCCCC" with your molecule of interest.

## 4.3. Running the script
In RStudio run:  
> source("4_predict_HC50.R")  

# 5. Dataset
The dataset used for model training is is sourced from the following publication:  
**doi:** [10.1021/acssuschemeng.0c03660](https://doi.org/10.1021/acssuschemeng.0c03660)  
The original dataset is available as **Supporting Information**.  
