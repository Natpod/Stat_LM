# Liquid Biopsy Data PDAC Logistic Model

*Author : Natalia García Sánchez*
Student in Master in Computational Biology

The code in this repository `DAT4LogRRoc.R`
1. Performs simple analysis with 2x2 contingency table for categorical variables describing presence of biomarkers/clinical groups and Mortality categorical variable `status`
2. Performs risk analysis on a single variable `Mutation`over `status`
3. Trains a Logistic Model with 3-Fold cross validation and plots the ROC curve

Package Dependencies : MASS, ROCR, pROC, caret, skimr, dplyr
