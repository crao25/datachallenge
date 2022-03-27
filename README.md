# Data Challenge
This repository contains the scripts to perform exploratory analysis on clinical and expression data of newly diagnosed Multiple Myeloma patients and to build a predictive model to predict the risk of death or relapse in these patients. 

1. data_challenge.R: R script to do exploratory and predictive analyis
2. data_challenge.Rmd: rmarkdown script which is a replicate of the R script, however with more text and comments.
3. data_challenge.html: Output of the rmarkdown script
4. Data_Challenge_Report.pdf: Short report to intereprt the model and the analysis results
5. pred_model.Rds: Final model trained on entire training dataset and which can be loaded and applied on an external dataset

Guideline on how to use the model pred_model.Rds on an external validation set:

* Features list: 
1. Clinical and cytogenetic features: "D_Age", "D_ISS", "CYTO_predicted_feature_01", "CYTO_predicted_feature_03", 
"CYTO_predicted_feature_14"
2. Expression data features: Entrez IDs: 71, 472, 595, 673, 1387, 1540, 1958, 2034, 2261, 3008, 3662, 3845, 4149, 4893, 7157, 7187, 7468, 8289, 9055, 10075, 22894, 25865, 26059, 26147, 51366, 54855, 55790, 58508

* Target variables: D_PFS, HR_FLAG 

Note: The expression dataset is reshaped so that the Entrez IDS or the gene names become column names and the patient IDs are values of a separate column and then the dataset is merged with the clinical dataset using the patient IDs. 


