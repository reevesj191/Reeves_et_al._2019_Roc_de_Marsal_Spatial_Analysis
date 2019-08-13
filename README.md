# Reeves et al. 2019 Roc de Marsal Spatial Analysis



This repository contains the code and markdown document (with supporting files) to carry out the moving window analysis previously published in the Journal of Archaeological and Anthropological Sciences (hyperlink fourth coming). The Rmarkdown document (Reeves_et_al._-_RDM_Time_Averaging.Rmd) is the paper itself and when knitted in R will produce a PDF of the paper. This document shows how the analysis can be used to start to understand time-averaged patterns of hominin space use in caves. The R document (RDM_Analysis.R) is the source script for the moving window functions used in the paper. There is currently no user-manual for these functions but the document provides comments on how the functions work. 

If you find problems with the code, I am happy to fix any issues as they arise. Please either contact me or you can issue a pull request. 

Please cite Reeves et al 2019 if the code is used in publications. 

### Using this code

Ensure that you have all of the packages listed in each of the markdown documents and scripts to execute the code used to run our analysis and reproduce our results.  If you wish to reproduce this paper, please note that the markdown documents make use of the 'rticles' package to format the outputted PDF file. This package calls out to latex to do alot of the heavily lifting when producing the PDF document. If you encounter errors please make sure you have a version of latex installed on your computer. 

### A note on parallel processing

The moving window function has the capacity to  run using parallelization. Take care when using these features as a misallocation of cores can quickly result maxing out the computers RAM or processing capacity. This ultimate results in crashes. 

### Bugs/Feature Fixes

I await your feedback 







