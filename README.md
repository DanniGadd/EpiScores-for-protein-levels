# Epigenetic scores for the circulatin proteome as tools for disease prediction

Published: eLife 2022;11:e71802 DOI: 10.7554/eLife.71802

EpiScores trained on 953 protein levels in the Scottish LBC1936 and German KORA cohorts. 109 EpiScores that explain between 1 and 58% of the variance in plasma protein levels are then projected into over 9000 individuals in the Generation Scotland cohort and related to the onset of 12 leading causes of morbidity.

The code follows the following naming order:
1) "prep_" - these are the preparations of the files for inputs into elastic net models and testing of models
2) "run_" - these are the elastic net penalised regressions run for LBC1936 and KORA protein datasets 
3) "project_" - these are projections of trained CpG weights for proteins into various testing cohorts 
4) "validation_" - this refers to the collation of which episcores passed robust thresholds in STRADL & LBC1921 test sets
5) "cox_" - this refers to the cox modelling between episcores and incident disease onset in Generation Scotland
6) "process_" - these scripts have been used to process and plot the results
7) "covid_" - these scripts detail the covid-19 specific analyses in the manuscript 

This project is licensed under the terms of the MIT license.
