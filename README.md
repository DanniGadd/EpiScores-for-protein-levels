# EpiScores for protein levels

EpiScores trained on 953 protein levels in LBC1936 and KORA. Robust EpiScores are then projected into Generation Scotland and related to the onset of 12 leading causes of morbidity.

The code follows the following naming order:
1) "prep_" - these are the preparations of the files for inputs into elastic net models and testing of models
2) "run_" - these are the elastic net penalised regressions run for LBC1936 and KORA protein datasets 
3) "project_" - these are projections of trained CpG weights for proteins into various testing cohorts 
4) "validation_" - this refers to the collation of which episcores passed robust thresholds in STRADL & LBC1921 test sets
5) "cox_" - this refers to the cox modelling between episcores and incident disease onset in Generation Scotland
6) "process_" - these scripts have been used to process and plot the results
