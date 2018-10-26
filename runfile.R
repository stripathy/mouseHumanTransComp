# run file

# first download datasets from AIBS portal
source('R/loadAllenExprData.R')

# we then need to map the mouse nuclei to clusters
source('R/mapCellsToClusters.R')

# go through analysis/mouse_v_human_analysis.R - this plots the mouse v human scatter plots and boxplots for indiv gene families
