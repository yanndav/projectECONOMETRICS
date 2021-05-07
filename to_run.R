########################################
# ECONOMETRICS PROJECT - MASTER SCRIPT #
#       YANN DAVID - ANDRE SANTOS      #
########################################



# 00. LOADING FUNCTIONS & PACKAGES ----------------------------------------
# sets wd as the path in which script is
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 
source("functions.R")


# 01. DATA GENERATING PROCESS ---------------------------------------------
data_homosk = dgp_homoskedastic(G=5)
data_heterosk = dgp_heteroskedastic(G=5)


# 02. RESULTS REPLICATION -------------------------------------------------


# 03. EXTENSION WITH SKEWED ESTIMATOR -------------------------------------

