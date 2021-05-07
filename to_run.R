########################################
# ECONOMETRICS PROJECT - MASTER SCRIPT #
#       YANN DAVID - ANDRE SANTOS      #
########################################



# 00. LOADING FUNCTIONS & PACKAGES ----------------------------------------
# sets wd as the path in which script is
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 
source("functions.R")


# 01. DATA GENERATING PROCESS ---------------------------------------------
data_homosk = data =  dgp_homoskedastic(G=150,N_g=100)
data_heterosk =data =  dgp_heteroskedastic(G=15,N_g=30)


# 02. RESULTS REPLICATION -------------------------------------------------
x_ig = data_heterosk$x
y_ig = data_heterosk$y

H_beta = solve(t(x_ig)%*%x_ig) %*%t(x_ig)%*%y_ig
lm(y~x-1,data=data_homosk)

# 03. EXTENSION WITH SKEWED ESTIMATOR -------------------------------------

