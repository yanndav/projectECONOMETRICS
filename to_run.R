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
data_homo_sk = dgp_homoskedastic_skewed(G=10,
                                              N_g=30,
                                              mu=0.5,
                                              sigma=1)
data_heterosk_sk = dgp_heteroskedastic_skewed(G=10,
                                              N_g=30,
                                              mu=0.5,
                                              sigma=1)

hist(data_homo_sk$x)

# 02. RESULTS REPLICATION -------------------------------------------------
k=1

reject_it = c()
reject_boot_it = c()

for (it in 1:num) {
  data =  dgp_homoskedastic(G=5,N_g=5)
  
    for (b in 1:B) {
    # Random sample of the data, with replacement
    index_b <- sample(N_g,N_g,replace=TRUE)
    data_b <- do.call(rbind,lapply(1:G,function(k){
      t = data[data$g==k ,]
      return(t[index_b,])}))

   
    
    # Wald 
    temp = theta_hat - theta_null
    Wald_t =  t(temp)%*%J1(y_b,x_b,theta_hat)%*%temp # Wald
    wald_vec = c(wald_vec,Wald_t)
    
  }
  
  wald_vec <- sort(wald_vec)
  cv_b <- wald_vec[floor(B*0.95)]
  
  if (is.nan(Wald) == 0) {
    if (Wald > cv_b) {
      reject_boot_it = c(reject_boot_it,1)
    }else{
      reject_boot_it = c(reject_boot_it,0)
    }
  }
}

mean(reject_boot_it)
mean(reject_it)


waldCR(data_b,1)>1.96

varCR(data)
sqrt(varJACK(data))

reg = lm(y~x-1,data=data)
library(lmtest)
library("clubSandwich")
reg2=coef_test(reg, df = data, vcov = vcovCR(reg, type = "CR3",cluster = data$g))
reg2$SE
# 03. EXTENSION WITH SKEWED ESTIMATOR -------------------------------------

