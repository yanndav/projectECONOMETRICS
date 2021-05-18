########################################
# ECONOMETRICS PROJECT - MASTER SCRIPT #
#       YANN DAVID - ANDRE SANTOS      #
########################################



# 00. LOADING FUNCTIONS & PACKAGES ----------------------------------------
# sets wd as the path in which script is
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 
source("functions.R")
if(!("clubSandwich"%in%installed.packages()[,'Package'])){
  install.packages('clubSandwich')
}
library("clubSandwich")


# 01. results_hoc REPLICATION -------------------------------------------------
# Setting parameters
mCarlo = 100 # Number of Monte Carlo Iteration
bStrap= 399 # Number of Boostrap iterations
clusters = c(5,10 ,15,20,25,30
             ) # Number of clusters to be tested
N_g = 30 # number of observation per cluster
beta_null = 1

# HOMOSKEDACTIC DGP
# Creating empty list if not yet in the environment
if(!exists("results_hoc")){
  results_hoc = list()
}

estimators_wald = c('ols','ols_Bt','crve','crve_Bt','cr3_Bt','cr3',
                     'resi_Bt','wild_Bt')

# results_hoc to be estimated:
for (G in clusters){
  
  # Initialising results_hoc storage if empty list
  if(length(results_hoc[[paste(G)]])==0){
    for(method in estimators_wald){
      results_hoc[[paste(G)]][[method]] = list("rejection" = c(),
                                              "wald" = c())
    }
    
  } # If not empty we want to add more iterations to improve our accuracy !
  
  
  for (it in 1:mCarlo) {# Start of Monte Carlo iteration
    print(it)
    set.seed(it)
    # Random data generation process
    data = dgp_homoskedastic(G=G,N_g=N_g)
    thresh = qchisq(0.95,df=1) # Significance level for rejection rate
    
    
    # Estimation of regression:
    reg = lm(y~x-1,data=data)
    beta_null_b = betaH(data)
    
    # Basic ols variance
    wald_OLS = wald(reg)
    results_hoc[[paste(G)]][["ols"]][["wald"]] = c(results_hoc[[paste(G)]][["ols"]][["wald"]],
                                               wald_OLS)
    results_hoc[[paste(G)]][["ols"]][["rejection"]] = c(results_hoc[[paste(G)]][["ols"]][["rejection"]],
                                               as.integer(wald_OLS>thresh))
    # CRVE variance estimator
    wald_CRVE = wald(reg, cluster = "CR1",data=data)
    results_hoc[[paste(G)]][["crve"]][["wald"]] = c(results_hoc[[paste(G)]][["crve"]][["wald"]],
                                                wald_CRVE)
    results_hoc[[paste(G)]][["crve"]][["rejection"]] = c(results_hoc[[paste(G)]][["crve"]][["rejection"]],
                                                as.integer(wald_CRVE>thresh))
    # CR3 variance estimator
    wald_CR3 = wald(reg, cluster = "CR3",data=data)
    results_hoc[[paste(G)]][["cr3"]][["wald"]] = c(results_hoc[[paste(G)]][["cr3"]][["wald"]],
                                               wald_CR3)
    results_hoc[[paste(G)]][["cr3"]][["rejection"]] = c(results_hoc[[paste(G)]][["cr3"]][["rejection"]],
                                               as.integer(wald_CR3>thresh))
    
    # Restricted beta and associated residuals
    restrict = restricted_OLS(data)
    data$ur = restrict[['residuals_r']]
    beta_r = restrict[['beta_r']]
    # Initialisation of storage for monte carlo
    boot_ols = c()
    boot_crve = c()
    boot_cr3 = c()
    boot_resi = c()
    boot_wild = c()
    
   
    
    for (b in 1:bStrap) {
      set.seed(it*b)
      # Random sample of the clusters, with replacement
      
      index_b <- sample(G,G,replace=TRUE)
      data_b <- do.call(rbind,lapply(index_b,function(k){
        return(data[data$g==k ,])}))
  
      while(length(unique(data_b$g))==1){
        index_b <- sample(G,G,replace=TRUE)
        data_b <- do.call(rbind,lapply(index_b,function(k){
          return(data[data$g==k ,])}))
      }
      # Estimation of regression:
      reg = lm(y~x-1,data=data_b)
      
      # Basic ols variance
      boot_ols = c(boot_ols, wald(reg, beta_null = beta_null_b,data=data_b))
      
      # CRVE variance estimator
      boot_crve = c(boot_crve,
                    wald(reg, beta_null = beta_null_b,cluster = "CR1",data=data_b))
      
            # CR3 variance estimator
      boot_cr3 = c(boot_cr3,
                     wald(reg, beta_null = beta_null_b,cluster = "CR3",data=data_b))
    
      
      # Residuals boostrap
      data_r = data.frame(x=data$x,
                          g=data$g,
                          y=beta_r*data$x+data_b$ur)
      
      reg_re = lm(y~x-1,data=data_r)
      boot_resi = c(boot_resi,
                    wald(reg_re, cluster = "CR1",data=data_r))
      
      
      # Wild boostrap
      wild_p = rep(sample(c(1,-1), G,
             prob = c(0.5,0.5), replace=T),
             each=N_g)
      
      data_r2 = data.frame(x=data$x,
                           g=data$g,
                           y=beta_r*data$x + wild_p * data$ur)
      reg_r2 = lm(y~x-1,data=data_r2)
      
      boot_wild = c(boot_wild,
                    wald(reg_r2, cluster = "CR1",data=data_r2))
   
  }
  # Asymptotic refinement moment
    # Ols variance:
    boot_ols <- sort(boot_ols)
    cv_boot_ols <- boot_ols[floor(bStrap*0.95)]
    results_hoc[[paste(G)]][["ols_Bt"]][["wald"]] = c(results_hoc[[paste(G)]][["ols_Bt"]][["wald"]],
                                                  cv_boot_ols)
    if (is.nan(wald_OLS) == 0) {
      results_hoc[[paste(G)]][["ols_Bt"]][["rejection"]] = c(results_hoc[[paste(G)]][["ols_Bt"]][["rejection"]],
                                                    as.integer(wald_OLS>cv_boot_ols))
    }
    
    # CRVE variance:
    boot_cvre <- sort(boot_crve)
    cv_boot_crve <- boot_cvre[floor(bStrap*0.95)]
    results_hoc[[paste(G)]][["crve_Bt"]][["wald"]] = c(results_hoc[[paste(G)]][["crve_Bt"]][["wald"]],
                                                  cv_boot_crve)
    if (is.nan(wald_CRVE) == 0) {
      results_hoc[[paste(G)]][["crve_Bt"]][["rejection"]] = c(results_hoc[[paste(G)]][["crve_Bt"]][["rejection"]],
                                                         as.integer(wald_CRVE>cv_boot_crve))
    }
    
    # CR3 variance:
    boot_cr3 <- sort(boot_cr3)
    cv_boot_cr3 <- boot_cr3[floor(bStrap*0.95)]
    results_hoc[[paste(G)]][["cr3_Bt"]][["wald"]] = c(results_hoc[[paste(G)]][["cr3_Bt"]][["wald"]],
                                                  cv_boot_cr3)
    if (is.nan(wald_CRVE) == 0) {
      results_hoc[[paste(G)]][["cr3_Bt"]][["rejection"]] = c(results_hoc[[paste(G)]][["cr3_Bt"]][["rejection"]],
                                                          as.integer(wald_CR3>cv_boot_cr3))
    }
    
    # Resi variance:
    boot_resi <- sort(boot_resi)
    cv_boot_resi <- boot_resi[floor(bStrap*0.95)]
    results_hoc[[paste(G)]][["resi_Bt"]][["wald"]] = c(results_hoc[[paste(G)]][["resi_Bt"]][["wald"]],
                                                       cv_boot_resi)
    if (is.nan(wald_CRVE) == 0) {
      results_hoc[[paste(G)]][["resi_Bt"]][["rejection"]] = c(results_hoc[[paste(G)]][["resi_Bt"]][["rejection"]],
                                                             as.integer(wald_CRVE>cv_boot_resi))
    }
    
    # Wild variance:
    boot_wild <- sort(boot_wild)
    cv_boot_wild<- boot_wild[floor(bStrap*0.95)]
    results_hoc[[paste(G)]][["wild_Bt"]][["wald"]] = c(results_hoc[[paste(G)]][["wild_Bt"]][["wald"]],
                                                       cv_boot_wild)
    if (is.nan(wald_CRVE) == 0) {
      results_hoc[[paste(G)]][["wild_Bt"]][["rejection"]] = c(results_hoc[[paste(G)]][["wild_Bt"]][["rejection"]],
                                                              as.integer(wald_CRVE>cv_boot_wild))
    }
    
  }
}

saveRDS(object = results_hoc, 'results_hoc_pairs_t.RDS')

# 03. EXTENSION WITH SKEWED ESTIMATOR -------------------------------------

