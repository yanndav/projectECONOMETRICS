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
clusters = c(5,10,15,20,25,30) # Number of clusters to be tested
thresh = 1.96 # Significance level for rejection rate
N_g = 30 # umber of observation per cluster
# HOMOSKEDACTIC DGP
# Creating empty list if not yet in the environment
if(!exists("results_hoc")){
  results_hoc = list()
}

# results_hoc to be estimated:
for (G in clusters){
  # Initialising results_hoc storage if empty list
  if(length(results_hoc[[paste(G)]])==0){
    results_hoc[[paste(G)]][["ols"]] = list("rejection" = c(),
                                        "wald" = c())
    results_hoc[[paste(G)]][["ols_Bt"]] = list("rejection" = c(),
                                           "wald" = c())
    results_hoc[[paste(G)]][["crve"]] = list("rejection" = c(),
                                        "wald" = c())
    results_hoc[[paste(G)]][["crve_Bt"]] = list("rejection" = c(),
                                            "wald" = c())
    results_hoc[[paste(G)]][["cr3"]] = list("rejection" = c(),
                                        "wald" = c())
    results_hoc[[paste(G)]][["cr3_Bt"]] = list("rejection" = c(),
                                           "wald" = c())
  } # If not empty we want to add more iterations to improve our accuracy !
  
  
  for (it in 1:mCarlo) {# Start of Monte Carlo iteration
    # Random data generation process
    data = dgp_homoskedastic(G=G,N_g=N_g)
  
    # Estimation of regression:
    reg = lm(y~x-1,data=data)

    # Basic ols variance
    wald_OLS = abs(waldOLS(reg))
    results_hoc[[paste(G)]][["ols"]][["wald"]] = c(results_hoc[[paste(G)]][["ols"]][["wald"]],
                                               wald_OLS)
    results_hoc[[paste(G)]][["ols"]][["rejection"]] = c(results_hoc[[paste(G)]][["ols"]][["rejection"]],
                                               as.integer(wald_OLS>thresh))
    # CRVE variance estimator
    wald_CRVE = abs(waldCR(data))
    results_hoc[[paste(G)]][["crve"]][["wald"]] = c(results_hoc[[paste(G)]][["crve"]][["wald"]],
                                                wald_CRVE)
    results_hoc[[paste(G)]][["crve"]][["rejection"]] = c(results_hoc[[paste(G)]][["crve"]][["rejection"]],
                                                as.integer(wald_CRVE>thresh))
    # CR3 variance estimator
    wald_CR3 = abs(waldOLS(coef_test(reg, df = data, vcov = vcovCR(reg, type = "CR3",cluster = data$g))))
    results_hoc[[paste(G)]][["cr3"]][["wald"]] = c(results_hoc[[paste(G)]][["cr3"]][["wald"]],
                                               wald_CR3)
    results_hoc[[paste(G)]][["cr3"]][["rejection"]] = c(results_hoc[[paste(G)]][["cr3"]][["rejection"]],
                                               as.integer(wald_CR3>thresh))
    
  
    # Initialisation of storage for monte carlo
    boot_ols = c()
    boot_crve = c()
    boot_cr3 = c()
    
    for (b in 1:bStrap) {
    # Random sample of the data, with replacement
    index_b <- sample(N_g,N_g,replace=TRUE)
    data_b <- do.call(rbind,lapply(1:G,function(k){
      t = data[data$g==k ,]
      return(t[index_b,])}))

    # Estimation of regression:
    reg = lm(y~x-1,data=data_b)
    
    # Basic ols variance
    boot_ols = c(boot_ols, abs(waldOLS(reg)))
    
    # CRVE variance estimator
    boot_crve = c(boot_crve,
                  abs(waldCR(data_b)))
    
    # CR3 variance estimator
    boot_cr3 = c(boot_cr3,
                 abs(waldOLS(coef_test(reg, df = data_b, 
                                       vcov = vcovCR(reg, 
                                                     type = "CR3",
                                                     cluster = data_b$g)))))
   
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
    
    
  }
}
saveRDS(object = results_hoc, 'results_hoc_pairs_t.RDS')

# Printing results_hoc
results_hoc = readRDS('results_hoc_pairs_t.RDS')
boot_results_hoc = c('ols','ols_Bt','crve','crve_Bt','cr3_Bt','cr3')
table_results_hoc = do.call(rbind,lapply(clusters, function(c){
  tempo = sapply(boot_results_hoc,function(res){
    return(mean(results_hoc[[paste(c)]][[res]][["rejection"]]))
  })
  tempo = t(as.data.frame(tempo))
  return(data.frame("n_cluster"=c,tempo,row.names = NULL))
}))
print(table_results_hoc)


# 03. EXTENSION WITH SKEWED ESTIMATOR -------------------------------------

