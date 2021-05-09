########################################
#    ECONOMETRICS PROJECT - RESULTS    #
#       YANN DAVID - ANDRE SANTOS      #
########################################


# 00. LOADING RESULTS -----------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 
estimators_wald = c('ols','ols_Bt','crve','crve_Bt','cr3_Bt','cr3',
                    'resi_Bt','wild_Bt')
results_hoc = readRDS('results_hoc_pairs_t.RDS')

# Printing results_hoc
table_results_hoc = do.call(rbind,lapply(clusters, function(c){
  tempo = sapply(estimators_wald,function(res){
    return(mean(results_hoc[[paste(c)]][[res]][["rejection"]]))
  })
  tempo = t(as.data.frame(tempo))
  return(data.frame("n_cluster"=c,tempo,row.names = NULL))
}))
print(table_results_hoc)



write.table(table_results_hoc)[2:ncol(table_results_hoc),], 
file = "output/blp_endline.tex",eol='\\\\', sep='&',quote = F,row.names = T,
            col.names = F,na="")
