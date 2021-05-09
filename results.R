########################################
#    ECONOMETRICS PROJECT - RESULTS    #
#       YANN DAVID - ANDRE SANTOS      #
########################################


# 00. LOADING RESULTS -----------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 
estimators_wald = c('ols','ols_Bt','crve','crve_Bt','cr3_Bt','cr3',
                    'resi_Bt','wild_Bt')

names <- c('Assume i.i.d.',
           'Cluster-robust','CR3 residual correction',
           'BDM bootstrap-t','Pairs cluster bootstrap-t',
           'Pairs CR3 bootstrap-t',
           'Residual cluster bootstrap-t','Wild cluster bootstrap-t')

results_hoc = readRDS('results_hoc_pairs_t.RDS')

# Printing results_hoc
table_results_hoc = do.call(rbind,lapply(clusters, function(c){
  tempo = sapply(estimators_wald,function(res){
    return(c(mean(results_hoc[[paste(c)]][[res]][["rejection"]]),
             sd(results_hoc[[paste(c)]][[res]][["rejection"]])))
  })
  tempo = t(as.data.frame(tempo))
  return(data.frame("n_cluster"=c,tempo,row.names = NULL))
}))
print(table_results_hoc)

hoc_tex = t(table_results_hoc)[2:ncol(table_results_hoc),]
hoc_tex = hoc_tex[c('ols','crve','cr3_Bt',
                    'ols_Bt','crve_Bt','cr3',
                    'resi_Bt','wild_Bt'),]

rownames(hoc_tex) <- names

write.table(hoc_tex,
            file = "output/homoske_norma.tex", eol='\\\\',
            sep='&',quote = F,row.names = F,
            col.names = F,na="")
