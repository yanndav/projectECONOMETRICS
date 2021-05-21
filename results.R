########################################
#    ECONOMETRICS PROJECT - RESULTS    #
#       YANN DAVID - ANDRE SANTOS      #
########################################


# 00. LOADING RESULTS -----------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 
library('tidyverse')
estimators_wald = c('ols','ols_Bt','crve','crve_Bt','cr3','cr3_Bt',
                    'resi_Bt','wild_Bt')

names <- c('Assume i.i.d.',' ',
           'Cluster-robust','','CR3 residual correction','',
           'BDM bootstrap-t','','Pairs cluster bootstrap-t','',
           'Pairs CR3 bootstrap-t','',
           'Residual cluster bootstrap-t','',
           'Wild cluster bootstrap-t','')

clusters = c(5,10 ,15,20,25,30             ) # Number of clusters to be tested


# HOMOSKEDASTICITY NORMAL RESULTS ------------------------------------------


if(!exists("results_hoc")){
  results_hoc = readRDS('results_hoc_pairs_t.RDS')
}


# Printing results_hoc
table_results_hoc = lapply(clusters, function(c){
  tempo = do.call(rbind,lapply(estimators_wald,function(res){
    dt = data.frame(method = c("mean","se"),
               "wald"=rep(res,2),
               row.names = NULL)
    rej = mean(results_hoc[[paste(c)]][[res]][["rejection"]])
    R = length(results_hoc[[paste(c)]][[res]][["rejection"]])
    dt[[paste(c)]] <- c(rej,
                        sqrt(rej*(1-rej)/(R-1)))
    
    return(dt)
  }))
  return(tempo)
}) %>% reduce(full_join, by=c("method","wald")) %>% 
  mutate(rank = factor(paste0(wald,method),
                       levels = paste0(rep(c('ols','crve','cr3',
                                         'ols_Bt','crve_Bt','cr3_Bt',
                                         'resi_Bt','wild_Bt'),each=2),
                                       c('mean','se'))))


hoc_tex = with(table_results_hoc, table_results_hoc[order(rank, rank),]) %>% 
  mutate_if(is.numeric, round, digits=3) %>% 
  mutate(across(`5`:`30`, 
         ~ifelse(method == "se",paste0("(",.x,")"),.x)))


  
hoc_tex["Method"] <- names

write.table(hoc_tex[,c('Method',seq(5,30,5))],
            file = "output/homoske_norma.tex", eol='\\\\',
            sep='&',quote = F,row.names = F,
            col.names = F,na="")
  


# HETEROSKEDASTICITY NORMAL RESULTS ---------------------------------------


if(!exists("results_het")){
  results_het = readRDS('results_het_pairs_t.RDS')
}


# Printing results_het
table_results_het = lapply(clusters, function(c){
  tempo = do.call(rbind,lapply(estimators_wald,function(res){
    dt = data.frame(method = c("mean","se"),
                    "wald"=rep(res,2),
                    row.names = NULL)
    rej = mean(results_het[[paste(c)]][[res]][["rejection"]])
    R = length(results_het[[paste(c)]][[res]][["rejection"]])
    dt[[paste(c)]] <- c(rej,
                        sqrt(rej*(1-rej)/(R-1)))
    
    return(dt)
  }))
  return(tempo)
}) %>% reduce(full_join, by=c("method","wald")) %>% 
  mutate(rank = factor(paste0(wald,method),
                       levels = paste0(rep(c('ols','crve','cr3',
                                             'ols_Bt','crve_Bt','cr3_Bt',
                                             'resi_Bt','wild_Bt'),each=2),
                                       c('mean','se'))))


het_tex = with(table_results_het, table_results_het[order(rank, rank),]) %>% 
  mutate_if(is.numeric, round, digits=3) %>% 
  mutate(across(`5`:`30`, 
                ~ifelse(method == "se",paste0("(",.x,")"),.x)))



het_tex["Method"] <- names

write.table(het_tex[,c('Method',seq(5,30,5))],
            file = "output/hetero_norma.tex", eol='\\\\',
            sep='&',quote = F,row.names = F,
            col.names = F,na="")


# HOMOSKEDASTICITY SKEWED RESULTS ---------------------------------------


if(!exists("results_skho")){
  results_skho = readRDS('results_skho_pairs_t.RDS')
}


# Printing results_skho
table_results_skho = lapply(clusters, function(c){
  tempo = do.call(rbind,lapply(estimators_wald,function(res){
    dt = data.frame(method = c("mean","se"),
                    "wald"=rep(res,2),
                    row.names = NULL)
    rej = mean(results_skho[[paste(c)]][[res]][["rejection"]])
    R = length(results_skho[[paste(c)]][[res]][["rejection"]])
    dt[[paste(c)]] <- c(rej,
                        sqrt(rej*(1-rej)/(R-1)))
    
    return(dt)
  }))
  return(tempo)
}) %>% reduce(full_join, by=c("method","wald")) %>% 
  mutate(rank = factor(paste0(wald,method),
                       levels = paste0(rep(c('ols','crve','cr3',
                                             'ols_Bt','crve_Bt','cr3_Bt',
                                             'resi_Bt','wild_Bt'),each=2),
                                       c('mean','se'))))


skho_tex = with(table_results_skho, table_results_skho[order(rank, rank),]) %>% 
  mutate_if(is.numeric, round, digits=3) %>% 
  mutate(across(`5`:`30`, 
                ~ifelse(method == "se",paste0("(",.x,")"),.x)))



skho_tex["Method"] <- names

write.table(skho_tex[,c('Method',seq(5,30,5))],
            file = "output/homo_ske.tex", eol='\\\\',
            sep='&',quote = F,row.names = F,
            col.names = F,na="")


# HETEROSKEDASTICITY SKEWED RESULTS ---------------------------------------


if(!exists("results_skhet")){
  results_skhet = readRDS('results_skhet_pairs_t.RDS')
}


# Printing results_skhet
table_results_skhet = lapply(clusters, function(c){
  tempo = do.call(rbind,lapply(estimators_wald,function(res){
    dt = data.frame(method = c("mean","se"),
                    "wald"=rep(res,2),
                    row.names = NULL)
    rej = mean(results_skhet[[paste(c)]][[res]][["rejection"]])
    R = length(results_skhet[[paste(c)]][[res]][["rejection"]])
    dt[[paste(c)]] <- c(rej,
                        sqrt(rej*(1-rej)/(R-1)))
    
    return(dt)
  }))
  return(tempo)
}) %>% reduce(full_join, by=c("method","wald")) %>% 
  mutate(rank = factor(paste0(wald,method),
                       levels = paste0(rep(c('ols','crve','cr3',
                                             'ols_Bt','crve_Bt','cr3_Bt',
                                             'resi_Bt','wild_Bt'),each=2),
                                       c('mean','se'))))


skhet_tex = with(table_results_skhet, table_results_skhet[order(rank, rank),]) %>% 
  mutate_if(is.numeric, round, digits=3) %>% 
  mutate(across(`5`:`30`, 
                ~ifelse(method == "se",paste0("(",.x,")"),.x)))



skhet_tex["Method"] <- names

write.table(skhet_tex[,c('Method',seq(5,30,5))],
            file = "output/hetero_sk.tex", eol='\\\\',
            sep='&',quote = F,row.names = F,
            col.names = F,na="")


# HOMOSKEDASTICITY SKEWED SIGMA 2 RESULTS ---------------------------------------


if(!exists("results_skhover")){
  results_skhover = readRDS('results_skhover_pairs_t.RDS')
}


# Printing results_skhover
table_results_skhover = lapply(clusters, function(c){
  tempo = do.call(rbind,lapply(estimators_wald,function(res){
    dt = data.frame(method = c("mean","se"),
                    "wald"=rep(res,2),
                    row.names = NULL)
    rej = mean(results_skhover[[paste(c)]][[res]][["rejection"]])
    R = length(results_skhover[[paste(c)]][[res]][["rejection"]])
    dt[[paste(c)]] <- c(rej,
                        sqrt(rej*(1-rej)/(R-1)))
    
    return(dt)
  }))
  return(tempo)
}) %>% reduce(full_join, by=c("method","wald")) %>% 
  mutate(rank = factor(paste0(wald,method),
                       levels = paste0(rep(c('ols','crve','cr3',
                                             'ols_Bt','crve_Bt','cr3_Bt',
                                             'resi_Bt','wild_Bt'),each=2),
                                       c('mean','se'))))


skhover_tex = with(table_results_skhover, table_results_skhover[order(rank, rank),]) %>% 
  mutate_if(is.numeric, round, digits=3) %>% 
  mutate(across(`5`:`30`, 
                ~ifelse(method == "se",paste0("(",.x,")"),.x)))



skhover_tex["Method"] <- names

write.table(skhover_tex[,c('Method',seq(5,30,5))],
            file = "output/homo_skevery.tex", eol='\\\\',
            sep='&',quote = F,row.names = F,
            col.names = F,na="")


# HETEROSKEDASTICITY SKEWED RESULTS ---------------------------------------


if(!exists("results_skhetver")){
  results_skhetver = readRDS('results_skhetver_pairs_t.RDS')
}


# Printing results_skhetver
table_results_skhetver = lapply(clusters, function(c){
  tempo = do.call(rbind,lapply(estimators_wald,function(res){
    dt = data.frame(method = c("mean","se"),
                    "wald"=rep(res,2),
                    row.names = NULL)
    rej = mean(results_skhetver[[paste(c)]][[res]][["rejection"]])
    R = length(results_skhetver[[paste(c)]][[res]][["rejection"]])
    dt[[paste(c)]] <- c(rej,
                        sqrt(rej*(1-rej)/(R-1)))
    
    return(dt)
  }))
  return(tempo)
}) %>% reduce(full_join, by=c("method","wald")) %>% 
  mutate(rank = factor(paste0(wald,method),
                       levels = paste0(rep(c('ols','crve','cr3',
                                             'ols_Bt','crve_Bt','cr3_Bt',
                                             'resi_Bt','wild_Bt'),each=2),
                                       c('mean','se'))))


skhetver_tex = with(table_results_skhetver, table_results_skhetver[order(rank, rank),]) %>% 
  mutate_if(is.numeric, round, digits=3) %>% 
  mutate(across(`5`:`30`, 
                ~ifelse(method == "se",paste0("(",.x,")"),.x)))



skhetver_tex["Method"] <- names

write.table(skhetver_tex[,c('Method',seq(5,30,5))],
            file = "output/hetero_skvery.tex", eol='\\\\',
            sep='&',quote = F,row.names = F,
            col.names = F,na="")

