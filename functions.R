########################################
#                                      #
#         ECONOMETRICS PROJECT         #
#         - FUNCTIONS SCRIPT -         #
#       YANN DAVID - ANDRE SANTOS      #
#                                      #
########################################


# 01. DGP FUNCTIONS -------------------------------------------------------

dgp_homoskedastic <- function(G,N_g=30,
                b_0 = 0, 
                b_1=1){
  # FUNCTION DESCRIPTION   - - - - - - - -
  
    # This function codes for the data-generating process as defined in Cameron, Gelbach and Miller (2008).
  # Errors are independent across clusters but correlated within clusters
  # Errors are normally distributed and homoskedastic.
  # We assume same number of units N_g per cluster
  
  # DRAWING DATA - - - - - - - - - - - - -
  
  ## Error part
  # Correlation of errors within each cluster
  e_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of error
  e_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  ## Covariates part 
  # Correlation of covariates within each cluster
  z_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of covariates
  z_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  # GENERATING OUTCOME VARIABLES - - - - -
  
  u_ig = e_g + e_ig
  x_ig = z_g + z_ig
  y_ig = b_0 + b_1*x_ig + u_ig
  
  # OUTPUT - - - - - - - - - - - - - - - -
  
  return(data.frame("y"=y_ig,
                    "x"=x_ig,
                    "g"=rep(seq(1:G), each=N_g)))
  
}

dgp_heteroskedastic <- function(G,N_g=30,
                              b_0 = 1, 
                              b_1=1){
  # FUNCTION DESCRIPTION   - - - - - - - -
  
  # This function codes for the data-generating process as defined in Cameron, Gelbach and Miller (2008).
  # Errors are independent across clusters but correlated within clusters
  # Errors are normally distributed and heteroskedastic.
  # We assume same number of units N_g per cluster
  
  # DRAWING DATA - - - - - - - - - - - - -
  
  ## Covariates part 
  # Correlation of covariates within each cluster
  z_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of covariates
  z_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  ## Error part
  # Correlation of errors within each cluster
  e_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of error: heteroskedastic
  e_ig = rnorm(n=N_g*G, mean=0, sd=9 * (z_g + z_ig)^2)
  
  # GENERATING OUTCOME VARIABLES - - - - -
  
  u_ig = e_g + e_ig
  x_ig = z_g + z_ig
  y_ig = b_0 + b_1*x_ig + u_ig
  
    # OUTPUT - - - - - - - - - - - - - - - -
  
  return(data.frame("y"=y_ig,
                    "x"=x_ig,
                    "g"=rep(seq(1:G), each=N_g)))
}



dgp_homoskedastic_skewed <- function(G,N_g=30,
                              b_0 = 0, 
                              b_1=1,
                              mu,sigma){
  # FUNCTION DESCRIPTION   - - - - - - - -
  
  # This function codes for the data-generating process as defined in Cameron, Gelbach and Miller (2008).
  # Errors are independent across clusters but correlated within clusters
  # Errors are normally distributed and homoskedastic.
  # We assume same number of units N_g per cluster
  
  # DRAWING DATA - - - - - - - - - - - - -
  
  ## Error part
  # Correlation of errors within each cluster
  e_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of error
  e_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  ## Covariates part 
  # Correlation of covariates within each cluster
  z_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of covariates
  z_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  # GENERATING OUTCOME VARIABLES - - - - -
  
  u_ig = e_g + e_ig
  x_ig = exp(mu + sigma*(z_g + z_ig))
  y_ig = b_0 + b_1*x_ig + u_ig
  
  # OUTPUT - - - - - - - - - - - - - - - -
  
  return(data.frame("y"=y_ig,
                    "x"=x_ig,
                    "g"=rep(seq(1:G), each=N_g)))
  
}

dgp_heteroskedastic_skewed <- function(G,N_g=30,
                                b_0 = 1, 
                                b_1=1,
                                mu,
                                sigma){
  # FUNCTION DESCRIPTION   - - - - - - - -
  
  # This function codes for the data-generating process as defined in Cameron, Gelbach and Miller (2008).
  # Errors are independent across clusters but correlated within clusters
  # Errors are normally distributed and heteroskedastic.
  # We assume same number of units N_g per cluster
  
  # DRAWING DATA - - - - - - - - - - - - -
  
  ## Covariates part 
  # Correlation of covariates within each cluster
  z_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of covariates
  z_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  ## Error part
  # Correlation of errors within each cluster
  e_g = as.vector(rep(rnorm(n=G, mean=0, sd=1),each=N_g))
  # Individual part of error: heteroskedastic
  e_ig = rnorm(n=N_g*G, mean=0, sd=9 * (z_g + z_ig)^2)
  
  # GENERATING OUTCOME VARIABLES - - - - -
  
  u_ig = e_g + e_ig
  x_ig = exp(mu+sigma*(z_g + z_ig))
  y_ig = b_0 + b_1*x_ig + u_ig
  
  # OUTPUT - - - - - - - - - - - - - - - -
  
  return(data.frame("y"=y_ig,
                    "x"=x_ig,
                    "g"=rep(seq(1:G), each=N_g)))
}

# 02. CLUSTER-ROBUST ESTIMATORS ----------------------------------------------------------
betaH = function(data){
  # FUNCTION DESCRIPTION - - - - - -
  
  # This function returns the beta_hat estimate
  
  # BETA ESTIMATION - - - - - - - -
  return(solve(t(data$x)%*%data$x) %*%t(data$x)%*%data$y)
}



wald <- function(reg,beta_null=1,cluster="",data=NULL){
  if(cluster==""){
    tempo = summary(reg)$coefficients
    beta = tempo["x","Estimate"]
    se = tempo["x","Std. Error"]
  }else{
    reg = coef_test(reg, df = data, 
              vcov = vcovCR(reg, 
                            type = cluster,
                            cluster = data$g))
    
    beta = reg[1,"beta"]
    se = reg[1,"SE"]
  }
  return(abs((beta-beta_null)/se))
}





