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
  e_g = as.vector(replicate(n = G,rnorm(n=N_g, mean=0, sd=1)))
  # Individual part of error
  e_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  ## Covariates part 
  # Correlation of covariates within each cluster
  z_g = as.vector(replicate(n = G,rnorm(n=N_g, mean=0, sd=1)))
  # Individual part of covariates
  z_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  # GENERATING OUTCOME VARIABLES - - - - -
  
  u_ig = e_g + e_ig
  x_ig = z_g + z_ig
  y_ig = b_0 + b_1*x_ig + u_ig
  
  # OUTPUT - - - - - - - - - - - - - - - -
  
  return(data.frame("y"=y_ig,
                    "x"=x_ig))
  
}

dgp_heteroskedastic <- function(G,N_g=30,
                              b_0 = 0, 
                              b_1=1){
  # FUNCTION DESCRIPTION   - - - - - - - -
  
  # This function codes for the data-generating process as defined in Cameron, Gelbach and Miller (2008).
  # Errors are independent across clusters but correlated within clusters
  # Errors are normally distributed and heteroskedastic.
  # We assume same number of units N_g per cluster
  
  # DRAWING DATA - - - - - - - - - - - - -
  
  ## Covariates part 
  # Correlation of covariates within each cluster
  z_g = as.vector(replicate(n = G,rnorm(n=N_g, mean=0, sd=1)))
  # Individual part of covariates
  z_ig = rnorm(n=N_g*G, mean=0, sd=1)
  
  ## Error part
  # Correlation of errors within each cluster
  e_g = as.vector(replicate(n = G,rnorm(n=N_g, mean=0, sd=1)))
  # Individual part of error: heteroskedastic
  e_ig = rnorm(n=N_g*G, mean=0, sd=9 * (z_g + z_ig)^2)
  
  # GENERATING OUTCOME VARIABLES - - - - -
  
  u_ig = e_g + e_ig
  x_ig = z_g + z_ig
  y_ig = b_0 + b_1*x_ig + u_ig
  
  # OUTPUT - - - - - - - - - - - - - - - -
  
  return(data.frame("y"=y_ig,
              "x"=x_ig))
}

# 02. CLUSTER-ROBUST ESTIMATORS ----------------------------------------------------------


