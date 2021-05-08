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
  x_ig = exp(mu + sigma(z_g + z_ig))
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

varCR <- function(data){
  # FUNCTION DESCRIPTION - - - - - -
  
  # This function returns the CRVE
  
  # SETTING PARAMETERS - - - - - - -
  
  G = max(data$g)
  x = data$x
  beta_hat = betaH(data)
  
  # VARIANCE ESTIMATION - - - - - -
  
  return(solve(t(x)%*% x) %*% 
    sum(sapply(1:G,function(g){
    x_g = data$x[data$g==g]
    u_g = (G/(G-1))*data$y[data$g==g] - x_g %*% beta_hat
    return(t(x_g)%*%u_g%*%t(u_g)%*%x_g)
  })) %*%  solve(t(x)%*% x) )

} 

waldOLS <- function(reg,beta_null=1){
  if(!is.data.frame(reg)){
    tempo = summary(reg)$coefficients
    beta = tempo["x","Estimate"]
    se = tempo["x","Std. Error"]
  }else{
    beta = reg[1,"beta"]
    se = reg[1,"SE"]
  }
  
  
  return((beta-beta_null)/se)
}


waldCR <- function(data,beta_null=1){
  return((betaH(data)-beta_null)%*%
           solve(sqrt(varCR(data))))
}


varJACK <- function(data){
  G=max(data$g)
  
  return(((G-1)/G)*
           sum(sapply(1:G, function(g){
    betahat = betaH(data)
    beta_ghat = betaH(data[data$g!=g,])
    return((beta_ghat-betahat)^2)
  })))

  
  }
