
"estimate_pin" <- function(pin_data)
  {
    # Create `fdata` object with dimension `r=1`.
    pin_fdata <- finmix::fdata(matrix(pin_data), type='discrete', r=1)
    
    # Create a finite mixture model of two Poisson distributions.
    pin_model <- finmix::model(dist='poisson', K=2, r=1)
    
    # Define the prior distribution as a conditional conjugate
    # prior (a Gamma distribtion) with a hierarchical prior 
    # (also Gamma).
    pin_prior <- finmix::priordefine(fdata=pin_fdata, model=pin_model)
    
    # Set up the MCMC algorithm (a Gibbs sampler) with 10,000 iterations
    # and a burnin of 1,000. Allow random permutations of labels and do 
    # not store the indicators `S`.
    pin_mcmc <- finmix::mcmc(burnin=1000, M=10000, ranperm=TRUE, storeS=FALSE, storepost=FALSE)
    
    # Generate starting parameters for the indicators `S`.
    (pin_fdata~pin_model~pin_mcmc) %=% finmix::mcmcstart(fdata=pin_fdata, model=pin_model, pin_mcmc)
    
    # Start the MCMC sampling.
    pin_results <- finmix::mixturemcmc(fdata=pin_fdata, model=pin_model, prior=pin_prior, mcmc=pin_mcmc)
    
    # Estimate the parameters. 
    pin_parestimates <- finmix::mcmcestimate(pin_results)
    
    # Calculate the PIN estimates.
    pin_estimates <- compute_bayespin(pin_parestimates)
    
    return(pin_estimates)
}

"compute_bayespin" <- function(pin_estimates, params=TRUE)
{
  # Calculate PIN from IEAVG
  ## Get the order of the Poisson parameters as the larger one
  ## will be (mu + epsilon). 
  ordered_map <- sort.int(pin_estimates@map$par$lambda, index.return=TRUE, decreasing=TRUE)
  ordered_bml <- sort.int(pin_estimates@bml$par$lambda, index.return=TRUE, decreasing=TRUE)
  ordered_ieavg <- sort.int(pin_estimates@ieavg$par$lambda, index.return=TRUE, decreasing=TRUE)
  
  ## Note that we use T = 390 minutes.
  ## Furthermore, note that lambda_2 = 2 * epsilon * T.
  epsilon_map <- ordered_map$x[2] / (2 * 390)
  epsilon_bml <- ordered_bml$x[2] / (2 * 390)
  epsilon_ieavg <- ordered_ieavg$x[2] / (2 * 390)

  ## Note that lambda_1 = (2* epsilon + mu) * T
  mu_map <- ordered_map$x[1] / 390 - 2 * epsilon_map
  mu_bml <- ordered_bml$x[1] / 390 - 2 * epsilon_bml
  mu_ieavg <- ordered_ieavg$x[1] / 390 - 2 * epsilon_ieavg
  
  ## Get the weight for the first component. 
  alpha_map <- pin_estimates@map$weight[ordered_map$ix[1]]
  alpha_bml <- pin_estimates@bml$weight[ordered_bml$ix[1]]
  alpha_ieavg <- pin_estimates@ieavg$weight[ordered_ieavg$ix[1]]
  
  ## Calculate the PIN
  pin_map <- compute_pin(alpha_map, epsilon_map, mu_bml)
  pin_bml <- compute_pin(alpha_bml, epsilon_bml, mu_bml)
  pin_ieavg <- compute_pin(alpha_ieavg, epsilon_ieavg, mu_ieavg)

  col_names <- c('alpha_map', 'epsilon_map', 'mu_map', 'pin_map', 
                 'alpha_bml', 'epsilon_bml', 'mu_bml', 'pin_bml',
                 'alpha_ieavg', 'epsilon_ieavg', 'mu_ieavg', 'pin_ieavg')
  estimates <- data.frame(alpha_map, epsilon_map, mu_map, pin_map,
                          alpha_bml, epsilon_bml, mu_bml, pin_bml,
                          alpha_ieavg, epsilon_ieavg, mu_ieavg, pin_ieavg)
  colnames(estimates) <- col_names
  
  return(estimates)
}

"compute_pin" <- function(alpha, epsilon, mu)
{
  pin <- alpha * mu / ( alpha * mu + 2 * epsilon )
  return(pin)
}

"estimate_mlekop" <- function(data, startpar, T = 390, methodLik = c("precise", "approx"),
                              fnLik = c("computeEKOPLik", "computeEKOPOrigLik"), 
                              sim = FALSE, mis = FALSE, mis.prob, misInd, fnscale=-0.05,
                              trace=0, grad_free=TRUE) 
{
  fnLik <- match.arg(fnLik)
  if (missing(startpar)) {
    if (trace > 0)
    {
      cat("Using default starting values...\n")
    }
    tmp         <- mean(data[, 5], na.rm = TRUE)/T
    startpar 	  <- c(0, tmp * 0.75/2, 0, tmp * 0.25/2)
  }
  
  ## optimization settings ##
  optim_fn 	    <- match.fun(fnLik)
  optim_Method 	<- "L-BFGS-B"
  optim_lower 	<- c(-1e+6, 1e-6, -1e+6, 1e-6)
  optim_upper 	<- c(1e+6, 1e+6, 1e+6, 1e+6)
  optim_fnscale <- fnscale 
  optim_maxit 	<- 200
  optim_ctrl	  <- list(fnscale=optim_fnscale, maxit=optim_maxit, trace=trace)
  
  optim_res <- optim(par=startpar, fn=optim_fn, data=data, T=T, 
                     methodLik=methodLik, method=optim_Method, 
                     lower=optim_lower, upper=optim_upper, 
                     control=optim_ctrl, hessian=TRUE)
  
  if (optim_res$convergence > 0 && grad_free)
  {
    dfoptim_ctrl  <- list(maximize=TRUE)
    # Bounded Nelder-Mead needs lower absolute bounds for the 
    # alpha and delta , otherwise the log-likelihood becomes NaN.
    dfoptim_lower <- c(-1e+1, 1e-6, -1e+1, 1e-6)
    dfoptim_upper <- c(1e+1, 1e+6, 1e+1, 1e+6)
    optim_res     <- dfoptim::nmkb(par=startpar, fn=optim_fn, data=data, T=T, 
                                   methodLik=methodLik, lower=dfoptim_lower,
                                   upper=dfoptim_upper, control=dfoptim_ctrl)
  }
  
  return(optim_res)	
}

"estimate_mlkokot" <- function(data, startpar, T = 390, methodLik = c("precise", "approx"),
                                sim = FALSE, mis = FALSE, mis.prob, misInd, fnscale=-0.05,
                                trace=0, grad_free=TRUE) 
{
  if (missing(startpar)) {
    if (trace > 0)
    {
      cat("Using default starting values...\n")
    }
    tmp         <- mean(data, na.rm = TRUE)/T
    startpar 	  <- c(0, tmp * 0.75/2, tmp * 0.25/2)
  }
  
  ## optimization settings ##
  optim_fn 	    <- computeKokotLik
  optim_Method 	<- "L-BFGS-B"
  optim_lower 	<- c(-1e+6, 1e-6, 1e-6)
  optim_upper 	<- c(1e+6, 1e+6, 1e+6)
  optim_fnscale <- fnscale 
  optim_maxit 	<- 200
  optim_ctrl	  <- list(fnscale=optim_fnscale, maxit=optim_maxit, trace=trace)
  
  optim_res <- optim(par=startpar, fn=optim_fn, data=data, T=T, 
                     methodLik=methodLik, method=optim_Method, 
                     lower=optim_lower, upper=optim_upper, 
                     control=optim_ctrl, hessian=TRUE)
  if (optim_res$convergence > 0 && grad_free)
  {
    dfoptim_ctrl  <- list(maximize=TRUE)
    dfoptim_lower <- c(-1e+1, 1e-6, 1e-6)
    dfoptim_upper <- c(1e+1, 1e+6, 1e+6)
    optim_res     <- dfoptim::nmkb(par=startpar, fn=optim_fn, data=data, T=T, 
                                   methodLik=methodLik, lower=dfoptim_lower,
                                   upper=dfoptim_upper, control=dfoptim_ctrl)
  }
  return(optim_res)	
}

