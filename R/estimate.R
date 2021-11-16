#' Estimates the Bayesian probability of informed trading
#' 
#' @description
#' Calling [estimate_pin()] on trades data estimates the probability of 
#' informed trading from the compressed EKOP model presented in Grammig, 
#' Theissen and Zehnder (2015). 
#' 
#' @details 
#' Internally this function uses the `finmix` package for estimation of the 
#' finite mixture parameters. The `finmix` package performs MCMC sampling for 
#' the parameters of the compressed EKOP model and computes three parameter 
#' estimates namely
#' * Maximum a posterior (MAP) the mode of the joint posterior of parameters 
#'   and indicators,
#' * Bayesian maximum likelihood (BML) the mode of the joint posterior of 
#'   parameters and indicators in case of a flat prior distribution, 
#' * Identified ergodic average (IEAVG) is the ergodic average over the MCMC
#'   parameter traces after relabeling the parameters.
#' 
#' ## Relabeling
#' The MCMC sampling uses a socalled random permutation algorithm to force the 
#' labels of the two distributions in the mixture to switch randomly between 
#' iterations. Label switching occurs naturally in estimation of finite mixture 
#' distributions and is documented in the literature. Random permutation tries 
#' not to eliminate label switching, but to bring it into a more controlled 
#' setting. As a result each MCMC parameter trace contains parameter values 
#' of the other trace. For this reason re-labeling is performed after sampling. 
#' 
#' The `finmix` package offers three different re-labeling algorithms: 
#' * `kmeans` re-labels parameters by applying K-means clustering to the 
#'   merged component parameters,
#' * `Stephens1997a` applies the algorithm suggested by Stephens (1997a) that 
#'   tries to make the posterior marginal distributions of parameters uni-modal,
#' * `Stephens1997b` uses the algorithm presented in Stephens (1997b) that 
#'   tries to re-label parameter values by minimizing the Kullback-Leibler 
#'   distance between a parameter clustering via assignment probabilities and 
#'   the classification probabilities defined by Poisson distributions. 
#' 
#' To offer the user more flexibility the function can also return the MCMC 
#' samples for further inspection. Setting `return_mcmc` to `TRUE` (default) 
#' returns also the MCMC traces. 
#' 
#' @param pin_data A `vector` of trades data (number of trades per day).
#' @param return_mcmc A logical indicating, if MCMC parameter samples should 
#'   be returned in addition to the final PIN estimates.
#' @param with_stephens A logical indicating, if in addition to `kmeans` re-
#' labeling also re-labeling via the algorithms from Stephens (1997a) and 
#' Stephens (1997b) should be applied.
#' @return A `data.frame()` with PIN estimates, if `return_mcmc` is `FALSE`. 
#'   Otherwise, a `list` containing the `data.frame` of PIN estimates together 
#'   with a `finmix::mcmcoutput` (see \code{\link[finmix]{mcmcoutput-class}}) 
#'   object containing the MCMC traces of the component paramaters.
#'   
#' @seealso 
#' * \code{\link[finmix]{mixturemcmc}} for performing MCMC sampling
#' * \code{\link[finmix]{mcmcestimate}} for estimating parameters from MCMC samples
#' * [compute_bayespin()] for computing the Bayesian PIN from Bayesian 
#'   estimates of the component parameters of the compressed EKOP model
#'
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the 
#'   Probability of Informed Trading. Conference on Financial Econometrics & 
#'   Empirical Asset Pricing 2016, Lancaster University
#' * Stephens, M., 1997a. Discussion 'On Bayesian analysis of mixtures with an
#'   unknown number of components' (by S. Richardson and P.J. Green), J. R. 
#'   Statist. Soc., B, 59, 768-769.
#' * Stephens, M., 1997b. Bayesian methods for mixtures of normal 
#'   distributions. DPhil. Thesis. University of Oxford, Oxford.
"estimate_pin" <- function(pin_data, return_mcmc=FALSE, with_stephens=FALSE)
  {
    # Create `fdata` object with dimension `r=1`.
    pin_fdata <- finmix::fdata(matrix(pin_data), type='discrete', r=1)
    
    # Create a finite mixture model of two Poisson distributions.
    pin_model <- finmix::model(dist='poisson', K=2, r=1)
    
    # Define the prior distribution as a conditional conjugate
    # prior (a Gamma distribution) with a hierarchical prior 
    # (also Gamma).
    pin_prior <- finmix::priordefine(fdata=pin_fdata, model=pin_model)
    
    # Set up the MCMC algorithm (a Gibbs sampler) with 10,000 iterations
    # and a burnin of 1,000. Allow random permutation of labels and do 
    # not store the indicators `S`.
    pin_mcmc <- finmix::mcmc(burnin=1000, M=10000, ranperm=TRUE, 
                             storeS=FALSE, storepost=FALSE)
    if (with_stephens) pin_mcmc@storepost <- TRUE
    
    # Generate starting parameters for the indicators `S`.
    (pin_fdata~pin_model~pin_mcmc) %=% finmix::mcmcstart(fdata=pin_fdata, 
                                                         model=pin_model, 
                                                         pin_mcmc)
    
    # Start the MCMC sampling.
    pin_results <- finmix::mixturemcmc(fdata=pin_fdata, model=pin_model, 
                                       prior=pin_prior, mcmc=pin_mcmc)
    
    # Estimate the parameters. 
    pin_parest    <- finmix::mcmcestimate(pin_results)
    # Calculate the PIN estimates.
    pin_estimates <- compute_bayespin(pin_parestimates)
    if (with_stephens) 
    {
      
      tryCatch(
        {
          pin_parest_st1997a <- finmix::mcmcestimate(pin_results, 
                                                     method='Stephens1997a')
          pin_est_st1997a    <- compute_bayespin(pin_parest_st1997a)[1,9:12]
          colnames(pin_est_st1997a) <- c('alpha_ieavg_stephens1997a', 
                                         'epsilon_ieavg_stephens1997a',
                                         'mu_ieavg_stephens1997a', 
                                         'pin_ieavg_stephens1997a')
          pin_estimates      <- cbind(pin_estimates, pin_est_st1997a)
        },
        error= function(err)
        {
          print(err)
          pin_estimates <- cbind(pin_estimates, rep(NA, 4))
        },
        warning=function(warn) {}
      )
      
      tryCatch(
        {
          # Stephens1997b needs the data to draw conclusion about the labels.
          pin_parest_st1997b <- finmix::mcmcestimate(pin_results, 
                                                     method='Stephens1997b', 
                                                     fdata=pin_fdata)
          pin_est_st1997b    <- compute_bayespin(pin_parest_st1997b)[1,9:12]
          colnames(pin_est_st1997b) <- c('alpha_ieavg_stephens1997b', 
                                         'epislon_ieavg_stephens1997b',
                                         'mu_ieavg_stephens1997b', 
                                         'pin_ieavg_stephens1997b')
          pin_estimates      <- cbind(pin_estimates, pin_est_st1997b)
        },
        error = function(err)
        {
          print(err)
          pin_estimates <- cbind(pin_estimates, rep(NA, 4))
        },
        warning=function(warn) {}
      )
    }
    
    if (return_mcmc) 
    {
      return(list("pin_estimates"=pin_estimates, "pin_results"=pin_results))
    } else 
    {
      return(pin_estimates)
    }
}

#' Calculates estimates of the probability of informed trading
#' 
#' @description 
#' Calling [compute_bayespin()] calculates the probability of informed trading
#' (PIN) from the paper of Easley et al. (1996). The input argument is an
#' `mcmcest` (see \code{\link[finmix]{mcmcest-class}}) object from the `finmix`
#' package containing all parameters estimates from the finite mixture
#' distribution of the compressed EKOP model in Grammig, Theissen and Zehnder
#' (2015).
#' 
#' @param pin_estimates An `mcmcest` object of the `finmix` package containing 
#'   all estimated parameters from the finite mixture distribution of the 
#'   compressed EKOP model.
#' @return A `data.frame` with estimated PINs from the maximum a posterior, the 
#'   Bayesian maximum likelihood and the identified ergodic average parameter
#'   estimates of the underlying finite mixture distribution of the compressed 
#'   EKOP model.
#' @export
#' 
#' @seealso 
#' * \code{\link[finmix]{mcmcest-class}} for the definition of the `mcmcest` class union
#' * [estimate_pin()] for estimating the PIN with the Bayesian approach 
#'   described in Grammig, Theissen and Zehnder (2015)
#' * [compute_pin()] for calculating the PIN for provided parameters
#' 
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the 
#'   Probability of Informed Trading. Conference on Financial Econometrics & 
#'   Empirical Asset Pricing 2016, Lancaster University
"compute_bayespin" <- function(pin_estimates)
{
  # Calculate PIN
  ## Get the order of the Poisson parameters as the larger one
  ## will be (mu + epsilon). 
  ordered_map   <- sort.int(pin_estimates@map$par$lambda, index.return=TRUE, 
                            decreasing=TRUE)
  ordered_bml   <- sort.int(pin_estimates@bml$par$lambda, index.return=TRUE, 
                            decreasing=TRUE)
  ordered_ieavg <- sort.int(pin_estimates@ieavg$par$lambda, index.return=TRUE, 
                            decreasing=TRUE)
  
  ## Note that we use T = 390 minutes.
  ## Furthermore, note that lambda_2 = 2 * epsilon * T.
  epsilon_map   <- ordered_map$x[2] / (2 * 390)
  epsilon_bml   <- ordered_bml$x[2] / (2 * 390)
  epsilon_ieavg <- ordered_ieavg$x[2] / (2 * 390)

  ## Note that lambda_1 = (2 * epsilon + mu) * T
  mu_map   <- ordered_map$x[1] / 390 - 2 * epsilon_map
  mu_bml   <- ordered_bml$x[1] / 390 - 2 * epsilon_bml
  mu_ieavg <- ordered_ieavg$x[1] / 390 - 2 * epsilon_ieavg
  
  ## Get the weight for the first component. 
  alpha_map   <- pin_estimates@map$weight[ordered_map$ix[1]]
  alpha_bml   <- pin_estimates@bml$weight[ordered_bml$ix[1]]
  alpha_ieavg <- pin_estimates@ieavg$weight[ordered_ieavg$ix[1]]
  
  ## Calculate the PIN
  pin_map   <- compute_pin(alpha_map, epsilon_map, mu_map)
  pin_bml   <- compute_pin(alpha_bml, epsilon_bml, mu_bml)
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

#' Computes the probability of informed trading
#' 
#' @description 
#' Calling [compute_pin()] computes the probability of informed trading (PIN)
#' from the EKOP model for provided parameters. See Easley et al. (1996) for
#' details.
#' 
#' @param alpha A double specifying the probability of an information event. 
#' @param epsilon A double specifying the arrival rate of uninformed traders.
#' @param mu A double specifying the arrival rate of informed traders.
#' @return The probability of informed trading.
#' @export
#' 
#' @examples 
#' compute_pin(.2, .25, .01)
#' 
#' @seealso
#' * [compute_bayespin()] for computing the PIN from bayesian parameter 
#'   estimates
#' 
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436. 
"compute_pin" <- function(alpha, epsilon, mu)
{
  pin <- alpha * mu / ( alpha * mu + 2 * epsilon )
  return(pin)
}

"estimate_mlekop" <- function(data, startpar, T = 390, methodLik = c("precise", "approx"),
                              fnLik = c("computeEKOPLik", "computeEKOPOrigLik"), 
                              sim = FALSE, mis = FALSE, mis.prob, misInd, fnscale=-1,
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
  optim_lower 	<- c(-1e+2, 1e-6, -1e+2, 1e-6)
  optim_upper 	<- c(1e+2, 1e+6, 1e+2, 1e+6)
  optim_fnscale <- fnscale 
  optim_maxit 	<- 200
  optim_ctrl	  <- list(fnscale=optim_fnscale, maxit=optim_maxit, trace=trace)
  
  optim_res <- optim(par=startpar, fn=optim_fn, data=data, T=T, 
                     methodLik=methodLik, method=optim_Method, 
                     lower=optim_lower, upper=optim_upper, 
                     control=optim_ctrl, hessian=TRUE)
  
  conv_msg  <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH" 
  
  if ((optim_res$convergence > 0 || conv_msg != optim_res$message) && grad_free)
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

"estimate_compml" <- function(data, startpar, T = 390, methodLik = c("precise", "approx"),
                              sim = FALSE, mis = FALSE, mis.prob, misInd, fnscale=-1,
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
  optim_fn 	    <- computeCompLik
  optim_Method 	<- "L-BFGS-B"
  optim_lower 	<- c(-1e+6, 1e-6, 1e-6)
  optim_upper 	<- c(1e+2, 1e+6, 1e+6)
  optim_fnscale <- fnscale 
  optim_maxit 	<- 200
  optim_ctrl	  <- list(fnscale=optim_fnscale, maxit=optim_maxit, trace=trace)
  
  optim_res <- optim(par=startpar, fn=optim_fn, data=data, T=T, 
                     methodLik=methodLik, method=optim_Method, 
                     lower=optim_lower, upper=optim_upper, 
                     control=optim_ctrl, hessian=TRUE)
  
  conv_msg  <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
  
  if ((optim_res$convergence > 0 || conv_msg != optim_res$message) && grad_free)
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

