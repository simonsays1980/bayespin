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
#' @export
#' 
#' @examples 
#' # Simulate trades data
#' trades_data <- simulateEKOP()
#' # Estimate the Bayesian PIN.
#' estimate_pin(trades_data$Trades)
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
    mcmc_fdata <- finmix::fdata(matrix(pin_data), type='discrete', r=1)
    
    # Create a finite mixture model of two Poisson distributions.
    mcmc_model <- finmix::model(dist='poisson', K=2, r=1)
    
    # Define the prior distribution as a conditional conjugate
    # prior (a Gamma distribution) with a hierarchical prior 
    # (also Gamma).
    mcmc_prior <- finmix::priordefine(fdata=mcmc_fdata, model=mcmc_model)
    
    # Set up the MCMC algorithm (a Gibbs sampler) with 10,000 iterations
    # and a burnin of 1,000. Allow random permutation of labels and do 
    # not store the indicators `S`.
    mcmc_mcmc <- finmix::mcmc(burnin=1000, M=10000, ranperm=TRUE, 
                             storeS=FALSE, storepost=FALSE)
    if (with_stephens) mcmc_mcmc@storepost <- TRUE
    
    # Generate starting parameters for the indicators `S`.
    (mcmc_fdata~mcmc_model~mcmc_mcmc) %=% finmix::mcmcstart(fdata=mcmc_fdata, 
                                                            model=mcmc_model, 
                                                            mcmc_mcmc)
    
    # Start the MCMC sampling.
    mcmc_results <- finmix::mixturemcmc(fdata=mcmc_fdata, model=mcmc_model, 
                                       prior=mcmc_prior, mcmc=mcmc_mcmc)
    
    # Estimate the parameters. 
    mcmc_parest    <- finmix::mcmcestimate(mcmc_results)
    # Calculate the PIN estimates.
    pin_estimates <- compute_bayespin(mcmc_parest)
    if (with_stephens) 
    {
      # Relabel with Stephens1997a.
      tryCatch(
        {
          mcmc_parest_st1997a <- finmix::mcmcestimate(mcmc_results, 
                                                      method='Stephens1997a')
          pin_est_st1997a    <- compute_bayespin(mcmc_parest_st1997a)[1,9:12]
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
      # Relabel with Stephens1997b.
      tryCatch(
        {
          # Stephens1997b needs the data to draw conclusion about the labels.
          mcmc_parest_st1997b <- finmix::mcmcestimate(mcmc_results, 
                                                      method='Stephens1997b', 
                                                      fdata=mcmc_fdata)
          pin_est_st1997b    <- compute_bayespin(mcmc_parest_st1997b)[1,9:12]
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
      return(list("pin_estimates"=pin_estimates, "mcmc_results"=mcmc_results))
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
"compute_bayespin" <- function(mcmcest)
{
  # Calculate PIN
  ## Get the order of the Poisson parameters as the larger one
  ## will be (mu + epsilon). 
  ordered_map   <- sort.int(mcmcest@map$par$lambda, index.return=TRUE, 
                            decreasing=TRUE)
  ordered_bml   <- sort.int(mcmcest@bml$par$lambda, index.return=TRUE, 
                            decreasing=TRUE)
  ordered_ieavg <- sort.int(mcmcest@ieavg$par$lambda, index.return=TRUE, 
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
  alpha_map   <- mcmcest@map$weight[ordered_map$ix[1]]
  alpha_bml   <- mcmcest@bml$weight[ordered_bml$ix[1]]
  alpha_ieavg <- mcmcest@ieavg$weight[ordered_ieavg$ix[1]]
  
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
  row.names(estimates) <- "Estimates"
  
  return(estimates)
}

#' Computes the probability of informed trading from ML estimates
#' 
#' @description
#' Calling [compute_mlpin()] returns the probability of informed trading as 
#' presented in the paper of Easley et al. (1996). 
#' 
#' The estimation procedure in [estimate_mlekop()] uses the logits of the 
#' rates `alpha` and `delta` to ensure a more stable optimization. To arrive
#' at the original parameter estimates the logistic transformation from 
#' [logistic()] is applied. 
#' 
#' @param par A vector containing the parameter estimates from the optimization 
#'   procedure. 
#' @return A double holding the PIN estimate.
#' @export
#' 
#' @examples 
#' # Simulate trades data.
#' trades_data <- simulateEKOP()
#' # Estimate the EKOP model.
#' opt_out <- estimate_mlekop(trades_data, methodLik="approx")
#' # Estimate the PIN from the parameter estimates.
#' compute_mlpin(opt_out$par)
#' 
#' @seealso
#' * [estimate_mlekop()] for the calling function.
#' 
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
"compute_mlpin" <- function(par)
{
  # Transform logits.
  alpha <- logistic(par[1])
  epsilon <- par[2]
  if (length(par) == 4)
  {
    # Original EKOP model.
    mu <- par[4]
  }
  else {
    mu <- par[3]
  }
  pin <- compute_pin(alpha, epsilon, mu)
  
  col_names            <- c('pin_ml')
  estimates            <- data.frame(pin)
  colnames(estimates)  <- col_names
  row.names(estimates) <- "Estimates"
  
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

#' Estimates the probability of informed trading with maximum likelihood
#' 
#' @description 
#' Calling [estimate_mlekop()] estimates the parameters from the model of Easley
#' et al. (1996, EKOP) together with the probability of informed trading (PIN).
#' As described in the EKOP model trade data for buy and sell trades are needed,
#' respectively. Note that in contrast the compressed EKOP model needs only data
#' for the number of trades per day (see the function references below).
#' 
#' @details 
#' ## Maximum Likelihood functions
#' There exist different likelihood variants that can be used in maximum 
#' likelihood estimation of the EKOP model.
#' * `compute_ekop_orig_lik` uses the likelihood function proposed in the original 
#'   paper of Easley et al. (1996).
#' * `compute_ekop_lik` uses a likelihood that is modified in regard to deal 
#'   better with computational overflow. This likelihood function had been 
#'   presented in slightly modified version by Easley et al. (2002).
#' Furthermore, to deal with undefined function values like `NaN` or `inf` an
#' approximation method can be chosen by the argument `methodLik`. Choosing
#' `"approx"` then approximates values of `NaN`, `-inf` or `inf` by values
#' `1e+6`, `-1e+6`, and `1e+6`, respectively and basically ignores the
#' occurrence of undefined values.
#' 
#' ## Optimization algorithms 
#' Optimization is performed gradient-based by using
#' the \code{\link[stats]{optim}} function. The algorithm used is the `L-BFGS-B`
#' that allows for parameter restrictions. This is necessary because we use for
#' the probabilities `alpha` and `delta` the logistic transformation
#' `exp()/(1+exp())` to produce valid probabilities in optimization.
#' 
#' In case that the gradient-based algorithm does not converge in between 100
#' steps, a gradient-free optimization is applied. For gradient-free
#' optimization \code{\link[dfoptim]{nmkb}}, a bounded `Nelder-Mead` algorithm
#' is used. Derivative-free optimization is only performed after gradient-based
#' optimization did not converge and the argument `grad_free` is `TRUE`.
#' 
#' The argument `fnscale` can be used to scale the likelihood function in case 
#' of very large values (e.g. very large volumes) that might lead to number 
#' overflow during computation. Note, `fnscale` must always be negative as the 
#' likelihood function should be maximized.  
#' 
#' @param data A `data.frame` containing the number of buyer- and
#'   seller-initiated trades. The data must be ordered in columns beginning with
#'   the number of mis-specified buys, mis-specified sells, number of buys,
#'   number of sells, and finally the sum of trades per day. See for an example
#'   [simulateEKOP()].
#' @param startpar A vector containing start parameters for maximum likelihood 
#'   estimation. These must be starting values for the logit of alpha, epsilon,
#'   the logit of delta, and mu. If no starting values are provided the function 
#'   chooses the values (0, mean(trades)x.75/2, 0, mean(trades)x.25/2). 
#' @param T A double specifying the minutes of a trading day.
#' @param methodLik A character specifying, if undefined function values in 
#'   optimization should be approximated by large defined values (`1e+6`).
#'   This can help to make maximum likelihood estimation more stable. 
#' @param fnLik A character specifying which likelihood function to use. Either
#'   the original function by Easley et al. (1996) or the slightly modified 
#'   variant of Easley et al. (2002) can be used. The latter one is known to 
#'   also work better with large trading volumes. 
#' @param fnscale A negatve double specifying a scaling factor for the
#'   likelihood function. This can in some cases help when the algorithm does
#'   not reach convergence or suffer from number overflow.
#' @param trace An integer specifying which level of tracing should be used. 
#'   see `?optim` for more details. 
#' @param grad_free A logical indicating if gradient-free optimization should 
#'   be used when gradient descent did not converge. If `TRUE` the optimization 
#'   procedure \code{\link[dfoptim]{nmkb}} is used. 
#' @param return_opt A logical indicating, if in addition to the PIN estimates 
#'   also the results from the optimization procedure should be returned. If 
#'   `TRUE` a list is returned with an element `pin_estimates` holding the 
#'   PIN estimates and an element `opt_results` holding the output of the 
#'   optimization procedure.
#' @param opt_out (Deprecated) A logical indicating if only the output of the 
#'   optimization procedure should be returned. Some older applications still 
#'   rely on this output. In the next version this feature will be removed. 
#'   Note that the default value is `TRUE`.  
#' @return A `list` with all components as returned by 
#'   \code{\link[stats]{optim}} or \code{\link[dfoptim]{nmkb}}.
#' @export
#' 
#' @examples 
#' # Simulate data from the EKOP model. 
#' trades_data <- simulateEKOP()
#' # Estimate the EKOP model by maximum likelihood.
#' pin_estml <- estimate_mlekop(trades_data, methodLik="approx", 
#'                              fnLik="computeEKOPOrigLik", opt_out=FALSE)
#'                    
#' @seealso
#' * [estimate_pin()] for estimating the PIN with a Bayesian approach that 
#'   needs only the total number of trades
#' * [estimate_compml()] for estimating the PIN with the compressed EKOP model 
#'   that needs only the total number of trades
#' * [compute_ekop_lik()] for the implementation of the likelihood function of 
#'   the paper of Easley et al. (2002)
#' * [compute_ekop_orig_lik()] for the implementation of the likelihood function 
#'   of the paper of Easley et al. (1996)
#'   
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Easley, David, Hvidkjaer, Soeren, and O’Hara, Maureen (2002). 
#'   “Is Information Risk a Determinant of Asset Returns?” In: The Journal of 
#'   Finance 57.5, pp. 2185–2221. DOI: 10.1111/1540-6261.00493.
"estimate_mlekop" <- function(data, startpar, T = 390, 
                              methodLik = c("precise", "approx"),
                              fnLik = c("compute_ekop_lik", 
                                        "compute_ekop_orig_lik"), 
                              fnscale=-1, trace=0, grad_free=TRUE,
                              return_opt=FALSE, opt_out = TRUE) 
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
  
  # If the old version of the function is expected 
  # return immediately. 
  if (opt_out) return(optim_res)
  
  # Compute the PIN from the parameter estimates.
  pin_estimates <- compute_mlpin(optim_res$par)
  if (return_opt) {
    return(list(pin_estimates=pin_estimates, opt_results=optim_res))
  } else
  {
    return(pin_estimates)
  }
}

"estimate_compml" <- function(data, startpar, T = 390, 
                              methodLik = c("precise", "approx"),
                              fnscale=-1, trace=0, grad_free=TRUE,
                              return_opt=FALSE, opt_out = TRUE) 
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
  optim_fn 	    <- compute_comp_lik
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
  
  # If the old version of the function is expected 
  # return immediately. 
  if (opt_out) return(optim_res)
  
  # Compute the PIN from the parameter estimates.
  pin_estimates <- compute_mlpin(optim_res$par)
  if (return_opt) {
    return(list(pin_estimates=pin_estimates, opt_results=optim_res))
  } else
  {
    return(pin_estimates)
  }

}