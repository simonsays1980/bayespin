"compute_comp_lik" <- function(data, par, T, methodLik = c("precise", "approx")) {
  
  ## Check if correct argument is given ##
  method <- match.arg(methodLik)
  
  ## Transform parameters ##
  lalpha 		<- logistic(par[1])
  lepsilon 	<- par[2]
  lmu 		  <- par[3]
  
  ## Allocate containers to hold values ##
  tradeter1 	<- matrix(0, NROW(data), 1)
  tradeter2 	<- matrix(0, NROW(data), 1)
  
  ## Compute Likelihood ##
  for (i in 1:NROW(data)) {
    tradeter1[i] <- data[i] * log(2 * lepsilon * T) - lgamma(data[i] + 1)
    tradeter2[i] <- data[i] * log((2 * lepsilon + lmu) * T) - lgamma(data[i] + 1)
  }
  
  epster1 	<- (-2 * lepsilon * T)
  epster2 	<- (-(2 * lepsilon + lmu) * T) 
  sum1 		<- (1 - lalpha) * exp(epster1 + tradeter1)
  sum2 		<- lalpha * exp(epster2 + tradeter2)
  
  likl <- sum(log((sum1 + sum2)))
  if(method == "approx") {
    likl[likl == -Inf] <- -1e+6
    likl[likl == Inf]  <- 1e+6
    likl[likl == NaN]  <- 1e+6
  }
  
  return(likl)
}

"compute_ekop_orig_lik" <- function(data, par, T, methodLik = c("precise", "approx"))
{
  ## Check if correct argument is given ##
  method <- match.arg(methodLik)
  
  ## Transform parameters ##
  lalpha      <- logistic(par[1]) 
  ldelta      <- logistic(par[3])
  lepsilon    <- par[2]
  lmu         <- par[4]
  
  ## Data variables ##
  B           <- data[, 3]
  S           <- data[, 4]
  
  ## Allocate containers to hold values ##
  tradeter1   <- matrix(0, NROW(data), 1)
  tradeter2   <- matrix(0, NROW(data), 1)
  tradeter3   <- matrix(0, NROW(data), 1)
  tradeter4   <- matrix(0, NROW(data), 1)
  
  ## Compute likelihood ##
  for (i in 1:NROW(data)) {
    tradeter1[i]    <- B[i] * log(lepsilon * T) - lgamma(B[i] + 1)
    tradeter2[i]    <- S[i] * log(lepsilon * T) - lgamma(S[i] + 1)
    tradeter3[i]    <- B[i] * log((lepsilon + lmu) * T) - lgamma(B[i] + 1)
    tradeter4[i]    <- S[i] * log((lepsilon + lmu) * T) - lgamma(S[i] + 1)
  }
  
  epster1     <- (-lepsilon * T) 
  epster2     <- (-(lepsilon + lmu) * T)
  sum1        <- (1 - lalpha) * exp(epster1 + tradeter1 + epster1 + tradeter2)
  sum2        <- lalpha * ldelta * exp(epster1 + tradeter1 + epster2 + tradeter4)
  sum3        <- lalpha * (1 - ldelta) * exp(epster2 + tradeter3 + epster1 + tradeter2)
  
  likl        <- sum(log(sum1 + sum2 + sum3))
  if (method == "approx") {
    likl[likl == -Inf] <- -1e+6
    likl[likl == Inf]  <- 1e+6
    likl[likl == NaN]  <- 1e+6
  }
  return(likl)
}

#' Computes the modified likelihood function from the EKOP model
#' 
#' @description
#' Calling [compute_ekop_lik()] computes the modified likelihood function from 
#' the model of Easley et al. (1996). The modifications are presented in the 
#' paper by Easley et al. (2002) and are applied to the scenario of equal 
#' buy and sell rates in the model. 
#' 
#' Optimization is performed over the logit transformation of the ratios `alpha`
#' and `delta` and therefore these parameter logits get transformed in the 
#' likelihood function via the logistic transformation `exp()/(1+exp())`.  
#' 
#' @param data A `data.frame` containing the trades data in the following order:
#'   mis-specified buys, mis-specified sells, buyer-initiated trades, 
#'   seller-initiated trades and finally the total number of trades on a 
#'   specific trading day.
#' @param par A vector specifying the parameter values at which the function 
#'   should be evaluated. The parameter order is `alpha`, `epsilon`, `delta`, 
#'   and `mu`.
#' @param T A double indicating the length of a trading day in minutes. 
#' @param A character specifying, if undefined function values in optimization
#'   should be approximated by large defined values (`1e+6`). This can help to
#'   make maximum likelihood estimation more stable.
#' @return A double holding the likelihood value of the data and parameter 
#'   values passed in.
#' @export
#' 
#' @seealso 
#' * [estimate_mlekop()] for the calling function.
#' 
#' @references
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Easley, David, Hvidkjaer, Soeren, and O’Hara, Maureen (2002). 
#'   “Is Information Risk a Determinant of Asset Returns?” In: The Journal of 
#'   Finance 57.5, pp. 2185–2221. DOI: 10.1111/1540-6261.00493.
"compute_ekop_lik" <- function(data, par, T, methodLik = c("precise", "approx"))
{
  ## Check if correct argument is given ##
  method <- match.arg(methodLik)
  
  ## Transform parameters ##
  lalpha      <- logistic(par[1]) 
  ldelta      <- logistic(par[3])
  lepsilon    <- par[2]
  lmu         <- par[4]
  x           <- lepsilon/(lepsilon + lmu)
  
  buys        <- data[,3]
  sells       <- data[,4]
  trades      <- data[,5]
  M           <- as.matrix(apply(data[,3:4], 1, min) 
                           + apply(data[, 3:4], 1, max)/2)
  
  ## Compute Likelihood ##
  ter1    <- (-2) * (lepsilon * T) + M * log(x) + 
    trades * log((lmu + lepsilon) * T)
  ter2    <- lalpha * (1 - ldelta) * exp((-1) * lmu * T) * x^(sells - M)
  ter2    <- ter2 + lalpha * ldelta * exp((-1) * lmu * T) * x^(buys - M) + 
    (1 - lalpha) * x^(trades - M)
  
  likl    <- sum(ter1, na.rm = TRUE) + sum(log(ter2), na.rm = TRUE)
    
  if (method == "approx") {
    likl[likl == -Inf] <- -1e+6
    likl[likl == Inf]  <- 1e+6
    likl[likl == NaN]  <- 1e+6
  }

  return(likl)
}

#' Transforms values by the logistic function
#' 
#' @description
#' Calling [logistic()] applies the logistic transformation `exp()/(1+exp())` 
#' to the passed in arguments. 
#' 
#' @param value A double or a vector of doubles to be transformed. 
#' @return A double or vector of doubles containing the logistic 
#'   transformations of the passed in arguments.
#' @export
#' 
#' @examples
#' logistic(-1.5)
#' 
#' @seealso 
#' * [computeEKOPLik()] for a function applying the transformation during 
#'   optimization of parameters.
"logistic" <- function(value)
{
  return(exp(value)/(1 + exp(value)))
}