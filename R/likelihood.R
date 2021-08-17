"computeKokotLik" <- function(data, par, T, methodLik = c("precise", "approx")) {
  
  ## Check if correct argument is given ##
  method <- match.arg(methodLik)
  
  ## Transform parameters ##
  lalpha 		<- exp(par[1])/(1 + exp(par[1]))
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

"computeEKOPOrigLik" <- function(data, par, T, methodLik = c("precise", "approx"))
{
  ## Check if correct argument is given ##
  method <- match.arg(methodLik)
  
  ## Transform parameters ##
  lalpha      <- exp(par[1])/(1 + exp(par[1]))
  ldelta      <- exp(par[3])/(1 + exp(par[3]))
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

"computeEKOPLik" <- function(data, par, T, methodLik = c("precise", "approx"))
{
  ## Check if correct argument is given ##
  method <- match.arg(methodLik)
  
  ## Transform parameters ##
  lalpha      <- exp(par[1])/(1 + exp(par[1]))
  ldelta      <- exp(par[3])/(1 + exp(par[3]))
  lepsilon    <- par[2]
  lmu         <- par[4]
  x           <- lepsilon/(lepsilon + lmu)
  
  buys        <- data[,3]
  sells       <- data[,4]
  trades      <- data[,5]
  M           <- as.matrix(apply(data[,3:4], 1, min) 
                           + apply(data[, 3:4], 1, max)/2)
  
  ## Compute Likelihood ##
  ter1    <- (-2) * (lepsilon * T) + M * log(x) + trades * log((lmu + lepsilon) * T)
  ter2    <- lalpha * (1 - ldelta) * exp((-1) * lmu * T) * x^(sells - M)
  ter2    <- ter2 + lalpha * ldelta * exp((-1) * lmu * T) * x^(buys - M) + (1 - lalpha) * x^(trades - M)
  
  likl    <- sum(ter1, na.rm = TRUE) + sum(log(ter2), na.rm = TRUE)
  # if (is.nan(likl) || is.infinite(likl))
  # {
  #   cat('\talpha:', lalpha)
  #   cat('\teps:',lepsilon)
  #   cat('\tdelta:',ldelta)
  #   cat('\tmu:',lmu)
  #   cat('\t likl before:', likl)
  # }
    
  if (method == "approx") {
    likl[likl == -Inf] <- -1e+6
    likl[likl == Inf]  <- 1e+6
    likl[likl == NaN]  <- 1e+6
  }

  return(likl)
}