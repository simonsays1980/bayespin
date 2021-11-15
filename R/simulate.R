## Load the dynamic library
#' @useDynLib bayespin
#' @importFrom Rcpp sourceCpp
NULL

#' Simulates trades data
#' 
#' @description
#' Simulates trades data from the model by Easley, Kiefer, O'Hara and Paperman 
#' (1996). The simulation function is implemented in C++ and allows therefore 
#' usage in specific performance-sensitive scenarios, like simulating huge 
#' amounts of trade data or calling this function repeatedly in large 
#' simulations. 
#' 
#' @details 
#' The functions returns a `data.frame` with named columns `MisBuy`, `MisSell`, 
#' `Buy`, `Sell`, and `Trades`. `Buy` and `Sell` are the number of buy and sell 
#' trades per trading day and `Trades` defines the sum of these amounts. 
#' `MisBuy` and `MisSell` are only added to allow for a standardization of 
#' input data to the estimation functions of the `bayespin` package. The 
#' corresponding simulation function [simulateEKOPMis()] simulates 
#' mis-specified trades data and returns in the fields `MisBuy` and `MisSell` 
#' the number of mis-specified buy and sell trades per trading day respectively.
#'  
#' @param size An integer specifying the number of trading days that should
#'   be simulated.
#' @param alpha A double specifying the probability of an information event.
#' @param epsilon A double specifying the arrival rate of uninformed traders.
#' @param delta A double specifying the probability of information being good 
#'   or bad.
#' @param mu A double specifying the arrival rate of informed traders.
#' @param T A double defining the length of a trading day in minutes.
#' @return A `data.frame` containing the trades per day.
#' @export
#' 
#' @examples
#' # Simulate trades data for 1000 trading dates with a trading day lasting 
#' # 6.5 hours.
#' trades_data <- simulateEKOP(size = 1000, alpha = .3, epsilon = .4,
#'                             delta = .5, mu = .04, T = 60*6.5)
#' head(trades_data)
#' 
#' @seealso 
#' * [simulateEKOPMis()] for simulating mis-specified trades data 
#' 
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the 
#'   Probability of Informed Trading. Conference on Financial Econometrics & 
#'   Empirical Asset Pricing 2016, Lancaster University
"simulateEKOP" <- function(size = 1000, alpha = 0.2, epsilon = 0.2, 
                           delta = 0.5, mu = 0.02, T = 60 * 6.5) 
{
  trade.sample <- simulateEKOP_cc(as.integer(size), as.double(alpha), 
                                  as.double(epsilon), as.double(delta), 
                                  as.double(mu), as.double(T))
  cnames <- c("MisBuy", "MisSell", "Buy", "Sell", "Trades")
  trade.sample <- data.frame(trade.sample)
  colnames(trade.sample) <- cnames
  
  return(trade.sample)
}

#' Simulates mis-specified trades data
#' 
#' @description
#' Simulates mis-specified trades data from the model by Easley, Kiefer, O'Hara
#' and Paperman (1996). The simulation function is implemented in C++ and allows
#' therefore usage in specific performance-sensitive scenarios, like simulating
#' huge amounts of trade data or calling this function repeatedly in large
#' simulations.
#' 
#' @details 
#' The functions returns a `data.frame` with named columns `MisBuy`, `MisSell`, 
#' `Buy`, `Sell`, and `Trades`. `Buy` and `Sell` are the number of buy and sell 
#' trades per trading day and `Trades` defines the sum of these amounts. 
#' `MisBuy` and `MisSell` define the number of mis-specified buy and sell 
#' trades per trading day respectively. 
#'  
#' @param size An integer specifying the number of trading days that should
#'   be simulated.
#' @param alpha A double specifying the probability of an information event.
#' @param epsilon A double specifying the arrival rate of uninformed traders.
#' @param delta A double specifying the probability of information being good 
#'   or bad.
#' @param mu A double specifying the arrival rate of informed traders.
#' @param T A double defining the length of a trading day in minutes.
#' @param mis A double specifying the mis-sepcification rate for trades. 
#' @return A `data.frame` containing the trades per day.
#' @export
#' 
#' @examples
#' # Simulate mis-specifed trades data for 1000 trading dates with a trading 
#' # day lasting 6.5 hours.
#' trades_data <- simulateEKOPMis(size = 1000, alpha = .3, epsilon = .4,
#'                                delta = .5, mu = .04, T = 60*6.5, mis = .1)
#' head(trades_data)
#' 
#' @seealso 
#' * [simulateEKOP()] for simulating trades data without mis-specification
#' * [correct_trades()] to reconstruct original trades data from mis-specified
#'   ones
#'   
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the 
#'   Probability of Informed Trading. Conference on Financial Econometrics & 
#'   Empirical Asset Pricing 2016, Lancaster University
"simulateEKOPMis" <- function(size = 1000, alpha = 0.2, epsilon = 0.2, 
                              delta = 0.5, mu = 0.02, T = 60 * 6.5, 
                              mis = 0.15) 
{
  trade.sample <- simulateEKOPMis_cc(as.integer(size), as.double(alpha), 
                                     as.double(epsilon), as.double(delta), 
                                     as.double(mu), as.double(mis), 
                                     as.double(T))		
  cnames <- c("MisBuy", "MisSell", "Buy", "Sell", "Trades")
  trade.sample <- data.frame(trade.sample)
  colnames(trade.sample) <- cnames
  
  return(trade.sample)
}

#' Corrects mis-specified trades data
#'
#' @description 
#' Calling [correct_trades()] on simulated mis-specified trades data corrects 
#' the mis-specifications. This becomes handy when using simulated data in 
#' comparison studies of e.g. estimators, etc.   
#' 
#' @details 
#' Mis-specified trades are recorded in the columns `MisBuy` and `MisSell` for 
#' buy and sell trades respectively. This information is used to re-construct 
#' the well-specified trades data. 
#' 
#' @param trades A `data.frame` containing the trades data. This must have the 
#'   columns `MisBuy`, `MisSell`, `Buy`, and `Sell` specified.
#' @return A `data.frame` containing the re-constructed well-specified trades 
#'   data.
#' @export
#' 
#' @examples
#' # Simulate mis-specifed trades data for 1000 trading dates with a trading 
#' # day lasting 6.5 hours.
#' trades_data <- simulateEKOPMis(size = 1000, alpha = .3, epsilon = .4,
#'                                delta = .5, mu = .04, T = 60*6.5, mis = .1)
#' trades_data_corr <- correct_trades(trades = trades_data)                        
#' head(trades_data_corr)
#' 
#' @seealso
#' * [simulateEKOPMis()] for simulating trades data with mis-specification
#' * [simulateEKOP()] for simulating trades data without mis-specification 
"correct_trades" <- function(trades)
{
  trades$Buy <- trades$Buy - trades$MisBuy + trades$MisSell
  trades$Sell <- trades$Sell - trades$MisSell + trades$MisBuy
  trades$MisBuy <- 0
  trades$MisSell <- 0
  
  return(trades)
}