## Load the dynamic library
#' @useDynLib bayespin
#' @importFrom Rcpp sourceCpp
NULL

#' @details 
#' `bayespin` implements the statistical methods for estimating the probability 
#' of informed trading (PIN) with a Bayesian approach as proposed by Grammig et al. 
#' (2015). This should simplify the usage of this rather complicated estimation 
#' procedure and offers researchers an API that is easy to integrate, stable, 
#' and fast in performance. 
#' 
#' The model by Grammig et al. (2015) comes along with some advantages in 
#' comparison to the original model of Easley et al. (1996) and other Bayesian 
#' approaches found in literature: 
#' 1. It uses only the number of trades per day instead of the number of 
#' seller- and buyer-initiated trades used by other approaches. This enables 
#' the researcher to collect data more easily - also for historical horizons 
#' and leads to less bias in case trade initiation must be estimated by 
#' using the Lee and Ready (1991) algorithm. 
#' 2. The Bayesian estimation of the PIN measure is found to be more stable 
#' especially when it comes to very large trading volumes as they occur on 
#' modern markets. 
#' 3. Especially in settings where the rates of informed trading, `mu` and/or 
#' the probability of information events are very small Bayesian estimation of 
#' the underlying finite mixture distribution leads to more robust parameter 
#' estimates. 
#' 
#' The package makes use of high-performance C++ algorithms for MCMC sampling 
#' for finite mixture distributions offered by the 
#' \href{https://github.com/simonsays1980/finmix}{`finmix`} 
#' package. Model estimation with a simple `K-means` relabeling takes around 
#' 4-6 seconds. 
#' 
#' ## Implementation of other estimation approaches
#' In addition to the Bayesian estimation approach from Grammig et al. (2015) 
#' the `bayespin` package also implements several other methods to compute 
#' the probability of informed trading:
#' * the original maximum likelihood procedure of the model by Easley et al.
#'   (1996),
#' * the maximum likelihood procedure of the model by Jackson (2007) that also 
#'   uses solely the number of trades per trading day (this is similar to 
#'   Grammig et al. (2015))
#' These models were implemented to ease their use by researchers and to enable 
#' comparisons between different models and estimation approaches. 
#' 
#' @references 
#' * Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the 
#'   Probability of Informed Trading. Conference on Financial Econometrics & 
#'   Empirical Asset Pricing 2016, Lancaster University
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Jackson, D., 2007. Infering trader behavior from transaction data: A trade 
#'   count model. Journal of Computational and Graphical Statistics 12, 55-79.
#' * Lee, C., Ready, M. J., 1991. Inferring trade direction from intraday data.
#'   The Journal of Finance 46, 733-746.
"_PACKAGE"  

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
#' corresponding simulation function [simulate_ekop_mis()] simulates 
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
#' trades_data <- simulate_ekop(size = 1000, alpha = .3, epsilon = .4,
#'                              delta = .5, mu = .04, T = 60*6.5)
#' head(trades_data)
#' 
#' @seealso 
#' * [simulate_ekop_mis()] for simulating mis-specified trades data 
#' 
#' @references 
#' * Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity, 
#'   information, and infrequently traded stocks. Journal of Finance 51, 
#'   1405–1436.
#' * Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the 
#'   Probability of Informed Trading. Conference on Financial Econometrics & 
#'   Empirical Asset Pricing 2016, Lancaster University
"simulate_ekop" <- function(size = 1000, alpha = 0.2, epsilon = 0.2, 
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
#' trades_data <- simulate_ekop_mis(size = 1000, alpha = .3, epsilon = .4,
#'                                  delta = .5, mu = .04, T = 60*6.5, mis = .1)
#' head(trades_data)
#' 
#' @seealso 
#' * [simulate_ekop()] for simulating trades data without mis-specification
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
"simulate_ekop_mis" <- function(size = 1000, alpha = 0.2, epsilon = 0.2, 
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
#' trades_data <- simulate_ekop_mis(size = 1000, alpha = .3, epsilon = .4,
#'                                  delta = .5, mu = .04, T = 60*6.5, mis = .1)
#' trades_data_corr <- correct_trades(trades = trades_data)                        
#' head(trades_data_corr)
#' 
#' @seealso
#' * [simulate_ekop_mis()] for simulating trades data with mis-specification
#' * [simulate_ekop()] for simulating trades data without mis-specification 
"correct_trades" <- function(trades)
{
  trades$Buy <- trades$Buy - trades$MisBuy + trades$MisSell
  trades$Sell <- trades$Sell - trades$MisSell + trades$MisBuy
  trades$MisBuy <- 0
  trades$MisSell <- 0
  
  return(trades)
}