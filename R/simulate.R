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