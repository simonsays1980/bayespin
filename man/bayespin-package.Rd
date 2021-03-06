% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\docType{package}
\name{bayespin-package}
\alias{bayespin}
\alias{bayespin-package}
\title{bayespin: An R package for Bayesian estimation of the probability of informed 
  trading}
\description{
An R package for Bayesian estimation of the probability of informed 
  trading from a finite mixture distribution. The original model by Easley et 
  al. (1996) can be converted into a compressed model that is a finite mixture 
  distribution as has been shown by Grammig et al. (1996). This package 
  implements the Bayesian estimation of the compressed model together with the 
  traditional approaches using maximum likelihood. The package uses C++ code 
  and performs a single estimation in around 4-5 seconds.
}
\details{
\code{bayespin} implements the statistical methods for estimating the probability
of informed trading (PIN) with a Bayesian approach as proposed by Grammig et al.
(2015). This should simplify the usage of this rather complicated estimation
procedure and offers researchers an API that is easy to integrate, stable,
and fast in performance.

The model by Grammig et al. (2015) comes along with some advantages in
comparison to the original model of Easley et al. (1996) and other Bayesian
approaches found in literature:
\enumerate{
\item It uses only the number of trades per day instead of the number of
seller- and buyer-initiated trades used by other approaches. This enables
the researcher to collect data more easily - also for historical horizons
and leads to less bias in case trade initiation must be estimated by
using the Lee and Ready (1991) algorithm.
\item The Bayesian estimation of the PIN measure is found to be more stable
especially when it comes to very large trading volumes as they occur on
modern markets.
\item Especially in settings where the rates of informed trading, \code{mu} and/or
the probability of information events are very small Bayesian estimation of
the underlying finite mixture distribution leads to more robust parameter
estimates.
}

The package makes use of high-performance C++ algorithms for MCMC sampling
for finite mixture distributions offered by the
\href{https://github.com/simonsays1980/finmix}{\code{finmix}}
package. Model estimation with a simple \code{K-means} relabeling takes around
4-6 seconds.
\subsection{Implementation of other estimation approaches}{

In addition to the Bayesian estimation approach from Grammig et al. (2015)
the \code{bayespin} package also implements several other methods to compute
the probability of informed trading:
\itemize{
\item the original maximum likelihood procedure of the model by Easley et al.
(1996),
\item the maximum likelihood procedure of the model by Jackson (2007) that also
uses solely the number of trades per trading day (this is similar to
Grammig et al. (2015))
These models were implemented to ease their use by researchers and to enable
comparisons between different models and estimation approaches.
}
}
}
\references{
\itemize{
\item Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the
Probability of Informed Trading. Conference on Financial Econometrics &
Empirical Asset Pricing 2016, Lancaster University
\item Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity,
information, and infrequently traded stocks. Journal of Finance 51,
1405–1436.
\item Jackson, D., 2007. Infering trader behavior from transaction data: A trade
count model. Journal of Computational and Graphical Statistics 12, 55-79.
\item Lee, C., Ready, M. J., 1991. Inferring trade direction from intraday data.
The Journal of Finance 46, 733-746.
}
}
