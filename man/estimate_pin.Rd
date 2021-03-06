% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{estimate_pin}
\alias{estimate_pin}
\title{Estimates the Bayesian probability of informed trading}
\usage{
estimate_pin(pin_data, return_mcmc = FALSE, with_stephens = FALSE)
}
\arguments{
\item{pin_data}{A \code{vector} of trades data (number of trades per day).}

\item{return_mcmc}{A logical indicating, if MCMC parameter samples should
be returned in addition to the final PIN estimates.}

\item{with_stephens}{A logical indicating, if in addition to \code{kmeans} re-
labeling also re-labeling via the algorithms from Stephens (1997a) and
Stephens (1997b) should be applied.}
}
\value{
A \code{data.frame()} with PIN estimates, if \code{return_mcmc} is \code{FALSE}.
Otherwise, a \code{list} containing the \code{data.frame} of PIN estimates together
with a \code{finmix::mcmcoutput} (see \code{\link[finmix]{mcmcoutput-class}})
object containing the MCMC traces of the component paramaters.
}
\description{
Calling \code{\link[=estimate_pin]{estimate_pin()}} on trades data estimates the probability of
informed trading from the compressed EKOP model presented in Grammig,
Theissen and Zehnder (2015).
}
\details{
Internally this function uses the \code{finmix} package for estimation of the
finite mixture parameters. The \code{finmix} package performs MCMC sampling for
the parameters of the compressed EKOP model and computes three parameter
estimates namely
\itemize{
\item Maximum a posterior (MAP) the mode of the joint posterior of parameters
and indicators,
\item Bayesian maximum likelihood (BML) the mode of the joint posterior of
parameters and indicators in case of a flat prior distribution,
\item Identified ergodic average (IEAVG) is the ergodic average over the MCMC
parameter traces after relabeling the parameters.
}
\subsection{Relabeling}{

The MCMC sampling uses a socalled random permutation algorithm to force the
labels of the two distributions in the mixture to switch randomly between
iterations. Label switching occurs naturally in estimation of finite mixture
distributions and is documented in the literature. Random permutation tries
not to eliminate label switching, but to bring it into a more controlled
setting. As a result each MCMC parameter trace contains parameter values
of the other trace. For this reason re-labeling is performed after sampling.

The \code{finmix} package offers three different re-labeling algorithms:
\itemize{
\item \code{kmeans} re-labels parameters by applying K-means clustering to the
merged component parameters,
\item \code{Stephens1997a} applies the algorithm suggested by Stephens (1997a) that
tries to make the posterior marginal distributions of parameters uni-modal,
\item \code{Stephens1997b} uses the algorithm presented in Stephens (1997b) that
tries to re-label parameter values by minimizing the Kullback-Leibler
distance between a parameter clustering via assignment probabilities and
the classification probabilities defined by Poisson distributions.
}

To offer the user more flexibility the function can also return the MCMC
samples for further inspection. Setting \code{return_mcmc} to \code{TRUE} (default)
returns also the MCMC traces.
}
}
\examples{
# Simulate trades data
trades_data <- simulate_ekop()
# Estimate the Bayesian PIN.
estimate_pin(trades_data$Trades)

}
\references{
\itemize{
\item Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity,
information, and infrequently traded stocks. Journal of Finance 51,
1405–1436.
\item Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the
Probability of Informed Trading. Conference on Financial Econometrics &
Empirical Asset Pricing 2016, Lancaster University
\item Stephens, M., 1997a. Discussion 'On Bayesian analysis of mixtures with an
unknown number of components' (by S. Richardson and P.J. Green), J. R.
Statist. Soc., B, 59, 768-769.
\item Stephens, M., 1997b. Bayesian methods for mixtures of normal
distributions. DPhil. Thesis. University of Oxford, Oxford.
}
}
\seealso{
\itemize{
\item \code{\link[finmix]{mixturemcmc}} for performing MCMC sampling
\item \code{\link[finmix]{mcmcestimate}} for estimating parameters from MCMC samples
\item \code{\link[=compute_bayespin]{compute_bayespin()}} for computing the Bayesian PIN from Bayesian
estimates of the component parameters of the compressed EKOP model
}
}
