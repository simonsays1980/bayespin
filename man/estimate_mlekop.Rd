% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{estimate_mlekop}
\alias{estimate_mlekop}
\title{Estimates the probability of informed trading with maximum likelihood}
\usage{
estimate_mlekop(
  data,
  startpar,
  T = 390,
  methodLik = c("precise", "approx"),
  fnLik = c("compute_ekop_lik", "compute_ekop_orig_lik"),
  fnscale = -1,
  trace = 0,
  grad_free = TRUE,
  return_opt = FALSE,
  opt_out = TRUE
)
}
\arguments{
\item{data}{A \code{data.frame} containing the number of buyer- and
seller-initiated trades. The data must be ordered in columns beginning with
the number of mis-specified buys, mis-specified sells, number of buys,
number of sells, and finally the sum of trades per day. See for an example
\code{\link[=simulate_ekop]{simulate_ekop()}}.}

\item{startpar}{A vector containing start parameters for maximum likelihood
estimation. These must be starting values for the logit of alpha, epsilon,
the logit of delta, and mu. If no starting values are provided the function
chooses the values (0, mean(trades)x.75/2, 0, mean(trades)x.25/2).}

\item{T}{A double specifying the minutes of a trading day.}

\item{methodLik}{A character specifying, if undefined function values in
optimization should be approximated by large defined values (\code{1e+6}).
This can help to make maximum likelihood estimation more stable.}

\item{fnLik}{A character specifying which likelihood function to use. Either
the original function by Easley et al. (1996) or the slightly modified
variant of Easley et al. (2002) can be used. The latter one is known to
also work better with large trading volumes.}

\item{fnscale}{A negatve double specifying a scaling factor for the
likelihood function. This can in some cases help when the algorithm does
not reach convergence or suffer from number overflow.}

\item{trace}{An integer specifying which level of tracing should be used.
see \code{?optim} for more details.}

\item{grad_free}{A logical indicating if gradient-free optimization should
be used when gradient descent did not converge. If \code{TRUE} the optimization
procedure \code{\link[dfoptim]{nmkb}} is used.}

\item{return_opt}{A logical indicating, if in addition to the PIN estimates
also the results from the optimization procedure should be returned. If
\code{TRUE} a list is returned with an element \code{pin_estimates} holding the
PIN estimates and an element \code{opt_results} holding the output of the
optimization procedure.}

\item{opt_out}{(Deprecated) A logical indicating if only the output of the
optimization procedure should be returned. Some older applications still
rely on this output. In the next version this feature will be removed.
Note that the default value is \code{TRUE}.}
}
\value{
A \code{list} with all components as returned by
\code{\link[stats]{optim}} or \code{\link[dfoptim]{nmkb}}.
}
\description{
Calling \code{\link[=estimate_mlekop]{estimate_mlekop()}} estimates the parameters from the model of Easley
et al. (1996, EKOP) together with the probability of informed trading (PIN).
As described in the EKOP model trade data for buy and sell trades are needed,
respectively. Note that in contrast the compressed EKOP model needs only data
for the number of trades per day (see the function references below).
}
\details{
\subsection{Maximum Likelihood functions}{

There exist different likelihood variants that can be used in maximum
likelihood estimation of the EKOP model.
\itemize{
\item \code{compute_ekop_orig_lik} uses the likelihood function proposed in the original
paper of Easley et al. (1996).
\item \code{compute_ekop_lik} uses a likelihood that is modified in regard to deal
better with computational overflow. This likelihood function had been
presented in slightly modified version by Easley et al. (2002).
Furthermore, to deal with undefined function values like \code{NaN} or \code{inf} an
approximation method can be chosen by the argument \code{methodLik}. Choosing
\code{"approx"} then approximates values of \code{NaN}, \code{-inf} or \code{inf} by values
\code{1e+6}, \code{-1e+6}, and \code{1e+6}, respectively and basically ignores the
occurrence of undefined values.
}
}

\subsection{Optimization algorithms}{

Optimization is performed gradient-based by using
the \code{\link[stats]{optim}} function. The algorithm used is the \code{L-BFGS-B}
that allows for parameter restrictions. This is necessary because we use for
the probabilities \code{alpha} and \code{delta} the logistic transformation
\code{exp()/(1+exp())} to produce valid probabilities in optimization.

In case that the gradient-based algorithm does not converge in between 100
steps, a gradient-free optimization is applied. For gradient-free
optimization \code{\link[dfoptim]{nmkb}}, a bounded \code{Nelder-Mead} algorithm
is used. Derivative-free optimization is only performed after gradient-based
optimization did not converge and the argument \code{grad_free} is \code{TRUE}.

The argument \code{fnscale} can be used to scale the likelihood function in case
of very large values (e.g. very large volumes) that might lead to number
overflow during computation. Note, \code{fnscale} must always be negative as the
likelihood function should be maximized.
}
}
\examples{
# Simulate data from the EKOP model. 
trades_data <- simulate_ekop()
# Estimate the EKOP model by maximum likelihood.
pin_estml <- estimate_mlekop(trades_data, methodLik="approx", 
                             fnLik="compute_ekop_orig_lik", opt_out=FALSE)
                   
}
\references{
\itemize{
\item Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity,
information, and infrequently traded stocks. Journal of Finance 51,
1405–1436.
\item Easley, David, Hvidkjaer, Soeren, and O’Hara, Maureen (2002).
“Is Information Risk a Determinant of Asset Returns?” In: The Journal of
Finance 57.5, pp. 2185–2221. DOI: 10.1111/1540-6261.00493.
}
}
\seealso{
\itemize{
\item \code{\link[=estimate_pin]{estimate_pin()}} for estimating the PIN with a Bayesian approach that
needs only the total number of trades
\item \code{\link[=estimate_compml]{estimate_compml()}} for estimating the PIN with the compressed EKOP model
that needs only the total number of trades
\item \code{\link[=compute_ekop_lik]{compute_ekop_lik()}} for the implementation of the likelihood function of
the paper of Easley et al. (2002)
\item \code{\link[=compute_ekop_orig_lik]{compute_ekop_orig_lik()}} for the implementation of the likelihood function
of the paper of Easley et al. (1996)
}
}
