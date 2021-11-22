[![R-CMD-check](https://github.com/simonsays1980/bayespin/actions/workflows/r-cmd-check.yml/badge.svg?branch=documentation)](https://github.com/simonsays1980/bayespin/actions/workflows/r-cmd-check.yml) [![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html) [![GitHub tag](https://img.shields.io/github/tag/simonsays1980/bayespin.svg)](https://GitHub.com/simonsays1980/bayespin/tags/)
[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg) [![macOS](https://svgshare.com/i/ZjP.svg)](https://svgshare.com/i/ZjP.svg) [![Windows](https://svgshare.com/i/ZhY.svg)](https://svgshare.com/i/ZhY.svg)



# bayespin
**An R package for Bayesian estimation of the probability of informed trading**

`bayespin` bayespin implements the statistical methods for estimating the
probability of informed trading (PIN) with a Bayesian approach as proposed by
Grammig et al. (2015). This should simplify the usage of this rather complicated
estimation procedure and offers researchers an API that is easy to integrate,
stable, and fast in performance.

The estimation method of Grammig et al. (2015) offers some advantages in
comparison to the original model of Easley et al. (1996) and other Bayesian
approaches found in literature:

1. It uses only the number of trades per day instead of the number of seller-
and buyer-initiated trades used by other approaches. This enables the researcher
to collect data more easily - also for historical periods and in turn leads to
less bias in case trade initiation had to be estimated by using the Lee and
Ready (1991) algorithm or similar procedures.

2. The Bayesian estimation of the PIN measure is found to be more stable,
especially when it comes to very large trading volumes as they occur regularly
on modern markets today.

3. Especially in settings where the rates of informed trading, $\mu$ and/or the
probability of information events, $\alpha$ are very small Bayesian estimation
of the underlying finite mixture distribution leads to more robust parameter
estimates.

The package makes use of high-performance C++ algorithms for MCMC sampling of
finite mixture distributions offered by the
[`finmix`](https://github.com/simonsays1980/finmix) package. Model estimation
with a simple `K-means` relabeling takes around 4-6 seconds.

## Implementation of other estimation approaches
In addition to the Bayesian estimation approach from Grammig et al. (2015) the
`bayespin` package also implements several other methods to estimate the
probability of informed trading:

* The original maximum likelihood procedure of the model by Easley et al. (1996).
* The maximum likelihood procedure of the model by Easley et al. (1996) using a
  variation of the likelihood function proposed in Easley et al. (2002).
* The maximum likelihood procedure of the model by Jackson (2007) that also uses
solely the   number of trades per trading day (this is similar to Grammig et al.
(2015)).

These models were implemented to ease their use for researchers and to enable
comparisons between different models and estimation approaches.

## Installation
The package can be installed directly from GitHub by using the function
`install_github()` in the `devtools` package. The package passed all checks from
`R CMD check` on all major platforms and hence, should be installable on MacOS
X, Windows, and Linux. Be sure that you installed appropriate developer tools
for your platform as a C++ compiler for the source code is needed.

Note that installation of the dependencies can take some time as `bayespin` 
depends on the `finmix` package and needs to compile the C++ code therein. 

### MacOS
For MacOS the XCode Command Line Tools are needed. You should have installed
these when installing `R`. See the
[MacOSX-FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
for more information on how to install source packages on MacOS.

### Windows
For Windows the [`rtools`](https://cran.r-project.org/bin/windows/Rtools/)
package is needed. Follow the link and install this package, if you have not
installed it, yet.

## Quick start
To start, simulate data and then estimate the model by Grammig et al. (2015) and
compare results to outputs of maximum likelihood estimation of the original
model of Easley et al. (1996):

```
# Set the random seed so results can be replicated.
set.seed(42)
# Simulate trades data from the model by Easley et al. (1996).
trades_data <- simulate_ekop(size = 1000, alpha = .3, epsilon = .3,
                             delta = .5, mu = .1, T = 60*6.5)
# Show first lines of data.
head(trades_data)
  MisBuy MisSell Buy Sell Trades
1      0       0 175  167    342
2      0       0 163  161    324
3      0       0 163  147    310
4      0       0 141  172    313
5      0       0 176  156    332
6      0       0 154  163    317

# Estimate the model of Grammig et al. (2015). 
bayesian_pin <- estimate_pin(trades_data$Trades)
# Show results.
bayesian_pin
          alpha   epsilon         mu        pin
MAP   0.3392890 0.2998388 0.09379760 0.05039492
BML   0.3356769 0.2999046 0.09406135 0.05000800
IEAVG 0.3402887 0.2998041 0.09369827 0.05049063

# Estimate the original model by Easley et al. (1996).
ml_pin <- estimate_mlekop(trades_data, methodLik="approx", 
                          fnLik = "compute_ekop_orig_lik", opt_out=FALSE)
# Show results.
ml_pin
       alpha   epsilon    delta         mu        pin
ML 0.3211642 0.3000834 0.486291 0.09695398 0.04932346
```
We can see that the original model by Easley et al. (1996) performs better 
parameter estimates. This is not surprising, as if we have more data available 
it helps to use it. Things become interesting, if the buyer- and 
seller-initiated trades suffer from mis-specification (see herefor the 
simulation function `simulate_ekop_mis()`). 

## References
* Grammig, J., Theissen, E., Zehnder, L.S., 2015. Bayesian Estimation of the
Probability of Informed Trading. Conference on Financial Econometrics &
Empirical Asset Pricing 2016, Lancaster University.
* Easley, D., Kiefer, N., O’Hara, M., Paperman, J., 1996. Liquidity,
information, and infrequently traded stocks. Journal of Finance 51, 1405–1436.
* Jackson, D., 2007. Infering trader behavior from transaction data: A trade
count model. Journal of Computational and Graphical Statistics 12, 55-79.
* Lee, C., Ready, M. J., 1991. Inferring trade direction from intraday data. The
Journal of Finance 46, 733-746.

## Some more information
This is a package worked on for years and still not fully implemented. As it is
still maintained by a single author, please by patient with issues.
