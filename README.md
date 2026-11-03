<!-- PROJECT LOGO -->
<br />
<h3 align="center">Local Whittle Estimation for Distributional Time Series: A Regularization Framework in Bayes Hilbert Space</h3>

</div>

<!-- ABOUT THE PROJECT -->
## Abstract
We study local Whittle estimation of long memory in density-valued time series modeled in Bayes Hilbert space. The centered log-ratio (clr) transformation maps densities to a linear Hilbert space, making standard Hilbert space methods available for their analysis. When densities are reconstructed from finite samples, however, small absolute estimation errors can remain large relative to the true density in low-density tails. Taking logarithms can amplify this noise and introduce systematic distortions that affect memory estimation. We develop a framework for regularizing the estimated densities before transformation to improve finite-sample performance while preserving the latent memory target. Using the clr representation, we apply local Whittle estimation to scores along an estimated leading direction of the long-run covariance operator. The framework is based on \emph{local Whittle admissibility}, defined by uniform asymptotic equivalence between the feasible and latent-score criteria. We separate density-estimation error from the population perturbation induced by regularization and give sufficient conditions for admissibility based on density-estimation rates and the temporal dependence of the error components. Under suitable regularity assumptions, these conditions yield consistency for the latent memory parameter. We study additional tail smoothing, a vanishing density floor, and their combination. Monte Carlo experiments show that these procedures can reduce finite-sample bias relative to unregularized kernel density reconstruction. An application to daily distributions of intraday EUR/USD returns illustrates the methodology.


### Code
1. sim_run.R --- R functions to generate required simulation results.
2. real_data.R --- R functions to generate required real data analysis results.

### Data
1. daily_returns.csv --- data for log return distribution data, with raw data obtained following <https://www.kaggle.com/datasets/arashnic/stock-data-intraday-minute-bar/data>.

### Procedure to reproduce simulation results
Step 1. Run simulation.R to reproduce simulation results.

### Procedure to reproduce real data results
Step 1. Run real_data.R to reproduce real data analysis results.




