# ocd_CI

This repository contains implementation of the ocd_CI algorithm described in the submitted manuscript. In this algorithm, we output a confidence interval for the changepoint location and estimate the set of indices of coordinates in which the mean changes, all under online changepoint monitoring. A real data example is also included in the repository.

## Main function
**`ocd_CI.R`** is the main function of this repository. The function requires the 'MASS' library installed. You can do it by running `library('MASS')`. The inputs of the main function are as follows:

Auguments
- `p`: the dimension of each observation;
-  `s`: the sparsity of the post-change vector;
-  `z`: the changepoint location;
- `n_max`: maximum monitoring time;
- `vartheta`: the $\ell_2$ norm of the vector of change;
- `signal_shape`: the shape of the post-change vector, can be one of four: random, uniform, inverse square root or harmonic;
- `spatial`: the covariance matrix structure, can be one of two: identity $\Sigma= I_p$ or Toeplitz $\Sigma = (\rho^{|j-k|})_{j,k}$;
- `rho`: the Toeplitz parameter $\rho$;
- `gamma`: patience level of the online changepoint monitoring;
- `beta`: a lower bound of the $\ell_2$ norm of the vector of change;
- `ocd_thresholds`: declaration thresholds for the ocd base procedure, can be one of two: theoretical choices or Monte Carlo thresholds;
- `l_choice`: extra sample size post declaration, can be zero, or given based on theory, or otherwise specified;
- `l`: user-specified extra sample size;
- `d1_choice`: the value of parameter d1, can be set to default, or two variant values, or otherwise specified;
- `d1`: user-specified d1 value;
- `alpha`: the desired confidence level is 1-alpha;
- `kaul_comp`: whether to compare with a confidence interval procedure based on Kaul et al., (2021); default option is *False*;
- `print_type`: options to print results during the running, can be one of four: CI only, support only, both, or muted.

Details
The function outputs a list of three (when `kaul_comp=FALSE`) elements. The first element of the output list is an ocd_CI confidence interval for the changepoint location; the second is an estimate of the set of indices of coordinates in which the mean changes; the third is an augmented support estimate (with addition of the anchor coordinate). When `kaul_comp=TRUE`, it also outputs a confidence interval constructed via a procedure based on Kaul et al. (2021).

Examples:
`ret <- ocd_CI(p=100, s=10, vartheta=1, signal_shape='random', spatial='identity', beta=1, ocd_thresholds='MC', l_choice='zero', d1_choice='default', alpha=0.05, kaul_comp = TRUE, print_type ='both')`
`ret$CI`
`ret$support`

## Other functions
Other functions are collected in **`ingredients.R`** and **`auxiliary.R`**. As the names suggest, the former contains many functions that are building blocks of the main `ocd_CI` function, while the latter contains other auxiliary functions.

## Tables and figures
The codes for generating all the tables and figures in the submitted manuscript are provided in **`tables.R`**. In order, the readers are able to reproduce Table 1, Table 2, Table S1 (in supplement), Figure 1, and Figure 2 of the paper by simply running the corresponding blocks of codes. Readers are warned that the running time can be very long, as the current implementation does not use parallel computing.

## Real dataset
**`US_weekly_deaths.csv`** is a dataset of weekly deaths in the United States between January 2017 and June 2020 (available at: [https://www.cdc.gov/nchs/nvss/vsrr/covid19/excess_deaths.htm](https://www.cdc.gov/nchs/nvss/vsrr/covid19/excess_deaths.htm)). There are 12506 rows (excl. headers) and 12 columns in this dataset. The columns are:

- **country**: This is always United States;
-  **region**: All 50 states plus Washington, D.C., as well as the entire United States;
- **region code**: The abbreviation of the state;
- **start date**: The start date of a reporting week;
- **end date**: The end date of a reporting week;
- **days**: Days in a week (always 7);
- **year**: The year this reporting week belongs to;
-  **week**: The week number in the year of this reporting week;
-  **population**: The population of the region;
-   **total_deaths**: The total number of deaths in this reporting week;
-  **covid_deaths**: The number of covid deaths in this reporting week;
-  **expected_deaths**: The number of expected deaths in this reporting week;

Among the three columns regarding deaths numbers, we only make use of the **total_deaths** for our analyses.

## Running `ocd_CI` on this dataset
  **`US_excess_deaths.R`** provides the complete codes for running the `ocd_CI` procedure on the above US excess deaths dataset. It also generates Figure 3 of the submitted manuscript. This figure is also available in the repository as **`US_excess_deaths_ocd_CI.pdf`**.
  
## Reference
Main paper
Inference in high-dimensional online changepoint detection. *Preprint*.

Other reference
Kaul, A., Fotopoulos, S. B., Jandhyala, V. K., and Safikhani, A. (2021) Inference on the change point under a high dimensional sparse mean shift. *Electron. J. Statist.*, **15**, 71--134.
