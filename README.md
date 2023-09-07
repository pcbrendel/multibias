# multibias
An R package for multi-bias analysis corresponding to the [article](https://doi.org/10.1093/ije/dyad001):

Paul Brendel and others, Simultaneous adjustment of uncontrolled confounding, selection bias and misclassification in multiple-bias modelling, *International Journal of Epidemiology*, Volume 52, Issue 4, Pages 1220–1230

## Overview

The `multibias` package includes a suite of functions that provide odds ratio estimates that are adjusted for any combination of uncontrolled confounding (uc), selection bias (sel), and exposure misclassification (emis):

  - `adjust_uc_sel()` adjusts for uncontrolled confounding and selection bias.
  - `adjust_emc_sel()` adjusts for exposure misclassification and selection bias.
  - `adjust_uc_emc()` adjusts for uncontrolled confounding and exposure misclassificaiton.
  - `adjust_multinom_uc_emc()` adjusts for uncontrolled confounding and exposure misclassificaiton (with the bias models for the uncontrolled confounder and true exposure jointly modeled via a multinomial regression).
  - `adjust_uc_emc_sel()` adjusts for all three biases.
  - `adjust_multinom_uc_emc_sel()` adjusts for all three biases (with the bias models for the uncontrolled confounder and true exposure jointly modeled via a multinomial regression).
 
If you are new to bias analysis, I'd suggest reading [Applying Quantitative Bias Analysis to Epidemiologic Data](https://www.springer.com/us/book/9780387879604) textbook or visiting my [website](https://www.paulbrendel.com/).

## Installation

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("pcbrendel/multibias")
```

## Example

We are interested in quantifying the effect of exposure X on outcome Y. The causal system can be represented in the following directed acyclic graph (DAG):

![uc_mc_sel_DAG](DAGs/uc_mc_sel_DAG.png)

The variables are defined:
 - X: true, unmeasured exposure
 - Y: outcome
 - C: measured confounder(s)
 - U: unmeasured confounder
 - X*: misclassified, measured exposure
 - S: study selection

It can be seen from this DAG that the data suffers from three sources of bias. There is uncontrolled confounding from (unobserved) variable U. The true exposure, X, is unobserved, and the misclassified exposure X* is dependent on both the exposure and outcome. Lastly, there is collider stratification at variable S since exposure and outcome both affect selection. The study naturally only examines those who were selected (i.e. those with S=1).

A simulated data set corresponding to this DAG, `df_uc_emc_sel` can be loaded from the multibias package. 

```{r, eval = TRUE}
library(multibias)
head(df_uc_emis_sel)
#>   Xstar Y C1 C2 C3
#> 1     0 1  1  0  1
#> 2     1 0  0  0  1
#> 3     1 1  0  0  1
#> 4     0 0  0  1  1
#> 5     0 0  0  0  0
#> 6     0 0  0  1  0
```

In this data set, the true, unbiased exposure-outcome odds ratio (OR<sub>YX</sub>) equals ~2. However, when we run a logistic regression of the outcome on the exposure and confounders, we do not observe an odds ratio of 2 due to the multiple bias sources.

```{r, eval = TRUE}
biased_model <- glm(Y ~ Xstar + C1 + C2 + C3, data = df_uc_emis_sel, family = binomial(link = "logit"))
biased_or_yx <- exp(coef(biased_model)[2])
round(biased_or_yx, 2)
#>    Xstar
#> 1.66
```

The `adjust` family of functions serves to "reconstruct" the unbiased data and return the exposure-outcome odds ratio that would be observed in the unbiased setting.

Models for the missing variables (U, X, S) are used to facilitate this data reconstruction. For the above DAG, the corresponding bias models are:
 - logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y
 - logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>
 - logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>j</sub>

where j indicates the number of measured confounders. 

To perform the bias adjustment, it is necessary to obtain values of these bias parameters. Potential sources of these bias parameters include internal validation data, estimates in the literature, and expert opinion. For purposes of demonstrating the methodology, we will obtain the exact values of these bias parameters. This is possible because we have access to the data of missing values that would otherwise be absent in real-world practice.

```{r, eval = TRUE}
u_model <- ''
x_model <- ''
s_model <- ''
```

We will run the analysis over 1,000 bootstrap samples to obtain a valid confidence interval. To improve performance we will run the for loop in parallel using the `foreach()` function in the `doParallel` package. First, we create a cluster, make a seed, and specify the desired number of bootstrap repitions.

```{r, eval = TRUE}
library(doParallel)

no_cores <- detectCores() - 1
registerDoParallel(cores = no_cores)
cl <- makeCluster(no_cores)

set.seed(1234)
nreps <- 1000
est <- vector(length = nreps)
```

Then we run the parallel for loop in which we apply the `adjust_uc_emc_sel()` function to bootstrap samples of the `df_uc_emc_sel` data. In this function we specify the following arguments: the data, the exposure variable, the outcome variable, the confounder(s), the parameters in the U model, the parameters in the X model, and the parameters in the S model. Since knowledge of the complete data was known, the correct bias parameters were known in advance. The bias parameters can be provided as fixed values, as seen in this example, or values from a probability distribution. This latter strategy is referred to as probabilistic bias analysis.

```{r, eval = TRUE}
# create triple-bias dataframe
df <- df_uc_emis_sel
# parallel for loop to obtain odds ratio estimates
or <- foreach(i = 1:nreps, .combine = c, .packages = 'dplyr') %dopar% {
  # bootstrap sample
  bdf <- df[sample(1:nrow(df), nrow(df), replace = TRUE),]
  # adjust for biases
  est[i] <- adjust_uc_emc_sel(
    bdf, 
    exposure = "Xstar", 
    outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    u_model_coefs = c(-0.40, 0.38, 0.46),
    x_model_coefs = c(-1.61, 2.71, 0.62, -0.41, -0.41, 0.40), 
    s_model_coefs = c(-0.39, 0.40, 0.75, -0.04, -0.04, 0.05)
  )[[1]]
}
```
Finally, we obtain the OR<sub>YX</sub> estimate and 95% confidence interval from the distribution of 1,000 bootstrap odds ratio estimates. As expected, OR<sub>YX</sub> ~ 2, indicating that we were able to obtain an unbiased odds ratio from the biased data.

```{r, eval = TRUE}
# odds ratio estimate
round(median(or), 2)
# confidence interval
round(quantile(or,c(.025, .975)), 2)
#> 2.03
#> 2.01 2.06
```


