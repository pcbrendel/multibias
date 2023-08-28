# multibias
An R package for multi-bias analysis corresponding to the article:

Paul Brendel and others, Simultaneous adjustment of uncontrolled confounding, selection bias and misclassification in multiple-bias modelling, *International Journal of Epidemiology*, Volume 52, Issue 4, Pages 1220–1230
[Link](https://doi.org/10.1093/ije/dyad001)

## Overview

multibias is a set of functions that provide odds ratio estimates that are adjusted for any combination of uncontrolled confounding (uc), selection bias (sel), and exposure misclassification (mc):

  - `adjust_uc_sel()` adjusts for uncontrolled confounding and selection bias.
  - `adjust_mc_sel()` adjusts for exposure misclassification and selection bias.
  - `adjust_uc_mc()` adjusts for uncontrolled confounding and exposure misclassificaiton.
  - `adjust_uc_mc2()` adjusts for uncontrolled confounding and exposure misclassificaiton.
  - `adjust_uc_mc_sel()` adjusts for all three biases.
  - `adjust_uc_mc_sel2()` adjusts for all three biases.
 
If you are new to bias analysis, I'd recommend checking out the [Applying Quantitative Bias Analysis to Epidemiologic Data](https://www.springer.com/us/book/9780387879604) textbook or visiting my [website](https://www.paulbrendel.com/).

## Installation

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("pcbrendel/multibias")
```

## Usage

Say you have a data set whose causal relationships can be represented in the following DAG:

![uc_mc_sel_DAG](DAGs/uc_mc_sel_DAG.png)

The variables are defined:
 - X: true, unmeasured exposure
 - Y: outcome
 - C: measured confounder(s)
 - U: unmeasured confounder
 - X*: misclassified, measured exposure
 - S: study selection

It can be seen from this DAG that the data suffers from three sources of bias. There is uncontrolled confounding from (unobserved) variable U. The true exposure, X, is unobserved, and the misclassified exposure X* is dependent on both the exposure and outcome. Lastly, there is collider stratification at variable S since exposure and outcome both affect selection. The study naturally only examines those who were selected (i.e. those with S=1).

A simulated data set corresponding to this DAG, `df_uc_mc_sel` can be loaded from the multibias package. 

```{r, eval = TRUE}
library(multibias)
head(df_uc_mc_sel)
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
biased_model <- glm(Y ~ Xstar + C1 + C2 + C3, data = df_uc_mc_sel, family = binomial(link = "logit"))
exp(coef(biased_model)[2])
#>    Xstar
#> 1.663475
```

The `adjust` family of functions serves to "reconstruct" the unbiased data and return the exposure-outcome odds ratio that would be observed in the unbiased setting.

Models for the missing variables (U, X, S) are used to facilitate this data reconstruction. For the above DAG, these bias models are:
 - logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y
 - logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>
 - logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>j</sub>

where j indicates the number of measured confounders. 

To perform the bias adjustment, it is necessary to obtain values of these bias parameters. Potential sources of these bias parameters include internal validation data, estimates in the literature, and expert opinion.

We will run the analysis over 1,000 bootstrap samples to obtain a valid confidence interval. To improve performance we will run the for loop in parallel using the `foreach()` function in the doParallel package. First, we create a cluster, make a seed, and specify the desired number of bootstrap repitions.

```{r, eval = TRUE}
library(doParallel)

no_cores <- detectCores() - 1
registerDoParallel(cores = no_cores)
cl <- makeCluster(no_cores)

set.seed(1234)
nreps <- 1000
est <- vector(length = nreps)
```

Then we run the parallel for loop in which we apply the `adjust_uc_mc_sel2()` function to bootstrap samples of the `df_uc_mc_sel` data. In this function we specify the following arguments: the data, the exposure variable, the outcome variable, the confounder(s), the parameters in the U model, the parameters in the X model, and the parameters in the S model. Since knowledge of the complete data was known, the correct bias parameters were known in advance. The bias parameters can be provided as fixed values, as seen in this example, or values from a probability distribution. This latter strategy is referred to as probabilistic bias analysis.

```{r, eval = TRUE}
or <- foreach(i = 1:nreps, .combine = c, .packages = 'dplyr') %dopar% {
  
  bdf <- df_uc_mc_sel[sample(1:nrow(df_uc_mc_sel), nrow(df_uc_mc_sel), replace = TRUE),]
  
  est[i] <- adjust_uc_mc_sel2(bdf, exposure = "Xstar", outcome = "Y",
                              confounders = c("C1", "C2", "C3"),
                              pu1_parameters = c(-.40, .38, .46),
                              px1_parameters = c(-1.61, 2.71, .62, -.41, -.41, .40), 
                              ps1_parameters = c(-.39, .40, .75, -.04, -.04, .05)
                              )[[1]]
}
```
Finally, we obtain the OR<sub>YX</sub> estimate and 95% confidence interval from the distribution of 1,000 bootstrap estimates. As expected, OR<sub>YX</sub> ~ 2, indicating that we were able to obtain an unbiased odds ratio from the biased data.

```{r, eval = TRUE}
median(or)
quantile(or,c(.025, .975))
#> 2.031248
#> 2.007671 2.055546
```


