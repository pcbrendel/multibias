# multibias
R package for multi-bias analysis

## Overview

multibias is a set of functions that provide odds ratio estimates that are adjusted for any combination of uncontrolled confounding, selection bias, and exposure misclassification:

  - `adjust_uc_sel()` adjusts for uncontrolled confounding and selection bias.
  - `adjust_mc_sel()` adjusts for exposure misclassification and selection bias.
  - `adjust_uc_mc()` adjusts for uncontrolled confounding and exposure misclassificaiton.
  - `adjust_uc_mc_sel()` adjusts for all three biases.
 
 
 
 If you are new to bias analysis, I'd recommend checking out the [Applying Quantitative Bias Analysis to Epidemiologic Data](https://www.springer.com/us/book/9780387879604) textbook or visiting my [website](https://pcbrendel.github.io/).

## Installation

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("pcbrendel/multibias")
```

## Usage

Say you have a data set whose causal relationships can be represented in the following DAG:

![uc_mc_sel_DAG](DAGs/uc_mc_sel_DAG.png)

The following variables are defined:

It can be seen from this DAG that the data suffers from three sources of bias. There is uncontrolled confounding from (unobserved) variable U. The true exposure, X, is unobserved, and the misclassified exposure X* is dependent on both the exposure and outcome. Lastly, there is collider stratification at variable S, representing selection into the study. The study naturally only examines those who were selected (i.e. those with S=1).

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
