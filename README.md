# multibias
R package for multi-bias analysis

## Overview

multibias is a set of functions that provide effect estimates that are adjusted for any combination of uncontrolled confounding, selection bias, and exposure misclassification:

  - `adjust_uc_sel()` adjusts for uncontrolled confounding and selection bias.
  - `adjust_mc_sel()` adjusts for exposure misclassification and selection bias.
  - `adjust_uc_mc()` adjusts for uncontrolled confounding and exposure misclassificaiton.
  - `adjust_uc_mc_sel()` adjusts for all three biases.
 
 
 
 If you are new to bias analysis, I'd recommend checking out the [Applying Quantitative Bias Analysis to Epidemiologic Data](https://www.springer.com/us/book/9780387879604) textbook or visit my [website](https://pcbrendel.github.io/).

## Installation

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("pcbrendel/multibias")
```

## Usage

Say you have a data set whose causal relationships can be represented in the following DAG:

![uc_mc_sel_DAG](DAGs/uc_mc_sel_DAG.png)

A simulated data set corresponding to this DAG, `df_uc_mc_sel` can be loaded from the multibias package. 
