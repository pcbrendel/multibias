# Represent bias parameters

`bias_params` is one of two different options to represent bias
assumptions for bias adjustment. The
[`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
function will apply the assumptions from these models and use them to
adjust for biases in the observed data. It takes one input, a list,
where each item in the list corresponds to the necessary models for bias
adjustment. See below for bias models.

For each of the following bias models, the variables are defined:

- X = True exposure

- X\* = Misclassified exposure

- Y = True outcome

- Y\* = Misclassified outcome

- C = Known confounder

- j = Number of known confounders

- U = Uncontrolled confounder

- S = Selection indicator

&nbsp;

- Uncontrolled confounding:

  logit(P(U=1)) = α₀ + α₁X + α₂Y + α_(2+j)C_(j)

- Exposure misclassification:

  logit(P(X=1)) = δ₀ + δ₁X\* + δ₂Y + δ_(2+j)C_(j)

- Outcome misclassification:

  logit(P(Y=1)) = δ₀ + δ₁X + δ₂Y\* + δ_(2+j)C_(j)

- Selection bias:

  logit(P(S=1)) = β₀ + β₁X + β₂Y

- Uncontrolled Confounding & Exposure Misclassification (Option 1):

  logit(P(U=1)) = α₀ + α₁X + α₂Y  
  logit(P(X=1)) = δ₀ + δ₁X\* + δ₂Y + δ_(2+j)C_(j)

- Uncontrolled Confounding & Exposure Misclassification (Option 2):

  log(P(X=1,U=0)/P(X=0,U=0)) = γ_(1,0) + γ_(1,1)X\* + γ_(1,2)Y +
  γ_(1,2+j)C_(j)  
  log(P(X=0,U=1)/P(X=0,U=0)) = γ_(2,0) + γ_(2,1)X\* + γ_(2,2)Y +
  γ_(2,2+j)C_(j)  
  log(P(X=1,U=1)/P(X=0,U=0)) = γ_(3,0) + γ_(3,1)X\* + γ_(3,2)Y +
  γ_(3,2+j)C_(j)

- Uncontrolled Confounding & Outcome Misclassification (Option 1):

  logit(P(U=1)) = α₀ + α₁X + α₂Y  
  logit(P(Y=1)) = δ₀ + δ₁X + δ₂Y\* + δ_(2+j)C_(j)

- Uncontrolled Confounding & Outcome Misclassification (Option 2):

  log(P(U=1,Y=0)/P(U=0,Y=0)) = γ_(1,0) + γ_(1,1)X + γ_(1,2)Y\* +
  γ_(1,2+j)C_(j)  
  log(P(U=0,Y=1)/P(U=0,Y=0)) = γ_(2,0) + γ_(2,1)X + γ_(2,2)Y\* +
  γ_(2,2+j)C_(j)  
  log(P(U=1,Y=1)/P(U=0,Y=0)) = γ_(3,0) + γ_(3,1)X + γ_(3,2)Y\* +
  γ_(3,2+j)C_(j)

- Uncontrolled Confounding & Selection Bias:

  logit(P(U=1)) = α₀ + α₁X + α₂Y + α_(2+j)C_(j)  
  logit(P(S=1)) = β₀ + β₁X + β₂Y

- Exposure Misclassification & Outcome Misclassification (Option 1):

  logit(P(X=1)) = δ₀ + δ₁X\* + δ₂Y\* + δ_(2+j)C_(j)  
  logit(P(Y=1)) = β₀ + β₁X + β₂Y\* + β_(2+j)C_(j)

- Exposure Misclassification & Outcome Misclassification (Option 2):

  log(P(X=1,Y=0) / P(X=0,Y=0)) = γ_(1,0) + γ_(1,1)X\* + γ_(1,2)Y\* +
  γ_(1,2+j)C_(j)  
  log(P(X=0,Y=1) / P(X=0,Y=0)) = γ_(2,0) + γ_(2,1)X\* + γ_(2,2)Y\* +
  γ_(2,2+j)C_(j)  
  log(P(X=1,Y=1) / P(X=0,Y=0)) = γ_(3,0) + γ_(3,1)X\* + γ_(3,2)Y\* +
  γ_(3,2+j)C_(j)

- Exposure Misclassification & Selection Bias:

  logit(P(X=1)) = δ₀ + δ₁X\* + δ₂Y + δ_(2+j)C_(j)  
  logit(P(S=1)) = β₀ + β₁X\* + β₂Y + β_(2+j)C_(j)

- Outcome Misclassification & Selection Bias:

  logit(P(Y=1)) = δ₀ + δ₁X + δ₂Y\* + δ_(2+j)C_(j)  
  logit(P(S=1)) = β₀ + β₁X + β₂Y\* + β_(2+j)C_(j)

- Uncontrolled Confounding, Exposure Misclassification, and Selection
  Bias (Option 1):

  logit(P(U=1)) = α₀ + α₁X + α₂Y  
  logit(P(X=1)) = δ₀ + δ₁X\* + δ₂Y + δ_(2+j)C_(j)  
  logit(P(S=1)) = β₀ + β₁X\* + β₂Y + β_(2+j)C_(j)

- Uncontrolled Confounding, Exposure Misclassification, and Selection
  Bias (Option 2):

  log(P(X=1,U=0)/P(X=0,U=0)) = γ_(1,0) + γ_(1,1)X\* + γ_(1,2)Y +
  γ_(1,2+j)C_(j)  
  log(P(X=0,U=1)/P(X=0,U=0)) = γ_(2,0) + γ_(2,1)X\* + γ_(2,2)Y +
  γ_(2,2+j)C_(j)  
  log(P(X=1,U=1)/P(X=0,U=0)) = γ_(3,0) + γ_(3,1)X\* + γ_(3,2)Y +
  γ_(3,2+j)C_(j)  
  logit(P(S=1)) = β₀ + β₁X\* + β₂Y + β_(2+j)C_(j)

- Uncontrolled Confounding, Outcome Misclassification, and Selection
  Bias (Option 1):

  logit(P(U=1)) = α₀ + α₁X + α₂Y  
  logit(P(Y=1)) = δ₀ + δ₁X + δ₂Y\* + δ_(2+j)C_(j)  
  logit(P(S=1)) = β₀ + β₁X + β₂Y\* + β_(2+j)C_(j)

- Uncontrolled Confounding, Outcome Misclassification, and Selection
  Bias (Option 2):

  log(P(U=1,Y=0)/P(U=0,Y=0)) = γ_(1,0) + γ_(1,1)X + γ_(1,2)Y\* +
  γ_(1,2+j)C_(j)  
  log(P(U=0,Y=1)/P(U=0,Y=0)) = γ_(2,0) + γ_(2,1)X + γ_(2,2)Y\* +
  γ_(2,2+j)C_(j)  
  log(P(U=1,Y=1)/P(U=0,Y=0)) = γ_(3,0) + γ_(3,1)X + γ_(3,2)Y\* +
  γ_(3,2+j)C_(j)  
  logit(P(S=1)) = β₀ + β₁X + β₂Y\* + β_(2+j)C_(j)

## Usage

``` r
bias_params(coef_list)
```

## Arguments

- coef_list:

  List of coefficient values from the above options of models. Each item
  of the list is an equation. The left side of the equation identifies
  the model (i.e., "u" for the model predicting the uncontrolled
  confounder). For the multinomial models, specify the value here based
  on the numerator (i.e., "x1u0", "x0u1", "x1u1" for the three
  multinomial models in Uncontrolled Confounding & Exposure
  Misclassification, Option 2) The right side of the equation is the
  vector of values corresponding to the model coefficients (from left to
  right).

## Examples

``` r
list_for_uc <- list(
  u = c(-0.19, 0.61, 0.70, -0.09, 0.10, -0.15)
)

bp_uc <- bias_params(coef_list = list_for_uc)

list_for_em_om <- list(
  x1y0 = c(-2.18, 1.63, 0.23, 0.36),
  x0y1 = c(-3.17, 0.22, 1.60, 0.40),
  x1y1 = c(-4.76, 1.82, 1.83, 0.72)
)

bp_em_om <- bias_params(coef_list = list_for_em_om)
```
