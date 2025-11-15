# Create a Forest Plot comparing observed and adjusted effect estimates

This function generates a forest plot comparing the observed effect
estimate with adjusted effect estimates from sensitivity analyses. The
plot includes point estimates and confidence intervals for each
analysis.

## Usage

``` r
multibias_plot(data_observed, multibias_result_list, log_scale = FALSE)
```

## Arguments

- data_observed:

  Object of class `data_observed` representing the observed causal data
  and effect of interest.

- multibias_result_list:

  A named list of sensitivity analysis results. Each element should be a
  result from
  [`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md).

- log_scale:

  Boolean indicating whether to display the x-axis on the log scale.
  Default is FALSE.

## Value

A ggplot object showing a forest plot with:

- Point estimates (blue dots)

- Confidence intervals (gray horizontal lines)

- A vertical reference line at x=1 (dashed)

- Appropriate labels and title

## Examples

``` r
if (FALSE) { # \dontrun{
df_observed <- data_observed(
  data = df_em,
  bias = "em",
  exposure = "Xstar",
  outcome = "Y",
  confounders = "C1"
)

bp1 <- bias_params(coef_list = list(x = c(-2.10, 1.62, 0.63, 0.35)))
bp2 <- bias_params(coef_list = list(x = c(-2.10 * 2, 1.62 * 2, 0.63 * 2, 0.35 * 2)))

result1 <- multibias_adjust(
  data_observed = df_observed,
  bias_params = bp1
)
result2 <- multibias_adjust(
  data_observed = df_observed,
  bias_params = bp2
)

multibias_plot(
  data_observed = df_observed,
  multibias_result_list = list(
    "Adjusted with bias params" = result1,
    "Adjusted with bias params doubled" = result2
  )
)
} # }
```
