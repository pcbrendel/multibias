# Print method for data_validation objects

Prints a formatted summary of a `data_validation` object, including:

- The true exposure and outcome variables

- Any confounders, misclassified variables, or selection indicators

- A preview of the first 5 rows of data

## Usage

``` r
# S3 method for class 'data_validation'
print(x, ...)
```

## Arguments

- x:

  A `data_validation` object

- ...:

  Additional arguments passed to print

## Value

The input object invisibly
