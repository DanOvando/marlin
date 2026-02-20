# Logistic link function mapping positive reals to (0, 1)

Applies the transformation `q = raw_q / (1 + raw_q)`, equivalent to
`plogis(log(raw_q))`. This smoothly maps positive real values to (0, 1),
approximately preserving values much less than 1 (when `raw_q << 1`,
`q ≈ raw_q`) while preventing q from ever reaching or exceeding 1. Used
during numerical optimisation of catchability to maintain smooth
gradients near the constraint boundary.

## Usage

``` r
q_link(raw_q)
```

## Arguments

- raw_q:

  Positive numeric; unconstrained catchability value.

## Value

Numeric in (0, 1).
