# Set catchability for a single fleet-species metier

Assigns scalar catchability to a metier, rescales its spatial
catchability vector to have the correct mean, and rebuilds `vul_p_a`.
Catchability is enforced to lie in (0, 1).

## Usage

``` r
set_metier_catchability(metier, q_value, use_link = FALSE)
```

## Arguments

- metier:

  A [`Metier`](https://danovando.github.io/marlin/reference/Metier.md)
  R6 object (modified in place via R6 reference semantics).

- q_value:

  Numeric. Target scalar catchability. Must be in `[0, 1]` when
  `use_link = FALSE`. When `use_link = TRUE`, any non-negative value is
  accepted and mapped to (0, 1) via
  [`q_link`](https://danovando.github.io/marlin/reference/q_link.md).

- use_link:

  Logical. If `TRUE`, applies
  [`q_link`](https://danovando.github.io/marlin/reference/q_link.md) to
  map `q_value` from the positive reals into (0, 1). Use during
  numerical optimisation to ensure smooth gradients. Default `FALSE`.

## Value

The modified metier, invisibly.
