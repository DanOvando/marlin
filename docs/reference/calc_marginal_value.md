# Calculate Marginal Revenue and Marginal Profit by Patch and Fleet

Uses forward finite differences on
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md) to
approximate \\\partial R / \partial E\_{p,f}\\ (marginal revenue) and
\\\partial \Pi / \partial E\_{p,f}\\ (marginal profit) for every
patch-fleet combination. The output matrices are designed to be merged
directly into the `buffet` list consumed by
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)
when a fleet's `spatial_allocation` is `"marginal_revenue"` or
`"marginal_profit"`, or when `fleet_model = "sole_owner"`.

## Usage

``` r
calc_marginal_value(
  e_p_fl,
  fauna,
  n_p_a,
  fleets,
  baseline = NULL,
  method = c("separable", "patch_loop"),
  epsilon = 0.001
)
```

## Arguments

- e_p_fl:

  Numeric matrix of effort by patch (rows) and fleet (columns).

- fauna:

  Named list of fauna objects (as in
  [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md)).

- n_p_a:

  Named list of numbers-at-age matrices by patch (one element per
  species), representing the current population state.

- fleets:

  Named list of fleet objects (as in
  [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md)).

- baseline:

  Optional. List returned by a prior call to
  [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md)`(output_format = "matrix")`
  with the same inputs. Supplying this avoids one redundant baseline
  evaluation and can improve performance when marginals are computed
  alongside an existing buffet.

- method:

  Character. Approximation strategy; see Details. One of `"separable"`
  (default) or `"patch_loop"`.

- epsilon:

  Numeric. Additive effort perturbation for the finite difference.
  Should be small relative to typical effort values but large enough to
  avoid floating-point noise. Default `1e-3`.

## Value

A named list with two elements, each a patches x fleets matrix with the
same dimension names as `e_p_fl`:

- `mr_p_fl`:

  Marginal revenue \\\partial R / \partial E\_{p,f}\\ by patch and
  fleet.

- `mp_p_fl`:

  Marginal profit \\\partial \Pi / \partial E\_{p,f}\\ by patch and
  fleet.

Append to the buffet from
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md) as
`buffet$mr_p_fl` and `buffet$mp_p_fl` for use by
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md).

## Details

Two approximation strategies are available:

- `"separable"`:

  Bumps effort in **all patches** of one fleet simultaneously by
  `epsilon`, then computes derivatives in a single extra `go_fish` call
  per fleet. Requires \\n\_{fleets}\\ extra evaluations. Fast, but
  assumes the marginal return in patch \\i\\ does not depend on the
  perturbation in patch \\j\\ (valid when per-patch depletion effects
  are negligible, which is typical for large grids).

- `"patch_loop"`:

  Bumps effort in **one patch at a time**, requiring \\n\_{patches}
  \times n\_{fleets}\\ extra evaluations. Fully accounts for cross-patch
  depletion effects; appropriate for small, dense grids where
  within-step movements are significant.

Called automatically by
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) when
any fleet uses marginal allocation or the sole-owner model. The results
are appended to the buffet as `mr_p_fl` and `mp_p_fl`.

## See also

[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute baseline buffet
buffet <- go_fish(e_p_fl, fauna, n_p_a, fleets, output_format = "matrix")

# Append marginal values (reuse baseline to avoid redundant evaluation)
marginals <- calc_marginal_value(
  e_p_fl   = e_p_fl,
  fauna    = fauna,
  n_p_a    = n_p_a,
  fleets   = fleets,
  baseline = buffet,
  method   = "separable"
)

buffet$mr_p_fl <- marginals$mr_p_fl
buffet$mp_p_fl <- marginals$mp_p_fl

# Now allocate effort using marginal profit
result <- allocate_effort(
  effort_by_patch = e_p_fl,
  buffet          = buffet,
  fleets          = fleets,
  open_patch      = open_patches
)
} # }
```
