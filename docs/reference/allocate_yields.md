# Allocate Catches, Revenue, and Profit by Fleet

Distributes the population-level catch matrix `c_p_a` (patch x age)
produced by a critter's `swim()` method across fleets, then calculates
revenue and profit for each fleet-patch combination.

## Usage

``` r
allocate_yields(
  f_p_a_fl,
  p_p_a_fl,
  e_p_fl,
  critter,
  pop,
  patches,
  ages,
  fleets,
  fauna
)
```

## Arguments

- f_p_a_fl:

  3-D numeric array `[patches, ages, fleets]` of proportional fishing
  mortality by fleet (normalised to sum to 1 across fleets for each
  patch-age cell).

- p_p_a_fl:

  3-D numeric array `[patches, ages, fleets]` of price per unit catch,
  by patch, age, and fleet.

- e_p_fl:

  Numeric matrix of effort by patch (rows) and fleet (columns).

- critter:

  Character. Name of the species being processed (used to index
  `fauna`).

- pop:

  List returned by `critter$swim()`, including `c_p_a` (catch by patch
  and age).

- patches:

  Integer. Total number of patches.

- ages:

  Integer. Number of age classes for this species.

- fleets:

  Named list of fleet objects from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md).

- fauna:

  Named list of fauna objects from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

## Value

A named list:

- `c_p_a_fl`:

  3-D array `[patches, ages, fleets]` of catch in numbers.

- `r_p_a_fl`:

  3-D array `[patches, ages, fleets]` of revenue.

- `r_p_fl`:

  Matrix `[patches, fleets]` of revenue summed across ages.

- `c_p_fl`:

  Matrix `[patches, fleets]` of catch summed across ages.

- `prof_p_fl`:

  Matrix `[patches, fleets]` of profit (revenue minus cost).

## Details

The fishing-mortality share of each fleet in each patch-age cell is
given by `f_p_a_fl / f_p_a` (already normalised). Total catch in each
cell is then multiplied by this share to give fleet-specific catch.
Revenue is calculated as catch times price (`p_p_a_fl`). Profit is
revenue minus the patch-level cost, where costs follow the formula:
\$\$C\_{l,p} = c\_{0,l} \\ E^{ref}\_l
\left\[\left(\frac{E\_{l,p}/n\_{sp}}{E^{ref}\_l}\right)^{\gamma_l} +
\theta_l \\ \tilde{d}\_{l,p} \\
\frac{E\_{l,p}/n\_{sp}}{E^{ref}\_l}\right\]\$\$ Note that effort is
divided by the number of species (\\n\_{sp}\\) to avoid double-counting
costs across the species loop in
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

This is an internal helper called once per critter per time step inside
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) and
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md).

## See also

[`aggregate_yields`](https://danovando.github.io/marlin/reference/aggregate_yields.md),
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
