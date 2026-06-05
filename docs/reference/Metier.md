# R6 class: fishing metier (fleet–species interaction)

The `Metier` R6 class stores the parameters that define how a single
fleet interacts with a single species (e.g., price, selectivity, and
catchability). A named list of metiers (one per species) is passed to
[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)
via the `metiers` argument.

## Details

This class is primarily used internally by `marlin`, but it can be
created directly with `Metier$new(...)`, then placed into a named list
(one element per species).

Common fields include `price`, `sel_form`, `sel_start`, `sel_delta`,
`catchability`, `spatial_catchability`, and derived quantities such as
selectivity-at-age and vulnerability matrices used by
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

## See also

[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

## Methods

### Public methods

- [`metier$new()`](#method-metier-initialize)

- [`metier$plot_selectivity()`](#method-metier-plot_selectivity)

- [`metier$plot_catchability()`](#method-metier-plot_catchability)

- [`metier$clone()`](#method-metier-clone)

------------------------------------------------------------------------

### `metier$new()`

Parameters for specific metier of a given fleet

#### Usage

    metier$new(
      critter = NA,
      price = 10,
      sel_form = "logistic",
      sel_unit = "p_of_mat",
      sel_start = 1,
      sel_delta = 0.1,
      sel05_anchor = NULL,
      sel_at_linf = NULL,
      catchability = 0.2,
      spatial_catchability = NULL,
      p_explt = 1,
      sel_at_age = NULL
    )

#### Arguments

- `critter`:

  the name of the critter in the fauna object this metier applies to

- `price`:

  the price per unit weight of the critter in question

- `sel_form`:

  the selectivity form, one of "logistic", "dome","uniform", or "manual"

- `sel_unit`:

  the unit of selectivity, one of "p_of_mat" which means selectivity is
  in proportion of age at maturity, or "length" where selectivity is in
  units of length

- `sel_start`:

  the value of sel_unit at which selectivity "starts"

- `sel_delta`:

  the delta parameter in the selectivity function

- `sel05_anchor`:

  lower anchor at which (contact) selectivity is 0.05, either length or
  p_of_mat

- `sel_at_linf`:

  contact selectivity at linf

- `catchability`:

  the catchability per uni effort pararmeter, generally overwritten by
  tune_fleet

- `spatial_catchability`:

  a matrix of spatial q

- `p_explt`:

  the proportion of total exploitation for a given critter coming from
  this metier. This value is relative to all other p_explt values for
  the critter in question. Set to 0 to have metier not catch critter at
  all plot selectivity \#'

- `sel_at_age`:

  a manual vector of gear (contact) selectivity at age, where values are
  between 0 and 1

------------------------------------------------------------------------

### `metier$plot_selectivity()`

#### Usage

    metier$plot_selectivity()

#### Returns

a plot of selectivity at age for the metier plot selectivity \#'

------------------------------------------------------------------------

### `metier$plot_catchability()`

#### Usage

    metier$plot_catchability()

#### Returns

a plot of selectivity at age for the metier

------------------------------------------------------------------------

### `metier$clone()`

The objects of this class are cloneable with this method.

#### Usage

    metier$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
