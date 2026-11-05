# Calculate Fishable Revenue by Patch and Fleet

Computes the revenue standing available to each fleet in each patch: the
intersection of biomass at age, the fleet's vulnerability at age, and
price. Used by
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)
when a fleet's `spatial_allocation` is `"fishable_revenue"`.

## Usage

``` r
calc_fishable_revenue(b_p_a, fauna, fleets)
```

## Arguments

- b_p_a:

  Named list of biomass-at-age matrices (patches x ages), one per
  species, named to match `fauna`. In
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
  this is the `b_p_a` stored for the previous step, which is the biomass
  the current step is about to fish.

- fauna:

  Named list of fauna objects (e.g. from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)).

- fleets:

  Named list of fleet objects (e.g. from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)),
  each with a `metiers` list.

## Value

Numeric matrix, patches x fleets, with fleet names as column names.

## Details

For fleet \\f\\ and patch \\p\\, summing over species \\s\\ and ages
\\a\\: \$\$FR\_{p,f} = \sum_s price\_{s,f} \sum_a B\_{p,a,s} \\
v\_{p,a,s,f}\$\$ where \\B\\ is biomass at age and \\v\\ is the metier's
`vul_p_a` (spatial catchability times selectivity at age).

### Why this exists

The realised spatial-allocation objectives (`rpue`, `revenue`, `ppue`,
`profit`, `cpue`, `catch`) are functions of realised fishing: they
record what a fleet *got* given the effort it actually deployed, so the
signal reaching
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)
at step \\t\\ necessarily describes step \\t-1\\. When habitat changes
on a timescale close to the time step – seasonal migration, say – a
fleet on one of them is always responding to where fish were rather than
where they are, and no value of `responsiveness` fixes it, because the
lag is structural rather than a matter of step size.

The marginal objectives (`marginal_revenue`, `marginal_profit`) avoid
that lag, because
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
evaluates them on the biomass the coming step will fish, but they cost
extra
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md)
evaluations every step and model a fleet that maximises its total
return. Fishable revenue is the cheap, economics-free alternative. It
takes no effort argument, so it can be evaluated against the biomass a
step is *about* to fish with no circularity, no fixed point to solve and
no Baranov solve.

It is additionally defined where effort is zero, unlike the
per-unit-effort objectives (`rpue`, `cpue`, `ppue`), which are `NA`
there. A fleet allocating on fishable revenue therefore sees unfished
patches honestly, rather than having to fish somewhere to learn anything
about it.

### What it gives up

- **No crowding.** The objective is standing stock, not return per unit
  effort, so within a step there is nothing to stop effort piling into
  the best patch. Depletion still pushes back, but only from the next
  step onward – delayed negative feedback, which can oscillate. Consider
  `memory_halflife > 0` to damp it.

- **No costs.** Cost depends on effort, so profit-side variants (`ppue`,
  `profit`) have no fishable-revenue analogue.

- **No fleet interaction within a step.** Two fleets on the same stock
  do not see each other's effort; they feel each other only through next
  step's depleted biomass.

- **No equalised returns.** The per-unit-effort rules (`rpue`, `ppue`,
  `cpue`) relax toward equal average returns across fished patches (the
  ideal free distribution), and the marginal rules toward equal marginal
  returns within each step (a myopic sole owner's allocation). This one
  relaxes toward effort tracking standing stock, which is a different
  equilibrium – "fish where the fish are". None of these is wrong; they
  are different fleets.

## See also

[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md),
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
[`calc_marginal_value`](https://danovando.github.io/marlin/reference/calc_marginal_value.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
