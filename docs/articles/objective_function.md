# The Objective Function of the Fleet

## Objective function

We define the objective for a single focal fleet $`f`$.

### Indices

- $`f`$: fleet  
- $`l`$: patch (location)  
- $`s`$: species  
- $`a`$: age

### Effort and fishing mortality

Let $`E_{l,f}`$ denote the effort of fleet $`f`$ in patch $`l`$.
Fleet-specific fishing mortality at age is

``` math
F_{l,s,a,f}
=
q_{l,s,f}\,\text{sel}_{a,s,f}\,E_{l,f},
```

where $`q_{l,s,f}`$ is catchability and $`\text{sel}_{a,s,f}`$ is
selectivity at age.

Let $`F^{\text{ext}}_{l,s,a}`$ denote exogenous fishing mortality from
all fleets other than $`f`$, and let $`m_{a,s}`$ denote natural
mortality. Total instantaneous mortality is

``` math
Z_{l,s,a}
=
F_{l,s,a,f}
+
F^{\text{ext}}_{l,s,a}
+
m_{a,s}.
```

------------------------------------------------------------------------

### Revenue (Baranov catch)

Revenue for fleet $`f`$ is given by Baranov catch multiplied by price:

``` math
\text{Rev}_f
=
\sum_{l}\sum_{s}\sum_{a}
\text{price}_{s,f}
\;
\frac{F_{l,s,a,f}}{Z_{l,s,a}}
\left(1 - e^{-\Delta_t Z_{l,s,a}}\right)
b_{l,s,a},
```

where $`b_{l,s,a}`$ is biomass and $`\Delta_t`$ is the time step.

------------------------------------------------------------------------

### Cost (Normalized Formulation)

#### Parameters

- $`c0_f`$: cost-per-unit-effort scalar (calibrated by `tune_fleets`)
- $`E^{\text{ref}}_f`$: reference effort level (mean per-patch effort at
  equilibrium)
- $`\gamma_f \ge 1`$: effort cost exponent (congestion parameter;
  default 1.2)
- $`\theta_f \ge 0`$: travel weight, controls relative importance of
  spatial costs
- $`\tilde{d}_{l,f}`$: normalized cost per patch (spatial shape with
  mean 1 across open patches)

#### User-facing parameters

Users specify the following interpretable parameters in `create_fleet`:

| Parameter | Meaning | Default |
|----|----|----|
| `cr_ratio` | Cost/revenue ratio at equilibrium | 1 |
| `effort_cost_exponent` ($`\gamma`$) | Convexity of congestion costs | 1.2 |
| `travel_fraction` | Share of costs from travel at equilibrium | 0 |
| `cost_per_patch` | Spatial cost pattern (normalized internally) | uniform |

The travel weight $`\theta_f`$ is derived from `travel_fraction` as:

``` math
\theta_f = \frac{\text{travel\_fraction}}{1 - \text{travel\_fraction}}
```

This parameterization ensures that travel costs constitute exactly
`travel_fraction` of total costs at equilibrium when effort is uniformly
distributed across open patches.

#### Cost formula

Total cost for fleet $`f`$ is:

``` math
\text{Cost}_f
=
c0_f \cdot E^{\text{ref}}_f
\left(
\sum_{l} \left(\frac{E_{l,f}}{E^{\text{ref}}_f}\right)^{\gamma_f}
+
\theta_f \sum_{l} \tilde{d}_{l,f} \cdot \frac{E_{l,f}}{E^{\text{ref}}_f}
\right)
```

**Key properties:**

1.  **Effort normalization:** Dividing by $`E^{\text{ref}}_f`$ makes
    costs scale-invariant. Changing grid resolution or `base_effort`
    doesn’t affect the relative importance of congestion vs. travel.

2.  **Separable components:**

    - Convex term $`\sum_l (E_l / E^{\text{ref}})^{\gamma}`$:
      congestion/diminishing returns
    - Travel term $`\theta \sum_l \tilde{d}_l (E_l / E^{\text{ref}})`$:
      spatial heterogeneity

3.  **Interpretable parameters:**

    - $`\gamma`$ controls congestion severity (1 = linear, \>1 = convex)
    - $`\theta`$ (via `travel_fraction`) controls how much spatial
      location matters
    - $`\tilde{d}_l`$ is a pure spatial shape (mean 1), making
      `travel_fraction` directly interpretable

#### Tuning

The `tune_fleets` function calibrates $`c0_f`$ to achieve the target
cost/revenue ratio:

``` math
c0_f = \frac{\text{cr\_ratio} \cdot \text{Rev}_f}{E^{\text{ref}}_f \cdot P_f \cdot (1 + \theta_f)}
```

where $`P_f`$ is the number of open fishing patches. The function also
computes and stores the equilibrium $`E^{\text{ref}}_f`$ as the mean
per-patch effort.

------------------------------------------------------------------------

### Final objective

The objective for fleet $`f`$ is

``` math
\begin{aligned}
\max_{\{E_{l,f}\}}
\;\; \text{obj}_f
&=
\sum_{l}\sum_{s}\sum_{a}
\text{price}_{s,f}
\;
\frac{q_{l,s,f}\,\text{sel}_{a,s,f}\,E_{l,f}}
{q_{l,s,f}\,\text{sel}_{a,s,f}\,E_{l,f} + F^{\text{ext}}_{l,s,a} + m_{a,s}}
\\[6pt]
&\quad \times
\left(1 - e^{-\Delta_t\left(
q_{l,s,f}\,\text{sel}_{a,s,f}\,E_{l,f}
+ F^{\text{ext}}_{l,s,a}
+ m_{a,s}
\right)}\right)
b_{l,s,a}
\\[10pt]
&\quad
-
c0_f \cdot E^{\text{ref}}_f
\left(
\sum_{l} \left(\frac{E_{l,f}}{E^{\text{ref}}_f}\right)^{\gamma_f}
+
\theta_f \sum_{l} \tilde{d}_{l,f} \cdot \frac{E_{l,f}}{E^{\text{ref}}_f}
\right).
\end{aligned}
```

------------------------------------------------------------------------

### Constraints

The optimization is subject to:

- Total effort constraint:
  ``` math
  \sum_l E_{l,f} = E^{\text{tot}}_f
  ```
- Non-negativity:
  ``` math
  E_{l,f} \ge 0
  ```
- Fleet-specific spatial closures, enforced by setting $`E_{l,f} = 0`$
  in closed patches.
