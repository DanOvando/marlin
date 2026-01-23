#include <Rcpp.h>
#include <Rmath.h>      // expm1, R_FINITE
#include <algorithm>    // std::fill, std::max
#include <cmath>        // std::exp, std::pow, std::fabs, std::sqrt
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

inline double harvest_term(const double x) {
  return -expm1(-x); // stable 1-exp(-x)
}

inline double cost_derivative_patch(const double E,
                                    const double cost_patch,
                                    const double c0,
                                    const double gamma,
                                    const bool include_costs) {
  if (!include_costs) return 0.0;
  if (gamma == 1.0) return c0 * (cost_patch + 1.0);
  const double Egm1 = (E > 0.0) ? std::pow(E, gamma - 1.0) : 0.0;
  return c0 * (cost_patch + gamma * Egm1);
}

inline double cost_total_patch(const double E,
                               const double cost_patch,
                               const double c0,
                               const double gamma) {
  const double Eg = (E > 0.0) ? std::pow(E, gamma) : 0.0;
  return c0 * (cost_patch * E + Eg);
}

// [[Rcpp::export]]
List cpp_allocate_ifd_kkt_fullsolve_fast(
    const double Etot_target,
    const List& alpha_mats,
    const List& other_mort_mats,
    const List& biomass_mats,
    const NumericVector& price_s,
    const NumericVector& cost_patch,
    const double c0,
    const double gamma,
    const IntegerVector& fishable_int,
    const double time_step,
    const bool include_costs,
    const int n_outer,
    const int n_inner,
    const double active_tol,
    const double flat_tol_sd,
    const double flat_tol_abs
) {
  const int P = fishable_int.size();
  const int n_species = alpha_mats.size();
  const double dt = time_step;

  NumericVector E_star(P);
  std::fill(E_star.begin(), E_star.end(), 0.0);

  // ---- ONE-TIME unwrap matrices + raw pointers ----
  std::vector< NumericMatrix > alphaV;
  std::vector< NumericMatrix > otherV;
  std::vector< NumericMatrix > biomV;
  std::vector< const double* > alphaPtr;
  std::vector< const double* > otherPtr;
  std::vector< const double* > biomPtr;
  std::vector<int> Acols;
  std::vector<double> priceV;

  alphaV.reserve(n_species);
  otherV.reserve(n_species);
  biomV.reserve(n_species);
  alphaPtr.reserve(n_species);
  otherPtr.reserve(n_species);
  biomPtr.reserve(n_species);
  Acols.reserve(n_species);
  priceV.reserve(n_species);

  for (int s = 0; s < n_species; ++s) {
    NumericMatrix a = as<NumericMatrix>(alpha_mats[s]);
    NumericMatrix o = as<NumericMatrix>(other_mort_mats[s]);
    NumericMatrix b = as<NumericMatrix>(biomass_mats[s]);

    // assume nrow == P (validated on R side)
    alphaV.push_back(a);
    otherV.push_back(o);
    biomV.push_back(b);

    alphaPtr.push_back(REAL(alphaV.back()));
    otherPtr.push_back(REAL(otherV.back()));
    biomPtr.push_back(REAL(biomV.back()));

    Acols.push_back(a.ncol());
    priceV.push_back(price_s[s]);
  }

  // ---- inline mp_rev(p,E) using raw pointers ----
  auto mp_rev_patch = [&](const int p, const double E) {
    double out = 0.0;
    for (int s = 0; s < n_species; ++s) {
      const int A = Acols[s];
      const double price = priceV[s];
      const double* aPtr = alphaPtr[s];
      const double* oPtr = otherPtr[s];
      const double* bPtr = biomPtr[s];

      // column-major index: p + P*a
      for (int a = 0; a < A; ++a) {
        const int idx = p + P * a;
        const double alpha_pa  = aPtr[idx];
        const double fish_mort = alpha_pa * E;
        const double other_m   = oPtr[idx];
        const double bio       = bPtr[idx];

        double z = other_m + fish_mort;
        if (!R_FINITE(z)) z = 0.0;

        const double z_safe = (z > 1e-15) ? z : 1e-15;
        const double x = dt * z_safe;

        const double harvest = harvest_term(x);
        const double surv    = std::exp(-x);

        double dcatch_dmort =
          ((other_m / (z_safe * z_safe)) * harvest +
          (fish_mort / z_safe) * (dt * surv)) * bio;

        if (!R_FINITE(dcatch_dmort)) dcatch_dmort = 0.0;
        out += alpha_pa * dcatch_dmort * price;
      }
    }
    if (!R_FINITE(out)) out = 0.0;
    return out;
  };

  auto revenue_patch = [&](const int p, const double E) {
    double rev = 0.0;
    for (int s = 0; s < n_species; ++s) {
      const int A = Acols[s];
      const double price = priceV[s];
      const double* aPtr = alphaPtr[s];
      const double* oPtr = otherPtr[s];
      const double* bPtr = biomPtr[s];

      for (int a = 0; a < A; ++a) {
        const int idx = p + P * a;
        const double alpha_pa  = aPtr[idx];
        const double fish_mort = alpha_pa * E;
        const double other_m   = oPtr[idx];
        const double bio       = bPtr[idx];

        double z = other_m + fish_mort;
        if (!R_FINITE(z)) z = 0.0;

        const double z_safe = (z > 1e-15) ? z : 1e-15;
        const double x = dt * z_safe;

        const double harvest = harvest_term(x);

        double catch_pa = (fish_mort / z_safe) * harvest * bio;
        if (!R_FINITE(catch_pa)) catch_pa = 0.0;

        rev += catch_pa * price;
      }
    }
    if (!R_FINITE(rev)) rev = 0.0;
    return rev;
  };

  // ---- precompute mp0 once ----
  std::vector<double> mp0(P, -1e300);
  std::vector<double> mp0_open;
  mp0_open.reserve(P);

  double mp0_max = -1e300;
  double mp0_maxabs = 0.0;

  for (int p = 0; p < P; ++p) {
    if (fishable_int[p] != 1) continue;

    const double mp_rev0  = mp_rev_patch(p, 0.0);
    const double mp_cost0 = cost_derivative_patch(0.0, cost_patch[p], c0, gamma, include_costs);
    const double mp       = mp_rev0 - mp_cost0;

    mp0[p] = mp;
    mp0_open.push_back(mp);

    if (mp > mp0_max) mp0_max = mp;
    const double ab = std::fabs(mp);
    if (ab > mp0_maxabs) mp0_maxabs = ab;
  }

  if (mp0_open.size() == 0) {
    NumericVector revenue_p(P); std::fill(revenue_p.begin(), revenue_p.end(), 0.0);
    NumericVector cost_p(P);    std::fill(cost_p.begin(), cost_p.end(), 0.0);
    NumericVector profit_p(P);  std::fill(profit_p.begin(), profit_p.end(), 0.0);
    NumericVector marginal_objective_p(P); std::fill(marginal_objective_p.begin(), marginal_objective_p.end(), NA_REAL);

    return List::create(
      _["E_target"] = E_star,
      _["lambda"] = NA_REAL,
      _["revenue_p"] = revenue_p,
      _["cost_p_total"] = cost_p,
      _["profit_p"] = profit_p,
      _["marginal_objective_p"] = marginal_objective_p,
      _["diagnostic"] = List::create(_["note"] = "No fishable patches")
    );
  }

  // sd(mp0_open)
  double mean = 0.0;
  for (size_t i = 0; i < mp0_open.size(); ++i) mean += mp0_open[i];
  mean /= (double) mp0_open.size();

  double var = 0.0;
  if (mp0_open.size() >= 2) {
    for (size_t i = 0; i < mp0_open.size(); ++i) {
      const double d = mp0_open[i] - mean;
      var += d * d;
    }
    var /= (double)(mp0_open.size() - 1);
  }
  const double mp_sd = std::sqrt(var);

  // ---- flat fallback (tiny-biomass fix) ----
  if (!include_costs && (mp_sd <= flat_tol_sd || mp0_maxabs <= flat_tol_abs)) {
    int n_open = 0;
    for (int p = 0; p < P; ++p) if (fishable_int[p] == 1) n_open++;
    const double share = (n_open > 0) ? (Etot_target / (double)n_open) : 0.0;

    for (int p = 0; p < P; ++p) E_star[p] = (fishable_int[p] == 1) ? share : 0.0;

    NumericVector revenue_p(P), cost_p(P), profit_p(P), marginal_objective_p(P);
    double tot_rev = 0.0, tot_cost = 0.0, tot_profit = 0.0;

    for (int p = 0; p < P; ++p) {
      if (fishable_int[p] != 1) {
        revenue_p[p] = 0.0;
        cost_p[p] = 0.0;
        profit_p[p] = 0.0;
        marginal_objective_p[p] = NA_REAL;
        continue;
      }

      const double E = E_star[p];
      const double rev = revenue_patch(p, E);
      const double cst = cost_total_patch(E, cost_patch[p], c0, gamma);

      revenue_p[p] = rev;
      cost_p[p] = cst;
      profit_p[p] = rev - cst;

      // marginal objective at solution effort
      const double mp_rev  = mp_rev_patch(p, E);
      const double mp_cost = cost_derivative_patch(E, cost_patch[p], c0, gamma, include_costs);
      marginal_objective_p[p] = mp_rev - mp_cost;

      tot_rev += rev;
      tot_cost += cst;
      tot_profit += (rev - cst);
    }

    return List::create(
      _["E_target"] = E_star,
      _["lambda"] = 0.0,
      _["revenue_p"] = revenue_p,
      _["cost_p_total"] = cost_p,
      _["profit_p"] = profit_p,
      _["marginal_objective_p"] = marginal_objective_p,
      _["total_revenue"] = tot_rev,
      _["total_cost"] = tot_cost,
      _["total_profit"] = tot_profit,
      _["diagnostic"] = List::create(
        _["include_costs"] = include_costs,
        _["flat_fallback_used"] = true,
        _["mp0_open_sd"] = mp_sd,
        _["mp0_open_maxabs"] = mp0_maxabs,
        _["flat_tol_sd"] = flat_tol_sd,
        _["flat_tol_abs"] = flat_tol_abs
      )
    );
  }

  // ---- bracket lambda ----
  const double bump = std::max(1e-12, 1e-6 * std::fabs(mp0_max));
  double lambda_hi = mp0_max + bump;

  // lower via mp at E=Etot
  double lambda_lo = 0.0;
  bool lo_set = false;
  for (int p = 0; p < P; ++p) {
    if (fishable_int[p] != 1) continue;
    const double mp_rev  = mp_rev_patch(p, Etot_target);
    const double mp_cost = cost_derivative_patch(Etot_target, cost_patch[p], c0, gamma, include_costs);
    const double mp      = mp_rev - mp_cost;
    if (!lo_set) { lambda_lo = mp; lo_set = true; }
    else if (mp < lambda_lo) lambda_lo = mp;
  }
  lambda_lo = lambda_lo - std::fabs(lambda_lo) - 1.0;

  // compute E(lambda) and sum
  auto compute_E_sum = [&](const double lambda, NumericVector& E_out) {
    double sumE = 0.0;
    for (int p = 0; p < P; ++p) {
      if (fishable_int[p] != 1) { E_out[p] = 0.0; continue; }
      if (!(mp0[p] > (lambda + active_tol))) { E_out[p] = 0.0; continue; }

      double lo = 0.0;
      double hi = Etot_target;

      for (int it = 0; it < n_inner; ++it) {
        const double mid = 0.5 * (lo + hi);

        const double mp_rev  = mp_rev_patch(p, mid);
        const double mp_cost = cost_derivative_patch(mid, cost_patch[p], c0, gamma, include_costs);
        const double mp      = mp_rev - mp_cost;

        if (mp > (lambda + active_tol)) lo = mid;
        else hi = mid;
      }

      E_out[p] = hi;
      sumE += hi;
    }
    return sumE;
  };

  // ---- outer bisection ----
  NumericVector E_tmp(P);
  for (int it = 0; it < n_outer; ++it) {
    const double lambda_mid = 0.5 * (lambda_lo + lambda_hi);
    const double sumE = compute_E_sum(lambda_mid, E_tmp);
    if (sumE > Etot_target) lambda_lo = lambda_mid;
    else lambda_hi = lambda_mid;
  }

  const double lambda_star = 0.5 * (lambda_lo + lambda_hi);
  const double sumE_final = compute_E_sum(lambda_star, E_star);

  // ---- accounting + marginal objective ----
  NumericVector revenue_p(P), cost_p(P), profit_p(P), marginal_objective_p(P);
  double tot_rev = 0.0, tot_cost = 0.0, tot_profit = 0.0;

  for (int p = 0; p < P; ++p) {
    if (fishable_int[p] != 1) {
      revenue_p[p] = 0.0;
      cost_p[p] = 0.0;
      profit_p[p] = 0.0;
      marginal_objective_p[p] = NA_REAL;
      continue;
    }

    const double E = E_star[p];
    const double rev = revenue_patch(p, E);
    const double cst = cost_total_patch(E, cost_patch[p], c0, gamma);

    revenue_p[p] = rev;
    cost_p[p] = cst;
    profit_p[p] = rev - cst;

    // marginal objective at solution effort
    const double mp_rev  = mp_rev_patch(p, E);
    const double mp_cost = cost_derivative_patch(E, cost_patch[p], c0, gamma, include_costs);
    marginal_objective_p[p] = mp_rev - mp_cost;

    tot_rev += rev;
    tot_cost += cst;
    tot_profit += (rev - cst);
  }

  int n_open = 0, n_active = 0;
  for (int p = 0; p < P; ++p) {
    if (fishable_int[p] == 1) {
      n_open++;
      if (E_star[p] > 1e-12) n_active++;
    }
  }

  return List::create(
    _["E_target"] = E_star,
    _["lambda"] = lambda_star,
    _["revenue_p"] = revenue_p,
    _["cost_p_total"] = cost_p,
    _["profit_p"] = profit_p,
    _["marginal_objective_p"] = marginal_objective_p,
    _["total_revenue"] = tot_rev,
    _["total_cost"] = tot_cost,
    _["total_profit"] = tot_profit,
    _["diagnostic"] = List::create(
      _["include_costs"] = include_costs,
      _["flat_fallback_used"] = false,
      _["mp0_open_sd"] = mp_sd,
      _["mp0_open_maxabs"] = mp0_maxabs,
      _["sum_E_final"] = sumE_final,
      _["n_fishable"] = n_open,
      _["n_active"] = n_active,
      _["n_outer"] = n_outer,
      _["n_inner"] = n_inner
    )
  );
}
