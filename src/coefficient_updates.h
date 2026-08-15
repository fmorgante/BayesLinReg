#ifndef BAYESLINREG_COEFFICIENT_UPDATES_H
#define BAYESLINREG_COEFFICIENT_UPDATES_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "sampler_types.h"

namespace bayeslinreg {

// Adapter that gives the serial coefficient kernel the same small RNG
// interface as the independent block streams while retaining R's RNG and,
// in particular, R's Bernoulli draw for fixed-seed compatibility.
class RSessionRng {
 public:
  double uniform() { return R::runif(0.0, 1.0); }
  double normal() { return R::rnorm(0.0, 1.0); }
  int bernoulli(const double probability) {
    return static_cast<int>(R::rbinom(1.0, probability));
  }
};

struct MixtureWorkspace {
  std::vector<double> weights;
  std::vector<double> conditional_vars;
  std::vector<double> conditional_means;

  void ensure_size(const int size) {
    if (static_cast<int>(weights.size()) >= size) return;
    weights.resize(size);
    conditional_vars.resize(size);
    conditional_means.resize(size);
  }
};

// Draw one coefficient conditional on the current partial residual. The
// serial and block-parallel sweeps share these formulas but supply different
// RNG implementations and retain their own matrix-state bookkeeping.
template <typename Rng>
double draw_coefficient_coordinate(
    const PriorModel model,
    const int prior_block,
    const int local_index,
    const double x_squared,
    const double residual_var,
    const double conditional_numerator,
    const std::vector<double>& normal_var,
    const std::vector<double>& pi,
    const std::vector<double>& slab_var,
    const std::vector<double>& tau_sq,
    const std::vector<double>& local_var,
    const std::vector< std::vector<double> >& multi_gamma,
    const std::vector< std::vector<double> >& multi_pi,
    const std::vector<double>& multi_var,
    std::vector<int>& inclusion,
    std::vector<int>& multi_component,
    Rng& generator,
    MixtureWorkspace& workspace) {
  if (model == PriorModel::SpikeMultiSlab) {
    const int component_count = multi_gamma[prior_block].size();
    workspace.ensure_size(component_count);
    std::vector<double>& weights = workspace.weights;
    std::vector<double>& conditional_vars = workspace.conditional_vars;
    std::vector<double>& conditional_means = workspace.conditional_means;
    double maximum_log_weight = -std::numeric_limits<double>::infinity();
    for (int component = 0; component < component_count; ++component) {
      weights[component] = std::log(std::max(
        multi_pi[prior_block][component],
        std::numeric_limits<double>::min()
      ));
      if (component > 0) {
        const double prior_var =
          multi_gamma[prior_block][component] * multi_var[prior_block];
        conditional_vars[component] = 1.0 / (
          x_squared / residual_var + 1.0 / prior_var
        );
        conditional_means[component] = conditional_vars[component] *
          conditional_numerator / residual_var;
        weights[component] +=
          0.5 * std::log(conditional_vars[component] / prior_var) +
          conditional_means[component] * conditional_means[component] /
            (2.0 * conditional_vars[component]);
      }
      maximum_log_weight = std::max(
        maximum_log_weight, weights[component]
      );
    }
    double weight_total = 0.0;
    for (int component = 0; component < component_count; ++component) {
      weights[component] = std::exp(
        weights[component] - maximum_log_weight
      );
      weight_total += weights[component];
    }
    const double threshold = generator.uniform() * weight_total;
    double cumulative_weight = 0.0;
    int selected_component = component_count - 1;
    for (int component = 0; component < component_count; ++component) {
      cumulative_weight += weights[component];
      if (threshold <= cumulative_weight) {
        selected_component = component;
        break;
      }
    }
    multi_component[local_index] = selected_component;
    return selected_component == 0
      ? 0.0
      : conditional_means[selected_component] +
          std::sqrt(conditional_vars[selected_component]) * generator.normal();
  }

  const double prior_precision = model == PriorModel::Fixed
    ? 0.0
    : (model == PriorModel::GlobalLocal
        ? 1.0 / tau_sq[prior_block] / local_var[local_index]
        : (model == PriorModel::SpikeSlab
            ? 1.0 / slab_var[prior_block]
            : 1.0 / normal_var[prior_block]));
  const double conditional_var = 1.0 / (
    x_squared / residual_var + prior_precision
  );
  const double conditional_mean =
    conditional_var * conditional_numerator / residual_var;
  if (model == PriorModel::SpikeSlab) {
    const double epsilon = std::numeric_limits<double>::epsilon();
    const double bounded_pi = std::min(
      std::max(pi[prior_block], epsilon), 1.0 - epsilon
    );
    const double log_inclusion_odds =
      std::log(bounded_pi) - std::log1p(-bounded_pi) +
      0.5 * std::log(conditional_var / slab_var[prior_block]) +
      conditional_mean * conditional_mean / (2.0 * conditional_var);
    const double inclusion_probability = log_inclusion_odds >= 0.0
      ? 1.0 / (1.0 + std::exp(-log_inclusion_odds))
      : std::exp(log_inclusion_odds) /
          (1.0 + std::exp(log_inclusion_odds));
    inclusion[local_index] = generator.bernoulli(inclusion_probability);
  }
  return model != PriorModel::SpikeSlab || inclusion[local_index] == 1
    ? conditional_mean + std::sqrt(conditional_var) * generator.normal()
    : 0.0;
}

}  // namespace bayeslinreg

#endif
