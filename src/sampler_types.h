#ifndef BAYESLINREG_SAMPLER_TYPES_H
#define BAYESLINREG_SAMPLER_TYPES_H

#include <Rcpp.h>

namespace bayeslinreg {

enum class PriorModel : int {
  Normal = 0,
  SpikeSlab = 1,
  GlobalLocal = 2,
  SpikeMultiSlab = 3,
  Fixed = 4
};

constexpr int prior_model_count = 5;

inline PriorModel prior_model_from_code(const int code) {
  switch (code) {
    case 0: return PriorModel::Normal;
    case 1: return PriorModel::SpikeSlab;
    case 2: return PriorModel::GlobalLocal;
    case 3: return PriorModel::SpikeMultiSlab;
    case 4: return PriorModel::Fixed;
    default: Rcpp::stop("Unknown coefficient-prior model code.");
  }
  return PriorModel::Normal;
}

inline int prior_model_index(const PriorModel model) {
  return static_cast<int>(model);
}

enum class GramStorage : int {
  Dense = 0,
  SparseGeneral = 1,
  SparseLowerTriangular = 2
};

inline GramStorage gram_storage_from_code(const int code) {
  switch (code) {
    case 0: return GramStorage::Dense;
    case 1: return GramStorage::SparseGeneral;
    case 2: return GramStorage::SparseLowerTriangular;
    default: Rcpp::stop("Unknown Gram-block storage type.");
  }
  return GramStorage::Dense;
}

enum class PveType : int {
  Standalone = 0,
  Allocated = 1
};

inline PveType pve_type_from_code(const int code) {
  switch (code) {
    case 0: return PveType::Standalone;
    case 1: return PveType::Allocated;
    default: Rcpp::stop("Unknown PVE type code.");
  }
  return PveType::Standalone;
}

}  // namespace bayeslinreg

#endif
