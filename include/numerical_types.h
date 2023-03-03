#pragma once

#include <array>

namespace numerical_types {

// Number of dimensions in the simulation.
constexpr int num_dimensions = 2;

// Machine precision type.
typedef double real;

using ndarray = std::array<real, num_dimensions>;

// Size of an ndarray in bytes.
constexpr int ndarray_bytes_size = sizeof(real) * num_dimensions;

// Newtons constant of gravitation.
constexpr real G = 6.6743015e-1;
// constexpr real G = 6.6743015e-11;

typedef uint32_t index_key_t;
typedef uint16_t depth_key_t;
struct NodeKey {
  depth_key_t depth;
  index_key_t index;
};

inline bool operator==(const NodeKey& lhs, const NodeKey& rhs) {
  return lhs.depth == rhs.depth && lhs.index == rhs.index;
}
inline bool operator!=(const NodeKey& lhs, const NodeKey& rhs) {
  return !(lhs == rhs);
}

constexpr NodeKey emptykey = { 0xffff, 0xffffffff };


}  // namespace numerical_types
