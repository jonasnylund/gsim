#pragma once

#include <array>
#include <vector>

#include "gsim/common/numerical_types.h"

namespace gsim {

// Class representing a particle in the simulation.
class Particle {
 public:
  static std::array<numerical_types::ndarray, 3> positionExtent(
    const std::vector<Particle>& particles);

  // Set the current accelleration and return the preferred update frequency.
  void setAccelleration(
    const numerical_types::ndarray& accelleration,
    numerical_types::real dt);

  void update(numerical_types::real previous_dt, numerical_types::real next_dt);

  numerical_types::ndarray position;
  numerical_types::ndarray velocity;
  numerical_types::ndarray current_accelleration;
  numerical_types::ndarray accelleration_change;

  numerical_types::real mass;
  numerical_types::NodeKey containing_node = numerical_types::emptykey;
  int update_frequency = 1;
};

}  // namespace gsim
