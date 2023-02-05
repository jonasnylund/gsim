#pragma once

#include <array>
#include <vector>

#include "numerical_types.h"

namespace model {

class Node;

// Class representing a particle in the simulation.
class Particle {
 public:
  static std::array<numerical_types::ndarray, 3> positionExtent(const std::vector<Particle>& particles);
  static void reset(std::vector<Particle>& particles);

  // Set the current accelleration and return the preferred update frequency.
  unsigned int setAccelleration(const numerical_types::ndarray& accelleration);

  void update(numerical_types::real dt,
              numerical_types::real substep_progress = 0.0,
              numerical_types::real substep_length = 1.0);

  numerical_types::ndarray position;
  numerical_types::ndarray velocity;
  numerical_types::ndarray accelleration;
  numerical_types::ndarray delta_accelleration;

  numerical_types::real mass;
  unsigned short update_frequency;
  Node* containing_node = nullptr;
};

}
