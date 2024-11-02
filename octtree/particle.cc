#include "gsim/octtree/particle.h"

#include <cmath>

#include <algorithm>
#include <array>
#include <vector>

#include "gsim/common/numerical_types.h"

namespace gsim {

std::array<numerical_types::ndarray, 3> Particle::positionExtent(
    const std::vector<Particle>& particles) {
  numerical_types::ndarray min = {0.0};
  numerical_types::ndarray max = {0.0};
  numerical_types::ndarray average = {0.0};

  for (const Particle& particle : particles) {
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      average[i] += particle.position[i];
      if (particle.position[i] < min[i])
        min[i] = particle.position[i];
      else if (particle.position[i] > max[i])
        max[i] = particle.position[i];
    }
  }

  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    average[i] /= particles.size();	 // Normalize center position.
  }

  return {average, min, max};
}

void Particle::setAccelleration(
    const numerical_types::ndarray& accelleration,
    numerical_types::real dt) {
  const numerical_types::real one_by_dt = 1.0 / dt;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    this->accelleration_change[i] = (
      accelleration[i] - this->current_accelleration[i]) * one_by_dt;
  }
}

void Particle::update(
    const numerical_types::real previous_dt,
    const numerical_types::real next_dt) {
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    this->current_accelleration[i] += this->accelleration_change[i] * next_dt;
    this->velocity[i] += this->current_accelleration[i] * 0.5 * previous_dt;

    this->velocity[i] += this->current_accelleration[i] * 0.5 * next_dt;
    this->position[i] += this->velocity[i] * next_dt;
  }
}

}  // namespace gsim