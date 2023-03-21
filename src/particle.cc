#include "particle.h"

#include <cmath>

#include <algorithm>
#include <array>
#include <vector>

namespace model {

constexpr int max_log_frequency = 8;

std::array<numerical_types::ndarray, 3> Particle::positionExtent(const std::vector<Particle>& particles) {
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

void Particle::reset(std::vector<Particle>& particles) {
  #pragma omp parallel for
  for (Particle& particle : particles) {
    particle.update_frequency = 1;
  }
}

unsigned int Particle::setAccelleration(
    const numerical_types::ndarray& accelleration,
    const numerical_types::real dtime) {
  numerical_types::real accelleration_magnitude = 0.0;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    this->accelleration[i] += this->delta_accelleration[i];
    this->delta_accelleration[i] = accelleration[i] - this->accelleration[i];
    accelleration_magnitude += accelleration[i] * accelleration[i];
  }

  accelleration_magnitude *= std::sqrt(dtime);
  // Calculate the frequency of the accelleration calculation of this particle.
  // The update frequency must be a multiple of 2.
  int log_acc = 2 + static_cast<int>(std::log2(accelleration_magnitude));
  log_acc = std::clamp(log_acc, 0, max_log_frequency);
  return 1U << log_acc;
}

void Particle::update(numerical_types::real dt,
                      numerical_types::real substep_progress,
                      numerical_types::real substep_length) {
  constexpr numerical_types::real one_half = 0.5;

  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    const numerical_types::real a1 = this->accelleration[i] + this->delta_accelleration[i] * substep_progress;
    const numerical_types::real a2 = a1 + this->delta_accelleration[i] * substep_length;

    this->velocity[i] += a1 * dt * one_half;
    this->position[i] += this->velocity[i] * dt;
    this->velocity[i] += a2 * dt * one_half;
  }
}

}  // namespace model