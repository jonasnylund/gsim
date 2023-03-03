#include "model.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <omp.h>

#include "numerical_types.h"
#include "octtree.h"
#include "particle.h"
#include "timers.h"

namespace model {
namespace {

float random() {
  return static_cast <numerical_types::real> (std::rand()) / static_cast <numerical_types::real> (RAND_MAX);
};

Particle randomParticle(numerical_types::real pseudo_mass) {
  Particle particle;

  numerical_types::real x = 5 + 10 * random() + std::pow(random(), 2) * 200.0;
  numerical_types::real angle = random() * 2 * M_PI;
  numerical_types::real v = 0.1 + std::sqrt(numerical_types::G * pseudo_mass / x);

  particle.position[0] = x * std::cos(angle);
  particle.position[1] = x * std::sin(angle);
  particle.velocity[0] = -v * std::sin(angle);
  particle.velocity[1] = v * std::cos(angle);

  for (int i = 2; i < numerical_types::num_dimensions; i++) {
    particle.velocity[i] = (2 * random() - 1) / 20;
    particle.position[i] = (2 * random() - 1) / 20;
  }

  particle.mass = std::pow(random(), 3);

  return particle;
}

}  // namespace

void Model::randomParticles(int num_particles) {
  assert(num_particles > 0);
  std::srand(12345);

  this->particles.resize(num_particles);

  Particle big;
  big.mass = 10.0 + static_cast<numerical_types::real>(num_particles) / 10.0;

  this->total_mass = big.mass;
  this->particles[0] = big;
  for (int i = 1; i < num_particles; i++) {
    this->particles[i] = randomParticle(big.mass);
    this->total_mass += this->particles[i].mass;
  }
}

void Model::initialize() {
  this->rebuildTree();
  this->substep_counter = 0;
  this->substep_frequency = 1;

  // Initialize the particles current accelleration.
  for (Particle& particle : this->particles) {
      numerical_types::ndarray accelleration = {0.0};
      num_interactions += this->tree.computeAccelleration(
          particle, this->theta, this->epsilon, accelleration);
      particle.setAccelleration(accelleration);
      particle.update_frequency = 1;
  }
}

void Model::rebuildTree() {
  Timer::byName("Tree: total")->set();

  if (!this->tree.empty()) {
    // Attempt to update the existing tree if possible.
    this->num_rebuilds += static_cast<unsigned int>(
      !this->tree.update(this->particles));
  }
  else {
    this->tree.rebuild(this->particles);
    this->num_rebuilds++;
  }

  Timer::byName("Tree: total")->reset();
}

void Model::updateTree() {
  Timer::byName("Tree: total")->set();
  this->tree.update();
  Timer::byName("Tree: total")->reset();
}

void Model::updateParticles() {
  assert(this->substep_counter >= 0);
  assert(this->substep_counter < this->substep_frequency);
  assert(this->substep_frequency > 0);

  Timer::byName("Particles: total")->set();

  unsigned int num_interactions = 0;
  unsigned int max_frequency = 1;

  // Update frequencies can only be lowered on even substeps to not loose steps.
  const bool substep_is_even = this->substep_counter % 2 == 0;

  Timer::byName("Particles: accelleration")->set();
  #pragma omp parallel for \
          schedule(static) \
          reduction(+: num_interactions) \
          reduction(max: max_frequency)
  for (Particle& particle : particles) {
    if (this->substep_counter % (
          this->substep_frequency / particle.update_frequency) != 0) {
      continue;
    }

    numerical_types::ndarray accelleration;
    accelleration.fill(0.0);
    num_interactions += this->tree.computeAccelleration(
        particle, this->theta, this->epsilon, accelleration);

    // Set the new accelleration and get the preferred new update frequency.
    unsigned int frequency = particle.setAccelleration(accelleration);

    // Check if the update frequency should and/or can be changed.
    if (particle.update_frequency < frequency) {
      particle.update_frequency = frequency;
    }
    else if (particle.update_frequency > frequency && substep_is_even) {
      particle.update_frequency /= 2;
    }
    // Since we only lower the substep_frequency by half on even substeps,
    // all particles that is not updated have a frequency that is lower or equal
    // to the possibly new one. Thus we only need to check the particles we
    // update.
    if (max_frequency < particle.update_frequency) {
      max_frequency = particle.update_frequency;
    }
  }
  Timer::byName("Particles: accelleration")->reset();

  // If the highest update frequency changed, we update the counter to reflect.
  if (this->substep_frequency < max_frequency) {
    int ratio = max_frequency / this->substep_frequency;
    this->substep_counter *= ratio;
    this->substep_frequency = max_frequency;
  }
  else if (this->substep_frequency > max_frequency && substep_is_even) {
    this->substep_frequency /= 2;
    this->substep_counter /= 2;
  }

  const numerical_types::real dtime = this->dtime / this->substep_frequency;

  Timer::byName("Particles: position")->set();
  // Step the particle position and velocity.
  #pragma omp parallel for schedule(static)
  for (Particle& particle : this->particles) {
    const int update_period = this->substep_frequency / particle.update_frequency;
    const int update_counter = this->substep_counter % update_period;

    const numerical_types::real substep_length = (
      1.0 / static_cast<numerical_types::real>(update_period));
    // substep_progress is this particles amount of progress made towards its next
    // accelleration computation, the the substep_length is the progress to be made
    // this iteration.
    const numerical_types::real substep_progress = update_counter * substep_length;

    particle.update(dtime, substep_progress, substep_length);
  }
  Timer::byName("Particles: position")->reset();

  this->num_interactions += num_interactions;
  Timer::byName("Particles: total")->reset();
}

void Model::step(numerical_types::real time) {
  Timer::byName("Timestepping")->set();
  
  numerical_types::real endtime = this->time + time;
  while (std::abs(this->time - endtime) > std::abs(this->dtime)) {
    this->substep_counter = 0;
    this->rebuildTree();

    while (this->substep_counter < this->substep_frequency) {
      // Only update the tree without calculating new particle
      // positions. But don't update the tree twice on the first
      // iteration when we already rebuilt the tree.
      if (this->substep_counter > 0)
        this->updateTree();

      this->updateParticles();

      this->num_iterations++;
      this->substep_counter++;
    }
    this->time += this->dtime;
  }
  // TODO: Reenable once new allocation is fixed.
  // this->tree.prune();

  Timer::byName("Timestepping")->reset();
}

void Model::setTimeStep(numerical_types::real dtime) {
  this->dtime = dtime;
}

void Model::setEpsilon(numerical_types::real epsilon) {
  assert(epsilon >= 0.0);
  this->epsilon = epsilon;
}

void Model::setTheta(numerical_types::real theta) {
  assert(theta >= 0.0);
  this->theta = theta;
}

numerical_types::real Model::getTime() const {
  return this->time;
}

void Model::writeParticles(std::ofstream& file) const {
  Timer::byName("IO")->set();

  file << this->time << ", ";
  for (const Particle& particle: this->particles) {
    file << particle.mass << ", ";
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      file << particle.position[i] << ", ";
    }
  }
  file << "\n";
  Timer::byName("IO")->reset();
}

void Model::writeTree(std::ofstream& file) const {
  Timer::byName("IO")->set();

  file << this->time << ", ";
  this->tree.write(file);
  file << "\n";

  Timer::byName("IO")->reset();
}

void Model::printStats() const {
  printf("\n--- Stats for run: ---\n");
  printf("Number of iterations:   %u\n", this->num_iterations);
  printf("Number of interactions: %lu\n", this->num_interactions);
  printf("Interactions/iteration: %lu\n", this->num_interactions / this->num_iterations);
  printf("Number of particles:    %lu\n", this->particles.size());
  printf("Number of rebuilds:     %u\n", this->num_rebuilds);
  printf("Fraction of rebuilds:   %.1f%%\n",
    100.0f * static_cast<float>(this->num_rebuilds) / static_cast<float>(this->num_iterations));
}

}  // namespace model
