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
  std::srand(12345);

  this->particles.resize(num_particles);

  Particle big;
  big.mass = 1.5 * std::sqrt(num_particles);

  this->total_mass = big.mass;
  this->particles[0] = big;
  for (int i = 1; i < num_particles; i++) {
    this->particles[i] = randomParticle(big.mass);
    this->total_mass += this->particles[i].mass;
  }
}

void Model::initialize() {
  this->buildTree();
  this->substep_counter = 0;
  this->substep_frequency = 1;

  // Initialize the particles current accelleration.
  for (Particle& particle : this->particles) {
      numerical_types::ndarray accelleration = {0.0};
      num_interactions += this->tree.getRoot()->computeAccelleration(particle, this->theta, this->epsilon, accelleration);
      particle.setAccelleration(accelleration);
      particle.update_frequency = 1;
  }
}

void Model::buildTree() {
  if (this->timer != nullptr)
    this->timer->start(Timers::TREE);

  std::array<numerical_types::ndarray, 3> stats = Particle::positionExtent(this->particles);
  
  if (this->tree.getRoot() != nullptr && 
      this->tree.getRoot()->contains(stats[1]) &&
      this->tree.getRoot()->contains(stats[2])) {
    // If min and max extent of particle positions still lies inside the
    // root node, reuse the tree.
    this->tree.update(this->particles);
  }
  else {
    // If particles lie outside the tree, we must rebuild the structure. 
    this->tree.rebuild(this->particles);
    this->num_rebuilds++;
  }

  if (this->timer != nullptr)
    this->timer->stop(Timers::TREE);
}

void Model::updateParticles() {
  if (this->timer != nullptr)
    this->timer->start(Timers::PART);

  unsigned int num_interactions = 0;
  unsigned int max_frequency = 1;

  const numerical_types::real dtime = this->dtime / this->substep_frequency;
  // Update frequencies can only be lowered on even substeps to not loose steps.
  const bool substep_is_even = this->substep_counter % 2 == 0;

  #pragma omp parallel for \
          schedule(static) \
          reduction(+: num_interactions) \
          reduction(max: max_frequency)
  for (Particle& particle : this->particles) {
    // Use dynamic timestepping for particles, and only recompute the accelleration for
    // particles with the highest accelleration each substep.
    if (this->substep_counter % (this->substep_frequency / particle.update_frequency) == 0) {
      numerical_types::ndarray accelleration = {0.0};
      num_interactions += this->tree.getRoot()->computeAccelleration(particle, this->theta, this->epsilon, accelleration);

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
    // Step the particle position and velocity.
    const numerical_types::real substep_length = 1.0 / static_cast<numerical_types::real>(particle.update_frequency);
    // substep_progress is this particles amount of progress made towards its next
    // accelleration computation, the the substep_length is the progress to be made
    // this iteration.
    const numerical_types::real substep_progress = 
      static_cast<numerical_types::real>(this->substep_counter % particle.update_frequency) * substep_length;

    particle.update(dtime, substep_progress, substep_length);
  }

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
  
  this->num_interactions += num_interactions;
  if (this->timer != nullptr)
    this->timer->stop(Timers::PART);
}

void Model::step(numerical_types::real time) {
  if (this->timer != nullptr)
    this->timer->start(Timers::ALL);

  if (this->tree.getRoot() != nullptr) {
    int balanced_depth = static_cast<int>(
      std::log2(this->particles.size()) / std::log2(num_subnodes) + 1);
    this->tree.getRoot()->prune(balanced_depth);
  }

  numerical_types::real endtime = this->time + time;
  while (std::abs(this->time - endtime) > std::abs(this->dtime)) {
    this->substep_counter = 0;

    while (this->substep_counter < this->substep_frequency) {
      this->buildTree();
      this->updateParticles();

      this->num_iterations++;
      this->substep_counter++;
    }
    this->time += this->dtime;
  }

  if (this->timer != nullptr)
    this->timer->stop(Timers::ALL);
}

void Model::setTimeStep(numerical_types::real dtime) {
  this->dtime = dtime;
}

void Model::setTimer(Timer* timer) {
  this->timer = timer;
}

void Model::setEpsilon(numerical_types::real epsilon) {
  this->epsilon = epsilon;
}

void Model::setTheta(numerical_types::real theta) {
  this->theta = theta;
}

numerical_types::real Model::getTime() const {
  return this->time;
}

void Model::writeParticles(std::ofstream& file) const {
  if (this->timer != nullptr)
    this->timer->start(Timers::IO);

  file << this->time << ", ";
  for (const Particle& particle: this->particles) {
    file << particle.mass << ", ";
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      file << particle.position[i] << ", ";
    }
  }
  file << "\n";
  if (this->timer != nullptr)
    this->timer->stop(Timers::IO);
}

void Model::printStats() const {
  printf("\n--- Stats for run ---\n");
  if (this->timer != nullptr)
    this->timer->write();

  printf("Number of iterations:   %u\n", this->num_iterations);
  printf("Number of particles:    %lu\n", this->particles.size());
  printf("Number of interactions: %lu\n", this->num_interactions);
  printf("Number of rebuilds:     %u\n", this->num_rebuilds);
  printf("Fraction of rebuilds:   %.1f%%\n",
    100.0f * static_cast<float>(this->num_rebuilds) / static_cast<float>(this->num_iterations));
}

}  // namespace model
