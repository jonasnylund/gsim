#include "gsim/model/model.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <omp.h>

#include "gsim/common/numerical_types.h"
#include "gsim/common/timers.h"
#include "gsim/octtree/octtree.h"
#include "gsim/octtree/particle.h"

namespace gsim {
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

numerical_types::real magnitude(const numerical_types::ndarray& vector) {
  numerical_types::real m = 0.0;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    m += vector[i] * vector[i];
  }
  return m;
} 

}  // namespace

void Model::addParticle(
    numerical_types::real mass,
    numerical_types::ndarray position,
    numerical_types::ndarray velocity) {
  Particle particle;
  particle.mass = mass;
  particle.position = position;
  particle.velocity = velocity;
  this->particles.push_back(particle);
}

void Model::randomParticles(int num_particles) {
  assert(num_particles > 0);
  std::srand(12345);

  this->particles.resize(num_particles);

  Particle big;
  big.mass = 10.0 + static_cast<numerical_types::real>(num_particles) / 10.0;

  this->particles[0] = big;
  for (int i = 1; i < num_particles; i++) {
    this->particles[i] = randomParticle(big.mass);
  }
}

void Model::initialize() {
  // Avoid hyperthreading.
  omp_set_num_threads(omp_get_num_procs() / 2);

  this->rebuildTree();

  // Initialize the particles current accelleration.
  for (Particle& particle : this->particles) {
      numerical_types::ndarray accelleration = {0.0};
      num_interactions += this->tree.computeAccelleration(
          particle, this->theta, this->epsilon, accelleration);
      particle.setAccelleration(accelleration, this->dtime);
  }
}

void Model::rebuildTree() {
  Timer::byName("Tree: total")->set();

  if (!this->tree.empty()) {
    // Attempt to update the existing tree if possible.
    if (!this->tree.relocate(this->particles)) {
      this->num_rebuilds++;
    }
  }
  else {
    this->tree.rebuild(this->particles);
    this->num_rebuilds++;
  }

  Timer::byName("Tree: total")->reset();
}

void Model::updateTree() {
  assert(this->substep_counter >= 0);
  assert(this->substep_counter < this->substep_frequency);
  assert(this->substep_frequency > 0);
  Timer::byName("Tree: total")->set();
  this->tree.update();
  Timer::byName("Tree: total")->reset();
}

int Model::getPreferredSubstepping(numerical_types::real accelleration_magnitude) const {
  assert(accelleration_magnitude >= 0.0);
  assert(this->substep_max_accelleration > 0.0);
  if (this->substep_max_accelleration <= accelleration_magnitude)
    return this->max_substep_frequency;

  const float ratio = accelleration_magnitude / this->substep_max_accelleration;
  return 1 << static_cast<int>(ratio * ratio * this->max_substep_frequency);
}

void Model::updateParticles() {
  assert(this->substep_counter >= 0);
  assert(this->substep_counter < this->substep_frequency);
  assert(this->substep_frequency > 0);

  Timer::byName("Particles: total")->set();

  unsigned int num_interactions = 0;
  int max_frequency = 1;

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
    const numerical_types::real accelleration_magnitude = magnitude(accelleration);

    int frequency = this->getPreferredSubstepping(accelleration_magnitude);
    if (frequency > this->max_substep_frequency)
      frequency = this->max_substep_frequency;

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

    particle.setAccelleration(accelleration, this->dtime / particle.update_frequency);
  }
  Timer::byName("Particles: accelleration")->reset();

  const numerical_types::real previous_dt = this->dtime / this->substep_frequency;

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

  const numerical_types::real next_dt = this->dtime / this->substep_frequency;
  numerical_types::real highest_accelleration = 0.0;

  Timer::byName("Particles: position")->set();
  // Step the particle position and velocity.
  #pragma omp parallel for schedule(static) \
          reduction(max: highest_accelleration)
  for (Particle& particle : this->particles) {
    particle.update(previous_dt, next_dt);

    const numerical_types::real accelleration_magnitude = magnitude(
      particle.current_accelleration);
    if (accelleration_magnitude > highest_accelleration) {
      highest_accelleration = accelleration_magnitude;
    }
  }
  Timer::byName("Particles: position")->reset();

  if (highest_accelleration > this->substep_max_accelleration) {
    this->substep_max_accelleration = highest_accelleration;
  }
  this->num_interactions += num_interactions;
  Timer::byName("Particles: total")->reset();
}

void Model::step(numerical_types::real time) {
  Timer::byName("Timestepping")->set();
  
  numerical_types::real t = 0.0;
  while (t < time) {
    this->substep_counter = 0;
    this->substep_max_accelleration *= static_cast<numerical_types::real>(0.95);

    while (this->substep_counter < this->substep_frequency) {
      this->rebuildTree();
      this->updateTree();

      this->updateParticles();
      
      this->num_iterations++;
      this->substep_counter++;
    }
    this->tree.prune();
    t += this->dtime;
    this->time += this->dtime;
  }

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

void Model::setSubstepMaxFrequency(int frequency) {
  assert(frequency > 0);
  assert(frequency ==  1U << static_cast<int>(std::log2(frequency)));
  this->max_substep_frequency = frequency;
}

numerical_types::real Model::getTime() const {
  return this->time;
}

void Model::writeParticles(const std::string& path) const {
  Timer::byName("IO")->set();
  std::ofstream file(path, std::ios_base::app);

  const size_t real_bytes = sizeof(numerical_types::real);
  const size_t n_particles = this->particles.size();
  const size_t bytes_to_write = real_bytes + n_particles * (
    numerical_types::ndarray_bytes_size + real_bytes);

  // First write out the number of bytes to write.
  file.write(reinterpret_cast<const char*>(&bytes_to_write), sizeof(bytes_to_write));
  file.write(reinterpret_cast<const char*>(&this->time), real_bytes);

  for (const Particle& particle: this->particles) {
    file.write(reinterpret_cast<const char*>(&particle.mass), real_bytes);
    file.write(reinterpret_cast<const char*>(particle.position.data()),
               numerical_types::ndarray_bytes_size);
  }

  file.close();
  Timer::byName("IO")->reset();
}

void Model::writeTree(const std::string& path) const {
  Timer::byName("IO")->set();
  std::ofstream file(path, std::ios_base::app);

  const size_t real_bytes = sizeof(numerical_types::real);
  const size_t num_nodes = this->tree.numNodes();
  const size_t bytes_to_write = real_bytes + num_nodes * (
    numerical_types::ndarray_bytes_size + real_bytes);

  // First write out the number of bytes to write.
  file.write(reinterpret_cast<const char*>(&bytes_to_write), sizeof(bytes_to_write));
  file.write(reinterpret_cast<const char*>(&this->time), real_bytes);
  this->tree.write(file);

  file.close();
  Timer::byName("IO")->reset();
}

void Model::printStats() const {
  printf("\n--- Final state: ---\n");
  printf("Time:          %.0f\n", this->time);
  printf("Num particles: %d\n", this->tree.countLeafNodes().second);
  printf("Tree height:   %d\n", this->tree.height());
  printf("Tree nodes:    %d\n", this->tree.numNodes());
  printf("Leaf nodes:    %d\n", this->tree.countLeafNodes().first);
  printf("Dimensions:    %d\n", numerical_types::num_dimensions);

  printf("\nParameters:\n");
  printf("Delta T:  %.3f\n", this->dtime);
  printf("Theta:    %.3f\n", this->theta);
  printf("Epsilon:  %.3f\n", this->epsilon);
  printf("Substeps: %u\n", this->max_substep_frequency);


  printf("\n--- Stats for run: ---\n");
  printf("Number of iterations:   %lu\n", this->num_iterations);
  printf("Number of interactions: %lu\n", this->num_interactions);
  printf("Interactions/iteration: %lu\n", this->num_interactions / this->num_iterations);
  printf("Number of particles:    %lu\n", this->particles.size());
  printf("Number of rebuilds:     %u\n", this->num_rebuilds);
  printf("Fraction of rebuilds:   %.1f%%\n",
    100.0f * static_cast<float>(this->num_rebuilds) / static_cast<float>(this->num_iterations));
}

}  // namespace gsim
