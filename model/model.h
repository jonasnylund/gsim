#pragma once

#include <fstream>
#include <string>
#include <memory>
#include <vector>

#include "gsim/common/numerical_types.h"
#include "gsim/common/timers.h"
#include "gsim/octtree/octtree.h"
#include "gsim/octtree/particle.h"

namespace gsim {

class Model {
 public:
  // Initialize the model with some random data.
  void randomParticles(int num_particles);

  void initialize();

  // Step the simulation forward one iteration.
  void step(numerical_types::real time);
  // Set the smallest time step.
  void setTimeStep(numerical_types::real dtime);

  // Set minimum distance to use in calculations. 
  void setEpsilon(numerical_types::real epsilon);

  void setTheta(numerical_types::real theta);

  void setSubstepMaxFrequency(int frequency);

  // Return the current time of the simulation.
  numerical_types::real getTime() const;

  // Write particle states to the output stream.
  void writeParticles(std::ofstream& file) const;

  // Write the current tree structure to file.
  void writeTree(std::ofstream& file) const;

  // Print timer and counter stats to stdout.
  void printStats() const;

 private:
  // Fully updates the tree.
  void rebuildTree();
  // Partially updates the tree.
  void updateTree();
  // Calculate the accelleration on each particle.
  void updateParticles();

  int getPreferredSubstepping(numerical_types::real accelleration_magnitude) const;

  std::vector<Particle> particles;
  Tree tree;

  numerical_types::real time = 0.0;
  numerical_types::real dtime = 0.1;
  numerical_types::real theta = 0.5;
  numerical_types::real epsilon = 1e-3;
  numerical_types::real substep_max_accelleration = 1.0;

  unsigned long num_iterations = 0;
  unsigned long num_interactions = 0;

  unsigned int num_rebuilds = 0;

  int substep_frequency = 1;
  int substep_counter = 0;
  int max_substep_frequency = 1;

  numerical_types::real total_mass = 0.0;
};

}  // namespace gsim
