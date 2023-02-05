#pragma once

#include <fstream>
#include <string>
#include <memory>
#include <vector>

#include "numerical_types.h"
#include "octtree.h"
#include "particle.h"
#include "timers.h"

namespace model {

class Model {
 public:
	// Initialize the model with some random data.
  void randomParticles(int num_particles);

	void initialize();

	void buildTree();
	// Calculate the accelleration on each particle.
	void updateParticles();
	// Step the simulation forward one iteration.
	void step(numerical_types::real time);
	// Set the smallest time step.
	void setTimeStep(numerical_types::real dtime);

	// Add a timer object to the model.
	void setTimer(Timer* timer);

	// Set minimum distance to use in calculations. 
	void setEpsilon(numerical_types::real epsilon);

	void setTheta(numerical_types::real theta);

	// Return the current time of the simulation.
	numerical_types::real getTime() const;

	// Wraite particle states to the output stream.
	void writeParticles(std::ofstream& file) const;

	// Print timer and counter stats to stdout.
	void printStats() const;

 private:
	std::vector<Particle> particles;
	Tree tree;

	numerical_types::real time = 0.0;
	numerical_types::real dtime = 0.1;
	numerical_types::real theta = 0.5;
	numerical_types::real epsilon = 1e-3;
	unsigned int num_iterations = 0;
	unsigned long num_interactions = 0;

	unsigned int num_rebuilds = 0;

	unsigned int substep_frequency = 1;
	unsigned int substep_counter = 0;
	Timer* timer;

};

}  // namespace model
