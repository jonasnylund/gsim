#include <cstring>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <string>

#include "gsim/common/timers.h"
#include "gsim/model/model.h"


int main(int argc, char *argv[]) {
  int num_particles = 300;
  int num_iterations = 300;
  float theta = 0.4;
  float epsilon = 0.1;
  float dt = 0.01;
  unsigned int substep_frequency = 4;
  bool verbose = true;

  std::string particles_path("particles.bin");
  std::string tree_path("tree.bin");

  int i = 1;
  while (i < argc) {
    if (strncmp(argv[i], "-p", 2) == 0 || strncmp(argv[i], "--particles", 11) == 0) {
      num_particles = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-n", 2) == 0 || strncmp(argv[i], "--iterations", 12) == 0) {
      num_iterations = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-o", 2) == 0 || strncmp(argv[i], "--output", 8) == 0) {
      particles_path = argv[++i];
    }
    else if (strncmp(argv[i], "-t", 2) == 0 || strncmp(argv[i], "--tree", 8) == 0) {
      tree_path = argv[++i];
    }
    else if (strncmp(argv[i], "--theta", 7) == 0) {
      theta = atof(argv[++i]);
    }
    else if (strncmp(argv[i], "--epsilon", 7) == 0) {
      epsilon = atof(argv[++i]);
    }
    else if (strncmp(argv[i], "--dtime", 7) == 0) {
      dt = atof(argv[++i]);
    }
    else if (strncmp(argv[i], "--substeps", 14) == 0) {
      substep_frequency = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-q", 2) == 0 || strncmp(argv[i], "--quiet", 8) == 0) {
      verbose = false;
    }
    else {
      printf("Unknown argument: %s\n", argv[i]);
    }

    i++;
  }
  if (verbose) {
    printf("Simulating %d particles for %d timesteps\n", num_particles, num_iterations);
  }

  gsim::Model model;

  model.setEpsilon(epsilon);
  model.setTheta(theta);
  model.setTimeStep(dt);
  model.setSubstepMaxFrequency(substep_frequency);

  model.randomParticles(num_particles);
  model.initialize();
  model.writeParticles(particles_path);
  model.writeTree(tree_path);

  for (int i = 1; i <= num_iterations; i++) {
    model.step(static_cast<numerical_types::real>(i) - model.getTime());

    model.writeParticles(particles_path);
    model.writeTree(tree_path);

    if (verbose) {
      printf("T: %.1f s\r", model.getTime());
      std::fflush(stdout);
    }
  }

  gsim::Timer::byName("Teardown")->set();

  if (verbose) {
    printf("Writing output to %s\n", particles_path.c_str());
  }

  model.printStats();

  return 0;
}
