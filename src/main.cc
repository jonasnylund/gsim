#include <cstring>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <typeinfo> 

#include <omp.h>

#include "numerical_types.h"
#include "model.h"
#include "timers.h"


int main(int argc, char *argv[]) {
  int num_particles = 300;
  int num_iterations = 300;
  numerical_types::real theta = 0.5;
  numerical_types::real epsilon = 0.25;
  numerical_types::real dt = 0.01;
  float substep_ratio = 0.25;
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
    else if (strncmp(argv[i], "--updateratio", 14) == 0) {
      substep_ratio = atof(argv[++i]);
    }
    else if (strncmp(argv[i], "-q", 2) == 0 || strncmp(argv[i], "--quiet", 8) == 0) {
      verbose = !atoi(argv[++i]);
    }
    else {
      printf("Unknown argument: %s\n", argv[i]);
    }

    i++;
  }
  if (verbose) {
    printf("Simulating %d particles for %d timesteps\n", num_particles, num_iterations);
  }

  std::ofstream particles_file(particles_path);
  std::ofstream tree_file(tree_path);

  model::Model model;

  model::Timer::byName("Setup")->set();

  model.setEpsilon(epsilon);
  model.setTheta(theta);
  model.setTimeStep(dt);
  model.setSubstepUpdateRatio(substep_ratio);

  model.randomParticles(num_particles);
  model.initialize();
  model.writeParticles(particles_file);
  model.writeTree(tree_file);

  model::Timer::byName("Setup")->reset();

  while (model.getTime() < num_iterations) {
    model.step(1.0);
    model.writeParticles(particles_file);
    model.writeTree(tree_file);

    printf("T: %.1f s\r", model.getTime());
    std::fflush(stdout);
  }

  printf("T: %.1f s\n", model.getTime());

  model::Timer::byName("Teardown")->set();

  particles_file.close();
  tree_file.close();
  if (verbose) {
    printf("Writing output to %s\n", particles_path.c_str());
  }

  model::Timer::byName("Teardown")->reset();

  if (verbose) {
    model.printStats();
    printf("\n");
    model::Timer::write();
  }

  return 0;
}
