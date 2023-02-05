#include <cstring>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <typeinfo> 

#include "numerical_types.h"
#include "model.h"
#include "timers.h"


int main(int argc, char *argv[]) {
	int num_particles = 10;
	int num_iterations = 1000;
	numerical_types::real theta = 0.5;
	numerical_types::real epsilon = 0.2;
	numerical_types::real dt = 0.01;
	bool verbose = true;

	std::string filename("output.csv");

	int i = 0;
	while (i < argc) {
		if (strncmp(argv[i], "-p", 2) == 0 || strncmp(argv[i], "--particles", 11) == 0) {
			num_particles = atoi(argv[++i]);
		}
		else if (strncmp(argv[i], "-n", 2) == 0 || strncmp(argv[i], "--iterations", 12) == 0) {
			num_iterations = atoi(argv[++i]);
		}
		else if (strncmp(argv[i], "-o", 2) == 0 || strncmp(argv[i], "--output", 8) == 0) {
			filename = argv[++i];
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
		else if (strncmp(argv[i], "-q", 2) == 0 || strncmp(argv[i], "--quiet", 8) == 0) {
			verbose = !atoi(argv[++i]);
		}

		i++;
	}
	if (verbose) {
		printf("Simulating %d particles for %d timesteps\n", num_particles, num_iterations);
		printf("Writing output to %s\n", filename.c_str());
	}

	std::ofstream file(filename);

	model::Timer timer;
	model::Model model;

	timer.start(model::Timers::SETUP);

	model.setTimer(&timer);
	model.setEpsilon(epsilon);
	model.setTheta(theta);
	model.setTimeStep(dt);

	model.randomParticles(num_particles);
	model.initialize();
	model.writeParticles(file);

	timer.stop(model::Timers::SETUP);

	for (int i = 0; i < num_iterations; i++) {
		model.step(1.0);
		model.writeParticles(file);

		printf("T: %.1f s\r", model.getTime());
		std::fflush(stdout);
	}

	printf("T: %.1f s\n", model.getTime());

	timer.start(model::Timers::TEARDOWN);
	file.close();

	timer.stop(model::Timers::TEARDOWN);

	if (verbose) {
		model.printStats();
	}

	return 0;
}
