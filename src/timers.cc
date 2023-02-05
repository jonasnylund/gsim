#include "timers.h"

#include <cstdio>

#include <array>

namespace model {

void Timer::stop(enum Timers timer) {
	const double now = omp_get_wtime();
	const double time = now - this->timers.at(timer);
	this->total_time.at(timer) += time;
	this->timers.at(timer) = now;
}

double Timer::get(enum Timers timer) const {
	return this->total_time.at(timer) * 1000;
}

void Timer::write() const {
	printf("Simulation time: %.1f ms\n", this->get(Timers::ALL));
	printf("Tree construct:  %.1f ms\n", this->get(Timers::TREE));
	printf("Particle motion: %.1f ms\n", this->get(Timers::PART));
	printf("Setup:           %.1f ms\n", this->get(Timers::SETUP));
	printf("Tear down:       %.1f ms\n", this->get(Timers::TEARDOWN));
	printf("File IO:         %.1f ms\n\n", this->get(Timers::IO));
}

}  // namespace model
