#pragma once

#include <array>
#include <omp.h>

namespace model {

enum Timers : int {
	ALL,
	PART,
	TREE,
	IO,
	SETUP,
	TEARDOWN,
	__COUNT,
};


class Timer {
 public:
	inline void start(enum Timers timer) {
		this->timers.at(timer) = omp_get_wtime();
	}

	void stop(enum Timers timer);
	double get(enum Timers timer) const;
	void write() const;

 private:
  std::array<double, Timers::__COUNT> total_time = {0.0};
	std::array<double, Timers::__COUNT> timers = {0.0};
};

}  // namespace model