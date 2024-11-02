#pragma once

#include <cstring>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <omp.h>

namespace model {

class Timer {
 public:
  static inline Timer* byName(const std::string& id) {
    if (id.size() > Timer::longest_key) {
      Timer::longest_key = id.size();
    }
    return &Timer::timers[id];
  }
  static void write();

  inline void set() { this->set_time = omp_get_wtime(); }
  inline void reset() {
    const double t = omp_get_wtime();
    this->total_time += t - this->set_time;
    this->set_time = t;
  }
  double timeMS() const { return this->total_time * 1000.0; }

 private:
  static std::map<std::string, Timer> timers;
  static int longest_key;
  double set_time = 0.0f;
  double total_time = 0.0f;
};

}  // namespace model