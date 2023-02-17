#include "timers.h"

namespace model {

std::map<std::string, Timer> Timer::timers;

void Timer::write() {
  printf("--- Timers: ---\n");
  auto it = timers.begin();
  while (it != timers.end()) {
    printf("%15s: %.1f ms\n", it->first.c_str(), it->second.timeMS());
    it++;
  }
}
}  // namespace model
