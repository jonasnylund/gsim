#include "gsim/common/timers.h"

namespace gsim {

std::map<std::string, Timer> Timer::timers;
int Timer::longest_key = 0;

void Timer::write() {
  std::string format;
  format.reserve(20);
  snprintf(format.data(), format.capacity(), "%%-%ds - %%.1f ms\n", Timer::longest_key);

  printf("--- Timers: ---\n");
  auto it = timers.begin();
  while (it != timers.end()) {
    printf(format.c_str(), it->first.c_str(), it->second.timeMS());
    it++;
  }
}
}  // namespace gsim
