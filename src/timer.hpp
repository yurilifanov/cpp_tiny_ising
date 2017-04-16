#ifndef _TIMER_
#define _TIMER_
#include <chrono>
template<typename Type> 
class Timer {
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = Clock::time_point;
    TimePoint tick_point;
  public:
    Timer() { tick_point = Clock::now(); }
    void tic() { tick_point = Clock::now(); }
    Type toc() { 
      using namespace std::chrono;
      return duration_cast<duration<Type>>(Clock::now() - tick_point).count(); 
    }
};
#endif