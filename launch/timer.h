#pragma once

#include <chrono>
#include <iomanip>
#include <cassert>
#include <iostream>

#include <skies/utils/mpi_wrapper.h>

namespace skies { namespace launch {

// milliseconds, microsecond and nanoseconds
constexpr double msec_per_sec = 1000.0;
constexpr double usec_per_sec = msec_per_sec * msec_per_sec;
constexpr double nsec_per_sec = msec_per_sec * msec_per_sec * msec_per_sec;


class Timer {
public:
  Timer() = default;

  void start(const std::string& text = "") {
      assert(!is_started);
      is_started = true;
      start_ = std::chrono::high_resolution_clock::now();

      int rank{ 0 };
#ifdef SKIES_MPI
      rank = skies::mpi::rank();
#endif
      if (!rank) std::cout << text << std::endl;
  }

  void stop(const std::string& text = "") {
      assert(is_started);
      is_started = false;
      stop_ = std::chrono::high_resolution_clock::now();
      
      int rank{ 0 };
#ifdef SKIES_MPI
      rank = skies::mpi::rank();
#endif
      if (!rank) std::cout << text << std::endl;
  }

  unsigned elapsed() {
      assert(!is_started);
      auto diff = stop_ - start_;
      auto diff_ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff);
      return diff_ms.count();
  }

  void print_start(const std::string& text) {
      const std::time_t t_c = std::chrono::system_clock::to_time_t(start_);
      std::cout << text << std::put_time(std::localtime(&t_c), "%F %T.\n") << std::endl;
  }

  auto print_stop(const std::string& text) {
      const std::time_t t_c = std::chrono::system_clock::to_time_t(stop_);
      std::cout << text << std::put_time(std::localtime(&t_c), "%F %T.\n") << std::endl;
  }

  void print_elapsed(const std::string& text) {
      int rank{ 0 };
#ifdef SKIES_MPI
      rank = skies::mpi::rank();
#endif
      if (!rank) std::cout << text << elapsed() << " ms" << std::endl;
  }


private:
    std::chrono::high_resolution_clock::time_point start_, stop_;
    bool is_started = false;
};

} // launch
} // skies