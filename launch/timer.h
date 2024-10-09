#pragma once

#include <chrono>
#include <iomanip>
#include <cassert>
#include <iostream>

namespace skies { namespace launch {

class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_, stop_;
    bool is_started = false;
public:
  Timer() = default;

  void start(const std::string& text = "")
  {
      is_started = true;
      start_ = std::chrono::high_resolution_clock::now();
      std::cout << text << std::endl;
  }

  void stop(const std::string& text = "")
  {
      is_started = false;
      stop_ = std::chrono::high_resolution_clock::now();
      std::cout << text << std::endl;
  }

  unsigned elapsed()
  {
      auto diff = stop_ - start_;
      auto diff_ms = std::chrono::duration_cast<std::chrono::seconds>(diff);
      return diff_ms.count();
  }

  void print_start(const std::string& text)
  {
      std::time_t t_c = std::chrono::system_clock::to_time_t(start_);
      std::cout << text << std::put_time(std::localtime(&t_c), "%F %T.\n") << std::endl;
  }

  auto print_stop(const std::string& text)
  {
      std::time_t t_c = std::chrono::system_clock::to_time_t(stop_);
      std::cout << text << std::put_time(std::localtime(&t_c), "%F %T.\n") << std::endl;
  }

  void print_elapsed(const std::string& text)
  {
      std::cout << text << elapsed() << " s" << std::endl;
  }
};

} // launch
} // skies