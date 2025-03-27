/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * 
    * (C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#pragma once

#include <chrono>
#include <iomanip>
#include <cassert>
#include <iostream>

namespace skies { namespace launch {

/**
 * @brief Class for measuring and reporting elapsed time.
 *
 * The `Timer` class provides utilities to measure elapsed time using high-resolution clocks.
 * It supports starting and stopping timers, calculating elapsed time, and printing formatted
 * timestamps or elapsed durations.
 */
class Timer {
private:
std::chrono::high_resolution_clock::time_point start_; ///< Time point when the timer starts.
std::chrono::high_resolution_clock::time_point stop_;  ///< Time point when the timer stops.
bool is_started = false; ///< Flag indicating whether the timer is currently running.

public:
/**
 * @brief Default constructor for the `Timer` class.
 */
Timer() = default;

/**
 * @brief Starts the timer and optionally prints a message.
 *
 * @param text Optional message to print when the timer starts.
 */
void start(const std::string& text = "")
{
    is_started = true;
    start_ = std::chrono::high_resolution_clock::now();
    if (!text.empty()) {
        std::cout << text << std::endl;
    }
}

/**
 * @brief Stops the timer and optionally prints a message.
 *
 * @param text Optional message to print when the timer stops.
 */
void stop(const std::string& text = "")
{
    is_started = false;
    stop_ = std::chrono::high_resolution_clock::now();
    if (!text.empty()) {
        std::cout << text << std::endl;
    }
}

/**
 * @brief Calculates the elapsed time in seconds.
 *
 * @return unsigned The elapsed time in seconds.
 */
unsigned elapsed()
{
    auto diff = stop_ - start_;
    auto diff_ms = std::chrono::duration_cast<std::chrono::seconds>(diff);
    return diff_ms.count();
}

/**
 * @brief Prints the start timestamp with an optional message.
 *
 * @param text Optional message to accompany the start timestamp.
 */
void print_start(const std::string& text)
{
    std::time_t t_c = std::chrono::system_clock::to_time_t(start_);
    std::cout << text << std::put_time(std::localtime(&t_c), "%F %T.\n") << std::endl;
}

/**
 * @brief Prints the stop timestamp with an optional message.
 *
 * @param text Optional message to accompany the stop timestamp.
 */
void print_stop(const std::string& text)
{
    std::time_t t_c = std::chrono::system_clock::to_time_t(stop_);
    std::cout << text << std::put_time(std::localtime(&t_c), "%F %T.\n") << std::endl;
}

/**
 * @brief Prints the elapsed time in seconds with an optional message.
 *
 * @param text Optional message to accompany the elapsed time.
 */
void print_elapsed(const std::string& text)
{
    std::cout << text << elapsed() << " s" << std::endl;
}
};

} // namespace launch
} // namespace skies
