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

#include <launch/types.h>

namespace skies { namespace launch {

/**
 * @brief Parses command-line options and stores them in the provided arguments and options structures.
 *
 * This function processes the command-line arguments (`argc` and `argv`) and populates the `args`
 * and `opts` structures with the parsed values. It supports extracting options and their associated values.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @param args Structure to store parsed arguments.
 * @param opts Structure to store parsed options.
 */
void parse_opts(int argc, char* argv[], TArgs& args, TOpts& opts);

/**
 * @brief Checks the length of a command-line option and updates the waiting state for additional options.
 *
 * This function validates whether a given option has the correct length and updates the `wait_for_opts`
 * flag to indicate whether the parser should expect additional values for the current option.
 *
 * @param opt The command-line option to check (e.g., "--option").
 * @param wait_for_opts Reference to a flag indicating whether the parser is waiting for additional options.
 * @param opts Structure containing parsed options.
 */
void check_opt_length(const std::string& opt, bool& wait_for_opts, TOpts& opts);

} // namespace launch
} // namespace skies
