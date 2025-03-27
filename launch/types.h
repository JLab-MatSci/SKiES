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

#include <vector>
#include <string>
#include <unordered_map>

namespace skies { namespace launch {

/**
 * @brief Namespace containing type aliases for command-line argument handling.
 *
 * This anonymous namespace defines type aliases for storing and managing command-line arguments and options.
 */
namespace {
    /**
     * @brief Type alias for a vector of command-line arguments.
     *
     * Represents a list of command-line arguments as strings. Each element corresponds to a single argument.
     */
    using TArgs = std::vector<std::string>;

    /**
     * @brief Type alias for an unordered map of command-line options.
     *
     * Represents a collection of command-line options, where each option is stored as a key-value pair:
     * - Key: The option name (e.g., "--option").
     * - Value: The associated value (if any).
     */
    using TOpts = std::unordered_map<std::string, std::string>;
}

} // namespace launch
} // namespace skies