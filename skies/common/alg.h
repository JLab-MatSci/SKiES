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
#include <functional>

namespace skies {

/**
 * \brief Split a string by a separator
 * @param str string to be splitted
 * @param separator given separator symbol
*/
std::vector<std::string>
custom_split(const std::string& str, const char separator = ' ');

std::vector<std::string>
custom_split(const std::string& str, const std::string& separator);

double find_root_bisect(std::function<double(double)> f, double a, double b, double crit);

} // skies
