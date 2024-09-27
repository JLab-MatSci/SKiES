/**
 @file
 @brief Some useful utils used in the project
 @author Galtsov Ilya
 */
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
