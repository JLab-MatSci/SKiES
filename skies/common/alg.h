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
 * \brief Divides the total number of elements into given number of chuncks. Used for MPI parallelization
 * @param total total number of elements
 * @param nchuncks number of chuncks to be divided into
*/
std::pair<std::vector<int>, std::vector<int>>
get_rcounts_displs(int total, int nchunks);

/**
 * \brief Split a string by a separator
 * @param str string to be splitted
 * @param separator given separator symbol
*/
std::vector<std::string>
custom_split(const std::string& str, const char separator);

double find_root_bisect(std::function<double(double)> f, double a, double b, double crit);

} // skies
