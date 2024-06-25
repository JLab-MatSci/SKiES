#include <cmath>
#include <cassert>
#include <regex>

#include <skies/common/alg.h>

#include <iostream>

namespace skies {

std::pair<std::vector<int>, std::vector<int>>
get_rcounts_displs(int total, int nchunks)
{
    if (total < nchunks)
        throw std::runtime_error("Number of chuncks must be less or equal than total number of elements to decompose\n");

    int ave = total / nchunks;
    int res = total % nchunks;

    std::vector<int> rcounts(nchunks);
    std::vector<int> displs(nchunks);

    for (int k = 0; k < nchunks; ++k) {
        if (k < res)
            rcounts[k] = ave + 1;
        else
            rcounts[k] = ave;
        if (k == 0)
            displs[k] = 0;
        else
            displs[k] = displs[k-1] + rcounts[k-1];
    }

    auto output = std::make_pair(rcounts, displs);
    return output;
}

std::vector<std::string>
custom_split_regex(const std::string& str, const std::string& separator)
{
    std::vector<std::string> strings;
    std::regex rgx(separator);
    std::sregex_token_iterator iter(str.begin(), str.end(), rgx, -1);
    std::sregex_token_iterator end;
    if (*iter == "") ++iter;
    for (; iter != end; ++iter)
        strings.push_back(*iter);
    return strings;
}

std::vector<std::string>
custom_split(const std::string& str, const char separator)
{
    // if white space, any number is assumed
    if (separator == ' ')
        return custom_split_regex(str, "\\s+");
    return custom_split_regex(str, std::string(1, separator));
}

std::vector<std::string>
custom_split(const std::string& str, const std::string& separator)
{
    if (separator.empty())
        throw std::runtime_error("A null-separator given!\n");
    if (separator.find(" ") != std::string::npos)
        throw std::runtime_error("A separator must not contain white spaces!\n");
    return custom_split_regex(str, separator);
}

double find_root_bisect(std::function<double(double)> f, double a, double b, double crit)
{
    int Nit{ 0 };
    while ((b - a) / std::pow(2, Nit) > crit)
    {
        assert((f(a) <= 0) && (f(b) >= 0));
        Nit++;
        f((a + b) / 2) > 0 ? b = (a + b) / 2 : a = (a + b) / 2;
    }
    return (a + b) / 2;
}

} // skies
