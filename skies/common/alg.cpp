#include <cmath>
#include <cassert>

#include <skies/common/alg.h>

#include <iostream>

namespace skies {

std::pair<std::vector<int>, std::vector<int>>
get_rcounts_displs(int total, int nchunks)
{
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
custom_split(const std::string& str, const char separator) {
    std::vector<std::string> strings;
    int startIndex = 0;
    int endIndex = 0;
    size_t i = 0;
    while (i < str.length() + 1) {
        // If we reached the end of the word or the end of the input
        int sep_count = 0;
        size_t j = i;
        while (str[j] == separator || j == str.length()) {
            endIndex = j;
            if (!sep_count) {
                std::string temp;
                temp.append(str, startIndex, endIndex - startIndex);
                strings.push_back(temp);
            }
            startIndex = endIndex + 1;
            sep_count++;
            j++;
        }
        i += sep_count + 1;
    }
    return strings;
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
