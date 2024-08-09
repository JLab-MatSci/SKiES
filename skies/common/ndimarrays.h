/**
 @file
 @brief Interface of vector-matrix operations used elsewhere in the project
 @author Galtsov Ilya
 */
#pragma once

#include <vector>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include <skies/common/alg.h>

namespace skies { namespace arrays {

using array1D = std::vector<double>;
using array2D = std::vector<array1D>;
using array3D = std::vector<array2D>;
using array4D = std::vector<array3D>;
using array5D = std::vector<array4D>;

/**
 * \brief Transposes a matrix given as a vector of vectors
 * @param w matrix to transpose
 * @result another matrix equal transposed of the given
*/
template <typename T>
std::vector<std::vector<T>>
transpose(const std::vector<std::vector<T>>& w)
{
    std::vector<std::vector<T>> tr;
    auto vsize = w.size();
    auto hsize = w[0].size();
    tr.resize(hsize,
              std::vector<T>(vsize));

    for (size_t i = 0; i < vsize; ++i)
        for (size_t j = 0; j < hsize; ++j)
            tr[j][i] = w[i][j];

    return tr;
}

/**
 * \brief Adds two vectors
 * @param v first vector
 * @param w second vector
 * @result another vector equal the result of the summation
*/
template <typename T>
std::vector<T> operator+(const std::vector<T>& v, 
                         const std::vector<T>& w)
    {
        if (v.size() != w.size())
            throw std::runtime_error("Arrays must have same sizes.");
        auto v_iter = v.begin();
        auto w_iter = w.begin();
        size_t i = 0;
        std::vector<T> sum;
        sum.resize(v.size());
        for (; v_iter != v.end();) {
            sum[i] = *v_iter + *w_iter;
            ++v_iter; ++w_iter;
            ++i;
        }
        return sum;
    }

/**
 * \brief Subtracts two vectors
 * @param v first vector
 * @param w second vector
 * @result another vector equal the result of the subtraction
*/
template <typename T>
std::vector<T> operator-(const std::vector<T>& v, 
                         const std::vector<T>& w)
    {
        if (v.size() != w.size())
            throw std::runtime_error("Arrays must have same sizes.");
        auto v_iter = v.cbegin();
        auto w_iter = w.cbegin();
        size_t i = 0;
        std::vector<T> res;
        res.resize(v.size());
        for (;v_iter != v.cend();) {
            res[i] = *v_iter - *w_iter;
            ++v_iter; ++w_iter;
            ++i;
        }
        return res;
    }

/**
 * \brief Evaluates the Euclidian norm of a vector
 * @param v vector to find norm
 * @result norm of the vector
*/
template <typename T>
T find_norm(const std::vector<T>& v) {
    double norm = 0;
    for (auto x : v) {
        norm += x*x;
    }
    return norm;
}

/**
 * \brief Multiplies a vector by a number
 * @param v vector
 * @param val number multiply to
 * @result another vector equal v * val
*/
template <typename T>
std::vector<T> operator* (const std::vector<T>& v, T val) {
    std::vector<T> res;
    for (auto&& x : v)
        res.push_back(x * val);
    return res;
}

/**
 * \brief The same as above
*/
template <typename T>
std::vector<T> operator* (T val, const std::vector<T>& v) {
    return v * val;
}

/**
 * \brief Multiplies a 2D-vector by a number
 * @param v 2D-vector
 * @param val number multiply to
 * @result another 2D-vector equal v * val
*/
template <typename T>
std::vector<std::vector<T>> operator* (const std::vector<std::vector<T>>& v, T val) {
    std::vector<std::vector<T>> res;
    for (auto&& l : v)
        res.push_back(l * val);
    return res;
}

/**
 * \brief The same as above
*/
template <typename T>
std::vector<std::vector<T>> operator* (T val, const std::vector<std::vector<T>>& v) {
    return v * val;
}

/**
 * \brief Evaluates the product of a matrix and a vector
 * matrix M is row-major\n
 * \t\t| v0 |\t{ {v_00, v_01, v_02},\n
 * M = | v1 | =\t{v_10, v_11, v_12},\n
 * \t\t| v2 |\t{v_20, v_21, v_22} }\n
 * @param M matrix
 * @param v vector to multiply by
 * @result another vector equal the result of the product
*/
template <typename T>
std::vector<T> matmul(const std::vector<std::vector<T>>& M,
                      const std::vector<T>& v)
{
    if (M.size() != v.size())
        throw std::runtime_error("Number of columns in matrix is not equal to size of vector.");
    std::vector<T> res;
    size_t N = v.size();
    for (auto&& l : M)
    {
        T sum{ 0 };
        for (size_t i = 0; i < N; ++i)
            sum += l[i] * v[i];
        res.push_back(sum);
    }
    return res;
}

/**
 * \brief Constructs a range
 * @param begin start point
 * @param end end point
 * @param nbins amount of bins in the range
 * @result range with nbins points starting at begin inclusive and ending at end inclusive
*/
inline array1D create_range(double begin, double end, int nbins)
{
    array1D range(nbins + 2);
    double d = (end - begin) / (nbins + 1.0);
    std::generate(range.begin(), range.end(), [d, pos = begin - d] () mutable { pos += d; return pos; });
    return range;
}

/**
 * \brief Makes a 1D-version of a 2D array
 * @param mat2d given 2D-array
 * @result new 1D vector which contains flattened elements
*/
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& mat2d)
{
    std::vector<T> flattened;
    auto n1 = mat2d.size();
    auto n2 = mat2d[0].size();
    flattened.resize(n1 * n2);
    for (size_t i1 = 0; i1 < n1; ++i1)
    {
        for (size_t i2 = 0; i2 < n2; ++i2)
        {
            flattened[i1 * n2 + i2] = mat2d[i1][i2];
        }
    }
    return flattened;
}

/**
 * \brief Makes a 1D-version of a 3D array
 * @param mat3d given 3D-array
 * @result new 1D vector which contains flattened elements
*/
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<std::vector<T>>>& mat3d)
{
    std::vector<T> flattened;
    auto n1 = mat3d.size();
    auto n2 = mat3d[0].size();
    auto n3 = mat3d[0][0].size();
    auto n23 = n2 * n3;
    flattened.resize(n1 * n2 * n3, 0);
    for (size_t i1 = 0; i1 < n1; ++i1)
    {
        for (size_t i2 = 0; i2 < n2; ++i2)
        {
            for (size_t i3 = 0; i3 < n3; ++i3)
            {
                flattened[i1 * n23 + i2 * n3 + i3] = mat3d[i1][i2][i3];
            }
        }
    }
    return flattened;
}

/**
 * \brief Creates 2D-array with given dimensions from a given 1D-array
 * @param v given 1D-array
 * @param n1 number of rows
 * @param n2 number of columns
 * @result new reshaped 2D array
*/
template <typename T>
std::vector<std::vector<T>>
reshape(const std::vector<T>& v, size_t n1, size_t n2)
{
    assert(v.size() == n1 * n2);
    std::vector<std::vector<T>> mat2d;
    mat2d.resize(n1, std::vector<T>(n2, 0));
    for (size_t i1 = 0; i1 < n1; ++i1)
    {
        std::vector<T> tmp(n2);
        for (size_t i2 = 0; i2 < n2; ++i2)
        {
            tmp[i2] = v[i2 + i1 * n2];
        }
        mat2d[i1] = tmp;
    }
    return mat2d;
}

/**
 * \brief Creates 3D-array with given dimensions from a given 1D-array
 * @param v given 1D-array
 * @param n1 number of layers
 * @param n2 number of rows in one layer
 * @param n3 number of columns in one layer
 * @result new reshaped 3D array
*/
template <typename T>
std::vector<std::vector<std::vector<T>>>
reshape(const std::vector<T>& v, size_t n1, size_t n2, size_t n3)
{
    assert(v.size() == n1 * n2 * n3);
    std::vector<std::vector<std::vector<T>>> mat3d;
    mat3d.resize(n1, std::vector<std::vector<T>>{n2, std::vector<T>(n3, 0.0)});
    auto n23 = n2 * n3;
    for (size_t i1 = 0; i1 < n1; ++i1)
    {
        auto start = i1 * n23;
        auto finish = start + n23;
        std::vector<T> v_part(v.begin() + start, v.begin() + finish);
        assert(v_part.size() == n23);
        mat3d[i1] = reshape<T>(v_part, n2, n3);
    }
    return mat3d;
}

inline void dump_array(const std::string& filename, const array2D& array, const char sep = ' ')
{
    std::ofstream ofs{filename};
    for (auto&& line : array)
    {
        for (auto&& v : line)
        {
            ofs << v << sep;
        }
        ofs << std::endl;
    }
    ofs.close();
}

inline void load_array(const std::string& filename, array2D& array, const char sep = ' ')
{
    assert(array.empty());
    std::ifstream ifs{filename};
    std::string line_str;

    auto parse_line = [&sep] (const std::string& line_str, array1D& line_arr) {
        for (auto&& v : custom_split(line_str, sep))
            line_arr.push_back(std::stod(v));
    };

    while (ifs.good())
    {
        std::getline(ifs, line_str);
        if (!line_str.empty())
        {
            array1D line_arr;
            parse_line(line_str, line_arr);
            array.push_back(line_arr);
        }
    }
    ifs.close();
}

} // arrays
} // skies
