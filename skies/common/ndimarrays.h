/**
 @file
 @brief Interface of vector-matrix operations used elsewhere in the project
 @author Galtsov Ilya
 */
#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>

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

} // arrays
} // skies
