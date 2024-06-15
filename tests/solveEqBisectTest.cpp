#include <gtest/gtest.h>

#include <cmath>

#include <skies/common/alg.h>

TEST(bisect, quad) {
    auto f = [] (double x) {
        return x*x - 9;
    };
    EXPECT_TRUE(std::abs(skies::find_root_bisect(f, 0, 10, 1e-10) - 3.0) <= 1e-5);
}

TEST(bisect, sin) {
    auto f = [] (double x) {
        return std::sin(x);
    };
    EXPECT_TRUE(std::abs(skies::find_root_bisect(f, -1.0, 1.0, 1e-10) - 0.0) <= 1e-5);
}