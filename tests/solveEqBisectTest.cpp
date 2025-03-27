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