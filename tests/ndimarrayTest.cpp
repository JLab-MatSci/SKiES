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

#include <skies/common/ndimarrays.h>

using namespace skies::arrays;

TEST(arraysND, binaryOps) {
    array1D v = { 1.0, -1.0, 2.5 };
    array1D w = { -1.0, 1.0, -2.5 };
    array1D s = v + w;
    array1D z = { 0.0, 0.0, 0.0 };
    EXPECT_EQ(s, z);
}

TEST(arraysND, flatten) {
    array2D m = {{1, 2, 3},
                 {4, 5, 6},
                 {7, 8, 9}};
    array1D f = flatten(m);
    EXPECT_EQ(f.size(), 9);
    EXPECT_EQ(f[0], 1);
    EXPECT_EQ(f[4], 5);
    EXPECT_EQ(f[8], 9);

    array3D mm = {{ {1, 2, 3},
                    {4, 5, 6} },
                  { {7, 8, 9},
                    {10, 11, 12} }};
    array1D ff = flatten(mm);
    EXPECT_EQ(ff.size(), 12);
    EXPECT_EQ(ff[0], 1);
    EXPECT_EQ(ff[4], 5);
    EXPECT_EQ(ff[7], 8);
}

TEST(arraysND, reshape) {
    array1D f = {1, 2, 3, 4, 5, 6, 7, 8};
    array2D m = reshape(f, 2, 4);
    EXPECT_EQ(m[0][0], 1);
    EXPECT_EQ(m[0][2], 3);
    EXPECT_EQ(m[1][3], 8);

    array1D ff = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    array3D mm = reshape(ff, 2, 2, 3);
    EXPECT_EQ(mm[0][0][0], 1);
    EXPECT_EQ(mm[0][0][2], 3);
    EXPECT_EQ(mm[0][1][1], 5);
    EXPECT_EQ(mm[1][0][2], 9);
    EXPECT_EQ(mm[1][1][1], 11);
}