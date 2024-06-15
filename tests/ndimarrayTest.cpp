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
