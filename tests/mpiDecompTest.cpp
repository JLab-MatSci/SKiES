#include <gtest/gtest.h>

#include <skies/common/alg.h>

TEST(manyProcs, mpiDecomp) {
    int total = 10;
    int nchuncks = 3;
    auto rcounts_displs = skies::get_rcounts_displs(total, nchuncks);
    auto rcounts = rcounts_displs.first;
    auto displs  = rcounts_displs.second;
    EXPECT_EQ(rcounts[0], 4);
    EXPECT_EQ(rcounts[1], 3);
    EXPECT_EQ(rcounts[2], 3);
    EXPECT_EQ(displs[0],  0);
    EXPECT_EQ(displs[1],  4);
    EXPECT_EQ(displs[2],  7);
}

TEST(equal, mpiDecomp) {
    int total = 10;
    int nchuncks = 10;
    auto rcounts_displs = skies::get_rcounts_displs(total, nchuncks);
    auto rcounts = rcounts_displs.first;
    auto displs  = rcounts_displs.second;
    EXPECT_EQ(rcounts[0], 1);
    EXPECT_EQ(rcounts[9], 1);
    EXPECT_EQ(displs[0],  0);
    EXPECT_EQ(displs[9],  9);
}

TEST(fewProcs, mpiDecomp) {
    int total = 3;
    int nchuncks = 10;
    try {
        auto rcounts_displs = skies::get_rcounts_displs(total, nchuncks);
    } catch (std::runtime_error& err) {
        EXPECT_STREQ(err.what(), "Number of chuncks must be less or equal than total number of elements to decompose\n");
    }
}