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
#include <algorithm>

#include <gtest/gtest.h>

#include <skies/lattices/kp_protocol.h>

TEST(checkMesh, kpProtocol) {
    skies::KPprotocol p{3, 4, 5};
    auto mesh = p.mesh();
    auto n1 = std::get<0>(mesh);
    auto n2 = std::get<1>(mesh);
    auto n3 = std::get<2>(mesh);
    EXPECT_EQ(n1, 3);
    EXPECT_EQ(n2, 4);
    EXPECT_EQ(n3, 5);
}

TEST(checkGrids, kpProtocol) {
    skies::KPprotocol p{3, 4, 5};

    EXPECT_EQ(p.igrid()[0][0], -2);
    EXPECT_EQ(p.igrid()[0][1], -3);
    EXPECT_EQ(p.igrid()[0][2], -4);
    EXPECT_EQ(p.grid()[0][0],  -1.0 / 3.0);
    EXPECT_EQ(p.grid()[0][1],  -0.375);
    EXPECT_EQ(p.grid()[0][2],  -0.4);

    EXPECT_EQ(p.igrid()[2][0], -2);
    EXPECT_EQ(p.igrid()[2][1],  1);
    EXPECT_EQ(p.igrid()[2][2], -4);
    EXPECT_EQ(p.grid()[2][0], -1.0 / 3.0);
    EXPECT_EQ(p.grid()[2][1],  0.125);
    EXPECT_EQ(p.grid()[2][2],  -0.4);
}

TEST(checkRange, kpProtocol) {
    skies::KPprotocol p{3, 4, 5};

    EXPECT_EQ(p.find_ind_from_vi(p.igrid()[0]), 0);
    EXPECT_EQ(p.find_ind_from_vi(p.igrid()[1]), 1);
    EXPECT_EQ(p.find_ind_from_vi(p.igrid()[2]), 2);

    auto range = p.range();
    EXPECT_EQ(range[0], 0);
    EXPECT_EQ(range[1], 1);
    EXPECT_EQ(range[2], 2);

    EXPECT_EQ(range.end(), std::unique(range.begin(), range.end()));
}