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
#include <string>
#include <stdexcept>

#include <gtest/gtest.h>

#include <skies/common/alg.h>

TEST(customSplitBlank, common) {
    std::string to_split = "\t  \t\t Hello, I am \t\t string    to be sepa rated!";
    auto separated = skies::custom_split(to_split, ' ');
    EXPECT_EQ(separated.size(), 8);
    EXPECT_EQ(separated[0], "Hello,");
    EXPECT_EQ(separated[1], "I");
    EXPECT_EQ(separated[2], "am");
    EXPECT_EQ(separated[3], "string");
    EXPECT_EQ(separated[4], "to");
    EXPECT_EQ(separated[5], "be");
    EXPECT_EQ(separated[6], "sepa");
    EXPECT_EQ(separated[7], "rated!");
}

TEST(customSplitX, common) {
    std::string to_split = "67x67x67 67x";
    auto separated = skies::custom_split(to_split, 'x');
    EXPECT_EQ(separated.size(), 3);
    EXPECT_EQ(separated[2], "67 67");
}

TEST(customSplitXX, common) {
    std::string to_split = "xxx67xx67x67 67x";
    auto separated = skies::custom_split(to_split, "xx");
    EXPECT_EQ(separated.size(), 2);
    EXPECT_EQ(separated[0], "x67");
    EXPECT_EQ(separated[1], "67x67 67x");
}

TEST(customSplitExcept, common) {
    std::string to_split = "I am string!";

    std::string empty;
    try {
        skies::custom_split(to_split, empty);
    } catch (std::runtime_error& err) {
        EXPECT_STREQ(err.what(), "A null-separator given!\n");
    }

    std::string wrong = "xxx xxx";
    try {
        skies::custom_split(to_split, wrong);
    } catch (std::runtime_error& err) {
        EXPECT_STREQ(err.what(), "A separator must not contain white spaces!\n");
    }
}