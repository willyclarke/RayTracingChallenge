#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>
//------------------------------------------------------------------------------
// Scenario: Creating a stripe pattern.
TEST(Ch10Patterns, CreatingStripePattern)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(Pattern.A == White, true);
  EXPECT_EQ(Pattern.B == Black, true);
}

