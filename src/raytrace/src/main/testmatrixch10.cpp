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

//------------------------------------------------------------------------------
// Scenario: A stripe pattern is constant in y.
TEST(Ch10Patterns, AStripePatternIsConstantInY)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 1.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 2.f, 0.f)) == White, true);
}

//------------------------------------------------------------------------------
// Scenario: A stripe pattern is constant in z.
TEST(Ch10Patterns, AStripePatternIsConstantInZ)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 1.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 2.f)) == White, true);
}

//------------------------------------------------------------------------------
// Scenario: A stripe pattern alternates in x.
TEST(Ch10Patterns, AStripePatternAlternatesInX)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.9f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(1.f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(-0.1f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(-1.f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(-1.1f, 0.f, 0.f)) == White, true);
}
