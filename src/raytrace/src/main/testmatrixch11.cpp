#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>
//------------------------------------------------------------------------------
// Scenario: Reflectivity for the default material
TEST(Ch11Patterns, ReflectivityForTheDefaultMaterial)
{
  ww::material const M = {};
  EXPECT_EQ(M.Reflective == 0.f, true);
}
