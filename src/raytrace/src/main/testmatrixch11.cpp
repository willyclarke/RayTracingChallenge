#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>
//------------------------------------------------------------------------------
// Scenario: Reflectivity for the default material
TEST(CH11ReflectionAndRefraction, ReflectivityForTheDefaultMaterial)
{
  ww::material const M = {};
  EXPECT_EQ(M.Reflective == 0.f, true);
}

//------------------------------------------------------------------------------
// Scenario: Precomputing the reflective vector
TEST(CH11ReflectionAndRefraction, PrecomputingTheReflectiveVector)
{
  ww::shared_ptr_plane const Shape = ww::PtrDefaultPlane();
  ww::ray const R = ww::Ray(ww::Point(0.f, 1.f, -1.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f));
  ww::intersection const I = ww::intersection{M_SQRT2, Shape};
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
  EXPECT_EQ(Comps.Reflect == ww::Vector(0.f, M_SQRT2 / 2.f, M_SQRT2 / 2.f), true);
}
