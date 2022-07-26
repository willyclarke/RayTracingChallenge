#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>

//------------------------------------------------------------------------------
// Scenario: A ray intersects a cube.
TEST(Ch12Cubes, ARayIntersectsACube)
{
  ww::shared_ptr_cube C = ww::PtrDefaultCube();

  auto const IntersectCube = [&](ww::ray const &Ray, float TMin, float TMax)
  {
    ww::intersections const XS = ww::LocalIntersect(C, Ray);

    EXPECT_EQ(XS.Count(), 2);
    if (XS.Count() > 1)
    {
      EXPECT_EQ(XS.vI[0].t, TMin);
      EXPECT_EQ(XS.vI[1].t, TMax);
    }
  };

  ww::ray const RPlusX = ww::Ray(ww::Point(5.f, .5f, 0.f), ww::Vector(-1.f, 0.f, 0.f));
  ww::ray const RMinuX = ww::Ray(ww::Point(-5.f, .5f, 0.f), ww::Vector(1.f, 0.f, 0.f));
  ww::ray const RPlusY = ww::Ray(ww::Point(.5f, 5.f, 0.f), ww::Vector(0.f, -1.f, 0.f));
  ww::ray const RMinuY = ww::Ray(ww::Point(.5f, -5.f, 0.f), ww::Vector(0.f, 1.f, 0.f));
  ww::ray const RPlusZ = ww::Ray(ww::Point(.5f, 0.f, 5.f), ww::Vector(0.f, 0.f, -1.f));
  ww::ray const RMinuZ = ww::Ray(ww::Point(.5f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::ray const RInside = ww::Ray(ww::Point(0.f, .5f, 0.f), ww::Vector(0.f, 0.f, 1.f));
  IntersectCube(RPlusX, 4.f, 6.f);
  IntersectCube(RMinuX, 4.f, 6.f);
  IntersectCube(RPlusY, 4.f, 6.f);
  IntersectCube(RMinuY, 4.f, 6.f);
  IntersectCube(RPlusZ, 4.f, 6.f);
  IntersectCube(RMinuZ, 4.f, 6.f);
  IntersectCube(RInside, -1.f, 1.f);
}

//------------------------------------------------------------------------------
// Scenario: A ray misses a cube.
TEST(Ch12Cubes, ARayMissesACube)
{
  ww::shared_ptr_cube C = ww::PtrDefaultCube();

  auto const MissesCube = [&](ww::ray const &Ray)
  {
    ww::intersections const XS = ww::LocalIntersect(C, Ray);
    EXPECT_EQ(XS.Count(), 0);
  };

  MissesCube(ww::Ray(ww::Point(-2.f, 0.f, 0.f), ww::Vector(0.2673f, 0.5345f, 0.8018f)));
  MissesCube(ww::Ray(ww::Point(0.f, -2.f, 0.f), ww::Vector(0.8018f, 0.2673f, 0.5345f)));
  MissesCube(ww::Ray(ww::Point(0.f, 0.f, -2.f), ww::Vector(0.5345f, 0.8018f, 0.2673f)));
  MissesCube(ww::Ray(ww::Point(2.f, 0.f, 2.f), ww::Vector(0.f, 0.f, -1.f)));
  MissesCube(ww::Ray(ww::Point(0.f, 2.f, 2.f), ww::Vector(0.f, -1.f, 0.f)));
  MissesCube(ww::Ray(ww::Point(2.f, 2.f, 0.f), ww::Vector(-1.f, 0.f, 0.f)));
}
