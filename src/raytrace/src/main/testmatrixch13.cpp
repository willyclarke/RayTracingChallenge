#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>

//------------------------------------------------------------------------------
// Scenario: A ray misses a cylinder.
TEST(Ch13Cylinders, ARayMissesACylinder)
{
  ww::shared_ptr_cylinder const Cyl = ww::PtrDefaultCylinder();
  auto CheckIntersect = [](ww::shared_ptr_cylinder Cyl, ww::ray const &Ray)
  {
    ww::intersections const XS = ww::LocalIntersect(Cyl, Ray);
    EXPECT_EQ(XS.Count(), 0);
  };

  CheckIntersect(Cyl, ww::Ray(ww::Point(1.f, 0.f, 0.f), ww::Vector(0.f, 1.f, 0.f)));
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 1.f, 0.f)));
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(1.f, 1.f, 1.f)));
}
