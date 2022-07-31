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

//------------------------------------------------------------------------------
// Scenario: A ray strikes  a cylinder.
TEST(Ch13Cylinders, ARayStrikesACylinder)
{
  ww::shared_ptr_cylinder const Cyl = ww::PtrDefaultCylinder();
  auto CheckIntersect = [](ww::shared_ptr_cylinder Cyl, ww::ray const &Ray, float t0, float t1)
  {
    ww::intersections const XS = ww::LocalIntersect(Cyl, Ray);
    EXPECT_EQ(XS.Count(), 2);
    if (XS.Count() == 2)
    {
      EXPECT_EQ(std::abs(XS.vI[0].t - t0) < ww::EPSILON, true);
      EXPECT_EQ(std::abs(XS.vI[1].t - t1) < ww::EPSILON, true);
    }
  };

  CheckIntersect(Cyl, ww::Ray(ww::Point(1.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f)), 5.f, 5.f);
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f)), 4.f, 6.f);
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.5f, 0.f, -5.f), ww::Vector(.1f, 1.f, 1.f)), 6.80798f, 7.08872f);
}
