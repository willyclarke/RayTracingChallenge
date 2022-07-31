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

//------------------------------------------------------------------------------
// Scenario: Normal vector on a cylinder.
TEST(Ch13Cylinders, NormalVectorOnACylinder)
{
  ww::shared_ptr_cylinder const Cyl = ww::PtrDefaultCylinder();
  EXPECT_EQ(ww::LocalNormalAtCylinder(*Cyl, ww::Point(1.f, 0.f, 0.f)) == ww::Vector(1.f, 0.f, 0.f), true);
  EXPECT_EQ(ww::LocalNormalAtCylinder(*Cyl, ww::Point(0.f, 5.f, -1.f)) == ww::Vector(0.f, 0.f, -1.f), true);
  EXPECT_EQ(ww::LocalNormalAtCylinder(*Cyl, ww::Point(0.f, -2.f, 1.f)) == ww::Vector(0.f, 0.f, 1.f), true);
  EXPECT_EQ(ww::LocalNormalAtCylinder(*Cyl, ww::Point(-1.f, 1.f, 0.f)) == ww::Vector(-1.f, 0.f, 0.f), true);
}

//------------------------------------------------------------------------------
// Scenario: Indefinite cylinders test for fun.
TEST(DISABLED_Ch13Cylinders, IndefiniteCylinders)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  ptrPlane->Material.Shininess = 10.f;
  ptrPlane->Material.Diffuse = 0.5f;
  // ptrPlane->Material.Reflective = 1.f;
  ptrPlane->Material.Pattern = ww::CheckersPattern(ww::Color(0.f, 0.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  ptrPlane->Transform = ww::RotateY(ww::Radians(45.f));
  ptrPlane->Material.Color = ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  World.vPtrObjects.push_back(ptrPlane);

  // Add a second plane - this will act as the backdrop.
  ww::shared_ptr_plane ptrPlane2 = ww::PtrDefaultPlane();
  ptrPlane2->Transform =
      ww::TranslateScaleRotate(0.f, 0.f, 10.f, 1.f, 1.f, 1.f, ww::Radians(90.f), 0.f, ww::Radians(45.f));
  ptrPlane2->Material.Shininess = 100.f;
  ptrPlane2->Material.Diffuse = 0.7f;
  ptrPlane2->Material.Color =
      ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  // ptrPlane2->Material.Pattern = ww::CheckersPattern(ww::Color(0.f, 0.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrObjects.push_back(ptrPlane2);

  // Add a third plane - this will act as the backdrop.
  ww::shared_ptr_plane ptrPlane3 = ww::PtrDefaultPlane();
  ptrPlane3->Transform =
      ww::TranslateScaleRotate(0.f, 0.f, 10.f, 1.f, 1.f, 1.f, ww::Radians(90.f), 0.f, ww::Radians(-45.f));
  ptrPlane3->Material.Shininess = 100.f;
  ptrPlane3->Material.Diffuse = 0.7f;
  ptrPlane3->Material.Ambient = 1.2f;
  ptrPlane3->Material.Pattern = ww::StripePattern(ColorBrown, ColorBrownLight);
  World.vPtrObjects.push_back(ptrPlane3);

  ww::shared_ptr_cube Mirror = ww::PtrDefaultCube();
  Mirror->Transform =
      ww::TranslateScaleRotate(6.f, 3.f, 6.f, 2.f, 2.f, 2.f, ww::Radians(0.f), ww::Radians(45.f), ww::Radians(0.f));
  Mirror->Material.Ambient = 2.0f;
  Mirror->Material.Color = ww::Color(.1f, .1f, .1f);
  Mirror->Material.Diffuse = 0.5f;
  Mirror->Material.Reflective = .7f;
  Mirror->Material.Specular = 1.0f;
  Mirror->Material.Shininess = 500.f;
  World.vPtrObjects.push_back(Mirror);

  ww::shared_ptr_cube CubeOnTable1 = ww::PtrDefaultCube();
  CubeOnTable1->Transform = ww::TranslateScaleRotate(-0.5f, 1.2f, 0.f, .2f, .2f, .2f, 0.f, ww::Radians(45.f), 0.f);
  CubeOnTable1->Material.Ambient = 0.4f;
  CubeOnTable1->Material.Color = ww::Color(.1f, .1f, .1f);
  CubeOnTable1->Material.Diffuse = 0.6f;
  CubeOnTable1->Material.Reflective = .3f;
  CubeOnTable1->Material.Specular = 1.0f;
  CubeOnTable1->Material.Shininess = 500.f;
  World.vPtrObjects.push_back(CubeOnTable1);

  ww::shared_ptr_cube CubeOnTable2 = ww::PtrDefaultCube();
  CubeOnTable2->Transform = ww::TranslateScaleRotate(0.5f, 1.2f, 0.f, .5f, .4f, .2f, 0.f, 0.f, 0.f);
  CubeOnTable2->Material.Ambient = 0.4f;
  CubeOnTable2->Material.Color = ww::Color(0.f, 0.f, 1.f);
  World.vPtrObjects.push_back(CubeOnTable2);

  // ---
  // Add Sphere
  // ---
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(1.f, 1.f, .0f);
  float const Alpha = M_PI * float(45.f / 180.f);
  ww::pattern P1 = ww::StripePattern(Blue, Yellow);
  P1.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, .1f, 1.f, 1.f, -Alpha, Alpha, 0.f);
  ww::pattern P1Perturbed = ww::PerturbPattern(P1);
  ww::shared_ptr_shape PtrMiddle = ww::PtrDefaultSphere();
  ww::shape &Middle = *PtrMiddle;
  Middle.Transform = ww::Translation(.5f, 1.f, 3.5f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  Middle.Material.Pattern = P1Perturbed;
  Middle.Material.Pattern.Transform = ww::Scaling(0.15f, 0.5f, 0.5f) * ww::RotateZ(0.78);
  World.vPtrObjects.push_back(PtrMiddle);

  // ---
  // Add a cylinder.
  // ---
  ww::shared_ptr_shape PtrCylinder1 = ww::PtrDefaultCylinder();
  World.vPtrObjects.push_back(PtrCylinder1);

  // ---
  // Add a cylinder.
  // ---
  ww::shared_ptr_shape PtrCylinder2 = ww::PtrDefaultCylinder();
  PtrCylinder2->Transform = ww::TranslateScaleRotate(2.f, 0.f, 3.f, 0.3f, 0.3f, 0.3f, 0.f, 0.f, ww::Radians(0.f));
  PtrCylinder2->Material.Color = Blue;
  World.vPtrObjects.push_back(PtrCylinder2);

  // ---
  // Add a cylinder.
  // ---
  ww::shared_ptr_shape PtrCylinder3 = ww::PtrDefaultCylinder();
  PtrCylinder3->Material.Color = Yellow;
  PtrCylinder3->Transform = ww::TranslateScaleRotate(-2.f, 0.f, -3.f, 0.3f, 0.3f, 0.3f, 0.f, 0.f, ww::Radians(-45.f));
  ww::pattern P2 = ww::StripePattern(Blue, Yellow);
  P2.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, .1f, 1.f, 1.f, -Alpha, Alpha, 0.f);
  ww::pattern P2Perturbed = ww::PerturbPattern(P2);
  PtrCylinder3->Material.Pattern = P2Perturbed;
  World.vPtrObjects.push_back(PtrCylinder3);

  // ---
  // We need some light, please.
  // ---
  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-5.f, 25.f, -5.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 4.5f, -7.f);
  ww::tup const ViewTo = ww::Point(2.0f, 1.f, 1.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch13IndefiniteCylinders.ppm");
}

//------------------------------------------------------------------------------
// Scenario: The default maximum and minimum for a cylinder.
TEST(Ch13Cylinders, TheDefaultMaximumAndMinimumForACylinder)
{
  ww::shared_ptr_cylinder const Cyl = ww::PtrDefaultCylinder();
  EXPECT_EQ(Cyl->Minimum == -ww::INIFINITY, true);
  EXPECT_EQ(Cyl->Maximum == ww::INIFINITY, true);
}

//------------------------------------------------------------------------------
// Scenario: Intersecting a constrained cylinder.
TEST(Ch13Cylinders, IntersectingAConstrainedCylinder)
{
  ww::shared_ptr_cylinder Cyl = ww::PtrDefaultCylinder();
  Cyl->Minimum = 1.f;
  Cyl->Maximum = 2.f;

  auto CheckIntersect = [](ww::shared_ptr_cylinder Cyl, ww::ray const &Ray, int NumIntersect)
  {
    ww::intersections const XS = ww::LocalIntersect(Cyl, Ray);
    EXPECT_EQ(XS.Count(), NumIntersect);
  };

  // ---
  // Cast a ray diagonally from inside the cylinder with the ray escaping without intersecting.
  // ---
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 1.5f, 0.f), ww::Normalize(ww::Vector(0.1f, 1.f, 0.f))), 0);

  // ---
  // Cast rays perpendicular to the y-axis, but from above and below the cylinder, should miss.
  // ---
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 3.f, 5.f), ww::Normalize(ww::Vector(0.f, 0.f, 1.f))), 0);
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Normalize(ww::Vector(0.f, 0.f, 1.f))), 0);

  // ---
  // Edge cases, shows that the minimum and maximum of the y values themselves are outside of the
  // bounds of the cylinder. Remember that the height of the cylinder is from +1 to +2.
  // ---
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 2.f, -5.f), ww::Normalize(ww::Vector(0.f, 0.f, 1.f))), 0);
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 1.f, -5.f), ww::Normalize(ww::Vector(0.f, 0.f, 1.f))), 0);

  // ---
  // The final example cast a ray perpendicularly through the middle of the cylinder and produces two
  // intersections.
  // ---
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 1.5f, -2.f), ww::Normalize(ww::Vector(0.f, 0.f, 1.f))), 2);
}

//------------------------------------------------------------------------------
// Scenario: The default closed value for a cylinder.
TEST(Ch13Cylinders, TheDefaultClosedValueForACylinder)
{
  ww::shared_ptr_cylinder Cyl = ww::PtrDefaultCylinder();
  EXPECT_EQ(Cyl->Closed, false);
}

//------------------------------------------------------------------------------
// Scenario: Intersecting the caps of a closed cylinder.
TEST(Ch13Cylinders, IntersectingTheCapsOfAClosedCylinder)
{
  ww::shared_ptr_cylinder Cyl = ww::PtrDefaultCylinder();
  Cyl->Minimum = 1.f;
  Cyl->Maximum = 2.f;
  Cyl->Closed = true;

  auto CheckIntersect = [](ww::shared_ptr_cylinder Cyl, ww::ray const &Ray, int NumIntersect)
  {
    ww::intersections const XS = ww::LocalIntersect(Cyl, Ray);
    EXPECT_EQ(XS.Count(), NumIntersect);
  };

  // ---
  // NOTE: Start above the cylinder and point a ray straight through downwards. Test #1.
  // ---
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 3.f, 0.f), ww::Normalize(ww::Vector(0.f, -1.f, 0.f))), 2);

  // ---
  // NOTE: Start above and below the cylinder and cast the ray diagonally. Tests #2 and #4.
  // ---
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 3.f, -2.f), ww::Normalize(ww::Vector(0.f, -1.f, 2.f))), 2);
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 0.f, -2.f), ww::Normalize(ww::Vector(0.f, 1.f, 2.f))), 2);

  // ---
  // NOTE: Check the corner cases, tests #3 and #5 in the book.
  // ---
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, 4.f, -2.f), ww::Normalize(ww::Vector(0.f, -1.f, 1.f))), 2);
  CheckIntersect(Cyl, ww::Ray(ww::Point(0.f, -1.f, -2.f), ww::Normalize(ww::Vector(0.f, 1.f, 1.f))), 2);
}
