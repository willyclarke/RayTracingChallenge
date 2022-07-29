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

//------------------------------------------------------------------------------
// Scenario: The normal on the surface of a cube.
TEST(Ch12Cubes, TheNormalOnTheSurfaceOfACube)
{
  ww::shared_ptr_cube C = ww::PtrDefaultCube();

  auto const CheckLocalNormalAt = [&](ww::tup const &Point, ww::tup const &NormalExpected)
  { EXPECT_EQ(ww::LocalNormalAtCube(*C, Point) == NormalExpected, true); };

  CheckLocalNormalAt(ww::Point(1.f, 0.5f, -0.8f), ww::Vector(1.f, 0.f, 0.f));
  CheckLocalNormalAt(ww::Point(-1.f, -0.2f, 0.9f), ww::Vector(-1.f, 0.f, 0.f));
  CheckLocalNormalAt(ww::Point(-0.4f, 1.0f, 0.1f), ww::Vector(0.f, 1.f, 0.f));
  CheckLocalNormalAt(ww::Point(0.3f, -1.0f, -0.7f), ww::Vector(0.f, -1.f, 0.f));
  CheckLocalNormalAt(ww::Point(-0.6f, 0.3f, 1.f), ww::Vector(0.f, 0.f, 1.f));
  CheckLocalNormalAt(ww::Point(0.4f, 0.4f, -1.f), ww::Vector(0.f, 0.f, -1.f));
  CheckLocalNormalAt(ww::Point(1.f, 1.0f, 1.f), ww::Vector(1.f, 0.f, 0.f));
  CheckLocalNormalAt(ww::Point(-1.f, -1.0f, -1.f), ww::Vector(-1.f, 0.f, 0.f));
}

//------------------------------------------------------------------------------
// Scenario: Putting it together.
TEST(DISABLED_Ch12Cubes, PuttingItTogetherTake1)
{
  ww::world W = ww::World();
  W.vPtrLights.clear();
  W.vPtrObjects.clear();

  ww::shared_ptr_cube C = ww::PtrDefaultCube();
  C->Transform = ww::Scaling(1.f, 1.f, 1.f);
  W.vPtrObjects.push_back(C);

  ww::shared_ptr_cube C1 = ww::PtrDefaultCube();
  C1->Transform = ww::TranslateScaleRotate(3.f, 1.f, 0.f, 1.f, 4.f, 4.f, 0.f, 0.f, 0.f);
  W.vPtrObjects.push_back(C1);

  ww::shared_ptr_cube C2 = ww::PtrDefaultCube();
  C2->Transform = ww::TranslateScaleRotate(0.f, -2.f, 0.f, 2.f, 1.f, 4.f, 0.f, 0.f, 0.f);
  C2->Material.Color = ww::Color(0.99f, 0.f, 0.f);
  W.vPtrObjects.push_back(C2);

  ww::shared_ptr_cube C3 = ww::PtrDefaultCube();
  C3->Transform = ww::TranslateScaleRotate(0.f, 1.f, 5.f, 2.f, 4.f, 1.f, 0.f, 0.f, 0.f);
  C3->Material.Color = ww::Color(0.f, 0.99f, 0.f);
  W.vPtrObjects.push_back(C3);

  ww::shared_ptr_cube C4 = ww::PtrDefaultCube();
  C4->Transform = ww::TranslateScaleRotate(-3.f, 1.f, 0.f, 1.f, 4.f, 4.f, 0.f, 0.f, 0.f);
  W.vPtrObjects.push_back(C4);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-5.f, 25.f, -5.f), ww::Color(1.f, 1.f, 1.f));
  W.vPtrLights.push_back(pLight);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 5.0f, -10.f);
  ww::tup const ViewTo = ww::Point(0.f, 0.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, W);
  ww::WriteToPPM(Canvas, "Ch12PuttingItTogether.ppm");
}

//------------------------------------------------------------------------------
// Scenario: Putting it together.
TEST(Ch12Cubes, PuttingItTogetherTake2)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  float constexpr TWidth = 1.5f;
  float constexpr TDepth = 1.5f;
  float constexpr TThick = 0.05f;

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

  ww::shared_ptr_cube Table = ww::PtrDefaultCube();
  Table->Transform = ww::TranslateScaleRotate(0.f, 1.f, 0.f, TWidth, TThick, TDepth, 0.f, 0.f, 0.f);
  Table->Material.Pattern = ww::StripePattern(ColorBrown, ColorBrownLight);
  Table->Material.Pattern.Transform =
      ww::TranslateScaleRotate(0.f, 0.f, 0.f, .1f, 1.f, 1.f, 0.f, ww::Radians(45.f), 0.f);
  Table->Material.Ambient = 0.7f;
  World.vPtrObjects.push_back(Table);

  ww::shared_ptr_cube Leg1 = ww::PtrDefaultCube();
  Leg1->Transform = ww::TranslateScaleRotate(-TWidth + 0.1f, 0.f, TDepth - 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg1->Material.Ambient = 2.0f;
  Leg1->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg1);

  ww::shared_ptr_cube Leg2 = ww::PtrDefaultCube();
  Leg2->Transform = ww::TranslateScaleRotate(-TWidth + 0.1f, 0.f, -TDepth + 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg2->Material.Ambient = 2.0f;
  Leg2->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg2);

  ww::shared_ptr_cube Leg3 = ww::PtrDefaultCube();
  Leg3->Transform = ww::TranslateScaleRotate(TWidth - 0.1f, 0.f, -TDepth + 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg3->Material.Ambient = 2.0f;
  Leg3->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg3);

  ww::shared_ptr_cube Leg4 = ww::PtrDefaultCube();
  Leg4->Transform = ww::TranslateScaleRotate(TWidth - 0.1f, 0.f, TDepth - 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg4->Material.Ambient = 2.0f;
  Leg4->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg4);

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
  ww::WriteToPPM(Canvas, "Ch12PuttingItTogether2.ppm");
}

//------------------------------------------------------------------------------
// Scenario: Test that single and multithread generates the same image.
TEST(DISABLED_Ch12Cubes, TestThatSingleAndMultithreadGeneratesTheSameImage)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  float constexpr TWidth = 1.5f;
  float constexpr TDepth = 1.5f;
  float constexpr TThick = 0.05f;

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

  ww::shared_ptr_cube Table = ww::PtrDefaultCube();
  Table->Transform = ww::TranslateScaleRotate(0.f, 1.f, 0.f, TWidth, TThick, TDepth, 0.f, 0.f, 0.f);
  Table->Material.Pattern = ww::StripePattern(ColorBrown, ColorBrownLight);
  Table->Material.Pattern.Transform =
      ww::TranslateScaleRotate(0.f, 0.f, 0.f, .1f, 1.f, 1.f, 0.f, ww::Radians(45.f), 0.f);
  Table->Material.Ambient = 0.7f;
  World.vPtrObjects.push_back(Table);

  ww::shared_ptr_cube Leg1 = ww::PtrDefaultCube();
  Leg1->Transform = ww::TranslateScaleRotate(-TWidth + 0.1f, 0.f, TDepth - 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg1->Material.Ambient = 2.0f;
  Leg1->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg1);

  ww::shared_ptr_cube Leg2 = ww::PtrDefaultCube();
  Leg2->Transform = ww::TranslateScaleRotate(-TWidth + 0.1f, 0.f, -TDepth + 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg2->Material.Ambient = 2.0f;
  Leg2->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg2);

  ww::shared_ptr_cube Leg3 = ww::PtrDefaultCube();
  Leg3->Transform = ww::TranslateScaleRotate(TWidth - 0.1f, 0.f, -TDepth + 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg3->Material.Ambient = 2.0f;
  Leg3->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg3);

  ww::shared_ptr_cube Leg4 = ww::PtrDefaultCube();
  Leg4->Transform = ww::TranslateScaleRotate(TWidth - 0.1f, 0.f, TDepth - 0.2f, .05f, 1.f, .1f, 0.f, 0.f, 0.f);
  Leg4->Material.Ambient = 2.0f;
  Leg4->Material.Color = ColorBrown;
  World.vPtrObjects.push_back(Leg4);

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

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-5.f, 25.f, -5.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256 >> 2, 256 >> 2, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 4.5f, -7.f);
  ww::tup const ViewTo = ww::Point(2.0f, 1.f, 1.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  Camera.RenderSingleThread = false;
  ww::canvas const CanvasMultiThread = ww::Render(Camera, World);

  Camera.RenderSingleThread = true;
  ww::canvas const CanvasSingleThread = ww::Render(Camera, World);

  EXPECT_EQ(CanvasMultiThread.vXY.size(), CanvasSingleThread.vXY.size());

  if (CanvasMultiThread.vXY.size() == CanvasSingleThread.vXY.size())
  {
    for (int Idx = 0;                         //!<
         Idx < CanvasMultiThread.vXY.size();  //!<
         ++Idx)
    {
      EXPECT_EQ(CanvasMultiThread.vXY[Idx] == CanvasSingleThread.vXY[Idx], true);
      ww::tup const Diff = CanvasMultiThread.vXY[Idx] - CanvasSingleThread.vXY[Idx];
    }
  }
  ww::WriteToPPM(CanvasMultiThread, "Ch12CanvasMultiThread.ppm");
  ww::WriteToPPM(CanvasSingleThread, "Ch12CanvasSingleThread.ppm");
}
