#include "gtest/gtest.h"
#include <datastructures.hpp>
#include <raymarch.hpp>

#include <cmath>
#include <memory>

//------------------------------------------------------------------------------
TEST(RayMarch, GetDistance)
{
  ww::shared_ptr_shape pDefaultSphere = ww::PtrDefaultSphere();
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -1.6f), ww::Vector(0.f, 0.f, 1.f));
  float const D = ww::rm::GetDistance(R.Origin, pDefaultSphere);
  EXPECT_FLOAT_EQ(D, 0.6f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, GetNormal)
{
  auto TestNormal = [](ww::ray const &R)
  {
    ww::shared_ptr_shape pDefaultSphere = ww::PtrDefaultSphere();
    ww::tup const N = ww::rm::GetNormal(R.Origin, pDefaultSphere);
    EXPECT_FLOAT_EQ(ww::Mag(N), 1.0f);
    EXPECT_LT(ww::Mag(N - R.Direction), 0.01f);
  };
  TestNormal(ww::Ray(ww::Point(0.f, 0.f, -1.f), ww::Vector(0.f, 0.f, -1.f)));
  TestNormal(ww::Ray(ww::Point(0.f, 0.f, 1.f), ww::Vector(0.f, 0.f, 1.f)));
  TestNormal(ww::Ray(ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, 1.f, 0.f)));
  TestNormal(ww::Ray(ww::Point(0.f, -1.f, 0.f), ww::Vector(0.f, -1.f, 0.f)));
  TestNormal(ww::Ray(ww::Point(1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f)));
  TestNormal(ww::Ray(ww::Point(-1.f, 0.f, 0.f), ww::Vector(-1.f, 0.f, 0.f)));
}

//------------------------------------------------------------------------------
TEST(RayMarch, RayMarch)
{
  auto TestDistance = [](ww::ray const &R, ww::shared_ptr_shape pDefaultSphere) -> float
  {
    float const Distance = ww::rm::RayMarch(R, pDefaultSphere);
    return Distance;
  };

  ww::shared_ptr_shape pDefaultSphere = ww::PtrDefaultSphere();
  // Hit the sphere straight on.
  ww::ray R = ww::Ray(ww::Point(0.f, 0.f, -2.f), ww::Vector(0.f, 0.f, 1.f));
  EXPECT_FLOAT_EQ(1.f, TestDistance(R, pDefaultSphere));

  // Hit the side of the sphere.
  // by moving sideways to X=-1 the distance is expected to be around value of 2.
  R.Origin = ww::Point(-1.f, 0.f, -2.f);
  EXPECT_LT(std::abs(2.f - TestDistance(R, pDefaultSphere)), 0.05);

  // or by moving sideways to X=1 the distance is expected to be around value of 2.
  R.Origin = ww::Point(1.f, 0.f, -2.f);
  EXPECT_LT(std::abs(2.f - TestDistance(R, pDefaultSphere)), 0.05);

  // or by moving up to Y=1 the distance is expected to be around value of 2.
  R.Origin = ww::Point(0.f, 1.f, -2.f);
  EXPECT_LT(std::abs(2.f - TestDistance(R, pDefaultSphere)), 0.05);

  // or by moving down to Y=-1 the distance is expected to be around value of 2.
  R.Origin = ww::Point(0.f, -1.f, -2.f);
  EXPECT_LT(std::abs(2.f - TestDistance(R, pDefaultSphere)), 0.05);

  // Miss the sphere.
  R.Origin = ww::Point(-2.f, 0.f, -2.f);
  EXPECT_GT(TestDistance(R, pDefaultSphere), 100.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, RayMarchMovedSphere)
{
  auto TestDistance = [](ww::ray const &R, ww::shared_ptr_shape pDefaultSphere) -> float
  {
    float const Distance = ww::rm::RayMarch(R, pDefaultSphere);
    return Distance;
  };

  ww::shared_ptr_shape pDefaultSphere = ww::PtrDefaultSphere();
  pDefaultSphere->Transform = ww::Translation(0.f, 0.f, 1.f);
  // Hit the sphere straight on.
  ww::ray R = ww::Ray(ww::Point(0.f, 0.f, -2.f), ww::Vector(0.f, 0.f, 1.f));
  EXPECT_FLOAT_EQ(2.f, TestDistance(R, pDefaultSphere));

  // Hit the side of the sphere.
  // by moving sideways to X=-1 the distance is expected to be around value of 3.
  R.Origin = ww::Point(-1.f, 0.f, -2.f);
  EXPECT_LT(std::abs(3.f - TestDistance(R, pDefaultSphere)), 0.05);

  // or by moving sideways to X=1 the distance is expected to be around value of 3.
  R.Origin = ww::Point(1.f, 0.f, -2.f);
  EXPECT_LT(std::abs(3.f - TestDistance(R, pDefaultSphere)), 0.05);

  // or by moving up to Y=1 the distance is expected to be around value of 3.
  R.Origin = ww::Point(0.f, 1.f, -2.f);
  EXPECT_LT(std::abs(3.f - TestDistance(R, pDefaultSphere)), 0.05);

  // or by moving down to Y=-1 the distance is expected to be around value of 3.
  R.Origin = ww::Point(0.f, -1.f, -2.f);
  EXPECT_LT(std::abs(3.f - TestDistance(R, pDefaultSphere)), 0.05);

  // Miss the sphere.
  R.Origin = ww::Point(-2.f, 0.f, -2.f);
  EXPECT_GT(TestDistance(R, pDefaultSphere), 100.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, MainImage)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  // ---
  // Add Sphere
  // ---
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(1.f, 1.f, .0f);

  ww::shared_ptr_shape pDefaultSphere = ww::PtrDefaultSphere();
  ww::shape &Middle = *pDefaultSphere;
  // Middle.Transform = ww::Translation(.5f, 1.f, 3.5f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(pDefaultSphere);

  ww::camera Camera{};

  for (int X = 0; X < 16; ++X)
    for (int Y = 0; Y < 16; ++Y) ww::rm::MainImage(Camera, World, X, Y, pDefaultSphere);
}

//------------------------------------------------------------------------------
// Scenario:
TEST(RayMarch, Test1)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  // ---
  // Add Sphere
  // ---
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(1.f, 1.f, .0f);

  ww::shared_ptr_shape pDefaultSphere = ww::PtrDefaultSphere();
  ww::shape &S = *pDefaultSphere;
  S.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.0f, 1.0f, 1.0f, 0.f, 0.f, 0.f);
  S.Material.Color = ww::Color(1.0f, 0.2f, 0.5f);
  S.Material.Diffuse = 0.7f;
  S.Material.Specular = 0.f;
  S.Material.Ambient = 0.f;
  World.vPtrObjects.push_back(pDefaultSphere);

  ww::shared_ptr_shape pDefaultSphere2 = ww::PtrDefaultSphere();
  ww::shape &S2 = *pDefaultSphere2;
  S2.Transform = ww::TranslateScaleRotate(-2.f, 0.f, 0.f, 2.f, 2.f, 2.f, 0.f, 0.f, 0.f);
  S2.Material.Color = ww::Color(0.0f, 0.2f, 0.5f) * 2.f;
  S2.Material.Diffuse = 0.0f;
  S2.Material.Specular = 0.0f;
  World.vPtrObjects.push_back(pDefaultSphere2);

  // ---
  // We need some light, please.
  // ---
  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  // *pLight = ww::PointLight(ww::Point(-5.f, 25.f, -5.f), ww::Color(1.f, 1.f, 1.f));
  *pLight = ww::PointLight(ww::Point(-3.f, 3.f,  -3.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(320, 320, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(-4.f, 2.0f, -7.f);
  ww::tup const ViewTo = ww::Point(0.0f, 0.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  std::cout << "vPtrObjects.size:" << World.Count() << std::endl;

  ww::canvas Canvas = ww::rm::Render(Camera, World);
  ww::WriteToPPM(Canvas, "RayMarchTest1.ppm");
}
