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

//------------------------------------------------------------------------------
// Scenario: The reflected color for a nonreflective material
TEST(CH11ReflectionAndRefraction, TheReflectedColorForANonReflectiveMaterial)
{
  ww::world W = ww::World();
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 0.f, 1.f));

  if (W.vPtrObjects.size() > 1)
  {
    ww::shared_ptr_shape Shape = W.vPtrObjects[1];
    Shape->Material.Ambient = 1.f;
    ww::intersection const I = ww::intersection{1.f, Shape};
    ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
    ww::tup const Color = ww::ReflectedColor(W, Comps);
    EXPECT_EQ(Color == ww::Color(0.f, 0.f, 0.f), true);
  }
}

//------------------------------------------------------------------------------
// Scenario: The reflected color for a reflective material
TEST(CH11ReflectionAndRefraction, TheReflectedColorForAReflectiveMaterial)
{
  ww::world W = ww::World();

  // ---
  // Add the first plane
  // ---
  ww::shared_ptr_plane Shape = ww::PtrDefaultPlane();
  Shape->Material.Reflective = 0.5f;
  Shape->Transform = ww::TranslateScaleRotate(0.f, -1.f, 0.f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  ww::WorldAddObject(W, Shape);

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -3.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f));
  ww::intersection const I = ww::intersection{M_SQRT2, Shape};
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
  ww::tup const Color = ww::ReflectedColor(W, Comps);
  EXPECT_EQ(Color == ww::Color(0.19032f, 0.2379f, 0.14274f), true);

  // ---
  // NOTE: Write the result to file.
  // ---
  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 50.f, 0.f), ww::Color(1.f, 1.f, 1.f) * 0.2f);
  W.vPtrLights.push_back(pLight);

  float const FieldOfView = 75.f / 180.f * M_PI;
  ww::camera Camera = ww::Camera(256, 256, FieldOfView);

  ww::tup const ViewFrom = ww::Point(0.f, 1.5f, -3.5f);
  ww::tup const ViewTo = ww::Point(0.f, 0.f, 25.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, W);
  ww::WriteToPPM(Canvas, "Ch11ReflectedColorReflectiveMaterial.ppm");
}

//------------------------------------------------------------------------------
// Scenario: shade_hit() with a reflective material
TEST(CH11ReflectionAndRefraction, ShadeHitWithAReflectiveMaterial)
{
  ww::world W = ww::World();

  // ---
  // Add the first plane
  // ---
  ww::shared_ptr_plane Shape = ww::PtrDefaultPlane();
  Shape->Material.Reflective = 0.5f;
  Shape->Transform = ww::TranslateScaleRotate(0.f, -1.f, 0.f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  ww::WorldAddObject(W, Shape);

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -3.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f));
  ww::intersection const I = ww::intersection{M_SQRT2, Shape};
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
  ww::tup const Color = ww::ShadeHit(W, Comps);
  EXPECT_EQ(Color == ww::Color(0.87677f, 0.92436f, 0.82918), true);

  // ---
  // NOTE: Write the result to file.
  // ---
  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 50.f, 0.f), ww::Color(1.f, 1.f, 1.f) * 0.2f);
  W.vPtrLights.push_back(pLight);

  float const FieldOfView = 75.f / 180.f * M_PI;
  ww::camera Camera = ww::Camera(256, 256, FieldOfView);

  ww::tup const ViewFrom = ww::Point(0.f, 1.5f, -3.5f);
  ww::tup const ViewTo = ww::Point(0.f, 0.f, 25.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, W);
  ww::WriteToPPM(Canvas, "Ch11ShadeHitWithAReflectiveMaterial.ppm");
}
