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

//------------------------------------------------------------------------------
// Scenario: color_at() with mutually reflective surface
TEST(CH11ReflectionAndRefraction, ColorAtWithMutuallyReflectiveSurface)
{
  ww::world W{};

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(0.f, 0.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  ww::WorldAddLight(W, pLight);

  // ---
  // Add the first plane
  // ---
  ww::shared_ptr_plane Lower = ww::PtrDefaultPlane();
  Lower->Material.Reflective = 1.0f;
  Lower->Transform = ww::Translation(0.f, -1.f, 0.f);
  // Lower->Print = true;
  ww::WorldAddObject(W, Lower);

  // ---
  // Add the second plane
  // ---
  ww::shared_ptr_plane Upper = ww::PtrDefaultPlane();
  Upper->Material.Reflective = 1.0f;
  Upper->Transform = ww::Translation(0.f, 1.f, 0.f);
  // Upper->Print = true;
  ww::WorldAddObject(W, Upper);

  EXPECT_EQ(W.Count(), 2);

  // W.Print = true;
  W.ColorAtCallCnt = 0;

  // ---
  // NOTE: Check that a ray pointing to the lower plane intersects.
  // ---
  {
    ww::ray const R{ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, -1.f, 0.f)};
    ww::intersections const XS = ww::LocalIntersect(Lower, R);

    EXPECT_EQ(XS.Count(), 1);
    if (XS.Count())
    {
      EXPECT_EQ(XS.vI[0].t, 1.f);
      EXPECT_EQ(XS.vI[0].pShape == Lower, true);
    }
  }

  // ---
  // NOTE: Check that a ray pointing to the upper plane intersects.
  // ---
  {
    ww::ray const R{ww::Point(0.f, -1.f, 0.f), ww::Vector(0.f, 1.f, 0.f)};
    ww::intersections const XS = ww::LocalIntersect(Upper, R);

    EXPECT_EQ(XS.Count(), 1);
    if (XS.Count())
    {
      EXPECT_EQ(XS.vI[0].t, 1.f);
      EXPECT_EQ(XS.vI[0].pShape == Upper, true);
    }
  }

  std::cout << __FUNCTION__ << ". Line: " << __LINE__ << ". \nUpper plane: " << Upper << "\nLower plane: " << Lower
            << "\n------------------------------------------------------------------------------" << std::endl;
  ww::ray const R{ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 1.f, 0.f)};
  ww::tup const Color = ww::ColorAt(W, R);
}
