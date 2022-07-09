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

  // std::cout << __FUNCTION__ << ". Line: " << __LINE__ << ". \nUpper plane: " << Upper << "\nLower plane: " << Lower
  //           << "\n------------------------------------------------------------------------------" << std::endl;

  // ---
  // NOTE: The test here checks that the recursion ends after a finite number of calls.
  //       If this test runs to finish it means success.
  // ---
  ww::ray const R{ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 1.f, 0.f)};
  ww::tup const Color = ww::ColorAt(W, R);

  // std::cout << __FUNCTION__ << ". Line: " << __LINE__ << ". \nUpper plane: " << Upper << "\nLower plane: " << Lower
  //           << "\n------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------
// Scenario: The reflected color at the maximum recursive depth.
TEST(CH11ReflectionAndRefraction, TheReflectedColorAtTheMaximumRecursiveDepth)
{
  ww::world W = ww::World();

  W.Print = true;

  // ---
  // Add the first plane
  // ---
  ww::shared_ptr_plane Shape = ww::PtrDefaultPlane();
  Shape->Material.Reflective = 0.5f;
  Shape->Transform = ww::Translation(0.f, -1.f, 0.f);
  Shape->Print = true;
  ww::WorldAddObject(W, Shape);

  ww::ray const R{ww::Point(0.f, 0.f, -3.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f)};
  ww::intersection const I = ww::intersection{M_SQRT2, Shape};
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
  ww::tup const Color = ww::ReflectedColor(W, Comps, 0);

  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  EXPECT_EQ(Color == Black, true);
}

//------------------------------------------------------------------------------
TEST(CH11ReflectionAndRefraction, PuttingItTogether)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  ww::shared_ptr_shape PtrRight = ww::PtrDefaultSphere();
  ww::shape &Right = *PtrRight;
  Right.Transform = ww::Translation(1.5f, 0.5f, -0.5f) *  //!<
                    ww::Scaling(0.5f, 0.5f, 0.5f);
  Right.Material.Color = ww::Color(0.5f, 1.0f, 0.1f);
  Right.Material.Diffuse = 0.7f;
  Right.Material.Specular = 0.3f;
  Right.Material.Reflective = 0.4f;
  World.vPtrObjects.push_back(PtrRight);

  ww::shared_ptr_shape PtrMiddle = ww::PtrDefaultSphere();
  ww::shape &Middle = *PtrMiddle;
  Middle.Transform = ww::Translation(-0.5f, 1.f, 0.5f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  Middle.Material.Reflective = 0.9f;
  // Middle.Material.Pattern = ww::CheckersPattern(ww::Color(0.f, 0.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  // Middle.Material.Pattern.Transform = ww::Scaling(0.4f, 0.4f, 0.4f);
  World.vPtrObjects.push_back(PtrMiddle);

  ww::shared_ptr_shape PtrLeft = ww::PtrDefaultSphere();
  ww::shape &Left = *PtrLeft;
  Left.Transform = ww::Translation(-1.5f, 0.33f, -0.75f) *  //!<
                   ww::Scaling(0.33f, 0.33f, 0.33f);
  Left.Material.Color = ww::Color(1.0f, 0.8f, 0.1f);
  Left.Material.Diffuse = 0.7f;
  Left.Material.Specular = 0.3f;
  Left.Material.Reflective = 0.4f;
  World.vPtrObjects.push_back(PtrLeft);

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::TranslateScaleRotate(0.f, 0.f, 00.f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
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
  ptrPlane2->Material.Pattern = ww::CheckersPattern(ww::Color(0.f, 0.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrObjects.push_back(ptrPlane2);

  // Add a third plane - this will act as the backdrop.
  ww::shared_ptr_plane ptrPlane3 = ww::PtrDefaultPlane();
  ptrPlane3->Transform =
      ww::TranslateScaleRotate(0.f, 0.f, 10.f, 1.f, 1.f, 1.f, ww::Radians(90.f), 0.f, ww::Radians(-45.f));
  ptrPlane3->Material.Shininess = 100.f;
  ptrPlane3->Material.Diffuse = 0.7f;
  ptrPlane3->Material.Color =
      ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  ptrPlane3->Material.Pattern = ww::CheckersPattern(ww::Color(0.f, 0.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrObjects.push_back(ptrPlane3);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 3.5f, -5.f);
  ww::tup const ViewTo = ww::Point(0.f, 1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch11MakingAScene.ppm");
}
