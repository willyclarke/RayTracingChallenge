#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>
//------------------------------------------------------------------------------
// Scenario: Reflectivity for the default material
TEST(DISABLED_CH11ReflectionAndRefraction, ReflectivityForTheDefaultMaterial)
{
  ww::material const M = {};
  EXPECT_EQ(M.Reflective == 0.f, true);
}

//------------------------------------------------------------------------------
// Scenario: Precomputing the reflective vector
TEST(DISABLED_CH11ReflectionAndRefraction, PrecomputingTheReflectiveVector)
{
  ww::shared_ptr_plane const Shape = ww::PtrDefaultPlane();
  ww::ray const R = ww::Ray(ww::Point(0.f, 1.f, -1.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f));
  ww::intersection const I = ww::intersection{M_SQRT2, Shape};
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
  EXPECT_EQ(Comps.vReflect == ww::Vector(0.f, M_SQRT2 / 2.f, M_SQRT2 / 2.f), true);
}

//------------------------------------------------------------------------------
// Scenario: The reflected color for a nonreflective material
TEST(DISABLED_CH11ReflectionAndRefraction, TheReflectedColorForANonReflectiveMaterial)
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
TEST(DISABLED_CH11ReflectionAndRefraction, TheReflectedColorForAReflectiveMaterial)
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
TEST(DISABLED_CH11ReflectionAndRefraction, ShadeHitWithAReflectiveMaterial)
{
  ww::world W = ww::World();

  // ---
  // NOTE: The floor should be element 2 in the default world.
  // ---
  ww::shared_ptr_shape Shape{};
  if (W.vPtrObjects.size() > 2)
  {
    Shape = W.vPtrObjects[2];
    Shape->Material.Reflective = 0.5f;
    Shape->Material.Transparency = 0.f;
    EXPECT_EQ(Shape->isA<ww::plane>(), true);
  }

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -3.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f));
  ww::intersection const I = ww::intersection{M_SQRT2, Shape};
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
  ww::tup const Color = ww::ShadeHit(W, Comps);
  EXPECT_EQ(Color == ww::Color(0.87677f, 0.92436f, 0.82918), true);

#if 0
  std::cout << "\n---\n" << std::endl;
  std::cout << "Computed Color:" << Color << std::endl;
  std::cout << "Expected Color:" << ww::Color(0.87677f, 0.92436f, 0.82918) << std::endl;
#endif

  // ---
  // NOTE: Write the result to file.
  // ---
  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 50.f, 0.f), ww::Color(1.f, 1.f, 1.f) * 0.2f);
  W.vPtrLights.push_back(pLight);

  float const FieldOfView = 75.f / 180.f * M_PI;
  ww::camera Camera = ww::Camera(256, 256, FieldOfView);

  ww::tup const ViewFrom = ww::Point(0.f, 3.5f, -10.f);
  ww::tup const ViewTo = ww::Point(0.f, 1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, W);
  ww::WriteToPPM(Canvas, "Ch11ShadeHitWithAReflectiveMaterial.ppm");
}

//------------------------------------------------------------------------------
// Scenario: color_at() with mutually reflective surface
TEST(DISABLED_CH11ReflectionAndRefraction, ColorAtWithMutuallyReflectiveSurface)
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
TEST(DISABLED_CH11ReflectionAndRefraction, TheReflectedColorAtTheMaximumRecursiveDepth)
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
TEST(DISABLED_CH11ReflectionAndRefraction, PuttingItTogether)
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

  ww::shared_ptr_shape PtrGlass = ww::PtrGlassSphere();
  ww::shape &Glass = *PtrGlass;
  Glass.Transform = ww::Translation(0.f, 0.33f, -1.75f) *  //!<
                    ww::Scaling(0.5f, 0.5f, 0.5f);
  Glass.Material.Color = ww::Color(0.09f, 0.09f, 0.09f);
  Glass.Material.Diffuse = 0.2f;
  Glass.Material.Ambient = 0.4f;
  Glass.Material.Specular = 0.3f;
  Glass.Material.Reflective = 1.0f;
  World.vPtrObjects.push_back(PtrGlass);

  // Add a second glass sphere
  ww::shared_ptr_shape PtrGlass2 = ww::PtrGlassSphere();
  ww::shape &Glass2 = *PtrGlass2;
  Glass2.Transform = ww::Translation(2.f, 1.0, 2.75f) *  //!<
                     ww::Scaling(1.f, 1.f, 1.f);
  Glass2.Material.Color = ww::Color(0.09f, 0.19f, 0.09f);
  Glass2.Material.Diffuse = 0.2f;
  Glass2.Material.Ambient = 0.4f;
  Glass2.Material.Specular = 0.3f;
  Glass2.Material.Reflective = 1.0f;
  World.vPtrObjects.push_back(PtrGlass2);

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
  ptrPlane3->Material.Reflective = 1.f;
  // ptrPlane3->Material.Pattern = ww::CheckersPattern(ww::Color(0.f, 0.f, 0.f), ww::Color(1.f, 1.f, 1.f));
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

//------------------------------------------------------------------------------
// Scenario: Transparency and Refractive index for the default material.
TEST(DISABLED_CH11ReflectionAndRefraction, TransparencyAndRefractiveIndexForTheDefaultMaterial)
{
  ww::material M = ww::material{};
  EXPECT_EQ(M.Transparency, 0.f);
  EXPECT_EQ(M.RefractiveIndex, 1.f);
}

//------------------------------------------------------------------------------
// Scenario: A helper for producing a sphere with a glassy material.
TEST(DISABLED_CH11ReflectionAndRefraction, AHelperForProducingASphereWithAGlassyMaterial)
{
  ww::shared_ptr_sphere Sphere = ww::PtrGlassSphere();
  EXPECT_EQ(Sphere->Transform == ww::I(), true);
  EXPECT_FLOAT_EQ(Sphere->Material.Transparency, 1.f);
  EXPECT_FLOAT_EQ(Sphere->Material.RefractiveIndex, 1.5f);
}

//------------------------------------------------------------------------------
// Scenario: Finding n1 and n2 at various intersections.
TEST(DISABLED_CH11ReflectionAndRefraction, FindingN1AndN2AtVariousIntersections)
{
  ww::shared_ptr_sphere A = ww::PtrGlassSphere();
  A->Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 2.f, 2.f, 2.f, 0.f, 0.f, 0.f);
  EXPECT_FLOAT_EQ(A->Material.RefractiveIndex, 1.5f);

  ww::shared_ptr_sphere B = ww::PtrGlassSphere();
  B->Transform = ww::TranslateScaleRotate(0.f, 0.f, -0.25f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  B->Material.RefractiveIndex = 2.f;
  EXPECT_FLOAT_EQ(B->Material.RefractiveIndex, 2.f);

  ww::shared_ptr_sphere C = ww::PtrGlassSphere();
  C->Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.25f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  C->Material.RefractiveIndex = 2.5f;
  EXPECT_FLOAT_EQ(C->Material.RefractiveIndex, 2.5f);

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -4.f), ww::Vector(0.f, 0.f, 1.f));

  ww::intersections XS{};
  XS = ww::Intersections(XS, {2.00f, A});
  XS = ww::Intersections(XS, {2.75f, B});
  XS = ww::Intersections(XS, {3.25f, C});
  XS = ww::Intersections(XS, {4.75f, B});
  XS = ww::Intersections(XS, {5.25f, C});
  XS = ww::Intersections(XS, {6.00f, A});

  struct refractive_index_for_test
  {
    float n1{1.f};
    float n2{1.f};
  };

  // ---
  // NOTE: Create a test vector of refractive indexes.
  // ---
  std::vector<refractive_index_for_test> vRefr{};
  vRefr.push_back(refractive_index_for_test{1.0f, 1.5f});
  vRefr.push_back(refractive_index_for_test{1.5f, 2.0f});
  vRefr.push_back(refractive_index_for_test{2.0f, 2.5f});
  vRefr.push_back(refractive_index_for_test{2.5f, 2.5f});
  vRefr.push_back(refractive_index_for_test{2.5f, 1.5f});
  vRefr.push_back(refractive_index_for_test{1.5f, 1.0f});

  // ---
  // NOTE: The vector need to be of the same lenght.
  // ---
  EXPECT_EQ(vRefr.size(), XS.Count());

  for (size_t Idx = 0;                          //!<
       Idx < vRefr.size() && Idx < XS.Count();  //!<
       ++Idx)
  {
    EXPECT_EQ(ww::PrepareComputations(XS.vI[Idx], R, &XS).n1, vRefr[Idx].n1);
    EXPECT_EQ(ww::PrepareComputations(XS.vI[Idx], R, &XS).n2, vRefr[Idx].n2);
  }
}

//------------------------------------------------------------------------------
// Scenario: The under point is offset below the surface.
TEST(DISABLED_CH11ReflectionAndRefraction, TheUnderPointIsOffsetBelowTheSurface)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_sphere Shape = ww::PtrGlassSphere();
  Shape->Transform = ww::TranslateScaleRotate(0.f, 0.f, 1.f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  EXPECT_FLOAT_EQ(Shape->Material.RefractiveIndex, 1.5f);

  ww::intersection const I{5.f, Shape};
  ww::intersections const XS = ww::Intersections(I);
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R, &XS);
  EXPECT_GT(Comps.UnderPoint.Z, ww::EPSILON / 2.f);
  EXPECT_LT(Comps.Point.Z, Comps.UnderPoint.Z);
}

//------------------------------------------------------------------------------
// Scenario: The refracted color with an opaque surface.
TEST(DISABLED_CH11ReflectionAndRefraction, TheRefractedColorWithAnOpaqueSurface)
{
  ww::world W = ww::World();
  ww::shared_ptr_shape Shape = W.vPtrObjects[0];

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));

  ww::intersections XS{};
  XS = ww::Intersections(XS, {4.f, Shape});
  XS = ww::Intersections(XS, {6.f, Shape});

  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[0], R, &XS);
  ww::tup const Color = ww::RefractedColor(W, Comps);
  EXPECT_EQ(Color == ww::Color(0.f, 0.f, 0.f), true);
}

//------------------------------------------------------------------------------
// Scenario: The refracted color at the maximum recursive depth.
TEST(DISABLED_CH11ReflectionAndRefraction, TheRefractedColorAtTheMaximumRecursiveDepth)
{
  ww::world W = ww::World();
  ww::shared_ptr_shape Shape = W.vPtrObjects[0];
  Shape->Material.Transparency = 1.f;
  Shape->Material.RefractiveIndex = 1.5;
  EXPECT_EQ(Shape->Material.Transparency, 1.f);
  EXPECT_EQ(Shape->Material.RefractiveIndex, 1.5f);

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));

  ww::intersections XS{};
  XS = ww::Intersections(XS, {4.f, Shape});
  XS = ww::Intersections(XS, {6.f, Shape});

  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[0], R, &XS);
  ww::tup const Color = ww::RefractedColor(W, Comps, 0);
  EXPECT_EQ(Color == ww::Color(0.f, 0.f, 0.f), true);
}

//------------------------------------------------------------------------------
// Scenario: The refracted color under total internal reflection.
TEST(DISABLED_CH11ReflectionAndRefraction, TheRefractedColorUnderTotalInternalReflection)
{
  ww::world W = ww::World();
  ww::shared_ptr_shape Shape = W.vPtrObjects[0];
  Shape->Material.Transparency = 1.f;
  Shape->Material.RefractiveIndex = 1.5;
  EXPECT_EQ(Shape->Material.Transparency, 1.f);
  EXPECT_EQ(Shape->Material.RefractiveIndex, 1.5f);

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, M_SQRT2 / 2.f), ww::Vector(0.f, 1.f, 0.f));

  ww::intersections XS{};
  XS = ww::Intersections(XS, {-M_SQRT2 / 2.f, Shape});
  XS = ww::Intersections(XS, {M_SQRT2 / 2.f, Shape});
  // NOTE: this time youre inside the sphere, so you need
  //       to look at the second intersection, XS[1] and not XS[0].

  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[1], R, &XS);
  ww::tup const Color = ww::RefractedColor(W, Comps, 5);
  EXPECT_EQ(Color == ww::Color(0.f, 0.f, 0.f), true);
}

//------------------------------------------------------------------------------
// Scenario: The refracted color with a refracted ray.
TEST(DISABLED_CH11ReflectionAndRefraction, TheRefractedColorWithARefractedRay)
{
  ww::world W = ww::World();

  ww::shared_ptr_shape A = W.vPtrObjects[0];
  A->Material.Ambient = 1.f;
  A->Material.Pattern = ww::TestPattern();
  A->Material.Pattern.Print = true;
  EXPECT_EQ(A->Material.Ambient, 1.f);

  ww::shared_ptr_shape B = W.vPtrObjects[1];
  B->Material.Transparency = 1.f;
  B->Material.RefractiveIndex = 1.5f;
  B->Material.Pattern = ww::TestPattern();
  B->Material.Pattern.Print = true;
  EXPECT_EQ(B->Material.Transparency, 1.f);
  EXPECT_EQ(B->Material.RefractiveIndex, 1.5f);

  EXPECT_EQ(W.vPtrObjects[0] == A, true);
  EXPECT_EQ(W.vPtrObjects[1] == B, true);

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 0.1f), ww::Vector(0.f, 1.f, 0.f));

  ww::intersections XS{};
  XS = ww::Intersections(XS, {-0.9899f, A});
  XS = ww::Intersections(XS, {-0.4899f, B});
  XS = ww::Intersections(XS, {0.4899f, B});
  XS = ww::Intersections(XS, {0.9899f, A});

  W.Print = true;

  ww::prepare_computation Comps{};
// Comps.PrintDebug = true;
#if 0
  std::cout << "Init::Comps.n1: " << Comps.n1 << std::endl;
  std::cout << "Init::Comps.n2: " << Comps.n2 << std::endl;
#endif

  Comps = ww::PrepareComputations(XS.vI[2], R, &XS);
  // Comps.PrintDebug = true;
  // std::cout << "Comps:\n" << Comps << std::endl;

#if 0
  std::cout << "\n---\n" << std::endl;
  std::cout << "\n---\n" << std::endl;
  std::cout << "Hit to examine: " << XS.vI[2].pShape << " with t=" << XS.vI[2].t << std::endl;
  std::cout << "Calc::Comps.n1: " << Comps.n1 << std::endl;
  std::cout << "Calc::Comps.n2: " << Comps.n2 << std::endl;
#endif

  ww::tup const Color = ww::RefractedColor(W, Comps, 5);

#if 0
  std::cout << "\n---\n" << std::endl;
  std::cout << "Input Intersections:\n" << XS << std::endl;
  std::cout << "Computed Color:" << Color << std::endl;
  std::cout << "Expected Color:" << ww::Color(0.f, 0.99888f, 0.04725f) << std::endl;
#endif

  EXPECT_EQ(Color == ww::Color(0.f, 0.99888f, 0.04725f), true);

  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 3.5f, -5.f);
  ww::tup const ViewTo = ww::Point(0.f, 1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, W);
  ww::WriteToPPM(Canvas, "Ch11RefractColorRefractRay.ppm");
}

//------------------------------------------------------------------------------
// Scenario: ShadeHit with a transparent material.
TEST(DISABLED_CH11ReflectionAndRefraction, ShadeHitWithATransparentMaterial)
{
  ww::world W = ww::World();

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -3.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f));
  ww::shared_ptr_shape Floor = W.vPtrObjects[2];

  ww::intersections XS{};
  XS = ww::Intersections(XS, {M_SQRT2, Floor});

  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[0], R, &XS);
  ww::tup Color = ww::ShadeHit(W, Comps);
  EXPECT_EQ(Color == ww::Color(0.93642f, 0.68642f, 0.68648f), true);
#if 0
  std::cout << "\n---\n" << std::endl;
  std::cout << "Input Intersections:\n" << XS << std::endl;
  std::cout << "Computed Color:" << Color << std::endl;
  std::cout << "Expected Color:" << ww::Color(0.93642f, 0.68642f, 0.68648f) << std::endl;
#endif

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 3.5f, -10.f);
  ww::tup const ViewTo = ww::Point(0.f, 1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, W);
  ww::WriteToPPM(Canvas, "Ch11ShadeHitWithATransparentMaterial.ppm");
}

//------------------------------------------------------------------------------
// Scenario: The Schlick approximation under total internal reflection.
TEST(DISABLED_CH11ReflectionAndRefraction, TheSchlickApproximationUnderTotalInternalReflection)
{
  ww::shared_ptr_sphere Shape = ww::PtrGlassSphere();
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, M_SQRT2 / 2.f), ww::Vector(0.f, 1.f, 0.f));
  ww::intersections XS{};
  XS = ww::Intersections(XS, {-M_SQRT2 / 2.f, Shape});
  XS = ww::Intersections(XS, {M_SQRT2 / 2.f, Shape});
  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[1], R, &XS);

  EXPECT_EQ(Comps.n1, 1.5f);
  EXPECT_EQ(Comps.n2, 1.f);
  float const Reflectance = ww::Schlick(Comps);

  EXPECT_FLOAT_EQ(Reflectance, 1.f);
}

//------------------------------------------------------------------------------
// Scenario: The Schlick approximation with a perpendicular viewing angle.
TEST(DISABLED_CH11ReflectionAndRefraction, TheSchlickApproximationWithAPerpendicularViewingAngle)
{
  ww::shared_ptr_sphere Shape = ww::PtrGlassSphere();
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 1.f, 0.f));
  ww::intersections XS{};
  XS = ww::Intersections(XS, {-1.f, Shape});
  XS = ww::Intersections(XS, {1.f, Shape});
  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[1], R, &XS);

  float const Reflectance = ww::Schlick(Comps);
  EXPECT_FLOAT_EQ(Reflectance, 0.04f);
}

//------------------------------------------------------------------------------
// Scenario: The Schlick approximation with small angle and n2 > n1 .
TEST(DISABLED_CH11ReflectionAndRefraction, TheSchlickApproximationWithSmallAngleAndn1GTn2)
{
  ww::shared_ptr_sphere Shape = ww::PtrGlassSphere();
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.99f, -2.f), ww::Vector(0.f, 0.f, 1.f));
  ww::intersections XS{};
  XS = ww::Intersections(XS, {1.8589f, Shape});
  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[0], R, &XS);

  float const Reflectance = ww::Schlick(Comps);
  EXPECT_EQ(Comps.n1 < Comps.n2, true);
  EXPECT_EQ(ww::Equal(Reflectance, 0.48873), true);
}

//------------------------------------------------------------------------------
// Scenario: shade_hit() with a reflective, transparent material.
TEST(DISABLED_CH11ReflectionAndRefraction, ShadeHitWithAReflectiveTranparentMaterial)
{
  ww::world W = ww::World();

  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -3.f), ww::Vector(0.f, -M_SQRT2 / 2.f, M_SQRT2 / 2.f));

  ww::shared_ptr_shape Floor{};
  if (W.vPtrObjects.size() > 2)
  {
    Floor = W.vPtrObjects[2];
    Floor->Material.Reflective = 0.5f;
  }

  ww::intersections XS{};
  XS = ww::Intersections(XS, {M_SQRT2, Floor});

  ww::prepare_computation const Comps = ww::PrepareComputations(XS.vI[0], R, &XS);
  ww::tup Color = ww::ShadeHit(W, Comps);
  EXPECT_EQ(Color == ww::Color(0.93391f, 0.69643f, 0.69243f), true);

#if 0
  std::cout << "\n---\n" << std::endl;
  std::cout << "Input Intersections:\n" << XS << std::endl;
  std::cout << "Computed Color:" << Color << std::endl;
  std::cout << "Expected Color:" << ww::Color(0.93391f, 0.69643f, 0.69243f) << std::endl;
#endif

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 3.5f, -10.f);
  ww::tup const ViewTo = ww::Point(0.f, 1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, W);
  ww::WriteToPPM(Canvas, "Ch11ShadeHitWithAReflectiveTransparentMaterial.ppm");
}
