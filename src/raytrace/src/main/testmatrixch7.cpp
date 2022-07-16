#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>

//------------------------------------------------------------------------------
// Scenario: Creating a world, an empty one...
TEST(Ch7MakingAScene, CreateAWorld)
{
  ww::world const World{};
  EXPECT_EQ(World.Count(), 0);
  EXPECT_EQ(World.vPtrLights.size(), 0);
  EXPECT_EQ(World.vPtrObjects.size(), World.Count());
}

//------------------------------------------------------------------------------
// Scenario: The default world
TEST(Ch7MakingAScene, DefaultWorld)
{
  ww::light const Light = ww::PointLight(ww::Point(-10.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));
  ww::sphere S1{};
  ww::sphere S2{};

  // NOTE: At this point the spheres are equal.
  EXPECT_EQ(S1 == S2, true);

  // NOTE: Set up material for sphere 1 so that it differs from sphere 2.
  S1.Material.Color = ww::Color(0.8f, 1.0f, 0.6f);
  S1.Material.Diffuse = 0.7f;
  S1.Material.Specular = 0.2f;

  S2.Transform = ww::Scaling(0.5f, 0.5f, 0.5f);
  S2.Radius = S2.Transform.R0.X;
  EXPECT_EQ(S2.Radius == 0.5f, true);

  // NOTE: And now the two spheres are different.
  EXPECT_EQ(!(S1 == S2), true);

  // NOTE: Assign our default world.
  ww::world const W = ww::World();

  EXPECT_EQ(W.vPtrLights.size(), 1);   //!< So now we expect there to be a light.
  EXPECT_EQ(W.vPtrObjects.size(), 4);  //!< and some objects.

  bool ContainsS1{};
  bool ContainsS2{};

  for (auto const &PtrObject : W.vPtrObjects)
  {
    if (PtrObject->isA<ww::sphere>())
    {
      ww::sphere *pSphere = dynamic_cast<ww::sphere *>(PtrObject.get());
      if (S1 == *pSphere)
      {
        ContainsS1 = true;
        // std::cerr << "Sphere1 found.\n" << *pSphere << std::endl;
      }
      if (S2 == *pSphere)
      {
        ContainsS2 = true;
        // std::cerr << "Sphere2 found.\n" << *pSphere << std::endl;
      }
    }
    else
      std::cerr << "Not a sphere" << std::endl;
  }
  EXPECT_EQ(ContainsS1, true);
  EXPECT_EQ(ContainsS2, true);
  if (W.vPtrLights.size() > 0)
  {
    ww::shared_ptr_light PtrLight = W.vPtrLights[0];
    ww::light &L = *PtrLight;
    bool const LightsAreEqual = L == Light;
    EXPECT_EQ(LightsAreEqual, true);
    EXPECT_EQ(bool(Light == *W.vPtrLights[0].get()), true);
  }
}

//------------------------------------------------------------------------------
// Scenario: Intersect a world with a ray
TEST(Ch7MakingAScene, IntersectWorldWithRay)
{
  ww::world const World = ww::World();
  ww::ray const Ray = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::intersections const XS = IntersectWorld(World, Ray);

  EXPECT_EQ(XS.Count(), 4);
  // NOTE: The expectation is that the hits are delivered in ascending order.
  if (XS.Count() >= 1) EXPECT_FLOAT_EQ(XS.vI[0].t, 4.f);
  if (XS.Count() >= 2) EXPECT_FLOAT_EQ(XS.vI[1].t, 4.5f);
  if (XS.Count() >= 3) EXPECT_FLOAT_EQ(XS.vI[2].t, 5.5f);
  if (XS.Count() >= 4) EXPECT_FLOAT_EQ(XS.vI[3].t, 6.f);
  // Assert(XS.Count() == -1, __FUNCTION__, __LINE__); // NOTE: for debugging
}

#if 1
//------------------------------------------------------------------------------
TEST(Ch7MakingAScene, PrecomputingStateOfIntersection)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape pSphere{};
  pSphere.reset(new ww::sphere);
  ww::intersection const I = ww::Intersection(4.f, pSphere);
  ww::prepare_computation Comps = ww::PrepareComputations(I, R);

  EXPECT_EQ(I.t, Comps.t);
  EXPECT_EQ(I.pShape == Comps.pShape, true);
  EXPECT_EQ(ww::Equal(Comps.Point, ww::Point(0.f, 0.f, -1.f)), true);
  EXPECT_EQ(ww::Equal(Comps.vEye, ww::Vector(0.f, 0.f, -1.f)), true);
  EXPECT_EQ(ww::Equal(Comps.vNormal, ww::Vector(0.f, 0.f, -1.f)), true);
}

//------------------------------------------------------------------------------
TEST(Ch7MakingAScene, HitFromInsideOrOutsideChangesTheNormal)
{
  {
    // NOTE: Scenario - The hit when an intersection occurs on the outside.
    ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
    ww::shared_ptr_shape pShape{};
    pShape.reset(new ww::sphere);
    ww::intersection const I = ww::Intersection(4.f, pShape);
    ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
    EXPECT_EQ(Comps.Inside, false);
  }

  {
    // NOTE: Scenario - The hit when an intersection occurs on the inside.
    ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 0.f, 1.f));
    ww::shared_ptr_shape pShape{};
    pShape.reset(new ww::sphere);
    ww::intersection const I = ww::Intersection(1.f, pShape);
    ww::prepare_computation const Comps = ww::PrepareComputations(I, R);

    EXPECT_EQ(ww::Equal(Comps.Point, ww::Point(0.f, 0.f, 1.f)), true);
    EXPECT_EQ(ww::Equal(Comps.vEye, ww::Vector(0.f, 0.f, -1.f)), true);
    EXPECT_EQ(Comps.Inside, true);
    // NOTE: Normal would have been 0,0,1 but is inverted.
    EXPECT_EQ(ww::Equal(Comps.vNormal, ww::Vector(0.f, 0.f, -1.f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch7MakingAScene, ShadingAnIntersection)
{
  // NOTE: the default world consist of two spheres.
  ww::world W = ww::World();

  // NOTE: the default world should have two objects
  EXPECT_EQ(W.vPtrObjects.size(), 4);
  ww::tup const C1 = ww::Color(0.8f, 1.0f, 0.6f);
  ww::tup const C2 = ww::Color(1.0f, 1.0f, 1.0f);
  int Idx{1};
  for (auto PtrShape : W.vPtrObjects)
  {
    if (1 == Idx) EXPECT_EQ(PtrShape->Material.Color == C1, true);
    if (2 == Idx) EXPECT_EQ(PtrShape->Material.Color == C2, true);
    ++Idx;
  }

  // Scenario - Shading an intersection
  if (W.vPtrObjects.size() > 0)
  {
    ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
    ww::shared_ptr_shape pShape = W.vPtrObjects[0];

    ww::intersection const I = ww::Intersection(4.f, pShape);
    ww::prepare_computation const Comps = ww::PrepareComputations(I, R);

    ww::tup const C = ww::ShadeHit(W, Comps);

    EXPECT_EQ(Comps.Inside, false);
    EXPECT_EQ(C == ww::Color(0.38066f, 0.47583f, 0.2855f), true);
  }

  // Scenario - Shading an intersection from the inside
  if ((W.vPtrObjects.size() > 1) && (W.vPtrLights.size() > 0))
  {
    ww::shared_ptr_light pLight{};
    pLight.reset(new ww::light);
    *pLight = ww::PointLight(ww::Point(0.f, 0.25f, 0.f), ww::Color(1.f, 1.f, 1.f));
    W.vPtrLights[0] = pLight;  // Assign the new light

    ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 0.f, 1.f));

    ww::shared_ptr_shape pShape = W.vPtrObjects[1];

    ww::intersection const I = ww::Intersection(0.5f, pShape);
    ww::prepare_computation const Comps = ww::PrepareComputations(I, R);
    ww::tup const C = ww::ShadeHit(W, Comps);

    // NOTE: When shading has been implemented the color changes
    //       The old/without shading color was:
    // EXPECT_EQ(ww::Equal(C, ww::Color(0.90498f, 0.90498f, 0.90498f)), true);
    // NOTE: The shaded color becomes:
    EXPECT_EQ(ww::Equal(C, ww::Color(0.1f, 0.1f, 0.1f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch7MakingAScene, TheColorWhenARayMisses)
{
  ww::world const W = ww::World();
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 1.f, 0.f));
  ww::tup const C = ww::ColorAt(W, R);
  EXPECT_EQ(C == ww::Color(0.f, 0.f, 0.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch7MakingAScene, TheColorWhenARayHits)
{
  ww::world const W = ww::World();
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::tup const C = ww::ColorAt(W, R);
  EXPECT_EQ(C == ww::Color(0.38066f, 0.47583f, 0.2855f), true);
}

//------------------------------------------------------------------------------
TEST(Ch7MakingAScene, TheColorWhenIntersectionBehindTheRay)
{
  ww::world W = ww::World();
  if (W.vPtrObjects.size() > 1)
  {
    int const OUTER{0};
    int const INNER{1};
    W.vPtrObjects[OUTER]->Material.Ambient = 1.f;  //!< The first object in the World. The outer sphere.
    W.vPtrObjects[INNER]->Material.Ambient = 1.f;  //!< The second object in the World. The inner sphere.

    // ---
    // NOTE: Here we put the ray inside the outer sphere, but outside the inner sphere.
    //       Also; we are pointing at the inner sphere (z=-1). We expect the hit to be
    //       on the inner sphere and return its color.
    // ---
    ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 0.75f), ww::Vector(0.f, 0.f, -1.f));

    ww::tup const C = ww::ColorAt(W, R);
    EXPECT_EQ(C == W.vPtrObjects[INNER]->Material.Color, true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch7DefiningAViewTransformation, TransformationMatrixForDefaultOrientation)
{
  ww::tup const From = ww::Point(0.f, 0.f, 0.f);
  ww::tup const To = ww::Point(0.f, 0.f, -1.f);
  ww::tup const Up = ww::Vector(0.f, 1.f, 0.f);
  ww::matrix const T = ww::ViewTransform(From, To, Up);
  EXPECT_EQ(ww::Equal(ww::I(), T), true);
}
//------------------------------------------------------------------------------
TEST(Ch7DefiningAViewTransformation, AViewTransformationMatrixInThePositiveZDirection)
{
  ww::tup const From = ww::Point(0.f, 0.f, 0.f);
  ww::tup const To = ww::Point(0.f, 0.f, 1.f);
  ww::tup const Up = ww::Vector(0.f, 1.f, 0.f);
  ww::matrix const T = ww::ViewTransform(From, To, Up);
  ww::matrix const CheckTransform = ww::Scaling(-1.f, 1.f, -1.f);

  // std::cout << "T:" << T << std::endl;
  // std::cout << "ViewTransform:\n" << ViewTransform << std::endl;

  EXPECT_EQ(ww::Equal(CheckTransform, T), true);
}

//------------------------------------------------------------------------------
TEST(Ch7DefiningAViewTransformation, TheViewTransformationMovesTheWorld)
{
  ww::tup const From = ww::Point(0.f, 0.f, 8.f);
  ww::tup const To = ww::Point(0.f, 0.f, 0.f);
  ww::tup const Up = ww::Vector(0.f, 1.f, 0.f);
  ww::matrix const T = ww::ViewTransform(From, To, Up);
  ww::matrix const CheckTransform = ww::Translation(0.f, 0.f, -8.f);
  EXPECT_EQ(ww::Equal(CheckTransform, T), true);
}

//------------------------------------------------------------------------------
TEST(Ch7DefiningAViewTransformation, AnArbitraryViewTransformation)
{
  ww::tup const From = ww::Point(1.f, 3.f, 2.f);
  ww::tup const To = ww::Point(4.f, -2.f, 8.f);
  ww::tup const Up = ww::Vector(1.f, 1.f, 0.f);
  ww::matrix const T = ww::ViewTransform(From, To, Up);
  ww::matrix const CheckTransform = ww::Matrix44(ww::tup{-0.50709, 0.50709, 0.67612, -2.36643},     //!<
                                                 ww::tup{0.76772f, 0.60609f, 0.12122f, -2.82843f},  //!<
                                                 ww::tup{-0.35857f, 0.59761f, -0.71714f, 0.f},      //!<
                                                 ww::tup{0.f, 0.f, 0.f, 1.f}                        //!<
  );
  // std::cout << "CheckTransform:\n" << CheckTransform << std::endl;

  EXPECT_EQ(ww::Equal(CheckTransform, T), true);
}

//------------------------------------------------------------------------------
TEST(Ch7ImplementingACamera, ConstructingACamera)
{
  int const HSize{160};
  int const VSize{120};
  int const Alfa{180};
  int const FieldOfView{Alfa / 2};

  ww::camera C = ww::Camera(HSize, VSize, ww::Radians(FieldOfView));
  EXPECT_EQ(HSize, C.HSize);
  EXPECT_EQ(VSize, C.VSize);
  EXPECT_EQ(ww::Radians(Alfa / 2), C.FieldOfView);
  EXPECT_EQ(ww::Equal(ww::I(), C.Transform), true);
  EXPECT_EQ(ww::I() == C.Transform, true);
}

//------------------------------------------------------------------------------
TEST(Ch7ImplementingACamera, PixelSizeForAHorisontalCanvas)
{
  int const Alfa{180};
  int const FieldOfView{Alfa / 2};
  {
    ww::camera const C = ww::Camera(200, 125, ww::Radians(FieldOfView));
    EXPECT_EQ(0.01f, C.PixelSize);
  }
  {
    ww::camera const C = ww::Camera(125, 200, ww::Radians(FieldOfView));
    EXPECT_EQ(0.01f, C.PixelSize);
  }
}

//------------------------------------------------------------------------------
TEST(Ch7ImplementingACamera, ConstructingARayThroughTheCenterOfTheCanvas)
{
  int const Alfa{180};
  int const FieldOfView{Alfa / 2};
  ww::camera const C = ww::Camera(201, 101, ww::Radians(FieldOfView));
  ww::ray const R = ww::RayForPixel(C, 100, 50);
  EXPECT_EQ(R.Origin == ww::Point(0.f, 0.f, 0.f), true);

  ww::tup const ExpectedVector = ww::Vector(0.f, 0.f, -1.f);
  EXPECT_EQ(R.Direction == ExpectedVector, true);
}

//------------------------------------------------------------------------------
TEST(Ch7ImplementingACamera, ConstructingARayThroughThroughACornerOfTheCanvas)
{
  int const Alfa{180};
  int const FieldOfView{Alfa / 2};
  ww::camera const C = ww::Camera(201, 101, ww::Radians(FieldOfView));
  ww::ray const R = ww::RayForPixel(C, 0, 0);

  EXPECT_EQ(R.Origin == ww::Point(0.f, 0.f, 0.f), true);

  ww::tup const ExpectedVector = ww::Vector(0.66519f, 0.33259f, -0.66851f);
  EXPECT_EQ(R.Direction == ExpectedVector, true);
}

//------------------------------------------------------------------------------
TEST(Ch7ImplementingACamera, ConstructingARayWhenTheCameraIsTransformed)
{
  int const Alfa{180};
  int const FieldOfView{Alfa / 2};
  ww::camera C = ww::Camera(201, 101, ww::Radians(FieldOfView));
  C.Transform = ww::RotateY(ww::Radians(Alfa / 4)) * ww::Translation(0.f, -2.f, 5.f);
  ww::ray const R = ww::RayForPixel(C, 100, 50);

  // NOTE: In this test the ray' origin winds up at (0,2,-5), despite the camera's
  //       transformation including a translation of (0,-2,5). That's not a typo!
  //       Remember that the camera's transformation describes how the world moves
  //       relative to the camera. Further, you're transforming everything by the
  //       inverse of that transformation. So moving the World is effectively the
  //       the same as moving the ray's origin in the opposite direction.
  EXPECT_EQ(R.Origin == ww::Point(0.f, 2.f, -5.f), true);
  float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
  EXPECT_EQ(R.Direction == ww::Vector(Sqrt2O2, 0.f, -Sqrt2O2), true);
}

//------------------------------------------------------------------------------
TEST(Ch7ImplementingACamera, RenderingAWorldWithACamera)
{
  ww::world const W = ww::World();

  int const n{180};
  int const FieldOfView{n / 2};
  ww::camera C = ww::Camera(11, 11, ww::Radians(FieldOfView));
  ww::tup const From = ww::Point(0.f, 0.f, -5.f);
  ww::tup const To = ww::Point(0.f, 0.f, 0.f);
  ww::tup const Up = ww::Vector(0.f, 1.f, 0.f);
  C.Transform = ww::ViewTransform(From, To, Up);
  ww::canvas const Image = ww::Render(C, W);
  EXPECT_EQ(ww::PixelAt(Image, 5, 5) == ww::Color(0.38066f, 0.47583f, 0.2855f), true);
}

//------------------------------------------------------------------------------
TEST(Ch7ImplementingACamera, PuttingItTogether)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  ww::shared_ptr_shape PtrFloor = ww::PtrDefaultSphere();
  ww::shape &Floor = *PtrFloor;
  Floor.Transform = ww::Scaling(10.f, 0.01f, 10.f);
  Floor.Material.Color = ww::Color(1.f, 0.9f, 0.9f);
  Floor.Material.Specular = 0.f;
  World.vPtrObjects.push_back(PtrFloor);

  ww::shared_ptr_shape PtrLeftWall = ww::PtrDefaultSphere();
  ww::shape &LeftWall = *PtrLeftWall;
  // NOTE: The order in which the transformations are multiplied;
  //       The wall need to be scaled, rotated in X, rotated in Y,
  //       and lastly translated, so the transformations are multiplied
  //       in reverse order.
  LeftWall.Transform = ww::Translation(0.f, 0.f, 5.f) *  //!<
                       ww::RotateY(-M_PI_4) *            //!<
                       ww::RotateX(M_PI_2) *             //!<
                       ww::Scaling(10.f, 0.01f, 10.f);
  LeftWall.Material = Floor.Material;
  World.vPtrObjects.push_back(PtrLeftWall);

  ww::shared_ptr_shape PtrRightWall = ww::PtrDefaultSphere();
  ww::shape &RightWall = *PtrRightWall;
  RightWall.Transform = ww::Translation(0.f, 0.f, 5.f) *  //!<
                        ww::RotateY(M_PI_4) *             //!<
                        ww::RotateX(M_PI_2) *             //!<
                        ww::Scaling(10.f, 0.01f, 10.f);
  RightWall.Material = Floor.Material;
  World.vPtrObjects.push_back(PtrRightWall);

  ww::shared_ptr_shape PtrMiddle = ww::PtrDefaultSphere();
  ww::shape &Middle = *PtrMiddle;
  Middle.Transform = ww::Translation(-0.5f, 1.f, 0.5f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrMiddle);

  ww::shared_ptr_shape PtrRight = ww::PtrDefaultSphere();
  ww::shape &Right = *PtrRight;
  Right.Transform = ww::Translation(1.5f, 0.5f, -0.5f) *  //!<
                    ww::Scaling(0.5f, 0.5f, 0.5f);
  Right.Material.Color = ww::Color(0.5f, 1.0f, 0.1f);
  Right.Material.Diffuse = 0.7f;
  Right.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrRight);

  ww::shared_ptr_shape PtrLeft = ww::PtrDefaultSphere();
  ww::shape &Left = *PtrLeft;
  Left.Transform = ww::Translation(-1.5f, 0.33f, -0.75f) * ww::Scaling(0.33f, 0.33f, 0.33f);
  Left.Material.Color = ww::Color(1.0f, 0.8f, 0.1f);
  Left.Material.Diffuse = 0.7f;
  Left.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrLeft);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  ww::camera Camera = ww::Camera(100, 50, M_PI / 3.f);
  // ww::camera Camera = ww::Camera(800, 400, M_PI / 3.f);

  Camera.Transform = ww::ViewTransform(ww::Point(0.f, 1.5f, -5.f), ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, 1.f, 0.f));

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch7MakingAScene.ppm");
}
#endif
