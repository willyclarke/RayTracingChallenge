#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>
//------------------------------------------------------------------------------
// Scenario: The default transformation
TEST(DISABLED_Ch9Planes, TheDefaultTransformation)
{
  ww::shape S = ww::TestShape();
  EXPECT_EQ(S.Transform == ww::I(), true);
}

//------------------------------------------------------------------------------
// Scenario: Assigning a transformation
TEST(DISABLED_Ch9Planes, AssigningATransformation)
{
  ww::shape S = ww::TestShape();

  // NOTE: Test that we can assign a translation.
  S.Transform = ww::Translation(2.f, 3.f, 4.f);
  EXPECT_EQ(S.Transform == ww::I(), false);
  EXPECT_EQ(S.Transform == ww::Translation(2.f, 3.f, 4.f), true);
}

//------------------------------------------------------------------------------
// Scenario: The default material
TEST(DISABLED_Ch9Planes, TheDefaultMaterial)
{
  ww::shape S = ww::TestShape();

  EXPECT_EQ(S.Transform == ww::I(), true);
  EXPECT_EQ(S.Material == ww::material{}, true);
}

//------------------------------------------------------------------------------
// Scenario: Assigning a material
TEST(DISABLED_Ch9Planes, AssignAMaterial)
{
  ww::shape S = ww::TestShape();
  ww::material M{};

  M.Ambient = 1.f;

  // NOTE: test that the materials are different before assignment
  EXPECT_EQ(S.Material == M, false);

  // NOTE: test that we can assign a material.
  S.Material = M;

  EXPECT_EQ(S.Material == M, true);
  EXPECT_EQ(S.Material.Ambient == 1.f, true);
}

//------------------------------------------------------------------------------
// Scenario: Intersecting a scaled shape with a ray
TEST(DISABLED_Ch9Planes, IntersectingAScaledAShapeWithRay)
{
  ww::ray R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shape S = ww::TestShape();
  S.Transform = ww::Scaling(2.f, 2.f, 2.f);

  ww::shared_ptr_shape PtrShape = ww::SharedPtrShape(S);

  ww::ray SavedRay{};
  ww::intersections XS = ww::Intersect(PtrShape, R, &SavedRay);

  EXPECT_EQ(SavedRay.Origin == ww::Point(0.f, 0.f, -2.5f), true);
  EXPECT_EQ(SavedRay.Direction == ww::Vector(0.f, 0.f, 0.5f), true);
  EXPECT_EQ(PtrShape->SavedRay.Origin == ww::Point(0.f, 0.f, -2.5f), true);
  EXPECT_EQ(PtrShape->SavedRay.Direction == ww::Vector(0.f, 0.f, 0.5f), true);
}

//------------------------------------------------------------------------------
// Scenario: Intersecting a translated shape with a ray
TEST(DISABLED_Ch9Planes, IntersectingATranslatedShapeWithARay)
{
  ww::ray R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shape S = ww::TestShape();
  S.Transform = ww::Translation(5.f, 0.f, 0.f);

  ww::shared_ptr_shape PtrShape = ww::SharedPtrShape(S);

  ww::ray SavedRay{};
  ww::intersections XS = ww::Intersect(PtrShape, R, &SavedRay);

  EXPECT_EQ(SavedRay.Origin == ww::Point(-5.f, 0.f, -5.f), true);
  EXPECT_EQ(SavedRay.Direction == ww::Vector(0.f, 0.f, 1.0f), true);
  EXPECT_EQ(PtrShape->SavedRay.Origin == ww::Point(-5.f, 0.f, -5.f), true);
  EXPECT_EQ(PtrShape->SavedRay.Direction == ww::Vector(0.f, 0.f, 1.0f), true);
}

//------------------------------------------------------------------------------
// Scenario: Computing the normal on a translated shape
TEST(DISABLED_Ch9Planes, ComputingTheNormalOnATranslatedSphere)
{
  ww::shape S = ww::TestShape();
  S.Transform = ww::Translation(0.f, 1.f, 0.f);

  ww::tup const N = ww::NormalAt(S, ww::Point(0.f, 1.70711f, -0.70711f));
  EXPECT_EQ(N == ww::Vector(0.f, 0.70711f, -0.70711f), true);
}

//------------------------------------------------------------------------------
// Scenario: Computing the normal on a transformed shape
TEST(DISABLED_Ch9Planes, ComputingTheNormalOnATransformedSphere)
{
  ww::shape S = ww::TestShape();
  S.Transform = ww::Scaling(1.f, 0.5f, 1.f) * ww::RotateZ(M_PI / 5.f);

  float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
  ww::tup const N = ww::NormalAt(S, ww::Point(0.f, Sqrt2O2, -Sqrt2O2));
  EXPECT_EQ(N == ww::Vector(0.f, 0.97014, -0.24254), true);
}

//------------------------------------------------------------------------------
// Scenario: Check that a sphere is a shape
//           This is part of the checklist of chapter 9
TEST(DISABLED_Ch9Planes, CheckThatASphereIsAShape)
{
  ww::sphere const S{};
  ww::shared_ptr_shape ptrShape = ww::SharedPtrShape(S);

  EXPECT_EQ(ptrShape->isA<ww::sphere>(), false);
  EXPECT_EQ(ptrShape->isA<ww::shape>(), true);

  std::shared_ptr<ww::sphere> ptrSphere2 = ww::SharedPtrSh<ww::sphere>(S);
  EXPECT_EQ(ptrSphere2->isA<ww::sphere>(), true);
  EXPECT_EQ(ptrSphere2->isA<ww::shape>(), true);
}

//------------------------------------------------------------------------------
// Scenario: The normal of a plane is constant everywhere
TEST(DISABLED_Ch9Planes, CheckThanNormalOnPlaneIsConstant)
{
  ww::plane P{};
  P.funcPtrLocalNormalAt = ww::LocalNormalAtPlane;
  ww::tup N1 = ww::LocalNormalAtPlane(P, ww::Point(0.f, 0.f, 0.f));
  ww::tup N2 = ww::LocalNormalAtPlane(P, ww::Point(10.f, 0.f, -10.f));
  ww::tup N3 = ww::LocalNormalAtPlane(P, ww::Point(-5.f, 0.f, 150.f));
  EXPECT_EQ(N1 == ww::Vector(0.f, 1.f, 0.f), true);
  EXPECT_EQ(N2 == ww::Vector(0.f, 1.f, 0.f), true);
  EXPECT_EQ(N3 == ww::Vector(0.f, 1.f, 0.f), true);
}

//------------------------------------------------------------------------------
// Scenario: Test that the shape's local NormalAt can be assigned.
TEST(DISABLED_Ch9Planes, CheckThatNormaFuncPointerOfShapeCanBeAssigned)
{
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ww::plane &P = *ptrPlane;
  ww::tup N1 = ww::NormalAt(P, ww::Point(0.f, 0.f, 0.f));
  ww::tup N2 = ww::NormalAt(P, ww::Point(10.f, 0.f, -10.f));
  ww::tup N3 = ww::NormalAt(P, ww::Point(-5.f, 0.f, 150.f));
  ww::tup N4 = ww::NormalAt(P, ww::Point(0.f, 0.f, 0.f));

  EXPECT_EQ(N1 == ww::Vector(0.f, 1.f, 0.f), true);
  EXPECT_EQ(N2 == ww::Vector(0.f, 1.f, 0.f), true);
  EXPECT_EQ(N3 == ww::Vector(0.f, 1.f, 0.f), true);
  EXPECT_EQ(N4 == ww::Vector(0.f, 1.f, 0.f), true);
}

//------------------------------------------------------------------------------
// Scenario: Intersect with a ray parallel to the plane.
TEST(DISABLED_Ch9Planes, IntersectRayParallelToPlane)
{
  std::shared_ptr<ww::plane> ptrPlane = ww::PtrDefaultPlane();
  ww::ray R = ww::Ray(ww::Point(0.f, 10.f, 0.f), ww::Vector(0.f, 0.f, 1.f));
  ww::intersections XS = ww::LocalIntersect(ptrPlane, R);

  EXPECT_EQ(ptrPlane->isA<ww::plane>(), true);
  EXPECT_EQ(ptrPlane->isA<ww::sphere>(), false);
  EXPECT_EQ(XS.Count(), 0);
}

//------------------------------------------------------------------------------
// Scenario: Intersect with a coplanar ray.
TEST(DISABLED_Ch9Planes, IntersectWithACoplanarRay)
{
  std::shared_ptr<ww::plane> ptrPlane = ww::PtrDefaultPlane();
  ww::ray R = ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 0.f, 1.f));
  ww::intersections XS = ww::LocalIntersect(ptrPlane, R);
  EXPECT_EQ(XS.Count(), 0);
}

//------------------------------------------------------------------------------
// Scenario: Intersecting a plane from above.
TEST(DISABLED_Ch9Planes, IntersectAPlaneFromAbove)
{
  std::shared_ptr<ww::plane> ptrPlane = ww::PtrDefaultPlane();
  ww::ray R = ww::Ray(ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, -1.f, 0.f));
  ww::intersections XS = ww::LocalIntersect(ptrPlane, R);

  EXPECT_EQ(XS.Count(), 1);
  if (XS.Count())
  {
    EXPECT_EQ(XS.vI[0].t, 1.f);
    EXPECT_EQ(XS.vI[0].pShape == ptrPlane, true);
  }
}

//------------------------------------------------------------------------------
// Scenario: A ray Intersecting a plane from below.
TEST(DISABLED_Ch9Planes, IntersectAPlaneFromBelow)
{
  std::shared_ptr<ww::plane> ptrPlane = ww::PtrDefaultPlane();
  ww::ray R = ww::Ray(ww::Point(0.f, -1.f, 0.f), ww::Vector(0.f, 1.f, 0.f));
  ww::intersections XS = ww::LocalIntersect(ptrPlane, R);

  EXPECT_EQ(XS.Count(), 1);
  if (XS.Count())
  {
    EXPECT_EQ(XS.vI[0].t, 1.f);
    EXPECT_EQ(XS.vI[0].pShape == ptrPlane, true);
  }
}

//------------------------------------------------------------------------------
TEST(DISABLED_Ch9Planes, PuttingItTogether)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  // ww::shared_ptr_shape PtrFloor = ww::PtrDefaultSphere();
  // ww::shape &Floor = *PtrFloor;
  // Floor.Transform = ww::Scaling(10.f, 0.01f, 10.f);
  // Floor.Material.Color = ww::Color(1.f, 0.9f, 0.9f);
  // Floor.Material.Specular = 0.f;
  // World.vPtrObjects.push_back(PtrFloor);

  // ww::shared_ptr_shape PtrLeftWall = ww::PtrDefaultSphere();
  // ww::shape &LeftWall = *PtrLeftWall;
  // // NOTE: The order in which the transformations are multiplied;
  // //       The wall need to be scaled, rotated in X, rotated in Y,
  // //       and lastly translated, so the transformations are multiplied
  // //       in reverse order.
  // LeftWall.Transform = ww::Translation(0.f, 0.f, 5.f) *  //!<
  //                      ww::RotateY(-M_PI_4) *            //!<
  //                      ww::RotateX(M_PI_2) *             //!<
  //                      ww::Scaling(10.f, 0.01f, 10.f);
  // LeftWall.Material = Floor.Material;
  // World.vPtrObjects.push_back(PtrLeftWall);

  // ww::shared_ptr_shape PtrRightWall = ww::PtrDefaultSphere();
  // ww::shape &RightWall = *PtrRightWall;
  // RightWall.Transform = ww::Translation(0.f, 0.f, 5.f) *  //!<
  //                       ww::RotateY(M_PI_4) *             //!<
  //                       ww::RotateX(M_PI_2) *             //!<
  //                       ww::Scaling(10.f, 0.01f, 10.f);
  // RightWall.Material = Floor.Material;
  // World.vPtrObjects.push_back(PtrRightWall);

  ww::shared_ptr_shape PtrRight = ww::PtrDefaultSphere();
  ww::shape &Right = *PtrRight;
  Right.Transform = ww::Translation(1.5f, 0.5f, -0.5f) *  //!<
                    ww::Scaling(0.5f, 0.5f, 0.5f);
  Right.Material.Color = ww::Color(0.5f, 1.0f, 0.1f);
  Right.Material.Diffuse = 0.7f;
  Right.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrRight);

  ww::shared_ptr_shape PtrMiddle = ww::PtrDefaultSphere();
  ww::shape &Middle = *PtrMiddle;
  Middle.Transform = ww::Translation(-0.5f, 1.f, 0.5f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrMiddle);

  ww::shared_ptr_shape PtrLeft = ww::PtrDefaultSphere();
  ww::shape &Left = *PtrLeft;
  Left.Transform = ww::Translation(-1.5f, 0.33f, -0.75f) *  //!<
                   ww::Scaling(0.33f, 0.33f, 0.33f);
  Left.Material.Color = ww::Color(1.0f, 0.8f, 0.1f);
  Left.Material.Diffuse = 0.7f;
  Left.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrLeft);

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::Translation(0.f, 0.f, 0.f)  //!<
                                                        // * ww::RotateX(0 * M_PI_2)       //!<
      ;                                                 //!<
  ptrPlane->Material.Shininess = 10.f;
  ptrPlane->Material.Diffuse = 0.5f;
  ptrPlane->Material.Color = ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  World.vPtrObjects.push_back(ptrPlane);

  // Add a second plane - this will act as the backdrop.
  ww::shared_ptr_plane ptrPlane2 = ww::PtrDefaultPlane();
  ptrPlane2->Transform = ww::Translation(0.f, 0.f, 10.f)  //!<
                                                          // ww::RotateY(M_PI_4) *             //!<
                         * ww::RotateX(1 * M_PI_2)        //!<
      ;                                                   //!<
  ptrPlane2->Material.Shininess = 100.f;
  ptrPlane2->Material.Diffuse = 0.7f;
  ptrPlane2->Material.Color =
      ww::Color(float(0x2f) / float(0xff), float(0xb5) / float(0xff), float(0xff) / float(0xff));
  World.vPtrObjects.push_back(ptrPlane2);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  ww::camera Camera = ww::Camera(256, 256, M_PI / 5.f);

  ww::tup const ViewFrom = ww::Point(0.f, 1.5f, -15.f);
  ww::tup const ViewTo = ww::Point(0.f, 1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch9MakingAScene.ppm");
}

//------------------------------------------------------------------------------
TEST(DISABLED_Ch9Planes, PuttingItTogetherHexagon)
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
  World.vPtrObjects.push_back(PtrRight);

  ww::shared_ptr_shape PtrMiddle = ww::PtrDefaultSphere();
  ww::shape &Middle = *PtrMiddle;
  Middle.Transform = ww::Translation(-0.5f, 1.f, 0.5f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrMiddle);

  ww::shared_ptr_shape PtrLeft = ww::PtrDefaultSphere();
  ww::shape &Left = *PtrLeft;
  Left.Transform = ww::Translation(-1.5f, 0.33f, -0.75f) *  //!<
                   ww::Scaling(0.33f, 0.33f, 0.33f);
  Left.Material.Color = ww::Color(1.0f, 0.8f, 0.1f);
  Left.Material.Diffuse = 0.7f;
  Left.Material.Specular = 0.3f;
  World.vPtrObjects.push_back(PtrLeft);

  // Add the first plane - the floor
  ww::shared_ptr_plane ptrFloor = ww::PtrDefaultPlane();
  ptrFloor->Transform = ww::Translation(0.f, 0.f, 0.f)  //!<
                                                        // * ww::RotateX(0 * M_PI_2)       //!<
      ;                                                 //!<
  ptrFloor->Material.Shininess = 10.f;
  ptrFloor->Material.Diffuse = 0.5f;
  ptrFloor->Material.Color = ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  World.vPtrObjects.push_back(ptrFloor);

  // Assert(0 == 1, __FUNCTION__, __LINE__);

  auto HexagonSide = [](float X, float Z, float Angle) -> ww::shared_ptr_plane
  {
    ww::shared_ptr_plane ptrPlaneHexagon = ww::PtrDefaultPlane();
    ptrPlaneHexagon->Transform = ww::Translation(X, 0.f, Z)  //!<
                                 * ww::RotateX(M_PI_2)       //!<
                                 * ww::RotateY(Angle)        //!<
        ;                                                    //!<
    return ptrPlaneHexagon;
  };

  // float const COS_PI_3 = std::cosf(M_PI / 3.f);
  // float const SIN_PI_3 = std::sinf(M_PI / 3.f);

  World.vPtrObjects.push_back(HexagonSide(0.f, 10.f, 0.f * M_PI / 2.f));

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  ww::camera Camera = ww::Camera(256, 256, M_PI / 2.f);

  ww::tup const ViewFrom = ww::Point(0.f, 10.f, -10.f);
  ww::tup const ViewTo = ww::Point(0.f, 0.0f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch9MakingASceneHexagon.ppm");
}
