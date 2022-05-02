#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
//------------------------------------------------------------------------------
TEST(Ch8Shadows, LightingWithTheSurfaceInShadow)
{
  ww::tup vEye = ww::Vector(0.f, 0.f, -1.f);
  ww::tup vNormal = ww::Vector(0.f, 0.f, -1.f);
  ww::light Light = ww::PointLight(ww::Point(0.f, 0.f, -10.f), ww::Color(1.f, 1.f, 1.f));
  ww::material M{};
  ww::tup Position = ww::Point(0.f, 0.f, 0.f);
  bool const InShadow{true};

  ww::tup Result = ww::Lighting(M, Light, Position, vEye, vNormal, InShadow);

  EXPECT_EQ(Result == ww::Color(0.1f, 0.1f, 0.1f), true);
}

//------------------------------------------------------------------------------
TEST(Ch8Shadows, ThereIsNoShadowWhenNothingIsCollinearWithPointAndLight)
{
  ww::world const W = ww::World();
  ww::tup const Point = ww::Point(0.f, 10.f, 0.f);
  EXPECT_EQ(ww::IsShadowed(W, Point), false);
}

//------------------------------------------------------------------------------
TEST(Ch8Shadows, TheShadowWhenAnObjectIsBetweenThePointAndLight)
{
  ww::world const W = ww::World();
  ww::tup const Point = ww::Point(10.f, -10.f, 10.f);
  EXPECT_EQ(ww::IsShadowed(W, Point), true);
}

//------------------------------------------------------------------------------
TEST(Ch8Shadows, ThereIsNoShadowWhenAnObjectIsBehindTheLight)
{
  ww::world const W = ww::World();
  ww::tup const Point = ww::Point(-20.f, 20.f, -20.f);
  EXPECT_EQ(ww::IsShadowed(W, Point), false);
}

//------------------------------------------------------------------------------
TEST(Ch8Shadows, ThereIsNoShadowWhenAnObjectIsBehindThePoint)
{
  ww::world const W = ww::World();
  ww::tup const Point = ww::Point(-2.f, 2.f, -2.f);
  EXPECT_EQ(ww::IsShadowed(W, Point), false);
}

//------------------------------------------------------------------------------
TEST(Ch8Shadows, ShadeHitIsGivenAnIntersectionInShadow)
{
  ww::world W = ww::World();
  ww::shared_ptr_light PtrLight = W.vPtrLights[0];
  ww::light &L = *PtrLight;
  L = ww::PointLight(ww::Point(0.f, 0.f, -10.f), ww::Color(1.f, 1.f, 1.f));

  ww::shared_ptr_shape PtrS1{};
  PtrS1.reset(new ww::sphere);
  W.vPtrObjects.push_back(PtrS1);

  ww::shared_ptr_shape PtrS2{};
  PtrS2.reset(new ww::sphere);
  PtrS2->Transform = ww::Translation(0.f, 0.f, 10.f);
  W.vPtrObjects.push_back(PtrS2);

  ww::ray R = ww::Ray(ww::Point(0.f, 0.f, 5.f), ww::Vector(0.f, 0.f, 1.f));

  ww::intersection I = ww::Intersection(4.f, PtrS2);

  ww::prepare_computation Comps = ww::PrepareComputations(I, R);

  ww::tup C = ww::ShadeHit(W, Comps);
  EXPECT_EQ(C == ww::Color(0.1f, 0.1f, 0.1f), true);
}

//------------------------------------------------------------------------------
TEST(Ch8Shadows, TheHitShouldOffsetThePoint)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape PtrShape{};
  PtrShape.reset(new ww::sphere);

  PtrShape->Transform = ww::Translation(0.f, 0.f, 1.f);
  ww::intersection I = ww::Intersection(5.f, PtrShape);
  ww::prepare_computation const Comps = ww::PrepareComputations(I, R);

  // NOTE: The point compares to z component to half of -EPSILON to make sure that the point
  //       has been adjusted in the correct direction.
  EXPECT_EQ(Comps.Point.Z < -ww::EPSILON / 2.f, true);
}

//------------------------------------------------------------------------------
TEST(Ch8Shadows, PuttingItTogether)
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

  ww::shared_ptr_shape PtrMiddleHigh = ww::PtrDefaultSphere();
  ww::shape &MiddleHigh = *PtrMiddleHigh;
  MiddleHigh.Transform = ww::Translation(0.0f, 1.75f, -2.0f) *  //!<
                         ww::RotateZ(5.f * M_PI_2 / 4.f) *      //!<
                         ww::RotateY(3.f * M_PI_2 / 4.f) *      //!<
                         ww::RotateX(1.f * M_PI_2 / 4.f) *      //!<
                         ww::Scaling(0.9f, 0.9f, 0.2f);
  MiddleHigh.Material.Color = ww::Color(0.8f, 0.999f, 0.f);
  MiddleHigh.Material.Diffuse = 0.99f;
  MiddleHigh.Material.Specular = 0.7f;
  World.vPtrObjects.push_back(PtrMiddleHigh);

  ww::shared_ptr_shape PtrRightHigh = ww::PtrDefaultSphere();
  ww::shape &RightHigh = *PtrRightHigh;
  RightHigh.Transform = ww::Translation(-1.0f, 1.0f, -3.5f) *  //!<
                        ww::RotateZ(5.f * M_PI_2 / 4.f) *      //!<
                        ww::RotateY(3.f * M_PI_2 / 4.f) *      //!<
                        ww::RotateX(1.f * M_PI_2 / 7.f) *      //!<
                        ww::Scaling(0.9f, 0.4f, 0.4f);
  RightHigh.Material.Color = ww::Color(0.1f, 0.999f, 8.f);
  RightHigh.Material.Diffuse = 0.99f;
  RightHigh.Material.Specular = 0.7f;
  World.vPtrObjects.push_back(PtrRightHigh);

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
  Camera.Transform = ww::ViewTransform(ww::Point(0.f, 1.5f, -5.f), ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, 1.f, 0.f));

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch8MakingAScene.ppm");
}

