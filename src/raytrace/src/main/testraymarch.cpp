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
TEST(RayMarch, GetDistanceBox)
{
  ww::shared_ptr_shape pDefaultBox = ww::PtrDefaultCube();

  {
    ww::tup const P1 = ww::Vector(0.f, 0.f, -2.f);
    float D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 1.5f);
  }

  {
    ww::tup const P1 = ww::Vector(0.f, 0.f, 2.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 1.5f);
  }

  {
    ww::tup const P1 = ww::Vector(0.f, 2.f, 0.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 1.5f);
  }

  {
    ww::tup const P1 = ww::Vector(0.f, -2.f, 0.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 1.5f);
  }

  {
    ww::tup const P1 = ww::Vector(-2.f, 0.f, 0.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 1.5f);
  }

  {
    ww::tup const P1 = ww::Vector(2.f, 0.f, 0.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 1.5f);
  }

  {
    ww::tup const P1 = ww::Vector(1.f, 1.f, 0.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(1.f / 2.f));
  }

  {
    ww::tup const P1 = ww::Vector(-1.f, 1.f, 0.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(1.f / 2.f));
  }

  // ---
  // NOTE: The length is Sqrt(1ˆ2 + 1ˆ2 + 1ˆ2) - Sqrt(0.5ˆ2 + 0.5ˆ2 +0.5ˆ2)
  //       since the length corresponds to a line going from 0,0,0 to 1,1,1
  //       via the point 0.5, 0.5, 0.5 (the corner of the box.)
  // ---
  {
    ww::tup const P1 = ww::Vector(1.f, 1.f, 1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(1.f, 1.f, -1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(1.f, -1.f, 1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(1.f, -1.f, -1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(-1.f, -1.f, -1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(-1.f, -1.f, 1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(-1.f, 1.f, 1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(-1.f, 1.f, -1.f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, std::sqrtf(3.f) - std::sqrtf(3.f / 4.f));
  }
  {
    ww::tup const P1 = ww::Vector(0.5f, 0.5f, 0.5f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 0.f);
  }
  {
    ww::tup const P1 = ww::Vector(-0.5f, -0.5f, -0.5f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, 0.f);
  }
  {  //!< Inside the box the value should be negative.
    ww::tup const P1 = ww::Vector(-0.25f, .0f, .0f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, -0.25f);
  }
  {  //!< Inside the box the value should be negative.
    ww::tup const P1 = ww::Vector(.0f, -.25f, .0f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, -0.25f);
  }
  {  //!< Inside the box the value should be negative.
    ww::tup const P1 = ww::Vector(.0f, .0f, .25f);
    float const D = ww::rm::GetDistance(P1, pDefaultBox);
    EXPECT_FLOAT_EQ(D, -0.25f);
  }
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
TEST(DISABLED_RayMarch, CalcNormal)
{
  auto TestNormal = [](ww::ray const &R, ww::funcPtrMap Map)
  {
    ww::tup const N = ww::rm::CalcNormal(R.Origin, Map);
    EXPECT_FLOAT_EQ(ww::Mag(N), 1.0f);
    EXPECT_LT(ww::Mag(N - R.Direction), 0.01f);
  };
  TestNormal(ww::Ray(ww::Point(0.f, 0.f, -1.f), ww::Vector(0.f, 0.f, -1.f)), ww::rm::MapBoxAndSphere);
  TestNormal(ww::Ray(ww::Point(0.f, 0.f, 1.f), ww::Vector(0.f, 0.f, 1.f)), ww::rm::MapBoxAndSphere);
  TestNormal(ww::Ray(ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, 1.f, 0.f)), ww::rm::MapBoxAndSphere);
  TestNormal(ww::Ray(ww::Point(0.f, -1.f, 0.f), ww::Vector(0.f, -1.f, 0.f)), ww::rm::MapBoxAndSphere);
  TestNormal(ww::Ray(ww::Point(1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f)), ww::rm::MapBoxAndSphere);
  TestNormal(ww::Ray(ww::Point(-1.f, 0.f, 0.f), ww::Vector(-1.f, 0.f, 0.f)), ww::rm::MapBoxAndSphere);
}

//------------------------------------------------------------------------------
TEST(RayMarch, RayMarch)
{
  auto TestDistance = [](ww::ray const &R, ww::shared_ptr_shape pDefaultSphere) -> float
  {
    ww::world World{};
    World.vPtrObjects.push_back(pDefaultSphere);

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
    ww::world World{};
    World.vPtrObjects.push_back(pDefaultSphere);
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
TEST(DISABLED_RayMarch, MainImage)
{
  return;
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
TEST(DISABLED_RayMarch, Test1)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();
  World.Map = ww::rm::MapDefault;

  // ---
  // Add Sphere
  // ---
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(1.f, 1.f, .0f);

  ww::shared_ptr_shape pDefaultSphere = ww::PtrDefaultSphere();
  ww::shape &S = *pDefaultSphere;
  S.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.0f, 1.0f, 1.0f, 0.f, 0.f, 0.f);
  S.Material.Color = ww::Color(0.0f, 0.9f, 0.0f);  //!< "normal white" material should be around 0.2 gray
  // S.Material.Diffuse = 0.7f;
  // S.Material.Specular = 0.f;
  S.Material.Ambient = 1.f;

  ww::shared_ptr_shape pDefaultSphere2 = ww::PtrDefaultSphere();
  ww::shape &S2 = *pDefaultSphere2;
  S2.Transform = ww::TranslateScaleRotate(-2.f, 0.f, 0.f, 1.1f, 1.1f, 1.1f, 0.f, 0.f, 0.f);
  // S2.Material.Color = ww::Color(0.0f, 0.2f, 0.5f) * 2.f;
  S2.Material.Color = ww::Color(0.f, 0.f, 1.f);
  S2.Material.Diffuse = 0.0f;
  S2.Material.Specular = 0.0f;
  // S2.Print = true;

  ww::shared_ptr_shape pDefaultBox = ww::PtrDefaultCube();
  pDefaultBox->Transform =
      ww::TranslateScaleRotate(2.f, 0.f, -1.f, 1.f, 1.f, 1.f, 1.f * ww::Radians(45.f), -1.f * ww::Radians(45.f), 0.f);
  pDefaultBox->Material.Color = ww::Color(0.7f, 0.4f, 0.8f);
  ww::cube *pCube = dynamic_cast<ww::cube *>(pDefaultBox.get());
  pCube->R = 0.1;

  ww::cube Box2{};
  Box2 = *pCube;
  Box2.R = 0.05;
  Box2.Transform =
      ww::TranslateScaleRotate(0.f, 0.f, 2.f, 1.f, 2.f, 1.f, 1.f * ww::Radians(45.f), -1.f * ww::Radians(45.f), 0.f);
  Box2.Material.Color = ww::Color(0.7f, 0.4f, 0.8f);
  ww::shared_ptr_shape pBox2 = ww::PtrDefaultCube();
  *pBox2 = Box2;

  ww::shared_ptr_shape pFloor = ww::PtrDefaultPlane();
  pFloor->Material.Color = ww::Color(0.4f, 0.4f, 0.4f);
  ww::plane *pFloorRaw = dynamic_cast<ww::plane *>(pFloor.get());
  pFloorRaw->H = 2.f;

  World.vPtrObjects.push_back(pDefaultBox);
  World.vPtrObjects.push_back(pDefaultSphere);
  World.vPtrObjects.push_back(pDefaultSphere2);
  World.vPtrObjects.push_back(pBox2);
  World.vPtrObjects.push_back(pFloor);
  // ---
  // We need some light, please.
  // ---
  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  // *pLight = ww::PointLight(ww::Point(-5.f, 25.f, 2.5f), ww::Color(1.f, 1.f, 1.f));
  *pLight = ww::PointLight(ww::Point(-6.f, 10.f, -0.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(512, 512, ww::Radians(50.f));
  // Camera.RenderSingleThread = true;

  ww::tup const ViewFrom = ww::Point(0.f, 7.0f, -7.f);
  ww::tup const ViewTo = ww::Point(0.0f, 0.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);

  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::rm::Render(Camera, World);
  ww::WriteToPPM(Canvas, "RayMarchTest1.ppm");
}

//------------------------------------------------------------------------------
// Scenario:
TEST(RayMarch, Test2)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();
  World.Map = ww::rm::MapBoxAndSphere;
  // World.Map = ww::rm::MapDefault;

  // ---
  // We need some light, please.
  // ---
  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  // *pLight = ww::PointLight(ww::Point(-5.f, 25.f, 2.5f), ww::Color(1.f, 1.f, 1.f));
  *pLight = ww::PointLight(ww::Point(-6.f, 10.f, -0.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(512, 512, ww::Radians(50.f));
  // Camera.RenderSingleThread = true;

  ww::tup const ViewFrom = ww::Point(0.f, 7.0f, -7.f);
  ww::tup const ViewTo = ww::Point(0.0f, 0.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);

  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::rm::Render(Camera, World);
  ww::WriteToPPM(Canvas, "RayMarchTest2.ppm");
}

//------------------------------------------------------------------------------
TEST(RayMarch, TestSignedDistanceFunctionsBox1)
{
  ww::tup BoxPos = ww::Vector(0.f, 0.f, 0.f);
  ww::tup BoxSize = ww::Vector(1.f, 1.f, 1.f);
  ww::tup BoxRot = ww::Vector(0.f, 0.f, 0.f);
  ww::tup P1 = ww::Point(2.f, 0.f, 0.f);
  ww::tup P2 = ww::Point(-2.f, 0.f, 0.f);
  ww::tup P3 = ww::Point(0.f, 2.f, 0.f);
  ww::tup P4 = ww::Point(0.f, -2.f, 0.f);
  ww::tup P5 = ww::Point(0.f, 0.f, 0.f);
  ww::tup P6 = ww::Point(0.f, 0.f, 0.f);

  ww::matrix MBox = ww::TranslateScaleRotate(BoxPos.X, BoxPos.Y, BoxPos.Z,     //!<
                                             BoxSize.X, BoxSize.Y, BoxSize.Z,  //!<
                                             BoxRot.X, BoxRot.Y, BoxRot.Z      //!<
  );
  ww::matrix MBoxI = ww::Inverse(MBox);
  float D1 = ww::rm::sdBox(ww::Vector(MBoxI * P1), BoxSize);
  float D2 = ww::rm::sdBox(ww::Vector(MBoxI * P2), BoxSize);
  float D3 = ww::rm::sdBox(ww::Vector(MBoxI * P3), BoxSize);
  float D4 = ww::rm::sdBox(ww::Vector(MBoxI * P4), BoxSize);
  float D5 = ww::rm::sdBox(ww::Vector(MBoxI * P5), BoxSize);

  // ---
  // NOTE: The Capsule is initially set up with A=1,0,0 and B=1,0,0.
  //       This means that the distance between A and B is of length 2.
  //       And that the overall length is 4 since the radius is 1.
  //       So when the capsule is located at 10,0,0 the distance to
  //       Origo will be 8.
  // ---
  ww::matrix MCapsule = ww::TranslateScaleRotate(10.f, 0.f, 0.f,  //!<
                                                 1.f, 1.f, 1.f,   //!<
                                                 0.f, 0.f, 0.f    //!<
  );
  ww::matrix MCapsuleI = ww::Inverse(MCapsule);
  float D6 = ww::rm::sdCapsule(Vector(MCapsuleI * P6), ww::Vector(-1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f), 1.f);
  EXPECT_FLOAT_EQ(D6, 8.f);
  // std::cout << "Distance to Capsule from point " << P6 << " is " << D6 << std::endl;

  // ---
  // NOTE: So now the Capsule is rotated around the Z-axis and the distance
  //       should increase to reflect that distance increases by the radius.
  // ---
  MCapsule = ww::TranslateScaleRotate(10.f, 0.f, 0.f,              //!<
                                      1.f, 1.f, 1.f,               //!<
                                      0.f, 0.f, ww::Radians(90.f)  //!<
  );
  MCapsuleI = ww::Inverse(MCapsule);
  float D7 = ww::rm::sdCapsule(Vector(MCapsuleI * P6), ww::Vector(-1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f), 1.f);
  EXPECT_FLOAT_EQ(D7, 9.f);

  auto CheckDist = [](ww::matrix const &M, ww::matrix const &MI, ww::tup const &P, ww::tup const &Box, float D,
                      float DExpected) -> void
  {
    // std::cout << "MBox :" << M << std::endl;
    // std::cout << "MBoxI:" << MI << std::endl;
    // std::cout << "P    :" << P << std::endl;
    // std::cout << "Pmov :" << MI * P << std::endl;
    // std::cout << "-> Distance: " << D << " to BOX:  " << Box << std::endl;
    EXPECT_FLOAT_EQ(D, DExpected);
  };

  CheckDist(MBox, MBoxI, P1, BoxPos, D1, 1.f);
  CheckDist(MBox, MBoxI, P2, BoxPos, D2, 1.f);
  CheckDist(MBox, MBoxI, P3, BoxPos, D3, 1.f);
  CheckDist(MBox, MBoxI, P4, BoxPos, D4, 1.f);
  CheckDist(MBox, MBoxI, P5, BoxPos, D5, -1.f);

  // For sanity check from the book.
  // ww::matrix Transform = ww::Translation(5.f, -3.f, 2.f);
  // ww::tup P = ww::Point(-3.f, 4.f, 5.f);
  // ww::tup PT = Transform * P;
  // ww::tup InvPT = ww::Inverse(Transform) * PT;
  // std::cout << "Translated point PT: " << PT << std::endl;
  // std::cout << "Translated point InvPT: " << InvPT << std::endl;
  // ww::matrix Inv = ww::Inverse(Transform);
  // std::cout << "Translated point InvPT2: " << Inv * P << std::endl;
}

//------------------------------------------------------------------------------
TEST(RayMarch, TestSignedDistanceFunctionsBox2)
{
  // char const *pTestDescription = "TestSignedDistanceFunctionsBox2";
  ww::tup BoxPos = ww::Vector(15.f, 10.f, 0.f);
  ww::tup BoxSiz = ww::Vector(5.f, 1.f, 1.f);
  ww::tup BoxRot = ww::Vector(0.f, 0.f, 0.f);

  ww::tup P1 = ww::Point(0.f, 0.f, 0.f);
  ww::tup P2 = ww::Point(30.f, 0.f, 0.f);
  ww::tup P3 = ww::Point(15.f, 0.f, 0.f);
  ww::tup P4 = ww::Point(0.f, 10.f, 0.f);
  ww::tup P5 = ww::Point(30.f, 10.f, 0.f);
  ww::tup P6 = ww::Point(0.f, 0.f, 0.f);

  ww::matrix MBox = ww::TranslateScaleRotate(BoxPos.X, BoxPos.Y, BoxPos.Z,  //!<
                                             1.f, 1.f, 1.f,                 //!<
                                             BoxRot.X, BoxRot.Y, BoxRot.Z   //!<
  );
  ww::matrix MBoxI = ww::Inverse(MBox);
  float D1 = ww::rm::sdBox(ww::Vector(MBoxI * P1), BoxSiz);
  float D2 = ww::rm::sdBox(ww::Vector(MBoxI * P2), BoxSiz);
  float D3 = ww::rm::sdBox(ww::Vector(MBoxI * P3), BoxSiz);
  float D4 = ww::rm::sdBox(ww::Vector(MBoxI * P4), BoxSiz);
  float D5 = ww::rm::sdBox(ww::Vector(MBoxI * P5), BoxSiz);

  // ---
  // NOTE: The Capsule is initially set up with A=1,0,0 and B=1,0,0.
  //       This means that the distance between A and B is of length 2.
  //       And that the overall length is 4 since the radius is 1.
  //       So when the capsule is located at 10,0,0 the distance to
  //       Origo will be 8.
  // ---
  ww::matrix MCapsule = ww::TranslateScaleRotate(10.f, 0.f, 0.f,  //!<
                                                 1.f, 1.f, 1.f,   //!<
                                                 0.f, 0.f, 0.f    //!<
  );
  ww::matrix MCapsuleI = ww::Inverse(MCapsule);
  float D6 = ww::rm::sdCapsule(Vector(MCapsuleI * P6), ww::Vector(-1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f), 1.f);
  EXPECT_FLOAT_EQ(D6, 8.f);
  // std::cout << pTestDescription << " -> Distance to Capsule from point " << P6 << " is " << D6 << std::endl;

  // ---
  // NOTE: So now the Capsule is rotated around the Z-axis and the distance
  //       should increase to reflect that distance increases by the radius.
  // ---
  MCapsule = ww::TranslateScaleRotate(10.f, 0.f, 0.f,              //!<
                                      1.f, 1.f, 1.f,               //!<
                                      0.f, 0.f, ww::Radians(90.f)  //!<
  );
  MCapsuleI = ww::Inverse(MCapsule);
  float D7 = ww::rm::sdCapsule(Vector(MCapsuleI * P6), ww::Vector(-1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f), 1.f);
  EXPECT_FLOAT_EQ(D7, 9.f);

  auto CheckDist = [&](ww::matrix const &M, ww::matrix const &MI, ww::tup const &P, ww::tup const &Box, float D,
                       float DExpected) -> void
  {
    // std::cout << "---------------------------------------------------" << std::endl;
    // std::cout << pTestDescription << ". Line:" << __LINE__ << std::endl;
    // std::cout << "MBox :" << M << std::endl;
    // std::cout << "MBoxI:" << MI << std::endl;
    //
    // std::cout << "Box Center: " << Box << std::endl;
    // std::cout << "P         : " << P << std::endl;
    // std::cout << "Pmov      : " << MI * P << std::endl;
    // std::cout << "Distance  : " << D << std::endl;

    EXPECT_FLOAT_EQ(D, DExpected);
  };

  CheckDist(MBox, MBoxI, P1, BoxPos, D1, std::sqrtf(10.f * 10.f + 9.f * 9.f));
  CheckDist(MBox, MBoxI, P2, BoxPos, D2, std::sqrtf(10.f * 10.f + 9.f * 9.f));
  CheckDist(MBox, MBoxI, P3, BoxPos, D3, 9.f);
  CheckDist(MBox, MBoxI, P4, BoxPos, D4, 10.f);
  CheckDist(MBox, MBoxI, P5, BoxPos, D5, 10.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, TestSignedDistanceFunctionsBox3)
{
  // char const *pTestDescription = "TestSignedDistanceFunctionsBox3";
  ww::tup BoxPos = ww::Vector(15.f, 10.f, 0.f);
  ww::tup BoxSiz = ww::Vector(5.f, 1.f, 2.f);
  ww::tup BoxRot = ww::Vector(ww::Radians(90.f), 0.f, 0.f);

  ww::tup P1 = ww::Point(0.f, 0.f, 0.f);
  ww::tup P2 = ww::Point(30.f, 0.f, 0.f);
  ww::tup P3 = ww::Point(15.f, 0.f, 0.f);
  ww::tup P4 = ww::Point(0.f, 10.f, 0.f);
  ww::tup P5 = ww::Point(30.f, 10.f, 0.f);
  ww::tup P6 = ww::Point(0.f, 0.f, 0.f);

  ww::matrix MBox = ww::TranslateScaleRotate(BoxPos.X, BoxPos.Y, BoxPos.Z,  //!<
                                             1.f, 1.f, 1.f,                 //!<
                                             BoxRot.X, BoxRot.Y, BoxRot.Z   //!<
  );
  ww::matrix MBoxI = ww::Inverse(MBox);
  float D1 = ww::rm::sdBox(ww::Vector(MBoxI * P1), BoxSiz);
  float D2 = ww::rm::sdBox(ww::Vector(MBoxI * P2), BoxSiz);
  float D3 = ww::rm::sdBox(ww::Vector(MBoxI * P3), BoxSiz);
  float D4 = ww::rm::sdBox(ww::Vector(MBoxI * P4), BoxSiz);
  float D5 = ww::rm::sdBox(ww::Vector(MBoxI * P5), BoxSiz);

  // ---
  // NOTE: The Capsule is initially set up with A=1,0,0 and B=1,0,0.
  //       This means that the distance between A and B is of length 2.
  //       And that the overall length is 4 since the radius is 1.
  //       So when the capsule is located at 10,0,0 the distance to
  //       Origo will be 8.
  // ---
  ww::matrix MCapsule = ww::TranslateScaleRotate(10.f, 0.f, 0.f,  //!<
                                                 1.f, 1.f, 1.f,   //!<
                                                 0.f, 0.f, 0.f    //!<
  );
  ww::matrix MCapsuleI = ww::Inverse(MCapsule);
  float D6 = ww::rm::sdCapsule(Vector(MCapsuleI * P6), ww::Vector(-1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f), 1.f);
  EXPECT_FLOAT_EQ(D6, 8.f);
  // std::cout << pTestDescription << " -> Distance to Capsule from point " << P6 << " is " << D6 << std::endl;

  // ---
  // NOTE: So now the Capsule is rotated around the Z-axis and the distance
  //       should increase to reflect that distance increases by the radius.
  // ---
  MCapsule = ww::TranslateScaleRotate(10.f, 0.f, 0.f,              //!<
                                      1.f, 1.f, 1.f,               //!<
                                      0.f, 0.f, ww::Radians(90.f)  //!<
  );
  MCapsuleI = ww::Inverse(MCapsule);
  float D7 = ww::rm::sdCapsule(Vector(MCapsuleI * P6), ww::Vector(-1.f, 0.f, 0.f), ww::Vector(1.f, 0.f, 0.f), 1.f);
  EXPECT_FLOAT_EQ(D7, 9.f);

  auto CheckDist = [&](ww::matrix const &M, ww::matrix const &MI, ww::tup const &P, ww::tup const &Box, float D,
                       float DExpected) -> void
  {
    // std::cout << "---------------------------------------------------" << std::endl;
    // std::cout << pTestDescription << ". Line:" << __LINE__ << std::endl;
    // std::cout << "MBox :" << M << std::endl;
    // std::cout << "MBoxI:" << MI << std::endl;
    //
    // std::cout << "Box Center: " << Box << std::endl;
    // std::cout << "P         : " << P << std::endl;
    // std::cout << "Pmov      : " << MI * P << std::endl;
    // std::cout << "Distance  : " << D << std::endl;

    EXPECT_FLOAT_EQ(D, DExpected);
  };

  CheckDist(MBox, MBoxI, P1, BoxPos, D1, std::sqrtf(10.f * 10.f + 8.f * 8.f));
  CheckDist(MBox, MBoxI, P2, BoxPos, D2, std::sqrtf(10.f * 10.f + 8.f * 8.f));
  CheckDist(MBox, MBoxI, P3, BoxPos, D3, 8.f);
  CheckDist(MBox, MBoxI, P4, BoxPos, D4, 10.f);
  CheckDist(MBox, MBoxI, P5, BoxPos, D5, 10.f);
}

//------------------------------------------------------------------------------
TEST(DISABLED_RayMarch, TestRaymarchPrimitives)
{
  ww::tup const Coordinates = ww::Point(0.f, 0.f, 0.f);
  ww::mainimage_config Cfg{};
  Cfg.Map = ww::rm::MapDefault;
  ww::tup FragColor = ww::rm::MainImage(Coordinates, Cfg);
  ww::tup Res = Cfg.Map(Coordinates);

  ww::tup const &Resolution = Cfg.Resolution;

  EXPECT_EQ(ww::IsPoint(Cfg.Resolution), ww::IsPoint(ww::Point(0.f, 0.f, 0.f)));
  ww::matrix MCoordXform = Matrix22(ww::tup{1.f / Resolution.X, 1.f / Resolution.Y},  //!<
                                    ww::tup{1.f / Resolution.X, 1.f / Resolution.X});

  ww::tup Coord1 = MCoordXform * ww::Vector(0.0f * Resolution.X, 0.0f * Resolution.Y, 0.f)  //!<
                   + ww::Vector(-1.f, -1.f, 0.f);
  ww::tup Coord2 = MCoordXform * ww::Vector(0.5f * Resolution.X, 0.5f * Resolution.Y, 0.f)  //!<
                   + ww::Vector(-1.f, -1.f, 0.f);
  ww::tup Coord3 = MCoordXform * ww::Vector(1.0f * Resolution.X, 1.0f * Resolution.Y, 0.f)  //!<
                   + ww::Vector(-1.f, -1.f, 0.f);

  EXPECT_EQ(Coord1.X, -1.f);
  EXPECT_EQ(Coord1.Y, -1.f);
  EXPECT_EQ(Coord2.X, 0.f);
  EXPECT_EQ(Coord2.Y, 0.f);
  EXPECT_EQ(Coord3.X, 1.f);
  EXPECT_EQ(Coord3.Y, 1.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, TestMap) {}

//------------------------------------------------------------------------------
TEST(RayMarch, TestSphereSphere1)
{
  auto CheckDist = [](ww::tup const &PosA, float RadA, ww::tup const &PosB, float RadB,  //!<
                      float ExpectD, float ExpectPd, bool ExpectC)                       //!<
      -> void
  {
    float const D = ww::rm::sdSphereSphereDistance(PosA, RadA, PosB, RadB);
    float const Pd = ww::rm::sdSphereSpherePenDist(PosA, RadA, PosB, RadB);
    bool const C = ww::rm::sdSphereSphereCollision(PosA, RadA, PosB, RadB);
    // std::cout << "Distance between \nSphA " << PosA << "\nSphB " << PosB << std::endl;
    // std::cout << "is " << D << ". Penetration is " << Pd << " and Collision is " << C << std::endl;
    EXPECT_FLOAT_EQ(D, ExpectD);
    EXPECT_FLOAT_EQ(Pd, ExpectPd);
    EXPECT_EQ(C, ExpectC);
  };

  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(3.f, 0.f, 0.f), 1.f, 1.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(2.f, 0.f, 0.f), 1.f, 0.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(1.f, 0.f, 0.f), 1.f, -1.f, -1.f, true);

  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(0.f, 3.f, 0.f), 1.f, 1.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(0.f, 2.f, 0.f), 1.f, 0.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(0.f, 1.f, 0.f), 1.f, -1.f, -1.f, true);

  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(0.f, 0.f, 3.f), 1.f, 1.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(0.f, 0.f, 2.f), 1.f, 0.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(0.f, 0.f, 1.f), 1.f, -1.f, -1.f, true);

  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(-3.f, 0.f, 0.f), 1.f, 1.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(-2.f, 0.f, 0.f), 1.f, 0.f, 0.f, false);
  CheckDist(ww::Point(0.f, 0.f, 0.f), 1.f, ww::Point(-1.f, 0.f, 0.f), 1.f, -1.f, -1.f, true);
}

//------------------------------------------------------------------------------
TEST(RayMarch, TestTriangle1)
{
  ww::tup V1 = ww::Point(-1.f, 0.f, 0.f);
  ww::tup V2 = ww::Point(0.f, 1.f, 0.f);
  ww::tup V3 = ww::Point(1.f, 0.f, 0.f);

  auto CheckDistTriangle = [&](ww::tup const &Pos, float Expect) -> void
  {
    float D = ww::rm::udTriangle(V1, V2, V3, Pos);
    EXPECT_FLOAT_EQ(D, Expect);
    // std::cout << "Distance from " << Pos << "\n                to triangle is " << D << std::endl;
  };

  CheckDistTriangle(ww::Point(-1.f, 0.f, 0.f), 0.f);
  CheckDistTriangle(ww::Point(0.f, 1.f, 0.f), 0.f);
  CheckDistTriangle(ww::Point(1.f, 0.f, 0.f), 0.f);
  CheckDistTriangle(ww::Point(0.f, 0.f, 0.f), 0.f);
  CheckDistTriangle(ww::Point(0.5f, 0.f, 0.f), 0.f);
  CheckDistTriangle(ww::Point(-0.5f, 0.f, 0.f), 0.f);
  CheckDistTriangle(ww::Point(0.f, 0.f, 1.f), 1.f);
  CheckDistTriangle(ww::Point(1.f, 0.f, 2.f), 2.f);
  CheckDistTriangle(ww::Point(1.f, 1.f, 0.f), 1.f / std::sqrtf(2.f));
  CheckDistTriangle(ww::Point(2.f, 2.f, 0.f), 1.f / std::sqrtf(2.f) + std::sqrtf(2.f));
}

//------------------------------------------------------------------------------
TEST(RayMarch, TestTriangle2)
{
  struct triangle
  {
    ww::tup V1 = ww::Point(-1.f, 0.f, 0.f);
    ww::tup V2 = ww::Point(0.f, 1.f, 0.f);
    ww::tup V3 = ww::Point(1.f, 0.f, 0.f);
  };

  triangle T{};
  // std::cout << "Triangle Pos\n" << T.V1 << "\n" << T.V2 << "\n" << T.V3 << "\n" << std::endl;

  auto CheckDistTriangleSphere = [](ww::tup const &PosSphere, float Radius, triangle const &Triangle,
                                    float Expect) -> void
  {
    auto const DSphereTriangle = ww::rm::sdSphereTriangle(Triangle.V1, Triangle.V2, Triangle.V3, PosSphere, Radius);
    EXPECT_FLOAT_EQ(DSphereTriangle, Expect);

    // auto const DPointTriangle = ww::rm::udTriangle(Triangle.V1, Triangle.V2, Triangle.V3, PosSphere);
    // std::cout << "Distance from Sphere at " << PosSphere << "\n                to triangle is " << DSphereTriangle
    //           << std::endl;
    // std::cout << "Distance from Point at " << PosSphere << "\n                to triangle is " << DPointTriangle
    //           << std::endl;
  };

  CheckDistTriangleSphere(ww::Point(1.f, 1.f, 0.f), 1.f, T, -(1.f - 1.f / std::sqrt(2.f)));

  // Calculate position of sphere:
  // h^2 = 1+1/sqrt(2)
  // h^2 = sqrt(2*x^2) // x and y are of same length.
  // => h^2 = sqrt(2)*x => x = h^2/sqrt(2) => x = (1+1/sqrt(2))/sqrt(2)
  float const SQRT2 = std::sqrt(2.f);
  float const X = (1.f + 1.f / SQRT2) / SQRT2;
  float const Y = X;
  CheckDistTriangleSphere(ww::Point(X, Y, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(-X, Y, 0.f), 1.f, T, 0.f);

  CheckDistTriangleSphere(ww::Point(0.f, 2.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(0.f, -1.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(0.f, -2.f, 0.f), 1.f, T, 1.f);
  CheckDistTriangleSphere(ww::Point(2.f, 0.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(-2.f, 0.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(-1.f, 0.f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(1.f, 0.f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(0.f, 1.f, 1.f), 1.f, T, 0.f);

  // Inside the extrusion of the triangle
  CheckDistTriangleSphere(ww::Point(0.f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(.5f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(-.5f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(0.f, .5f, -1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(.5f, .5f, -1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(ww::Point(-.5f, .5f, -1.f), 1.f, T, 0.f);
}

//------------------------------------------------------------------------------
// Test to confirm a triangle that has moved - see the translation.
TEST(RayMarch, TestTriangle3)
{
  struct triangle
  {
    ww::tup V1 = ww::Point(-1.f, 0.f, 0.f);
    ww::tup V2 = ww::Point(0.f, 1.f, 0.f);
    ww::tup V3 = ww::Point(1.f, 0.f, 0.f);
  };

  triangle T{};
  // std::cout << "Triangle Pos\n" << T.V1 << "\n" << T.V2 << "\n" << T.V3 << "\n" << std::endl;

  ww::matrix M = ww::Translation(10.f, 0.f, 0.f);
  ww::matrix MI = ww::Inverse(M);

  // std::cout << "M\n" << M << std::endl;
  // std::cout << "MI\n" << MI << std::endl;
  // std::cout << "M*MI\n" << M * MI << std::endl;  // sanity check, should become identity matrix.

  T.V1 = M * T.V1;
  T.V2 = M * T.V2;
  T.V3 = M * T.V3;

  auto CheckDistTriangleSphere = [](ww::tup const &PosSphere, float Radius, triangle const &T, float Expect) -> void
  {
    auto const DSphereTriangle = ww::rm::sdSphereTriangle(T.V1, T.V2, T.V3, PosSphere, Radius);
    EXPECT_NEAR(DSphereTriangle, Expect, 1.e-6);

    // auto const DPointTriangle = ww::rm::udTriangle(T.V1, T.V2, T.V3, PosSphere);
    // std::cout << "Triangle Pos\n" << T.V1 << "\n" << T.V2 << "\n" << T.V3 << "\n" << std::endl;
    // std::cout << "Distance from Sphere at " << PosSphere << "\n                to triangle is " << DSphereTriangle
    //           << std::endl;
    // std::cout << "Distance from Point at " << PosSphere << "\n                to triangle is " << DPointTriangle
    //           << std::endl;
  };

  CheckDistTriangleSphere(M * ww::Point(1.f, 1.f, 0.f), 1.f, T, -(1.f - 1.f / std::sqrt(2.f)));
  CheckDistTriangleSphere(ww::Point(11.f, 1.f, 0.f), 1.f, T, -(1.f - 1.f / std::sqrt(2.f)));

  // Calculate position of sphere:
  // h^2 = 1+1/sqrt(2)
  // h^2 = sqrt(2*x^2) // x and y are of same length.
  // => h^2 = sqrt(2)*x => x = h^2/sqrt(2) => x = (1+1/sqrt(2))/sqrt(2)
  float const SQRT2 = std::sqrt(2.f);
  float const X = (1.f + 1.f / SQRT2) / SQRT2;
  float const Y = X;
  CheckDistTriangleSphere(M * ww::Point(X, Y, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-X, Y, 0.f), 1.f, T, 0.f);
  //
  CheckDistTriangleSphere(M * ww::Point(0.f, 2.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, -1.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, -2.f, 0.f), 1.f, T, 1.f);
  CheckDistTriangleSphere(M * ww::Point(2.f, 0.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-2.f, 0.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-1.f, 0.f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(1.f, 0.f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, 1.f, 1.f), 1.f, T, 0.f);

  // Inside the extrusion of the triangle
  CheckDistTriangleSphere(M * ww::Point(0.f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(.5f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-.5f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, .5f, -1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(.5f, .5f, -1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-.5f, .5f, -1.f), 1.f, T, 0.f);
}

//------------------------------------------------------------------------------
// Test to confirm a triangle that has moved - see the translation.
TEST(RayMarch, TestTriangle4)
{
  struct triangle
  {
    ww::tup V1 = ww::Point(-1.f, 0.f, 0.f);
    ww::tup V2 = ww::Point(0.f, 1.f, 0.f);
    ww::tup V3 = ww::Point(1.f, 0.f, 0.f);
  };

  triangle T{};
  // std::cout << "Triangle Pos\n" << T.V1 << "\n" << T.V2 << "\n" << T.V3 << "\n" << std::endl;

  ww::matrix M = ww::Translation(10.f, 10.f, 0.f);
  ww::matrix MI = ww::Inverse(M);

  // std::cout << "M\n" << M << std::endl;
  // std::cout << "MI\n" << MI << std::endl;
  // std::cout << "M*MI\n" << M * MI << std::endl;  // sanity check, should become identity matrix.

  T.V1 = M * T.V1;
  T.V2 = M * T.V2;
  T.V3 = M * T.V3;

  auto CheckDistTriangleSphere = [](ww::tup const &PosSphere, float Radius, triangle const &T, float Expect) -> void
  {
    auto const DSphereTriangle = ww::rm::sdSphereTriangle(T.V1, T.V2, T.V3, PosSphere, Radius);
    EXPECT_NEAR(DSphereTriangle, Expect, 1.e-6);

    // auto const DPointTriangle = ww::rm::udTriangle(T.V1, T.V2, T.V3, PosSphere);
    // std::cout << "Triangle Pos\n" << T.V1 << "\n" << T.V2 << "\n" << T.V3 << "\n" << std::endl;
    // std::cout << "Distance from Sphere at " << PosSphere << "\n                to triangle is " << DSphereTriangle
    //           << std::endl;
    // std::cout << "Distance from Point at " << PosSphere << "\n                to triangle is " << DPointTriangle
    //           << std::endl;
  };

  CheckDistTriangleSphere(M * ww::Point(1.f, 1.f, 0.f), 1.f, T, -(1.f - 1.f / std::sqrt(2.f)));
  CheckDistTriangleSphere(ww::Point(1.f + ww::Get(M, 0, 3), 1.f + ww::Get(M, 1, 3), ww::Get(M, 2, 3)), 1.f, T,
                          -(1.f - 1.f / std::sqrt(2.f)));

  // Calculate position of sphere:
  // h^2 = 1+1/sqrt(2)
  // h^2 = sqrt(2*x^2) // x and y are of same length.
  // => h^2 = sqrt(2)*x => x = h^2/sqrt(2) => x = (1+1/sqrt(2))/sqrt(2)
  float const SQRT2 = std::sqrt(2.f);
  float const X = (1.f + 1.f / SQRT2) / SQRT2;
  float const Y = X;
  CheckDistTriangleSphere(M * ww::Point(X, Y, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-X, Y, 0.f), 1.f, T, 0.f);
  //
  CheckDistTriangleSphere(M * ww::Point(0.f, 2.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, -1.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, -2.f, 0.f), 1.f, T, 1.f);
  CheckDistTriangleSphere(M * ww::Point(2.f, 0.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-2.f, 0.f, 0.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-1.f, 0.f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(1.f, 0.f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, 1.f, 1.f), 1.f, T, 0.f);

  // Inside the extrusion of the triangle
  CheckDistTriangleSphere(M * ww::Point(0.f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(.5f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-.5f, .5f, 1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(0.f, .5f, -1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(.5f, .5f, -1.f), 1.f, T, 0.f);
  CheckDistTriangleSphere(M * ww::Point(-.5f, .5f, -1.f), 1.f, T, 0.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, CreateCapsule)
{
  ww::tup A = ww::Point(0.f, 0.f, 0.f);
  ww::tup B = ww::Point(1.f, 0.f, 0.f);

  auto CheckDist = [](ww::tup const &A, ww::tup const &B, ww::tup const &Pos, ww::tup const &Expect,
                      bool Debug = false) -> ww::tup
  {
    auto const D = ww::rm::ClosestPointOnLineSegment(A, B, Pos);
    EXPECT_NEAR(D.X, Expect.X, 1e-6);
    EXPECT_NEAR(D.Y, Expect.Y, 1e-6);
    EXPECT_NEAR(D.Z, Expect.Z, 1e-6);

    return D;
  };

  CheckDist(A, B, ww::Point(0.f, 1.f, 0.f), ww::Point(0.f, 0.f, 0.f));
  CheckDist(A, B, ww::Point(0.f, 0.f, 1.f), ww::Point(0.f, 0.f, 0.f));

  CheckDist(A, B, ww::Point(0.f, -1.f, 0.f), ww::Point(0.f, 0.f, 0.f));
  CheckDist(A, B, ww::Point(0.f, 0.f, -1.f), ww::Point(0.f, 0.f, 0.f));

  CheckDist(A, B, ww::Point(0.5f, 1.f, 0.f), ww::Point(0.5f, 0.f, 0.f));
  CheckDist(A, B, ww::Point(0.5f, 0.f, 1.f), ww::Point(0.5f, 0.f, 0.f));

  CheckDist(A, B, ww::Point(0.5f, 0.f, 0.f), ww::Point(0.5f, 0.f, 0.f));
  CheckDist(A, B, ww::Point(1.5f, 0.f, 0.f), ww::Point(1.0f, 0.f, 0.f));

  CheckDist(A, B, ww::Point(-0.5f, 0.f, 0.f), ww::Point(0.0f, 0.f, 0.f));
  CheckDist(A, B, ww::Point(-1.5f, 0.f, 0.f), ww::Point(0.0f, 0.f, 0.f));

  CheckDist(A, B, ww::Point(1.f, 1.f, 0.f), ww::Point(1.0f, 0.f, 0.f));
  CheckDist(A, B, ww::Point(1.f, 1.f, 1.f), ww::Point(1.0f, 0.f, 0.f));
  CheckDist(A, B, ww::Point(1.f, 1.f, -1.f), ww::Point(1.0f, 0.f, 0.f));

  /*
   * Create a capsule struct for test.
   */

  auto CreateCapsule = [](ww::tup const &A, ww::tup const &B, float R) -> ww::capsule
  {
    Assert(ww::IsPoint(A), __FUNCTION__, __LINE__);
    Assert(ww::IsPoint(B), __FUNCTION__, __LINE__);

    ww::capsule C{};
    C.R = R;
    C.A = A;
    C.B = B;

    ww::tup const AB = B - A;  // Vector created from A to B.
    Assert(ww::IsVector(AB), __FUNCTION__, __LINE__);

    C.Base = A - ww::Normalize(AB) * C.R;  // Subtract from point B to get to Base.
    C.Tip = B + ww::Normalize(AB) * C.R;   // Add to point A to get to Tip.

    std::cout << "Capsule ---" << std::endl;
    std::cout << "From A: " << C.A << "\n  to B: " << C.B << "\n"
              << "Base  : " << C.Base << "\nTip   : " << C.Tip << "\nRadius : " << C.R << std::endl;

    return C;
  };

  {
    ww::capsule const Ca = CreateCapsule(ww::Point(0.f, 0.f, 0.f), ww::Point(1.f, 0.f, 0.f), 0.1f);
    ww::capsule const Cb = CreateCapsule(ww::Point(0.f, 1.f, 0.f), ww::Point(1.f, 1.f, 0.f), 0.1f);
    auto const D = ww::rm::sdCapsuleCapsule(Ca, Cb);
    EXPECT_FLOAT_EQ(D, 0.8f);
  }
  {
    ww::capsule const Ca = CreateCapsule(ww::Point(0.f, 0.f, 0.f), ww::Point(1.f, 0.f, 0.f), 0.1f);
    ww::capsule const Cb = CreateCapsule(ww::Point(0.f, 1.f, 0.f), ww::Point(0.f, 2.f, 0.f), 0.1f);
    auto const D = ww::rm::sdCapsuleCapsule(Ca, Cb);
    EXPECT_FLOAT_EQ(D, 0.8f);
  }
  {
    ww::capsule const Ca = CreateCapsule(ww::Point(0.f, 0.f, 0.f), ww::Point(1.f, 0.f, 0.f), 0.1f);
    ww::capsule const Cb =
        CreateCapsule(ww::Point(0.722183f, 1.177817f, 0.f), ww::Point(1.429289f, 0.470711f, 0.f), 0.1f);
    auto const D = ww::rm::sdCapsuleCapsule(Ca, Cb);
    // ---
    // NOTE: Handle floating point rounding error and compare against EPSILON.
    // ---
    EXPECT_LT(std::abs(D - 0.437195), ww::EPSILON);
  }
  {
    ww::capsule const Ca = CreateCapsule(ww::Point(0.f, 0.f, 0.f), ww::Point(1.f, 0.f, 0.f), 0.1f);
    ww::capsule const Cb = CreateCapsule(ww::Point(0.5f, -0.1f, 0.f), ww::Point(0.5f, -1.1f, 0.f), 0.1f);
    auto const D = ww::rm::sdCapsuleCapsule(Ca, Cb);
    EXPECT_FLOAT_EQ(D, -0.1f);
  }
  {
    ww::capsule const Ca = CreateCapsule(ww::Point(0.f, 0.f, 0.f), ww::Point(1.f, 0.f, 0.f), 0.1f);
    ww::capsule const Cb = CreateCapsule(ww::Point(0.5f, 2.1f, 0.f), ww::Point(0.5f, 0.3f, 0.f), 0.1f);
    auto const D = ww::rm::sdCapsuleCapsule(Ca, Cb);
    // ---
    // NOTE: Handle floating point rounding error and compare against EPSILON.
    // ---
    EXPECT_LT(std::abs(D - 0.1f), ww::EPSILON);
  }

  {
    ww::capsule const Ca = CreateCapsule(ww::Point(1.f, 0.f, 0.f), ww::Point(0.f, 0.f, 0.f), 0.1f);
    ww::capsule const Cb = CreateCapsule(ww::Point(0.f, 2.f, 0.f), ww::Point(0.f, 1.f, 0.f), 0.1f);
    auto const D = ww::rm::sdCapsuleCapsule(Ca, Cb);
    // ---
    // NOTE: Handle floating point rounding error and compare against EPSILON.
    // ---
    EXPECT_LT(std::abs(D - 0.8f), ww::EPSILON);
  }

  {
    ww::capsule const Ca = CreateCapsule(ww::Point(3.23f, 25.16f, 1.84f), ww::Point(2.32f, 24.61f, 4.77f), 0.1f);
    ww::capsule const Cb = CreateCapsule(ww::Point(18.06f, 23.f, -2.53f), ww::Point(0.94f, 25.f, 4.84f), 0.01f);
    auto const D = ww::rm::sdCapsuleCapsule(Ca, Cb);
    EXPECT_FLOAT_EQ(D, 0.23186623);
  }

  {
    ww::capsule const C1 =
        CreateCapsule(ww::Point(2.82621f, 25.44422, 1.914072f), ww::Point(2.322477f, 24.60991f, 4.769926), 0.1f);
    ww::capsule const C2 =
        CreateCapsule(ww::Point(21.34589f, 24.f, -4.173028f), ww::Point(-6.346084, 26.f, 5.390827f), 0.01f);
    auto const D = ww::rm::sdCapsuleCapsule(C1, C2);
    EXPECT_FLOAT_EQ(D, 0.23186623);
  }
}

//------------------------------------------------------------------------------
TEST(RayMarch, udSquare)
{
  ww::quad Q{};
  EXPECT_FLOAT_EQ(Q.A.X, -1.f);
  EXPECT_FLOAT_EQ(Q.A.Y, 0.f);
  EXPECT_FLOAT_EQ(Q.A.Z, -1.f);
  EXPECT_FLOAT_EQ(Q.A.W, 1.f);

  EXPECT_FLOAT_EQ(Q.B.X, 1.f);
  EXPECT_FLOAT_EQ(Q.B.Y, 0.f);
  EXPECT_FLOAT_EQ(Q.B.Z, -1.f);
  EXPECT_FLOAT_EQ(Q.B.W, 1.f);

  EXPECT_FLOAT_EQ(Q.C.X, 1.f);
  EXPECT_FLOAT_EQ(Q.C.Y, 0.f);
  EXPECT_FLOAT_EQ(Q.C.Z, 1.f);
  EXPECT_FLOAT_EQ(Q.C.W, 1.f);

  EXPECT_FLOAT_EQ(Q.D.X, -1.f);
  EXPECT_FLOAT_EQ(Q.D.Y, 0.f);
  EXPECT_FLOAT_EQ(Q.D.Z, 1.f);
  EXPECT_FLOAT_EQ(Q.D.W, 1.f);

  auto CheckDist = [&Q](ww::tup const &P, float const Expect)
  {
    float const Dist = ww::rm::udQuad(P, Q.A, Q.B, Q.C, Q.D);
    EXPECT_FLOAT_EQ(Dist, Expect);
  };

  CheckDist(ww::Point(0.f, 1.f, 0.f), 1.f);
  CheckDist(ww::Point(1.f, 1.f, 0.f), 1.f);
  CheckDist(ww::Point(-1.f, 1.f, -1.f), 1.f);
  CheckDist(ww::Point(-1.f, 1.f, +1.f), 1.f);
  CheckDist(ww::Point(+1.f, 1.f, -1.f), 1.f);
  CheckDist(ww::Point(+1.f, 1.f, +1.f), 1.f);
  CheckDist(ww::Point(+2.f, 0.f, +0.f), 1.f);
  CheckDist(ww::Point(-2.f, 0.f, +0.f), 1.f);
  CheckDist(ww::Point(0.f, 0.f, +2.f), 1.f);
  CheckDist(ww::Point(0.f, 0.f, -2.f), 1.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, udSquareTranslate)
{
  ww::quad QOrigin{};
  ww::matrix M =
      ww::TranslateScaleRotate(10.f, 0.f, 0.f, 1.f, 1.f, 1.f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));
  ww::quad QTranslated{};

  // Move the Quad as specified by the matrix to the desired world coordinates.
  QTranslated.A = M * QOrigin.A;
  QTranslated.B = M * QOrigin.B;
  QTranslated.C = M * QOrigin.C;
  QTranslated.D = M * QOrigin.D;

  // Compute the inverse of the matrix in order to be able to use local coordinates.
  ww::matrix MI = ww::Inverse(M);

  // Compute the local position of the quad.
  ww::quad QLocal{};
  QLocal.A = MI * QTranslated.A;
  QLocal.B = MI * QTranslated.B;
  QLocal.C = MI * QTranslated.C;
  QLocal.D = MI * QTranslated.D;

  auto CheckDist = [&QLocal](ww::tup const &P, ww::matrix const &MI, float const Expect)
  {
    float const Dist = ww::rm::udQuad(MI * P, QLocal.A, QLocal.B, QLocal.C, QLocal.D);
    EXPECT_FLOAT_EQ(Dist, Expect);
  };

  // Check distances from points given in world coordinates.
  CheckDist(ww::Point(0.f, 0.f, 0.f), MI, 9.f);
  CheckDist(ww::Point(9.f, 0.f, 0.f), MI, 0.f);
  CheckDist(ww::Point(11.f, 0.f, 0.f), MI, 0.f);
  CheckDist(ww::Point(20.f, 0.f, 0.f), MI, 9.f);
  CheckDist(ww::Point(10.f, 1.f, 0.f), MI, 1.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, CapsuleVsTriangle)
{
  // ---
  // NOTE: Create a triangle with one vertex at origin.
  // ---
  ww::tup V0 = ww::Point(0.f, 0.f, 0.f);
  ww::tup V1 = ww::Point(1.f, 0.f, 0.f);
  ww::tup V2 = ww::Point(0.f, 0.f, 1.f);

  /*
   * Create a capsule struct for test.
   */

  // auto CreateCapsule = [](ww::tup const &A, ww::tup const &B, float R, bool Print) -> ww::capsule
  // {
  //   Assert(ww::IsPoint(A), __FUNCTION__, __LINE__);
  //   Assert(ww::IsPoint(B), __FUNCTION__, __LINE__);
  //
  //   ww::capsule C{};
  //   C.R = R;
  //   C.A = A;
  //   C.B = B;
  //
  //   ww::tup const AB = B - A;  // Vector created from A to B.
  //   Assert(ww::IsVector(AB), __FUNCTION__, __LINE__);
  //
  //   C.Base = A - ww::Normalize(AB) * C.R;  // Subtract from point B to get to Base.
  //   C.Tip = B + ww::Normalize(AB) * C.R;   // Add to point A to get to Tip.
  //
  //   if (Print)
  //   {
  //     std::cout << "Capsule ---" << std::endl;
  //     std::cout << "From A: " << C.A << "\n  to B: " << C.B << "\n"
  //               << "Base  : " << C.Base << "\nTip   : " << C.Tip << "\nRadius : " << C.R << std::endl;
  //   }
  //   return C;
  // };

  auto TestIntersect = [&](ww::tup const &CapA, ww::tup const &CapB, float Radius, bool Expect, bool Print) -> void
  {
    ww::capsule const C = CreateCapsule(CapA, CapB, Radius, Print);

    float Distance{};
    ww::tup HitPoint{};
    ww::tup HitNormal{};

    bool const Intersect =
        ww::rm::CapsuleTriangleIntersect(C.A, C.B, C.R, V0, V1, V2, &Distance, &HitPoint, &HitNormal, Print);
    EXPECT_EQ(Intersect, Expect);
    if (Print)
    {
      std::cout << "\nDistance :" << Distance;
      std::cout << "\nHitPoint :" << HitPoint;
      std::cout << "\nHitNormal:" << HitNormal;
      std::cout << "\nIntersect:" << Intersect;
      std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    }
  };
  auto Print = [](bool V) -> bool { return V; };
  auto Expect = [](bool V) -> bool { return V; };
  // auto Elevate = [](float h) -> float { return h; };

  TestIntersect(ww::Point(-1.f, -0.7f, 0.f), ww::Point(1.f, 1.f, 0.f), 0.1, Expect(true), Print(false));
  TestIntersect(ww::Point(-1.f, -1.f, 0.f), ww::Point(-1.f, 1.f, 0.f), 0.1, Expect(false), Print(false));
  TestIntersect(ww::Point(-.1f, -1.f, 0.f), ww::Point(-.1f, 1.f, 0.f), 0.11, Expect(true), Print(false));
  TestIntersect(ww::Point(0.f, 0.f, 0.f), ww::Point(0.f, 1.f, 0.f), 0.1, Expect(true), Print(false));
  TestIntersect(ww::Point(0.25f, 0.1f, 0.25f), ww::Point(0.25f, 1.f, 0.25f), 0.15, Expect(true), Print(true));
}

//------------------------------------------------------------------------------
TEST(RayMarch, DistanceBetweenPoints)
{
  /**
   *            c
   *            [0.5,1]
   *            x
   *            ^ \
   *          V2|  \
   *            |   \
   *      x<----<----x b[1,0]
   *      a[0,0]  V1
   *
   */
  ww::tup Pa = ww::Point(0.f, 0.f, 0.f);
  ww::tup Pb = ww::Point(1.f, 0.f, 0.f);
  ww::tup Pc = ww::Point(0.5f, 1.f, 0.f);

  ww::tup Vba = Pa - Pb;
  ww::tup Vbc = Pc - Pb;

  // Find the projection of vector from b->c onto vector from b->a.
  float l1 = ww::Dot(Vba, Vbc);
  EXPECT_FLOAT_EQ(l1, 0.5f);

  // V1 + V2 = Vbc
  // => V2 = Vbc - V1
  ww::tup V1 = l1 * Vba;
  ww::tup V2 = Vbc - V1;
  EXPECT_FLOAT_EQ(ww::Mag(V2), 1.f);
}

//------------------------------------------------------------------------------
TEST(RayMarch, RayIntersectTriangle)
{
  auto Expect = [](bool V) -> bool { return V; };

  ww::triangle T{};
  EXPECT_EQ(T.VA == ww::Point(0.f, 0.f, 0.f), true);
  EXPECT_EQ(T.VB == ww::Point(1.f, 0.f, 0.f), true);
  EXPECT_EQ(T.VC == ww::Point(0.f, 0.f, 1.f), true);

  auto CheckRayIntersectTriangle = [](ww::ray const &R, ww::triangle const &T, bool TestValue)
  { EXPECT_EQ(ww::rm::RayTriangleIntersect(R, T).Hit, TestValue); };

  // Ray pointing down
  CheckRayIntersectTriangle(ww::Ray(ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, -1.f, 0.f)), T, Expect(true));
  CheckRayIntersectTriangle(ww::Ray(ww::Point(1.f, 1.f, 0.f), ww::Vector(0.f, -1.f, 0.f)), T, Expect(true));
  CheckRayIntersectTriangle(ww::Ray(ww::Point(0.f, 1.f, 1.f), ww::Vector(0.f, -1.f, 0.f)), T, Expect(true));
  CheckRayIntersectTriangle(ww::Ray(ww::Point(1.f, 1.f, 1.f), ww::Vector(0.f, -1.f, 0.f)), T, Expect(false));
  CheckRayIntersectTriangle(ww::Ray(ww::Point(.1f, 1.f, .1f), ww::Vector(0.f, -1.f, 0.f)), T, Expect(true));

  // Ray pointing up
  CheckRayIntersectTriangle(ww::Ray(ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, 1.f, 0.f)), T, Expect(false));

  // Move the triangle up to +y=2:
  T.VA.Y = 2.f;
  T.VB.Y = 2.f;
  T.VC.Y = 2.f;
  CheckRayIntersectTriangle(ww::Ray(ww::Point(0.f, 1.f, 0.f), ww::Vector(0.f, 1.f, 0.f)), T, Expect(true));
}

//------------------------------------------------------------------------------
TEST(RayMarch, LineSegmentIntersectTriangle)
{
  // ---
  // NOTE: Convenience lambda for setting flags.
  // ---
  auto Print = [](bool V) -> bool { return V; };
  auto Expect = [](bool V) -> bool { return V; };

  ww::triangle T{};
  ww::capsule C = ww::CreateCapsule(ww::Point(.1f, 1.f, .2f), ww::Point(.1f, -1.f, -.0f), 0.1f, Print(true));

  auto CheckRayIntersectTriangle = [&Print](ww::ray const &R, ww::triangle const &T, bool TestValue)
  {
    auto const Result = ww::rm::RayTriangleIntersect(R, T, Print(false));
    std::cout << "\n"
              << R << "\nHit:" << Result.Hit << ". t:" << Result.t << ". u:" << Result.u << ". v:" << Result.v
              << std::endl;
    if (Result.Hit) std::cout << "\nPoint of Hit:" << Result.P << std::endl;

    EXPECT_EQ(Result.Hit, TestValue);
    return Result;
  };

  // ---
  // NOTE: Cast a ray from the base in the direction of the tip to see if there is a triangle hit.
  // ---
  ww::rm::rti_result HitTipBase{};
  ww::rm::rti_result HitBaseTip{};

  {
    ww::ray const R = ww::Ray(C.Base, C.Tip - C.Base);
    HitTipBase = CheckRayIntersectTriangle(R, T, Expect(true));
  }

  // ---
  // NOTE: Cast a ray from the tip in the direction of the base to see if there is a triangle hit.
  // ---
  {
    ww::ray const R = ww::Ray(C.Tip, C.Base - C.Tip);
    HitBaseTip = CheckRayIntersectTriangle(R, T, Expect(true));
  }

  // ---
  // NOTE: So there is an intersection when both rays hits the triangle.
  // ---
  EXPECT_EQ(HitBaseTip.Hit && HitTipBase.Hit, true);
}
