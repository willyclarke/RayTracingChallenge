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
TEST(RayMarch, Test1)
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
  World.Map = ww::rm::MapDefault;

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
  std::cout << "Distance to Capsule from point " << P6 << " is " << D6 << std::endl;

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
  char const *pTestDescription = "TestSignedDistanceFunctionsBox2";
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
  std::cout << pTestDescription << " -> Distance to Capsule from point " << P6 << " is " << D6 << std::endl;

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
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << pTestDescription << ". Line:" << __LINE__ << std::endl;
    std::cout << "MBox :" << M << std::endl;
    std::cout << "MBoxI:" << MI << std::endl;

    std::cout << "Box Center: " << Box << std::endl;
    std::cout << "P         : " << P << std::endl;
    std::cout << "Pmov      : " << MI * P << std::endl;
    std::cout << "Distance  : " << D << std::endl;

    EXPECT_FLOAT_EQ(D, DExpected);
  };

  CheckDist(MBox, MBoxI, P1, BoxPos, D1, std::sqrtf(10.f * 10.f + 9.f * 9.f));
  CheckDist(MBox, MBoxI, P2, BoxPos, D2, std::sqrtf(10.f * 10.f + 9.f * 9.f));
  CheckDist(MBox, MBoxI, P3, BoxPos, D3, 9.f);
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
