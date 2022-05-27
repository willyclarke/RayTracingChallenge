#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>
//------------------------------------------------------------------------------
// Scenario: Creating a stripe pattern.
TEST(Ch10Patterns, CreatingStripePattern)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(Pattern.A == White, true);
  EXPECT_EQ(Pattern.B == Black, true);
}

//------------------------------------------------------------------------------
// Scenario: A stripe pattern is constant in y.
TEST(Ch10Patterns, AStripePatternIsConstantInY)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 1.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 2.f, 0.f)) == White, true);
}

//------------------------------------------------------------------------------
// Scenario: A stripe pattern is constant in z.
TEST(Ch10Patterns, AStripePatternIsConstantInZ)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 1.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 2.f)) == White, true);
}

//------------------------------------------------------------------------------
// Scenario: A stripe pattern alternates in x.
TEST(Ch10Patterns, AStripePatternAlternatesInX)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern const Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(0.9f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(1.f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(-0.1f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(-1.f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(ww::StripeAt(Pattern, ww::Point(-1.1f, 0.f, 0.f)) == White, true);
}

//------------------------------------------------------------------------------
// Scenario: Lighting with a pattern applied.
TEST(Ch10Patterns, LightingWithAPatternApplied)
{
  ww::material Material{};
  Material.Pattern = ww::StripePattern(ww::Color(1.f, 1.f, 1.f), ww::Color(0.f, 0.f, 0.f));
  Material.Ambient = 1.f;
  Material.Diffuse = 0.f;
  Material.Specular = 0.f;

  ww::tup vEye = ww::Vector(0.f, 0.f, -1.f);
  ww::tup vNormal = ww::Vector(0.f, 0.f, -1.f);
  ww::light const Light = ww::PointLight(ww::Point(0.f, 0.f, -10.f), ww::Color(1.f, 1.f, 1.f));

  bool const NotInShadow{};
  ww::tup C1 = ww::Lighting(Material, ww::sphere{}, Light, ww::Point(0.9f, 0.f, 0.f), vEye, vNormal, NotInShadow);
  ww::tup C2 = ww::Lighting(Material, ww::sphere{}, Light, ww::Point(1.1f, 0.f, 0.f), vEye, vNormal, NotInShadow);

  EXPECT_EQ(C1 == ww::Color(1.f, 1.f, 1.f), true);
  EXPECT_EQ(C2 == ww::Color(0.f, 0.f, 0.f), true);
}

//------------------------------------------------------------------------------
// Scenario: Stripes with an object transformation
TEST(Ch10Patterns, StripesWithAnObjectTransformation)
{
  ww::sphere Object{};
  Object.Transform = ww::Scaling(2.f, 2.f, 2.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::pattern const Pattern = ww::StripePattern(White, Black);
  ww::tup const C = ww::StripeAtObject(Pattern, Object, ww::Point(1.5f, 0.f, 0.f));

  EXPECT_EQ(C == White, true);
}

//------------------------------------------------------------------------------
// Scenario: Stripes with a pattern transformation
TEST(Ch10Patterns, StripesWithAPatternTransformation)
{
  ww::sphere Object{};
  Object.Transform = ww::Scaling(2.f, 2.f, 2.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);

  ww::pattern Pattern = ww::StripePattern(White, Black);
  Pattern.Transform = ww::Scaling(2.f, 2.f, 2.f);

  ww::tup const C = ww::StripeAtObject(Pattern, Object, ww::Point(1.5f, 0.f, 0.f));

  EXPECT_EQ(C == White, true);
}

//------------------------------------------------------------------------------
// Scenario: Stripes with both an object and a pattern transformation
TEST(Ch10Patterns, StripesWithBothAnObjectAndAPatternTransformation)
{
  ww::sphere Object{};
  Object.Transform = ww::Scaling(2.f, 2.f, 2.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::pattern Pattern = ww::StripePattern(White, Black);
  Pattern.Transform = ww::Translation(0.5f, 0.f, 0.f);

  ww::tup const C = ww::StripeAtObject(Pattern, Object, ww::Point(2.5f, 0.f, 0.f));

  EXPECT_EQ(C == White, true);
}

//------------------------------------------------------------------------------
TEST(Ch10Patterns, AlmostThere)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(5.5f, 0.0f, 0.f);

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
  Middle.Material.Pattern = ww::StripePattern(White, Black);
  Middle.Material.Pattern.Transform = ww::Scaling(0.15f, 0.5f, 0.5f) * ww::RotateZ(0.78);
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

  ww::tup const ViewFrom = ww::Point(0.f, 1.5f, -5.f);
  ww::tup const ViewTo = ww::Point(0.f, 1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch10AlmostThere.ppm");
}

//------------------------------------------------------------------------------
// Scenario: The default pattern transformation.
TEST(Ch10Patterns, TheDefaultPatternTransformation)
{
  ww::pattern const Pattern{};
  EXPECT_EQ(Pattern.Transform == ww::I(), true);
}

//------------------------------------------------------------------------------
// Scenario: Assigning a transformation.
TEST(Ch10Patterns, AssigningATransformation)
{
  ww::pattern const Pattern = ww::TestPattern();
  EXPECT_EQ(Pattern.Transform == ww::Translation(1.f, 2.f, 3.f), true);
}
