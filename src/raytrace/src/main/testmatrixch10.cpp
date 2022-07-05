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

  ww::pattern Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.f, 1.f, 0.f)) == White, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.f, 2.f, 0.f)) == White, true);
}

//------------------------------------------------------------------------------
// Scenario: A stripe pattern is constant in z.
TEST(Ch10Patterns, AStripePatternIsConstantInZ)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.f, 0.f, 1.f)) == White, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.f, 0.f, 2.f)) == White, true);
}

//------------------------------------------------------------------------------
// Scenario: A stripe pattern alternates in x.
TEST(Ch10Patterns, AStripePatternAlternatesInX)
{
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);

  ww::pattern Pattern = ww::StripePattern(White, Black);

  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(0.9f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(1.f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(-0.1f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(-1.f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(Pattern.funcPtrPatternAt(Pattern, ww::Point(-1.1f, 0.f, 0.f)) == White, true);
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
  ww::pattern Pattern = ww::StripePattern(White, Black);
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

//------------------------------------------------------------------------------
// Scenario: A pattern with an object transformation.
// NOTE: This test will only test that the scaling works for the
//       PatternAtShape. Tweaks in PatternAtShape is neccessary
//       for this to work. Test is described on page 178 of the Challenge.
TEST(Ch10Patterns, APatternWithAnObjectTransformation)
{
  return;
  ww::sphere Shape = *ww::PtrDefaultSphere();
  Shape.Transform = ww::Scaling(2.f, 2.f, 2.f);
  ww::pattern Pattern = ww::TestPattern();
  ww::tup const C = ww::PatternAtShape(Pattern, Shape, ww::Point(2.f, 3.f, 4.f));

  EXPECT_EQ(C == ww::Color(1.f, 1.5f, 2.f), true);
  EXPECT_EQ(Shape.isA<ww::sphere>() == true, true);
}

//------------------------------------------------------------------------------
// Scenario: A pattern with a pattern transformation.
TEST(Ch10Patterns, APatternWithAPatternTransformation)
{
  ww::sphere Shape = *ww::PtrDefaultSphere();
  ww::pattern Pattern = ww::TestPattern();
  Pattern.Transform = ww::Scaling(2.f, 2.f, 2.f);
  ww::tup const C = ww::PatternAtShape(Pattern, Shape, ww::Point(2.f, 3.f, 4.f));

  EXPECT_EQ(C == ww::Color(1.f, 1.5f, 2.f), true);
  EXPECT_EQ(Shape.isA<ww::sphere>() == true, true);
}

//------------------------------------------------------------------------------
// Scenario: A pattern with both an object and a pattern transformation.
TEST(Ch10Patterns, APatternWithBothAnObjectAndAPatternTransformation)
{
  ww::sphere Shape = *ww::PtrDefaultSphere();
  Shape.Transform = ww::Scaling(2.f, 2.f, 2.f);
  ww::pattern Pattern = ww::TestPattern();
  Pattern.Transform = ww::Translation(0.5f, 1.f, 1.5f);
  ww::tup const C = ww::PatternAtShape(Pattern, Shape, ww::Point(2.5f, 3.f, 3.5f));

  EXPECT_EQ(C == ww::Color(0.75f, 0.5f, 0.25f), true);
  EXPECT_EQ(Shape.isA<ww::sphere>() == true, true);
}

//------------------------------------------------------------------------------
// Scenario: A gradient lineaerly interpolates between colors.
TEST(Ch10Patterns, AGradientLinearlyInterpolatesBetweenColors)
{
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::pattern Pattern{White, Black};
  Pattern.funcPtrPatternAt = &ww::GradientPatternAt;
  EXPECT_EQ(ww::PatternAt(Pattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::PatternAt(Pattern, ww::Point(0.25f, 0.f, 0.f)) == ww::Color(0.75f, 0.75f, 0.75f), true);
  EXPECT_EQ(ww::PatternAt(Pattern, ww::Point(0.5f, 0.f, 0.f)) == ww::Color(0.5f, 0.5f, 0.5f), true);
  EXPECT_EQ(ww::PatternAt(Pattern, ww::Point(0.75f, 0.f, 0.f)) == ww::Color(0.25f, 0.25f, 0.25f), true);
}

//------------------------------------------------------------------------------
// Scenario: A ring should extend in both x and z.
TEST(Ch10Patterns, ARingShouldExtendInBothXandZ)
{
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);

  ww::pattern RingPattern{White, Black};
  RingPattern.funcPtrPatternAt = &ww::RingPatternAt;

  EXPECT_EQ(ww::PatternAt(RingPattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::PatternAt(RingPattern, ww::Point(1.f, 0.f, 0.f)) == Black, true);
  EXPECT_EQ(ww::PatternAt(RingPattern, ww::Point(0.f, 0.f, 1.f)) == Black, true);

  // ---
  // NOTE: 0.708 is slightly more than sqrt(2)/2.
  // ---
  EXPECT_EQ(ww::PatternAt(RingPattern, ww::Point(0.708f, 0.f, 0.708f)) == Black, true);
}

//------------------------------------------------------------------------------
// Scenario: Checkers should repeat in x.
TEST(Ch10Patterns, CheckersShouldRepeatInX)
{
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);

  // ---
  // NOTE: Use default CTOR to set up pointer to function.
  // ---
  ww::pattern CheckersPattern{White, Black, ww::I(), true, &ww::CheckersPatternAt};
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.99f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(1.01f, 0.f, 0.f)) == Black, true);
}

//------------------------------------------------------------------------------
// Scenario: Checkers should repeat in y.
TEST(Ch10Patterns, CheckersShouldRepeatInY)
{
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);

  // ---
  // NOTE: Use default CTOR to set up pointer to function.
  // ---
  ww::pattern CheckersPattern{White, Black, ww::I(), true, &ww::CheckersPatternAt};
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.f, 0.99f, 0.f)) == White, true);
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.f, 1.01f, 0.f)) == Black, true);
}

//------------------------------------------------------------------------------
// Scenario: Checkers should repeat in z.
TEST(Ch10Patterns, CheckersShouldRepeatInZ)
{
  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);

  // ---
  // NOTE: Use default CTOR to set up pointer to function.
  // ---
  ww::pattern CheckersPattern{White, Black, ww::I(), true, &ww::CheckersPatternAt};
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.f, 0.f, 0.f)) == White, true);
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.f, 0.f, 0.99f)) == White, true);
  EXPECT_EQ(ww::PatternAt(CheckersPattern, ww::Point(0.f, 0.f, 1.01f)) == Black, true);
}

//------------------------------------------------------------------------------
TEST(Ch10Patterns, PuttingItTogether)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const Red = ww::Color(1.f, 0.f, 0.f);
  ww::tup const Green = ww::Color(0.f, 1.f, 0.f);
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(0.f, .5f, .5f);

  ww::shared_ptr_shape PtrRight = ww::PtrDefaultSphere();
  ww::shape &Right = *PtrRight;
  Right.Transform = ww::Translation(1.5f, 2.5f, -1.0f) *  //!<
                    ww::Scaling(0.5f, 0.5f, 0.5f);
  Right.Material.Color = ww::Color(0.5f, 1.0f, 0.1f);
  Right.Material.Diffuse = 0.7f;
  Right.Material.Specular = 0.3f;
  Right.Material.Pattern = ww::GradientPattern(Blue, Red);
  World.vPtrObjects.push_back(PtrRight);

  ww::shared_ptr_shape PtrMiddle = ww::PtrDefaultSphere();
  ww::shape &Middle = *PtrMiddle;
  Middle.Transform = ww::Translation(0.f, 1.f, 0.0f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  Middle.Material.Pattern = ww::CheckersPattern(White, Green);
  Middle.Material.Pattern.Transform = ww::Scaling(0.15f, 0.5f, 0.5f) * ww::RotateZ(0.78);
  World.vPtrObjects.push_back(PtrMiddle);

  ww::shared_ptr_shape PtrLeft = ww::PtrDefaultSphere();
  ww::shape &Left = *PtrLeft;
  Left.Transform = ww::Translation(-2.66f, 1.33f, 0.f) *  //!<
                   ww::Scaling(1.33f, 1.33f, 1.33f);
  Left.Material.Color = ww::Color(0.0f, 0.8f, 0.1f);
  Left.Material.Diffuse = 0.7f;
  Left.Material.Specular = 0.3f;
  Left.Material.Pattern = ww::RingPattern(White, Blue);
  World.vPtrObjects.push_back(PtrLeft);

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::Translation(0.f, 0.f, 0.f)  //!<
                                                        // * ww::RotateX(0 * M_PI_2)       //!<
      ;                                                 //!<
  ptrPlane->Material.Shininess = 10.f;
  ptrPlane->Material.Diffuse = 0.5f;
  ptrPlane->Material.Color = ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  ptrPlane->Material.Pattern = ww::RadialGradientPattern(Red, Green);
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
  ptrPlane2->Material.Pattern = ww::CheckersGradientPattern(White, Yellow);
  World.vPtrObjects.push_back(ptrPlane2);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(0.f, 5.f, -5.f), ww::Color(2.f, 2.f, 2.f));
  World.vPtrLights.push_back(pLight);

  float const FieldOfView = 75.f / 180.f * M_PI;
  ww::camera Camera = ww::Camera(256, 256, FieldOfView);

  ww::tup const ViewFrom = ww::Point(0.f, 1.5f, -5.f);
  ww::tup const ViewTo = ww::Point(0.f, 0.f, 2.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch10PuttingItTogether.ppm");
}

//------------------------------------------------------------------------------
TEST(Ch10Patterns, PuttingItTogetherNestedPattern)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const Grey = ww::Color(float(84.f / 255.f), float(84.f / 255.f), float(84.f / 255.f));
  ww::tup const Grey2 = ww::Color(float(44.f / 255.f), float(44.f / 255.f), float(44.f / 255.f));
  ww::tup const Red = ww::Color(1.f, 0.f, 0.f);
  ww::tup const Red2 = ww::Color(1.f, .4f, 0.f);
  ww::tup const Green = ww::Color(0.f, 1.f, 0.f);
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(1.f, 1.f, .0f);

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::Translation(0.f, 0.f, 0.f)  //!<
                                                        // * ww::RotateY(0.5 * M_PI_2)       //!<
      ;                                                 //!<

  float const Alpha = M_PI * float(45.f / 180.f);
  ptrPlane->Material.Shininess = 1.f;
  ptrPlane->Material.Diffuse = 0.5f;
  ptrPlane->Material.Color = ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  ww::pattern PM = ww::CheckersPattern(Blue, Yellow);
  ww::pattern P1 = ww::GradientPattern(Grey, Grey2);
  ww::pattern P2 = ww::GradientPattern(Red, Red2);
  // P1.Continue = false;
  PM.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 4.f, 4.f, 4.f, 0.f, Alpha, 0.f);
  P1.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, .25f, 1.f, .25f, 0.f, .5f * Alpha, 0.f);
  P2.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 0.1f, .1f, 1.0f, 0.f, Alpha, 0.f);
  ptrPlane->Material.Pattern = ww::NestedPattern(PM, P1, P2);
  World.vPtrObjects.push_back(ptrPlane);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-2.f, 5.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  float const FieldOfView = 75.f / 180.f * M_PI;
  ww::camera Camera = ww::Camera(256, 256, FieldOfView);

  ww::tup const ViewFrom = ww::Point(0.f, 5.f, -5.f);
  ww::tup const ViewTo = ww::Point(0.f, 0.f, 2.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  // Assert(false, __FUNCTION__, __LINE__);
  ww::WriteToPPM(Canvas, "Ch10PuttingItTogetherNestedPattern.ppm");
}

//------------------------------------------------------------------------------
TEST(Ch10Patterns, PuttingItTogetherBlendedPattern)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const Grey = ww::Color(float(84.f / 255.f), float(84.f / 255.f), float(84.f / 255.f));
  ww::tup const Grey2 = ww::Color(float(44.f / 255.f), float(44.f / 255.f), float(44.f / 255.f));
  ww::tup const Red = ww::Color(1.f, 0.f, 0.f);
  ww::tup const Red2 = ww::Color(1.f, .4f, 0.f);
  ww::tup const Green = ww::Color(0.2982f, 0.6039f, 0.2348f);
  ww::tup const Green2 = ww::Color(.2471, 0.4399f, 0.f);
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(1.f, 1.f, .0f);

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::Translation(0.f, 0.f, 0.f)  //!<
                        * ww::RotateY(0.25f * M_PI_2)   //!<
      ;                                                 //!<

  float const Alpha = M_PI * float(45.f / 180.f);
  ptrPlane->Material.Shininess = 1.f;
  ptrPlane->Material.Diffuse = 0.5f;
  ptrPlane->Material.Color = ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  ww::pattern P1 = ww::StripePattern(Green, White);
  ww::pattern P2 = ww::StripePattern(Green, Green2);
  // P1.Continue = false;
  P1.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f, 1.f * Alpha, 0.f);
  P2.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f, -1.f * Alpha, 0.f);
  ptrPlane->Material.Pattern = ww::BlendedPattern(P1, P2);
  World.vPtrObjects.push_back(ptrPlane);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(0.f, 50.f, 0.f), ww::Color(.7f, .7f, .7f));
  World.vPtrLights.push_back(pLight);

  float const FieldOfView = 75.f / 180.f * M_PI;
  ww::camera Camera = ww::Camera(256, 256, FieldOfView);

  ww::tup const ViewFrom = ww::Point(0.f, 1.5f, -5.f);
  ww::tup const ViewTo = ww::Point(0.f, 0.f, 25.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  // Assert(false, __FUNCTION__, __LINE__);
  ww::WriteToPPM(Canvas, "Ch10PuttingItTogetherBlendedPattern.ppm");
}

//------------------------------------------------------------------------------
TEST(Ch10Patterns, PuttingItTogetherPerturbedPattern)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  ww::tup const White = ww::Color(1.f, 1.f, 1.f);
  ww::tup const Black = ww::Color(0.f, 0.f, 0.f);
  ww::tup const Grey = ww::Color(float(84.f / 255.f), float(84.f / 255.f), float(84.f / 255.f));
  ww::tup const Grey2 = ww::Color(float(44.f / 255.f), float(44.f / 255.f), float(44.f / 255.f));
  ww::tup const Red = ww::Color(1.f, 0.f, 0.f);
  ww::tup const Red2 = ww::Color(1.f, .4f, 0.f);
  ww::tup const Green = ww::Color(0.2982f, 0.6039f, 0.2348f);
  ww::tup const Green2 = ww::Color(.2471, 0.4399f, 0.f);
  ww::tup const Blue = ww::Color(0.f, 0.f, 1.f);
  ww::tup const Yellow = ww::Color(1.f, 1.f, .0f);
  float const Alpha = M_PI * float(45.f / 180.f);

  ww::pattern P1 = ww::StripePattern(Blue, Yellow);
  P1.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, .3f, 1.f, 1.f, -Alpha, Alpha, 0.f);
  ww::pattern P1Perturbed = ww::PerturbPattern(P1);
  EXPECT_EQ(P1.Transform == P1Perturbed.ptrNext->Transform, true);
  EXPECT_EQ(P1.funcPtrPatternAt == P1Perturbed.ptrNext->funcPtrPatternAt, true);

  ww::pattern const P2 = ww::RingPattern(Red, Red2);
  ww::pattern const P2Perturbed = ww::PerturbPattern(P2);

  // ---
  // Add the first plane
  // ---
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::Translation(0.f, 0.f, 0.f)  //!<
                        * ww::RotateY(0.25f * M_PI_2)   //!<
      ;                                                 //!<

  ptrPlane->Material.Shininess = 1.f;
  ptrPlane->Material.Diffuse = 0.5f;
  ptrPlane->Material.Color = ww::Color(float(0xff) / float(0xff), float(0xe9) / float(0xff), float(0xca) / float(0xff));
  // ww::pattern PlanePatternP1 = ww::StripePattern(Green, White);
  // ww::pattern PlanePatternP2 = ww::StripePattern(Green, Green2);
  // // P1.Continue = false;
  // PlanePatternP1.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f, 1.f * Alpha, 0.f);
  // PlanePatternP2.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f, -1.f * Alpha, 0.f);
  // ptrPlane->Material.Pattern = ww::BlendedPattern(PlanePatternP1, PlanePatternP2);
  // ptrPlane->Material.Pattern = ww::CheckersPattern(Green, Blue);
  ptrPlane->Material.Pattern = P2Perturbed;
  World.vPtrObjects.push_back(ptrPlane);

  // ---
  // Add Sphere
  // ---
  ww::shared_ptr_shape PtrMiddle = ww::PtrDefaultSphere();
  ww::shape &Middle = *PtrMiddle;
  Middle.Transform = ww::Translation(0.f, 1.f, 0.0f);
  Middle.Material.Color = ww::Color(0.1f, 1.0f, 0.5f);
  Middle.Material.Diffuse = 0.7f;
  Middle.Material.Specular = 0.3f;
  Middle.Material.Pattern = P1Perturbed;
  Middle.Material.Pattern.Transform = ww::Scaling(0.15f, 0.5f, 0.5f) * ww::RotateZ(0.78);
  World.vPtrObjects.push_back(PtrMiddle);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-10.f, 25.f, 0.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  float const FieldOfView = 75.f / 180.f * M_PI;
  ww::camera Camera = ww::Camera(256, 256, FieldOfView);

  ww::tup const ViewFrom = ww::Point(0.f, 1.5f, -3.5f);
  ww::tup const ViewTo = ww::Point(0.f, 0.f, 25.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  // Assert(false, __FUNCTION__, __LINE__);
  ww::WriteToPPM(Canvas, "Ch10PuttingItTogetherPerturbedPattern.ppm");
}
