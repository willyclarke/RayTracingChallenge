#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, SphereNormals)
{
  ww::sphere S{};
  {
    ww::tup const N = ww::NormalAt(S, ww::Point(1.f, 0.f, 0.f));
    EXPECT_EQ(ww::Equal(N, ww::Vector(1.f, 0.f, 0.f)), true);
  }
  {
    ww::tup const N = ww::NormalAt(S, ww::Point(0.f, 1.f, 0.f));
    EXPECT_EQ(ww::Equal(N, ww::Vector(0.f, 1.f, 0.f)), true);
  }
  {
    ww::tup const N = ww::NormalAt(S, ww::Point(0.f, 0.f, 1.f));
    EXPECT_EQ(ww::Equal(N, ww::Vector(0.f, 0.f, 1.f)), true);
  }
  {
    float const Sqrt3O3 = std::sqrt(3.f) / 3.f;
    ww::tup const N = ww::NormalAt(S, ww::Point(Sqrt3O3, Sqrt3O3, Sqrt3O3));
    EXPECT_EQ(ww::Equal(N, ww::Vector(Sqrt3O3, Sqrt3O3, Sqrt3O3)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, SphereNormalsTranslated)
{
  ww::sphere S{};

  {
    // NOTE: Move the sphere up to +y=1
    S.Transform = ww::Translation(0.f, 1.f, 0.f);

    // ---
    // Scenario: Computing the nomal on a translated sphere.
    //           The normal should now be at y= 0.707 for its y and
    //           the x normal should be a negative number, -.707.
    // ---
    ww::tup const N1 = ww::NormalAt(S, ww::Point(0.f, 1.70711f, -0.70711f));
    EXPECT_EQ(ww::Equal(N1, ww::Vector(0.f, 0.70711f, -0.70711f)), true);

    // NOTE: Set up the scaling.
    //       First rotate the normal by some radians
    ww::tup Rot = ww::RotateZ(1.f / 5.f) * N1;

    EXPECT_EQ(ww::IsVector(Rot), true);

    // thereafter apply scaling and rotation.
    S.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.0f, 0.5f, 1.0f, Rot.X, Rot.Y, Rot.Z);

    // ---
    // Scenario: Computing the nomal on a translated sphere.
    //           The normal should now be at y= 0.707 for its y and
    //           the x normal should be a negative number, -.707.
    // ---
    float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
    ww::tup const N2 = ww::NormalAt(S, ww::Point(0.f, Sqrt2O2, -Sqrt2O2));
    EXPECT_EQ(ww::Equal(N2, ww::Vector(0.f, 0.97014f, -0.24254f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, ReflectingVector45Degrees)
{
  {
    //
    // \   A   /
    //  \  |  /
    //   \ | /
    //    \|/
    // NOTE: Think of a ball that travels at 45°, hits the ground and is reflected.
    //       The ball is travelling in the positive x and downward with a negative z.
    //       The normal will then point straight up from the ground, hence y=1.
    ww::tup const V = ww::Vector(1.f, -1.f, 0.f);
    ww::tup const N = ww::Vector(0.f, 1.f, 0.f);
    ww::tup const R = ww::Reflect(V, N);
    EXPECT_EQ(ww::Equal(R, ww::Vector(1.f, 1.f, 0.f)), true);
  }
  {
    // The surface is tilting 45°. Ball is falling straight down...
    // The result is that the ball travels straight in the X direction.
    ww::tup const V = ww::Vector(0.f, -1.f, 0.f);
    float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
    ww::tup const N = ww::Vector(Sqrt2O2, Sqrt2O2, 0.f);
    ww::tup const R = ww::Reflect(V, N);
    EXPECT_EQ(ww::Equal(R, ww::Vector(1.f, 0.f, 0.f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, PointLightPositionIntensity)
{
  ww::tup const Position = ww::Point(0.f, 0.f, 0.f);
  ww::tup const Intensity = ww::Color(1.f, 1.f, 1.f);
  ww::light const Light = ww::PointLight(Position, Intensity);
  EXPECT_EQ(ww::Equal(Position, Light.Position), true);
  EXPECT_EQ(ww::Equal(Intensity, Light.Intensity), true);
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, MaterialInitialization)
{
  ww::material const Material{};
  float const Ambient{0.1f};     //!< Typical value between 0 and 1. Non-negative.
  float const Diffuse{0.9f};     //!< Typical value between 0 and 1. Non-negative.
  float const Specular{0.9f};    //!< Typical value between 0 and 1. Non-negative.
  float const Shininess{200.f};  //!< Typical value between 10 and 200. Non-negative.

  EXPECT_EQ(Material.Ambient, Ambient);
  EXPECT_EQ(Material.Diffuse, Diffuse);
  EXPECT_EQ(Material.Specular, Specular);
  EXPECT_EQ(Material.Shininess, Shininess);
  EXPECT_EQ(ww::Equal(Material.Color, ww::Color(1.f, 1.f, 1.f)), true);
  EXPECT_EQ(ww::IsPoint(Material.Color), false);
  EXPECT_EQ(ww::IsVector(Material.Color), true);
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, SphereMaterialInitialization)
{
  ww::sphere const S{};
  ww::material const M{};

  // NOTE: Thest that the sphere default material is the same as the default material.
  EXPECT_EQ(M.Ambient, S.Material.Ambient);
  EXPECT_EQ(M.Diffuse, S.Material.Diffuse);
  EXPECT_EQ(M.Specular, S.Material.Specular);
  EXPECT_EQ(M.Shininess, S.Material.Shininess);
  EXPECT_EQ(ww::Equal(M.Color, ww::Color(1.f, 1.f, 1.f)), true);
  EXPECT_EQ(ww::Equal(M.Color, S.Material.Color), true);
  EXPECT_EQ(ww::IsPoint(M.Color), false);
  EXPECT_EQ(ww::IsVector(M.Color), true);
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, SphereMaterialAssignment)
{
  ww::sphere S{};
  ww::material M{};

  // NOTE: Test assignment of material.
  M.Ambient = 1.f;
  S.Material = M;

  EXPECT_EQ(M.Ambient, S.Material.Ambient);
  EXPECT_EQ(M.Diffuse, S.Material.Diffuse);
  EXPECT_EQ(M.Specular, S.Material.Specular);
  EXPECT_EQ(M.Shininess, S.Material.Shininess);
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, SphereLighting)
{
  // Background: Given the default material and a position.
  ww::material const M{};
  ww::tup const Position = ww::Point(0.f, 0.f, 0.f);

  {
    // Scenario: Lighting with the eye between the light and the surface.
    //                              |
    // sun/light O ->  eye <o   <---| Normal vector
    //                              |
    ww::tup const vEye = ww::Vector(0.f, 0.f, -1.f);
    ww::tup const vN = ww::Vector(0.f, 0.f, -1.f);

    ww::light const Light = ww::PointLight(ww::Point(0.f, 0.f, -10.f), ww::Color(1.f, 1.f, 1.f));

    ww::tup Result = ww::Lighting(M, ww::sphere{}, Light, Position, vEye, vN);
    EXPECT_EQ(ww::Equal(Result, ww::Color(1.9f, 1.9f, 1.9f)), true);
  }

  {
    // Scenario: Lighting with the eye between the light and the surface, eye offset 45°.
    //                       eye <o |
    //                             \|
    // sun/light O ->           <---| Normal vector
    //                              |
    float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
    ww::tup const vEye = ww::Vector(0.f, Sqrt2O2, -Sqrt2O2);
    ww::tup const vN = ww::Vector(0.f, 0.f, -1.f);

    ww::light const Light = ww::PointLight(ww::Point(0.f, 0.f, -10.f), ww::Color(1.f, 1.f, 1.f));

    ww::tup Result = ww::Lighting(M, ww::sphere{}, Light, Position, vEye, vN);
    EXPECT_EQ(ww::Equal(Result, ww::Color(1.0f, 1.0f, 1.0f)), true);
  }

  {
    // Scenario: Lighting with the eye opposite surface, eye offset 45° each.
    //                  Sun/light O |
    //                             \|
    //              eye <o      <---| Normal vector
    //                              |
    //                              |
    float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
    float const ColorValue = 0.1f + 0.9f * Sqrt2O2 + 0.f;  // 0.7364
    ww::tup const vEye = ww::Vector(0.f, 0.f, -1.f);
    ww::tup const vN = ww::Vector(0.f, 0.f, -1.f);

    ww::light const Light = ww::PointLight(ww::Point(0.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));

    ww::tup Result = ww::Lighting(M, ww::sphere{}, Light, Position, vEye, vN);

    // NOTE: From the book:
    // Because the angle between the light and the normal vectors has changed the diffuse
    // component becomes 0.9xSqrt2O2. The specular component again falls off to zero, so the total
    // intensity becomes
    // 0.1 + 0.9xSqrt2O2 + 0 = 0.7364
    //
    EXPECT_EQ(ww::Equal(Result, ww::Color(ColorValue, ColorValue, ColorValue)), true);
  }

  {
    // Scenario: Lighting with the eye opposite surface, eye/light offset 45° each.
    //                       eye <o |
    //                             \|
    //                          <---| Normal vector
    //                             /|
    //              Sun/light     O |
    float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
    float const ColorValue = 0.1f + 0.9f * Sqrt2O2 + 0.f;
    ww::tup const vEye = ww::Vector(0.f, 0.f, -1.f);
    ww::tup const vN = ww::Vector(0.f, 0.f, -1.f);

    ww::light const Light = ww::PointLight(ww::Point(0.f, -10.f, -10.f), ww::Color(1.f, 1.f, 1.f));

    ww::tup Result = ww::Lighting(M, ww::sphere{}, Light, Position, vEye, vN);
    EXPECT_EQ(ww::Equal(Result, ww::Color(ColorValue, ColorValue, ColorValue)), true);
  }

  {
    // Scenario: Lighting with the eye in the path of the reflection vector.
    //               Sun/light    O |
    //                             \|
    //                          <---| Normal vector
    //                             /|
    //                 eye       <o |
    float const Sqrt2O2 = std::sqrt(2.f) / 2.f;
    float const ColorValue = 0.1f + 0.9f * Sqrt2O2 + 0.9f;  // 1.6364
    ww::tup const vEye = ww::Vector(0.f, -Sqrt2O2, -Sqrt2O2);
    ww::tup const vN = ww::Vector(0.f, 0.f, -1.f);

    ww::light const Light = ww::PointLight(ww::Point(0.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));

    ww::tup Result = ww::Lighting(M, ww::sphere{}, Light, Position, vEye, vN);
    EXPECT_EQ(ww::Equal(Result, ww::Color(ColorValue, ColorValue, ColorValue)), true);
  }

  {
    // Lighting with the light behind the surface.
    //                              |
    //                             \|
    //                 eye <o   <---| Normal vector   <---- O Sun/light
    //                              |
    //                              |
    // The light does not illuminate the surface. So the diffuse and specular components goes to zero.
    // The total intensity should therefore be the same as the ambient color.
    float const ColorValue = 0.1f;
    ww::tup const vEye = ww::Vector(0.f, 0.f, -1.f);
    ww::tup const vN = ww::Vector(0.f, 0.f, -1.f);

    ww::light const Light = ww::PointLight(ww::Point(0.f, 0.f, 10.f), ww::Color(1.f, 1.f, 1.f));

    ww::tup Result = ww::Lighting(M, ww::sphere{}, Light, Position, vEye, vN);
    EXPECT_EQ(ww::Equal(Result, ww::Color(ColorValue, ColorValue, ColorValue)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch6LightAndShading, SpherePuttingItTogether)
{
  // ---
  // NOTE: Create a lambda for the ray/sphere intersection.
  //       This allows us to change the transform of the sphere
  //       to create interesting effects, like scaling, skewing etc.
  // ---
  auto CreateSphere = [&](ww::shared_ptr_shape S, std::string const &FileName) -> void {
    ww::tup const RayOrigin = ww::Point(0.f, 0.f, -5.f);
    ww::tup Color = ww::Color(1.f, 0.f, 0.f);
    ww::tup const LightPosition = ww::Point(-10.f, 10.f, -10.f);
    ww::tup const LightColor = ww::Color(1.f, 1.f, 1.f);
    ww::light const Light = ww::PointLight(LightPosition, LightColor);

    float const WallZ{10.f};
    float const WallSize{7.f};  // Assume a square wall
    float const Half{WallSize / 2.f};

    int const CanvasPixels{100};
    float const PixelSize{WallSize / CanvasPixels};
    ww::canvas Canvas(CanvasPixels, CanvasPixels);

    for (size_t IdxY = 0;          ///<!
         IdxY < CanvasPixels - 1;  ///<!
         ++IdxY)
    {
      // NOTE: Compute the world y coordinate : top = +half, bottom = -half
      float const WorldY = Half - PixelSize * IdxY;

      for (size_t IdxX = 0;          ///<!
           IdxX < CanvasPixels - 1;  ///<!
           ++IdxX)
      {
        // NOTE: Compute the world x coordinate : left = -half, right = half
        float const WorldX = -Half + PixelSize * IdxX;

        ww::tup const Position = ww::Point(WorldX, WorldY, WallZ);
        ww::ray const R = ww::Ray(Position, ww::Normalize(Position - RayOrigin));

        ww::intersections const XS = ww::Intersect(S, R);
        if (XS.Count())
        {
          ww::tup const Point = ww::PositionAt(R, XS.vI[0].t);
          ww::tup const vNormal = ww::NormalAt(*S, Point);
          ww::tup const &vEye = -R.Direction;

          Color = ww::Lighting(S->Material, ww::sphere{}, Light, Point, vEye, vNormal);

          // std::cout << "Got hit at " << WorldX << "," << WorldY << ". PixelPos:" << IdxX << "," << IdxY << std::endl;
          ww::WritePixel(Canvas, IdxX, IdxY, Color);
        }
      }
    }
    ww::WriteToPPM(Canvas, FileName);
    // NOTE: Just a dummy test.
    EXPECT_EQ(WallZ, 10.f);
  };

  ww::shared_ptr_shape S = ww::PtrDefaultSphere();
  EXPECT_EQ(ww::Equal(S->Material.Color, ww::material().Color), true);

  // Is this one really necessary
  S->Material = ww::material();
  S->Material.Color = ww::Color(1.f, 0.2f, 1.f);

  CreateSphere(S, "SpherePuttingItTogetherCh6.ppm");

#if (1)

  // NOTE: Scale the sphere
  S->Transform = ww::Scaling(1.f, 0.5f, 1.f);
  CreateSphere(S, "SphereScaledYCh6.ppm");

  S->Transform = ww::Scaling(0.5f, 1.0f, 1.f);
  CreateSphere(S, "SphereScaledXCh6.ppm");

  // NOTE: Scaling in the Z direction is not visible from this viewpoint.
  //       So the end result should still be a circle.
  S->Transform = ww::Scaling(1.0f, 1.0f, 0.5f);
  CreateSphere(S, "SphereScaledZCh6.ppm");

  S->Transform = ww::RotateZ(3.1415 / 4.) * ww::Scaling(0.5f, 1.f, 1.f);
  CreateSphere(S, "SphereScaledX1Ch6.ppm");
#endif
}

