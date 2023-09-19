#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>

//------------------------------------------------------------------------------
TEST(Ch14Test1, CreatingANewGroup)
{
  auto Group = ww::Group();
  EXPECT_EQ(Group.Transform == ww::I(), true);
  EXPECT_EQ(Group.isA<ww::shape>(), true);
}

TEST(Ch14Test2, ShapeHasAParentAttribute)
{
  auto const S = ww::TestShape();
  EXPECT_EQ(S.PtrParent == nullptr, true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test3, AddingChildToAGroup)
{
  auto G = ww::Group();
  auto S = ww::TestShape();
  auto pS = ww::SharedPtrShape(S);
  auto pG = ww::SharedPtrGroup(G);
  ww::AddChild(pG, pS);
  EXPECT_EQ(pG->vChildren.empty(), false);  // The group is not empty.
  if (!pG->vChildren.empty())
  {
    EXPECT_EQ(pG->vChildren.front() == pS, true);  // The group includes S (a copy of).
    EXPECT_EQ(pS->PtrParent == pG, true);          // The shape Parent is the group.
  }
}

//------------------------------------------------------------------------------
TEST(Ch14Test4, IntersectingARayWithAnEmptyGroup)
{
  auto G = ww::Group();
  auto pG = ww::SharedPtrGroup(G);
  auto R = ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 0.f, 1.f));
  auto XS = ww::LocalIntersect(pG, R);
  EXPECT_EQ(XS.Count(), 0);
}

//------------------------------------------------------------------------------
TEST(Ch14Test5, IntersectingARayWithAnNonEmptyGroup)
{
  auto G = ww::Group();
  auto pG = ww::SharedPtrGroup(G);
  auto pS1 = ww::PtrDefaultSphere();
  auto pS2 = ww::PtrDefaultSphere();
  pS2->Transform = ww::Translation(0.f, 0.f, -3.f);
  auto pS3 = ww::PtrDefaultSphere();
  pS3->Transform = ww::Translation(5.f, 0.f, 0.f);

  ww::AddChild(pG, pS1);
  ww::AddChild(pG, pS2);
  ww::AddChild(pG, pS3);

  auto R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));

  auto const XS = ww::LocalIntersect(pG, R);

  EXPECT_EQ(XS.Count(), 4);
  if (XS.Count() == 4)
  {
    EXPECT_EQ(XS.vI[0].pShape == pS2, true);
    EXPECT_EQ(XS.vI[1].pShape == pS2, true);
    EXPECT_EQ(XS.vI[2].pShape == pS1, true);
    EXPECT_EQ(XS.vI[3].pShape == pS1, true);
  }
}

//------------------------------------------------------------------------------
TEST(Ch14Test6, IntersectingATransformedGroup)
{
  auto G = ww::Group();
  G.Transform = ww::Scaling(2.f, 2.f, 2.f);
  auto pG = ww::SharedPtrGroup(G);

  auto pS1 = ww::PtrDefaultSphere();
  pS1->Transform = ww::Translation(5.f, 0.f, 0.f);
  ww::AddChild(pG, pS1);

  auto R = ww::Ray(ww::Point(10.f, 0.f, -10.f), ww::Vector(0.f, 0.f, 1.f));
  auto const XS = ww::LocalIntersect(pG, R);
  EXPECT_EQ(XS.Count(), 2);
}

//------------------------------------------------------------------------------
TEST(Ch14Test7, ConvertingAPointFromWorldToObjectSpace)
{
  auto Ptest = ww::Point(-2.f, 0.f, -10.f);
  auto PExpect = ww::Point(0.f, 0.f, -1.f);

  auto Tm = ww::Translation(5.f, 0.f, 0.f);
  auto Sm = ww::Scaling(2.f, 2.f, 2.f);
  auto Rm = ww::RotateY(M_PI_2);
  auto M = Rm * Sm * Tm;
  auto LocalPoint1 = ww::Inverse(M) * Ptest;
  EXPECT_EQ(LocalPoint1 == PExpect, true);

  auto G = ww::Group();

  G.Transform = Rm;
  auto pG1 = ww::SharedPtrGroup(G);

  G.Transform = Sm;

  auto pG2 = ww::SharedPtrGroup(G);
  ww::AddChild(pG1, pG2);

  auto pS = ww::PtrDefaultSphere();
  pS->Transform = Tm;

  ww::AddChild(pG2, pS);

  auto const P = ww::WorldToObject(pS, Ptest);

  EXPECT_EQ(P == PExpect, true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test8, ConvertingANormalFromObjectToWorldSpace)
{
  auto G = ww::Group();
  G.PtrShapeInfo = std::make_shared<ww::shape_info>();

  auto pG1 = ww::SharedPtrGroup(G);
  pG1->Transform = ww::RotateY(M_PI_2);

  auto pG2 = ww::SharedPtrGroup(G);
  pG2->Transform = ww::Scaling(1.f, 2.f, 3.f);
  ww::AddChild(pG1, pG2);

  auto pS = ww::PtrDefaultSphere();
  pS->Transform = ww::Translation(5.f, 0.f, 0.f);
  ww::AddChild(pG2, pS);

  auto const SQRT3_OVER_3 = std::sqrt(3.f) / 3.f;
  auto const I = ww::Vector(SQRT3_OVER_3, SQRT3_OVER_3, SQRT3_OVER_3);

  auto const N = ww::NormalToWorld(pS, I);
  auto const E = ww::Vector(0.2857f, 0.4286f, -0.8571f);

  // std::cout << "I test    : " << I << ". Len: " << ww::Mag(I) << std::endl;
  // std::cout << "N computed: " << N << ". Len: " << ww::Mag(N) << std::endl;
  // std::cout << "N expected: " << E << ". Len: " << ww::Mag(E) << std::endl;

  EXPECT_EQ(N == E, true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test9, FindTheNormalOnAnObjectInAGroup)
{
  auto G = ww::Group();

  auto pG1 = ww::SharedPtrGroup(G);
  pG1->Transform = ww::RotateY(M_PI_2);

  auto pG2 = ww::SharedPtrGroup(G);
  pG2->Transform = ww::Scaling(1.f, 2.f, 3.f);

  ww::AddChild(pG1, pG2);

  auto pS = ww::PtrDefaultSphere();
  pS->Transform = ww::Translation(5.f, 0.f, 0.f);

  ww::AddChild(pG2, pS);
  EXPECT_EQ(pG1->isA<ww::group>(), true);
  EXPECT_EQ(pG2->isA<ww::group>(), true);
  EXPECT_EQ(pS->isA<ww::group>(), false);

  auto const I = ww::Point(1.7321f, 1.1547f, -5.5774f);
  auto const N = ww::NormalAt(*pS, I);
  auto const E = ww::Vector(0.2857f, 0.4286f, -0.8571f);
  // std::cout << "Computed normal N: " << N << std::endl;
  // std::cout << "Expected normal E: " << E << std::endl;

  EXPECT_EQ(N == E, true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test10, PuttingItAllTogether)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  // float constexpr TWidth = 1.5f;
  // float constexpr TDepth = 1.5f;
  // float constexpr TThick = 0.05f;

  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  ptrPlane->Material.Shininess = 10.f;
  ptrPlane->Material.Diffuse = 0.2f;
  ptrPlane->Material.Ambient = 2.f;
  ptrPlane->Material.Reflective = 1.f;
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
  ptrPlane3->Material.Ambient = 1.2f;
  ptrPlane3->Material.Pattern = ww::StripePattern(ColorBrown, ColorBrownLight);
  World.vPtrObjects.push_back(ptrPlane3);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-5.f, 7.5f, -5.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  /**
   * Test case for checking normals.
   */
  {
    for (int N = 0; N < 6; ++N)
    {
      auto MRot = ww::RotateY(float(N) * M_PI / 3.f);
      auto Corner = ww::PtrDefaultSphere();
      Corner->Material.Color = ColorBrown;
      Corner->Transform = MRot * ww::Translation(0.f, 2.f, -1.f) * ww::Scaling(0.25f, 0.25f, 0.25f);

      auto Edge = ww::PtrDefaultCylinder();
      Edge->Material.Color = ColorBrownLight;
      Edge->Closed = true;
      Edge->Minimum = 0.f;
      Edge->Maximum = 1.f;
      /**
       * Matrix multiplication ensures that the shape is first scaled, rotated and finally translated.
       */
      Edge->Transform = MRot * ww::Translation(0.f, 2.f, -1.f) * ww::RotateY(-M_PI / 6.f) * ww::RotateZ(-M_PI_2) *
                        ww::Scaling(0.25f, 1.f, 0.25f);

      World.vPtrObjects.push_back(Corner);
      World.vPtrObjects.push_back(Edge);
    }
  }

  auto HexagonCorner = []() -> ww::shared_ptr_sphere
  {
    auto Corner = ww::PtrDefaultSphere();
    /**
     * Matrix multiplication ensures that the shape is first scaled, then translated.
     */
    Corner->Transform = ww::Translation(0.f, 0.f, -1.f) * ww::Scaling(0.25f, 0.25f, 0.25f);
    return Corner;
  };

  auto HexagonEdge = []() -> ww::shared_ptr_cylinder
  {
    auto Edge = ww::PtrCappedCylinder(0.f, 1.f);
    Edge->Transform = ww::Translation(0.f, 0.f, -1.f) * ww::RotateY(-M_PI / 6.f) * ww::RotateZ(-M_PI_2) *
                      ww::Scaling(0.25f, 1.f, 0.25f);
    return Edge;
  };

  auto HexagonSide = [&]() -> ww::shared_ptr_group
  {
    auto Side = ww::SharedPtrGroup({});
    ww::AddChild(Side, HexagonCorner());
    ww::AddChild(Side, HexagonEdge());

    /**
     * After using shared pointers to create everything we can return the actual object(aka group).
     */
    return Side;
  };

  /**
   * Create a group of six sides to make up the Hexagon.
   */
  auto Hexagon = [&]() -> ww::shared_ptr_group
  {
    auto Hex = ww::SharedPtrGroup({});
    for (int N{}; N < 6; ++N)
    {
      auto Side = HexagonSide();
      // Side->Transform = ww::RotateY(float(N) * M_PI / 3.f);
      Side->Transform = ww::Translation(0.f, 2.f, 0.f) * ww::RotateZ(45.f) * ww::RotateY(float(N) * M_PI / 3.f);
      ww::AddChild(Hex, Side);
    }
    return Hex;
  };

  auto MyHexagon = Hexagon();

  for (auto Side : MyHexagon->vChildren)
  {
    ww::group *SideGroup = dynamic_cast<ww::group *>(Side.get());
    for (int Jdx{}; Jdx < SideGroup->vChildren.size(); ++Jdx)
    {
      auto PtrChild = SideGroup->vChildren[Jdx];
      auto PtrParent = PtrChild->PtrParent;
      auto T = PtrChild->Transform;

      /**
       * Move up until there are no more parents.
       */
      while (PtrParent)
      {
        T = PtrParent->Transform * T;
        PtrParent = PtrParent->PtrParent;
      }

      PtrChild->Transform = T;

      if (PtrChild->isA<ww::sphere>())
      {
        auto Corner = ww::PtrDefaultSphere();
        Corner->Transform = T;
        World.vPtrObjects.push_back(Corner);
      }

      else if (PtrChild->isA<ww::cylinder>())
      {
        auto Side = ww::PtrCappedCylinder(0.f, 1.f);
        Side->Transform = T;
        World.vPtrObjects.push_back(Side);
      }
    }
  }

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(-0.f, 4.5f, -5.f);
  ww::tup const ViewTo = ww::Point(0.0f, -1.f, 5.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch14PuttingItTogether.ppm");
}

//------------------------------------------------------------------------------
TEST(Ch14Test11, AddGroupToWorld)
{
  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();
  EXPECT_EQ(World.Count(), 0);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-5.f, 7.5f, -5.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  auto G = ww::Group();

  auto pG1 = ww::SharedPtrGroup(G);
  pG1->Transform = ww::Translation(0.f, 0.f, 5.f) * ww::RotateZ(M_PI_4) * ww::RotateY(M_PI_2);

  auto pG2 = ww::SharedPtrGroup(G);
  pG2->Transform = ww::Scaling(1.f, 2.f, 3.f);

  ww::AddChild(pG1, pG2);

  auto pS = ww::PtrDefaultSphere();
  pS->Transform = ww::Translation(5.f, 0.f, 0.f);
  auto pS1 = ww::PtrDefaultSphere();
  pS->Transform = ww::Translation(1.f, 1.f, 0.f);

  ww::AddChild(pG2, pS);
  ww::AddChild(pG2, pS1);

  ww::AddGroupToWorld(World, pG1);

  EXPECT_EQ(pG1->isA<ww::group>(), true);
  EXPECT_EQ(pG2->isA<ww::group>(), true);
  EXPECT_EQ(pS->isA<ww::group>(), false);
  EXPECT_EQ(pS1->isA<ww::group>(), false);

  auto HexagonCorner = []() -> ww::shared_ptr_sphere
  {
    auto Corner = ww::PtrDefaultSphere();
    /**
     * Matrix multiplication ensures that the shape is first scaled, then translated.
     */
    Corner->Transform = ww::Translation(0.f, 0.f, -1.f) * ww::Scaling(0.25f, 0.25f, 0.25f);
    return Corner;
  };

  auto HexagonEdge = []() -> ww::shared_ptr_cylinder
  {
    auto Edge = ww::PtrCappedCylinder(0.f, 1.f);
    Edge->Transform = ww::Translation(0.f, 0.f, -1.f) * ww::RotateY(-M_PI / 6.f) * ww::RotateZ(-M_PI_2) *
                      ww::Scaling(0.25f, 1.f, 0.25f);
    return Edge;
  };

  auto HexagonSide = [&]() -> ww::shared_ptr_group
  {
    auto Side = ww::SharedPtrGroup({});
    ww::AddChild(Side, HexagonCorner());
    ww::AddChild(Side, HexagonEdge());

    /**
     * After using shared pointers to create everything we can return the actual object(aka group).
     */
    return Side;
  };

  /**
   * Create a group of six sides to make up the Hexagon.
   */
  auto Hexagon = [&]() -> ww::shared_ptr_group
  {
    auto Hex = ww::SharedPtrGroup({});
    for (int N{}; N < 6; ++N)
    {
      auto Side = HexagonSide();
      Side->Transform = ww::Translation(5.f, 0.f, 0.f) * ww::RotateZ(0.f) * ww::RotateY(float(N) * M_PI / 3.f);
      ww::AddChild(Hex, Side);
    }
    return Hex;
  };

  auto MyHexagon = Hexagon();
  ww::AddGroupToWorld(World, MyHexagon);

  auto MyHexagon2 = Hexagon();
  MyHexagon2->Transform = ww::Translation(0.f, 0.f, 0.f) * ww::RotateY(M_PI_4);
  ww::AddGroupToWorld(World, MyHexagon2);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(-0.f, 4.5f, -15.f);
  ww::tup const ViewTo = ww::Point(0.0f, -1.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch14AddGroupToWorld.ppm");
}

//------------------------------------------------------------------------------
TEST(Ch14Test12, TransformingPatterns)
{
  ww::tup const ColorBrown = ww::Color(float(0x87) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0xba) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  ww::pattern Stripes = ww::StripePattern(ColorBrown, ColorBrownLight);

  ww::material Material{};
  Material.Pattern = ww::StripePattern(ww::Color(1.f, 1.f, 1.f), ww::Color(0.f, 0.f, 0.f));
  Material.Ambient = 1.f;
  Material.Diffuse = 0.f;
  Material.Specular = 0.f;

  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();
  EXPECT_EQ(World.Count(), 0);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-5.f, 7.5f, -5.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  auto G = ww::Group();

  auto pG1 = ww::SharedPtrGroup(G);
  pG1->Transform = ww::Translation(-4.f, -2.f, 5.f) * ww::RotateX(2.f * M_PI_4) * ww::RotateY(2.0f * M_PI_2);
  // pG1->Material.Pattern = ww::CheckersPattern(ww::Color(1.f, 0.5f, 0.f), ww::Color(0.f, 1.f, 0.2f));
  // pG1->Material.Pattern.Transform =
  //     ww::TranslateScaleRotate(0.f, 0.f, 0.f, .4f, .4f, .4f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));
  //
  // auto pG2 = ww::SharedPtrGroup(G);
  // pG2->Transform = ww::Scaling(1.f, 2.f, 3.f);
  //
  // ww::AddChild(pG1, pG2);
  //
  // auto pS = ww::PtrDefaultSphere();
  // pS->Transform = ww::Translation(5.f, 0.f, 0.f);
  // auto pS1 = ww::PtrDefaultSphere();
  // pS->Transform = ww::Translation(1.f, 1.f, 0.f);
  //
  // ww::AddChild(pG2, pS);
  // ww::AddChild(pG2, pS1);

  // ww::AddGroupToWorld(World, pG1);

  // EXPECT_EQ(pG1->isA<ww::group>(), true);
  // EXPECT_EQ(pG2->isA<ww::group>(), true);
  // EXPECT_EQ(pS->isA<ww::group>(), false);
  // EXPECT_EQ(pS1->isA<ww::group>(), false);

  auto HexagonCorner = [&]() -> ww::shared_ptr_sphere
  {
    auto Corner = ww::PtrDefaultSphere();
    /**
     * Matrix multiplication ensures that the shape is first scaled, then translated.
     */
    Corner->Transform = ww::Translation(0.f, 0.f, -1.f) * ww::Scaling(0.25f, 0.25f, 0.25f);
    // Corner->Material.Pattern = ww::CheckersPattern(ColorBrown, ColorBrownLight);
    return Corner;
  };

  auto HexagonEdge = [&]() -> ww::shared_ptr_cylinder
  {
    auto Edge = ww::PtrCappedCylinder(0.f, 1.f);
    Edge->Transform = ww::Translation(0.f, 0.f, -1.f) * ww::RotateY(-M_PI / 6.f) * ww::RotateZ(-M_PI_2) *
                      ww::Scaling(0.25f, 1.f, 0.25f);
    // Edge->Material.Pattern = ww::CheckersPattern(ColorBrown, ColorBrownLight);
    return Edge;
  };

  auto HexagonSide = [&]() -> ww::shared_ptr_group
  {
    auto Side = ww::SharedPtrGroup({});
    ww::AddChild(Side, HexagonCorner());
    ww::AddChild(Side, HexagonEdge());

    /**
     * After using shared pointers to create everything we can return the actual object(aka group).
     */
    return Side;
  };

  /**
   * Create a group of six sides to make up the Hexagon.
   */
  auto Hexagon = [&]() -> ww::shared_ptr_group
  {
    auto Hex = ww::SharedPtrGroup({});
    for (int N{}; N < 6; ++N)
    {
      auto Side = HexagonSide();
      Side->Transform = ww::Translation(0.f, 0.f, 0.f) * ww::RotateZ(0.f) * ww::RotateY(float(N) * M_PI / 3.f);
      ww::AddChild(Hex, Side);
    }
    return Hex;
  };

  auto MyHexagon = Hexagon();
  MyHexagon->Material.Pattern = ww::StripePattern(ColorBrown, ColorBrownLight);
  MyHexagon->Material.Pattern.Transform = ww::Scaling(0.32f, 0.32f, 0.32f) * ww::RotateZ(ww::Radians(95.f));
  ww::AddGroupToWorld(World, MyHexagon);

  auto MyHexagon2 = Hexagon();
  MyHexagon2->Transform =
      ww::TranslateScaleRotate(2.f, 1.5f, 1.f, 1.f, 1.f, 1.f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));
  MyHexagon2->Material.Pattern = ww::CheckersGradientPattern(ww::Color(1.f, 0.5f, 0.f), ww::Color(0.f, 1.f, 0.2f));
  MyHexagon2->Material.Pattern.Transform = ww::Scaling(1.5f, 1.5f, 1.5f) * ww::RotateZ(ww::Radians(45.f));
  ww::AddGroupToWorld(World, MyHexagon2);

  // ---
  // NOTE: Add a primitive object for comparison.
  // ---
  auto pObject1 = ww::PtrCappedCylinder(-1.f, 3.f);
  auto pObject2 = ww::PtrDefaultSphere();
  pObject2->Transform =
      ww::TranslateScaleRotate(-2.f, 1.5f, 1.f, 1.f, 1.f, 1.f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));
  auto pObject3 = ww::PtrDefaultCube();
  pObject3->Transform =
      ww::TranslateScaleRotate(-1.5f, 0.f, 0.f, 1.f, 1.f, 1.f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));

  ww::AddChild(pG1, pObject1);
  ww::AddChild(pG1, pObject2);
  ww::AddChild(pG1, pObject3);
  // ww::AddGroupToWorld(World, pG1);

  auto pG3 = ww::SharedPtrGroup(G);
  ww::AddChild(pG3, pG1);

  auto pObject4 = ww::PtrCappedCylinder(-1.f, 3.f);
  auto pObject5 = ww::PtrDefaultSphere();
  pObject5->Transform =
      ww::TranslateScaleRotate(-2.f, 1.5f, 1.f, 1.f, 1.f, 1.f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));
  auto pObject6 = ww::PtrDefaultCube();
  pObject6->Transform =
      ww::TranslateScaleRotate(-1.5f, 0.f, 0.f, 1.f, 1.f, 1.f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));

  auto pG2 = ww::SharedPtrGroup(G);
  pG2->Transform =
      ww::TranslateScaleRotate(3.f, 0.f, 0.f, .25f, .25f, .25f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));
  ww::AddChild(pG2, pObject4);
  ww::AddChild(pG2, pObject5);
  ww::AddChild(pG2, pObject6);

  ww::AddChild(pG3, pG2);
  pG3->Material.Pattern = ww::CheckersPattern(ww::Color(1.f, 0.5f, 0.f), ww::Color(0.f, 1.f, 0.2f));
  pG3->Material.Pattern.Transform =
      ww::TranslateScaleRotate(0.f, 0.f, 0.f, .4f, .4f, .4f, ww::Radians(0.f), ww::Radians(0.f), ww::Radians(0.f));

  ww::AddGroupToWorld(World, pG3);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  ww::camera Camera = ww::Camera(256, 256, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 5.0f, -5.f);
  ww::tup const ViewTo = ww::Point(0.0f, 0.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch14AddPatternToGroup.ppm");
}

//------------------------------------------------------------------------------
TEST(DISABLED_Ch14Test13, UseBoxes)
{
  ww::tup const ColorBrown = ww::Color(float(0x07) / float(0xff), float(0x63) / float(0xff), float(0x3b) / float(0xff));
  ww::tup const ColorBrownLight =
      ww::Color(float(0x1a) / float(0xff), float(0x8f) / float(0xff), float(0x5d) / float(0xff));

  ww::pattern Stripes = ww::StripePattern(ColorBrown, ColorBrownLight);

  ww::material Material{};
  Material.Pattern = ww::CheckersPattern(ww::Color(1.f, 0.4f, 1.f), ww::Color(1.f, 0.f, 0.f));
  Material.Pattern.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 0.5f, 0.5f, 0.5f, 0.f, ww::Radians(40.f), 0.f);
  Material.Ambient = 1.f;
  Material.Diffuse = 1.f;
  Material.Specular = 0.5f;
  ww::material Material2 = Material;
  Material2.Pattern = ww::GradientPattern(ww::Color(.9f, 0.5f, 1.f), ww::Color(.8f, .4f, 0.f));
  Material2.Pattern.Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 0.5f, 0.5f, 0.5f, 0.f, ww::Radians(40.f), 0.f);

  ww::world World = ww::World();
  World.vPtrLights.clear();
  World.vPtrObjects.clear();
  EXPECT_EQ(World.Count(), 0);

  // Add the first plane
  ww::shared_ptr_plane ptrPlane = ww::PtrDefaultPlane();
  ptrPlane->Transform = ww::TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f, 0.f, 0.f);
  ptrPlane->Material.Shininess = 10.f;
  ptrPlane->Material.Diffuse = 0.2f;
  ptrPlane->Material.Ambient = 2.f;
  ptrPlane->Material.Reflective = 1.f;
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
  ptrPlane3->Material.Ambient = 1.2f;
  ptrPlane3->Material.Pattern = ww::StripePattern(ColorBrown, ColorBrownLight);
  World.vPtrObjects.push_back(ptrPlane3);

  ww::shared_ptr_light pLight{};
  pLight.reset(new ww::light);
  *pLight = ww::PointLight(ww::Point(-5.f, 20.f, -15.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight);

  ww::shared_ptr_light pLight2{};
  pLight2.reset(new ww::light);
  *pLight2 = ww::PointLight(ww::Point(5.f, 40.f, -25.f), ww::Color(1.f, 1.f, 1.f));
  World.vPtrLights.push_back(pLight2);

  auto G = ww::Group();

  auto pG1 = ww::SharedPtrGroup(G);
  pG1->Transform = ww::Translation(0.f, 6.f, -20.f) * ww::RotateX(0.f * M_PI_4) * ww::RotateY(0.f * M_PI_2);
  pG1->Material = Material;

  // ---
  // NOTE: Add a primitive object for comparison.
  // ---
  constexpr int NumObjects = 6;
  auto CreateGroup = [](int NumObjects, float X, float Y, float Z, ww::material const &Material) -> ww::shared_ptr_group
  {
    auto G = ww::Group();
    auto pG = ww::SharedPtrGroup(G);

    pG->Transform = ww::Translation(X, Y, Z) * ww::RotateX(0.f * M_PI_4) * ww::RotateY(0.f * M_PI_2);
    pG->Material = Material;

    for (int Idx = -NumObjects / 2; Idx < NumObjects / 2; ++Idx)
    {
      auto pObject = ww::PtrDefaultCube();
      pObject->Transform = ww::TranslateScaleRotate(float(Idx) * 3.f, 0.f, 0.f, 1.f, 1.f, 1.f, ww::Radians(0.f),
                                                    ww::Radians(0.f), ww::Radians(0.f));
      pObject->Material.Shininess = 10.f;
      pObject->Material.Diffuse = 0.2f;
      pObject->Material.Ambient = 2.f;
      pObject->Material.Reflective = 1.f;

      ww::AddChild(pG, pObject);
    }
    return pG;
  };

  for (int Idy = -3; Idy < 3; ++Idy)
  {
    for (int Idz = -3; Idz < 3; ++Idz)
    {
      auto pG = CreateGroup(NumObjects, 0.f, float(Idy) * 3.f, float(Idz) * 3.f, Material2);
      ww::AddChild(pG1, pG);
    }
  }
  ww::AddGroupToWorld(World, pG1);

  // ---
  // NOTE: Write out the result so that it is possible to see whats going on.
  // ---
  // ww::camera Camera = ww::Camera(2256, 2256, ww::Radians(80.f));
  ww::camera Camera = ww::Camera(556, 556, ww::Radians(50.f));

  ww::tup const ViewFrom = ww::Point(0.f, 20.0f, -49.f);
  ww::tup const ViewTo = ww::Point(0.0f, 0.f, 0.f);
  ww::tup const UpIsY = ww::Vector(0.f, 1.f, 0.f);
  Camera.Transform = ww::ViewTransform(ViewFrom, ViewTo, UpIsY);

  ww::canvas Canvas = ww::Render(Camera, World);
  ww::WriteToPPM(Canvas, "Ch14UseBoxes.ppm");
}

//------------------------------------------------------------------------------
TEST(Ch14Test14, CreateBounds)
{
  auto pBounds = ww::PtrDefaultBounds();
  auto &Bounds = *pBounds.get();
  EXPECT_EQ(Bounds.Minimum == ww::Point(-1.f, -1.f, -1.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(1.f, 1.f, 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test14, VerifyBoundsOfSphere)
{
  auto pSphere = ww::PtrDefaultSphere();
  auto Bounds = ww::Bounds(pSphere);

  EXPECT_EQ(Bounds.Minimum == ww::Point(-1.f, -1.f, -1.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(1.f, 1.f, 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test14, VerifyBoundsOfCube)
{
  auto pCube = ww::PtrDefaultCube();
  auto Bounds = ww::Bounds(pCube);
  EXPECT_EQ(Bounds.Minimum == ww::Point(-1.f, -1.f, -1.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(1.f, 1.f, 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test15, VerifyBoundsOfCappedCylinder)
{
  constexpr float CYLMINIMUM = 0.f;
  constexpr float CYLMAXIMUM = 2.f;
  auto pCylinder = ww::PtrCappedCylinder(CYLMINIMUM, CYLMAXIMUM);
  auto Bounds = ww::Bounds(pCylinder);
  EXPECT_EQ(Bounds.Minimum == ww::Point(-1.f, CYLMINIMUM, -1.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(1.f, CYLMAXIMUM, 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test16, VerifyBoundsOfCylinder)
{
  auto pCylinder = ww::PtrDefaultCylinder();
  auto Bounds = ww::Bounds(pCylinder);
  EXPECT_EQ(Bounds.Minimum == ww::Point(-1.f, -ww::INIFINITY, -1.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(1.f, ww::INIFINITY, 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test17, VerifyBoundsOfCone)
{
  auto pCone = ww::PtrDefaultCone();
  auto Bounds = ww::Bounds(pCone);
  EXPECT_EQ(Bounds.Minimum == ww::Point(-1.f, -ww::INIFINITY, -1.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(1.f, ww::INIFINITY, 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test18, VerifyBoundsOfCappedCone)
{
  constexpr float CONEMINIMUM = 0.f;
  constexpr float CONEMAXIMUM = 2.f;
  auto pCone = ww::PtrCappedCone(CONEMINIMUM, CONEMAXIMUM);
  auto Bounds = ww::Bounds(pCone);
  EXPECT_EQ(Bounds.Minimum == ww::Point(-1.f, CONEMINIMUM, -1.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(1.f, CONEMAXIMUM, 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Ch14Test19, VerifyBoundsOfPlane)
{
  auto pPlane = ww::PtrDefaultPlane();
  auto Bounds = ww::Bounds(pPlane);
  EXPECT_EQ(Bounds.Minimum == ww::Point(-ww::INIFINITY, -ww::INIFINITY, 0.f), true);
  EXPECT_EQ(Bounds.Maximum == ww::Point(ww::INIFINITY, ww::INIFINITY, 0.f), true);
}
