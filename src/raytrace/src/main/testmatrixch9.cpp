#include "gtest/gtest.h"
#include <datastructures.hpp>

#include <cmath>
#include <memory>
//------------------------------------------------------------------------------
// Scenario: The default transformation
TEST(Ch9Planes, TheDefaultTransformation)
{
  ww::shape S = ww::TestShape();
  EXPECT_EQ(S.Transform == ww::I(), true);
}

//------------------------------------------------------------------------------
// Scenario: Assigning a transformation
TEST(Ch9Planes, AssigningATransformation)
{
  ww::shape S = ww::TestShape();

  // NOTE: Test that we can assign a translation.
  S.Transform = ww::Translation(2.f, 3.f, 4.f);
  EXPECT_EQ(S.Transform == ww::I(), false);
  EXPECT_EQ(S.Transform == ww::Translation(2.f, 3.f, 4.f), true);
}

//------------------------------------------------------------------------------
// Scenario: The default material
TEST(Ch9Planes, TheDefaultMaterial)
{
  ww::shape S = ww::TestShape();

  EXPECT_EQ(S.Transform == ww::I(), true);
  EXPECT_EQ(S.Material == ww::material{}, true);
}

//------------------------------------------------------------------------------
// Scenario: Assigning a material
TEST(Ch9Planes, AssignAMaterial)
{
  ww::shape S = ww::TestShape();
  ww::material M{};

  M.Ambient = 1.f;

  // NOTE: test that the materials are different before assignment
  EXPECT_EQ(S.Material == M, false);

  // NOTE: test that we can assign a material.
  S.Material = M;

  EXPECT_EQ(S.Material == M, true);
}

//------------------------------------------------------------------------------
// Scenario: Intersecting a scaled shape with a ray
TEST(Ch9Planes, IntersectingAScaledAShapeWithRay)
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
TEST(Ch9Planes, IntersectingATranslatedShapeWithARay)
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
TEST(Ch9Planes, ComputingTheNormalOnATranslatedSphere)
{
  ww::shape S = ww::TestShape();
  S.Transform = ww::Translation(0.f, 1.f, 0.f);

  ww::tup const N = ww::NormalAt(S, ww::Point(0.f, 1.70711f, -0.70711f));
  EXPECT_EQ(N == ww::Vector(0.f, 0.70711f, -0.70711f), true);
}

//------------------------------------------------------------------------------
// Scenario: Computing the normal on a transformed shape
TEST(Ch9Planes, ComputingTheNormalOnATransformedSphere)
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
TEST(Ch9Planes, CheckThatASphereIsAShape)
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
TEST(Ch9Planes, CheckThanNormalOnPlaneIsConstant)
{
  ww::plane P{};
  ww::tup N1 = ww::LocalNormalAt(P, ww::Point(0.f, 0.f, 0.f));
  ww::tup N2 = ww::LocalNormalAt(P, ww::Point(10.f, 0.f, -10.f));
  ww::tup N3 = ww::LocalNormalAt(P, ww::Point(-5.f, 0.f, 150.f));
  EXPECT_EQ(N1 == ww::Vector(0.f, 1.f, 0.f), true);
  EXPECT_EQ(N2 == ww::Vector(0.f, 1.f, 0.f), true);
  EXPECT_EQ(N3 == ww::Vector(0.f, 1.f, 0.f), true);
}
