/******************************************************************************
 * Filename : featuretest.hpp
 * Date     : 2018 Aug 29
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Test routines for the Raytracer project.
 ******************************************************************************/
#ifndef SRC_RAYTRACE_SRC_MAIN_FEATURETEST_HPP
#define SRC_RAYTRACE_SRC_MAIN_FEATURETEST_HPP

#include <datastructures.hpp>

#include <cmath>
#include <iostream>

#include "gtest/gtest.h"

//------------------------------------------------------------------------------
TEST(Tuples, Assign)
{
   ww::tup A{};
   ww::tup B{};
   EXPECT_EQ(ww::Equal(A, B), true);
}
//------------------------------------------------------------------------------

TEST(Tuples, IsPointFalse)
{
   ww::tup A = ww::Vector(0.f, 0.f, 0.f);
   EXPECT_EQ(ww::IsPoint(A), false);
}

//------------------------------------------------------------------------------
TEST(Tuples, IsPointTrue)
{
   ww::tup B = ww::Point(4.3, -4.2, 3.1);
   EXPECT_EQ(ww::IsPoint(B), true);
   EXPECT_EQ(ww::IsVector(B), false);
}

//------------------------------------------------------------------------------
TEST(Tuples, IsVectorTrue)
{
   ww::tup B = ww::Vector(4, -4, 3);
   EXPECT_EQ(B.W == 0, true);
}

TEST(Tuples, Equality1)
{
   float const A{1.0};
   float const B = A + 0.5f * ww::EPSILON;
   EXPECT_EQ(ww::Equal(A, B), true);
}

TEST(Tuples, Equality2)
{
   ww::tup const A = ww::Vector(1.f, 2.f, 3.f);
   ww::tup const B = ww::Vector(1.f, 2.f, 3.f);
   EXPECT_EQ(ww::Equal(A, B), true);
}

TEST(Tuples, Equality3)
{
   ww::tup const A = ww::Vector(1.f, 2.f, 3.f);
   ww::tup const B = ww::Vector(1.f, 2.000000001f, 3.f);
   EXPECT_EQ(ww::Equal(A, B), true);
}

TEST(Tuples, Addition)
{
   ww::tup const A = {3.f, -2.f, 5.f, 1.f};
   ww::tup const B = {-2.f, 3.f, 1.f, 0.f};
   ww::tup const Expect = {1.f, 1.f, 6.f, 1.f};
   ww::tup const C = ww::Add(A, B);
   EXPECT_EQ(ww::Equal(C, Expect), true);
}

TEST(Tuples, SubtractionPointFromPoint)
{
   ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
   ww::tup const P2 = ww::Point(5.f, 6.f, 7.f);
   ww::tup const Expect = {-2.f, -4.f, -6.f, 0.f};
   ww::tup const V = ww::Sub(P1, P2);
   EXPECT_EQ(ww::Equal(V, Expect) && IsVector(V), true);
}

TEST(Tuples, SubtractionVectorFromPoint)
{
   ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
   ww::tup const V1 = ww::Vector(5.f, 6.f, 7.f);
   ww::tup const Expect = {-2.f, -4.f, -6.f, 1.f};
   ww::tup const P = ww::Sub(P1, V1);
   EXPECT_EQ(ww::Equal(P, Expect) && IsPoint(P), true);
}

TEST(Tuples, SubtractionVectorFromZeroVector)
{
   ww::tup const Z = ww::Vector(0.f, 0.f, 0.f);
   ww::tup const V1 = ww::Vector(-1.f, -2.f, -3.f);
   ww::tup const Expect = ww::Vector(1.f, 2.f, 3.f);
   ww::tup const P = ww::Sub(Z, V1);
   EXPECT_EQ(ww::Equal(P, Expect), true);
   EXPECT_EQ(ww::IsPoint(P), false);
}

TEST(Tuples, NegateTuple)
{
   ww::tup const T{-1.f, -2.f, 3.f, -4.f};
   ww::tup const N = ww::Negate(T);
   ww::tup const Expect{1.f, 2.f, -3.f, 4.f};
   EXPECT_EQ(ww::Equal(N, Expect), true);
}

TEST(Tuples, NegateTupleOperator)
{
   ww::tup const T{-1.f, -2.f, 3.f, -4.f};
   ww::tup const N = -T;
   ww::tup const Expect{1.f, 2.f, -3.f, 4.f};
   EXPECT_EQ(ww::Equal(N, Expect), true);
}

TEST(Tuples, Multiply)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R{ww::Multiply(3.5f, A)};
   ww::tup const Expect{3.5f, -7.f, 10.5f, -14.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

TEST(Tuples, MultiplyOperatorMul1)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R = A * 3.5f;
   ww::tup const Expect{3.5f, -7.f, 10.5f, -14.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

TEST(Tuples, MultiplyOperatorMul2)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R = 3.5f * A;
   ww::tup const Expect{3.5f, -7.f, 10.5f, -14.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

TEST(Tuples, DivideOperatorDiv1)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R = A / 2.f;
   ww::tup const Expect{0.5f, -1.f, 1.5f, -2.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

TEST(Tuples, MagnitudeSquared)
{
   EXPECT_EQ(ww::Equal(ww::MagSquared(ww::Vector(1.f, 0.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::MagSquared(ww::Vector(0.f, 1.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::MagSquared(ww::Vector(0.f, 0.f, 1.f)), 1.f), true);
}

TEST(Tuples, MagnitudeVectors)
{
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(1.f, 0.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(0.f, 1.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(0.f, 0.f, 1.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(1.f, 2.f, 3.f)), std::sqrt(14.f)), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(-1.f, -2.f, -3.f)), std::sqrt(14.f)), true);
}

TEST(Tuples, Normal)
{
   EXPECT_EQ(ww::Equal(ww::Normal(ww::Vector(4.f, 0.f, 0.f)), ww::Vector(1.f, 0.f, 0.f)), true);
   // ---
   // NOTE: The length/magnitude of the vector 1,2,3 is Sqrt(14).
   // ---
   float const Sq14 = std::sqrt(14.f);
   EXPECT_EQ(ww::Equal(ww::Normal(ww::Vector(1.f, 2.f, 3.f)), ww::Vector(1.f / Sq14, 2.f / Sq14, 3.f / Sq14)), true);

   // ---
   // NOTE: Test that the length/magnitude of the normal vector of the vector 1,2,3 is 1.0.
   // ---
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Normal(ww::Vector(1.f, 2.f, 3.f))), 1.f), true);
}

TEST(Tuples, DotProduct)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(2.f, 3.f, 4.f)};
   EXPECT_EQ(ww::Equal(ww::Dot(A, B), 20.f), true);
}

TEST(Tuples, CrossProduct)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(2.f, 3.f, 4.f)};
   ww::tup const Result = ww::Cross(A, B);
   EXPECT_EQ(ww::Equal(Result, ww::Vector(-1.f, 2.f, -1.f)), true);
}

TEST(Tuples, AddOpertor)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const Expect{ww::Add(A, B)};
   ww::tup const Result = A + B;
   EXPECT_EQ(ww::Equal(Expect, Result), true);
}

TEST(Tuples, SubOpertor)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const Expect{ww::Sub(A, B)};
   ww::tup const Result = A - B;
   EXPECT_EQ(ww::Equal(Expect, Result), true);
}

TEST(Colors, AssignColor)
{
   ww::tup const C{1.f, 0.f, 0.f, 1.f};
   EXPECT_EQ(ww::Equal(1.f, C.R) + ww::Equal(0.f, C.G) + ww::Equal(0.f, C.B) + ww::Equal(1.f, C.I),  //<!
             4.f);
}

// ---
// NOTE: Test function for the tuples.
//     : Prints out expected results for some simple tests define in the
//     : raytracing book.
// ---
void RunTupleTest(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   int const Result = RUN_ALL_TESTS();
   std::cout << "Testing result was " << (Result ? "Failure." : "Success.") << std::endl;
}
#endif
