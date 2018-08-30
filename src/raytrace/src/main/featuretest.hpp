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

TEST(Tuples, Subtraction1)
{
   ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
   ww::tup const P2 = ww::Point(5.f, 6.f, 7.f);
   ww::tup const Expect = {-2.f, -4.f, -6.f, 0.f};
   ww::tup const V = ww::Sub(P1, P2);
   EXPECT_EQ(ww::Equal(V, Expect) && IsVector(V), true);
}

TEST(Tuples, Subtraction2)
{
   ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
   ww::tup const V1 = ww::Vector(5.f, 6.f, 7.f);
   ww::tup const Expect = {-2.f, -4.f, -6.f, 1.f};
   ww::tup const P = ww::Sub(P1, V1);
   EXPECT_EQ(ww::Equal(P, Expect) && IsPoint(P), true);
}

TEST(Tuples, Subtraction3)
{
   ww::tup const Z = ww::Vector(0.f, 0.f, 0.f);
   ww::tup const V1 = ww::Vector(-1.f, -2.f, -3.f);
   ww::tup const Expect = ww::Vector(1.f, 2.f, 3.f);
   ww::tup const P = ww::Sub(Z, V1);
   EXPECT_EQ(ww::Equal(P, Expect), true);
   EXPECT_EQ(ww::IsPoint(P), false);
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
