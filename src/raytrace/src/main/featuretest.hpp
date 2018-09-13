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
#include <memory>  // for shared pointer.

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

//------------------------------------------------------------------------------
TEST(Tuples, Equality1)
{
   float const A{1.0};
   float const B = A + 0.5f * ww::EPSILON;
   EXPECT_EQ(ww::Equal(A, B), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, Equality2)
{
   ww::tup const A = ww::Vector(1.f, 2.f, 3.f);
   ww::tup const B = ww::Vector(1.f, 2.f, 3.f);
   EXPECT_EQ(ww::Equal(A, B), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, Equality3)
{
   ww::tup const A = ww::Vector(1.f, 2.f, 3.f);
   ww::tup const B = ww::Vector(1.f, 2.000000001f, 3.f);
   EXPECT_EQ(ww::Equal(A, B), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, Addition)
{
   ww::tup const A = {3.f, -2.f, 5.f, 1.f};
   ww::tup const B = {-2.f, 3.f, 1.f, 0.f};
   ww::tup const Expect = {1.f, 1.f, 6.f, 1.f};
   ww::tup const C = ww::Add(A, B);
   EXPECT_EQ(ww::Equal(C, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, SubtractionPointFromPoint)
{
   ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
   ww::tup const P2 = ww::Point(5.f, 6.f, 7.f);
   ww::tup const Expect = {-2.f, -4.f, -6.f, 0.f};
   ww::tup const V = ww::Sub(P1, P2);
   EXPECT_EQ(ww::Equal(V, Expect) && IsVector(V), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, SubtractionVectorFromPoint)
{
   ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
   ww::tup const V1 = ww::Vector(5.f, 6.f, 7.f);
   ww::tup const Expect = {-2.f, -4.f, -6.f, 1.f};
   ww::tup const P = ww::Sub(P1, V1);
   EXPECT_EQ(ww::Equal(P, Expect) && IsPoint(P), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, SubtractionVectorFromZeroVector)
{
   ww::tup const Z = ww::Vector(0.f, 0.f, 0.f);
   ww::tup const V1 = ww::Vector(-1.f, -2.f, -3.f);
   ww::tup const Expect = ww::Vector(1.f, 2.f, 3.f);
   ww::tup const P = ww::Sub(Z, V1);
   EXPECT_EQ(ww::Equal(P, Expect), true);
   EXPECT_EQ(ww::IsPoint(P), false);
}

//------------------------------------------------------------------------------
TEST(Tuples, NegateTuple)
{
   ww::tup const T{-1.f, -2.f, 3.f, -4.f};
   ww::tup const N = ww::Negate(T);
   ww::tup const Expect{1.f, 2.f, -3.f, 4.f};
   EXPECT_EQ(ww::Equal(N, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, NegateTupleOperator)
{
   ww::tup const T{-1.f, -2.f, 3.f, -4.f};
   ww::tup const N = -T;
   ww::tup const Expect{1.f, 2.f, -3.f, 4.f};
   EXPECT_EQ(ww::Equal(N, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, Multiply)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R{ww::Multiply(3.5f, A)};
   ww::tup const Expect{3.5f, -7.f, 10.5f, -14.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, MultiplyOperatorMul1)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R = A * 3.5f;
   ww::tup const Expect{3.5f, -7.f, 10.5f, -14.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, MultiplyOperatorMul2)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R = 3.5f * A;
   ww::tup const Expect{3.5f, -7.f, 10.5f, -14.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, DivideOperatorDiv1)
{
   ww::tup const A{1.f, -2.f, 3.f, -4.f};
   ww::tup const R = A / 2.f;
   ww::tup const Expect{0.5f, -1.f, 1.5f, -2.0f};
   EXPECT_EQ(ww::Equal(R, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, MagnitudeSquared)
{
   EXPECT_EQ(ww::Equal(ww::MagSquared(ww::Vector(1.f, 0.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::MagSquared(ww::Vector(0.f, 1.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::MagSquared(ww::Vector(0.f, 0.f, 1.f)), 1.f), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, MagnitudeVectors)
{
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(1.f, 0.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(0.f, 1.f, 0.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(0.f, 0.f, 1.f)), 1.f), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(1.f, 2.f, 3.f)), std::sqrt(14.f)), true);
   EXPECT_EQ(ww::Equal(ww::Mag(ww::Vector(-1.f, -2.f, -3.f)), std::sqrt(14.f)), true);
}

//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
TEST(Tuples, DotProduct)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(2.f, 3.f, 4.f)};
   EXPECT_EQ(ww::Equal(ww::Dot(A, B), 20.f), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, CrossProduct)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(2.f, 3.f, 4.f)};
   ww::tup const Result = ww::Cross(A, B);
   EXPECT_EQ(ww::Equal(Result, ww::Vector(-1.f, 2.f, -1.f)), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, AddOpertor)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const Expect{ww::Add(A, B)};
   ww::tup const Result = A + B;
   EXPECT_EQ(ww::Equal(Expect, Result), true);
}

//------------------------------------------------------------------------------
TEST(Tuples, SubOpertor)
{
   ww::tup const A{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const B{ww::Vector(1.f, 2.f, 3.f)};
   ww::tup const Expect{ww::Sub(A, B)};
   ww::tup const Result = A - B;
   EXPECT_EQ(ww::Equal(Expect, Result), true);
}

//------------------------------------------------------------------------------
TEST(Colors, AssignColor)
{
   {
      ww::tup const C{1.f, 0.f, 0.f, 1.f};
      EXPECT_EQ(ww::Equal(1.f, C.R) + ww::Equal(0.f, C.G) + ww::Equal(0.f, C.B) + ww::Equal(1.f, C.I),  //<!
                4.f);
   }
   {
      ww::tup const C = ww::Color(-0.5f, 0.4f, 1.7f);
      EXPECT_EQ(ww::Equal(-0.5f, C.R) && ww::Equal(0.4f, C.G) &&  //<!
                    ww::Equal(1.7f, C.B) && ww::Equal(0.f, C.I),  //<!
                true);
   }
}

//------------------------------------------------------------------------------
TEST(Colors, AddColors)
{
   ww::tup const C1 = ww::Color(0.9f, 0.6f, 0.75f);
   ww::tup const C2 = ww::Color(0.7f, 0.1f, 0.25f);
   ww::tup const C = C1 + C2;
   ww::tup const Expect = ww::Color(0.9f + 0.7f, 0.6f + 0.1f, 0.75f + 0.25f);
   EXPECT_EQ(ww::Equal(C, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Colors, SubColors)
{
   ww::tup const C1 = ww::Color(0.9f, 0.6f, 0.75f);
   ww::tup const C2 = ww::Color(0.7f, 0.1f, 0.25f);
   ww::tup const C = C1 - C2;
   ww::tup const Expect = ww::Color(0.9f - 0.7f, 0.6f - 0.1f, 0.75f - 0.25f);
   EXPECT_EQ(ww::Equal(C, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Colors, MulColors)
{
   ww::tup const C1 = ww::Color(1.f, 0.2f, 0.4f);
   ww::tup const C2 = ww::Color(0.9f, 1.f, 0.1f);
   ww::tup const C = C1 * C2;
   ww::tup const Expect = ww::Color(0.9f, 0.2f, 0.04f);
   EXPECT_EQ(ww::Equal(C, Expect), true);
}

//------------------------------------------------------------------------------
TEST(Colors, CreateCanvas)
{
   ww::canvas Canvas{};
   EXPECT_EQ(Canvas.vXY.size() > 0, true);

   // NOTE: Check that the inital colors are all 0.f.
   float Sum{};
   for (size_t Idx = 0;           //<!
        Idx < Canvas.vXY.size();  //<!
        ++Idx)
   {
      Sum += (Canvas.vXY[Idx].R + Canvas.vXY[Idx].G + Canvas.vXY[Idx].B);
   }
   EXPECT_EQ(Sum, 0.f);
}

//------------------------------------------------------------------------------
TEST(Colors, WritePixel)
{
   ww::canvas Canvas{};

   ww::tup const Red{ww::Color(1.f, 0.f, 0.f)};
   WritePixel(Canvas, 2, 3, Red);
   ww::tup const ReadRed = ww::ReadPixel(Canvas, 2, 3);
   EXPECT_EQ(ww::Equal(Red, ReadRed), true);
}

//------------------------------------------------------------------------------
TEST(Canvas, WritePPM)
{
   char const *ptrFilename = "WritePPM.ppm";
   ww::canvas Canvas(5, 3);
   ww::tup const C1{ww::Color(1.f, 0.f, 0.6f)};
   ww::tup const C2{ww::Color(0.f, 0.5f, 0.f)};
   ww::tup const C3{ww::Color(-0.5f, 0.f, 1.f)};
   WritePixel(Canvas, 0, 0, C1);
   WritePixel(Canvas, 2, 1, C2);
   WritePixel(Canvas, 4, 2, C3);
   ww::WriteToPPM(Canvas, ptrFilename);

   // ---
   // NOTE: Now test that things are actually equal.
   // ---
   {
      std::shared_ptr<ww::canvas> ptrCanvas = ww::ReadFromPPM(ptrFilename);
      if (ptrCanvas)
      {
         ww::canvas &Cv = *ptrCanvas;
         EXPECT_EQ(ptrCanvas != 0, true);
         EXPECT_EQ(Cv.W, Canvas.W);
         EXPECT_EQ(Cv.H, Canvas.H);
         EXPECT_EQ(Cv.vXY.size(), Canvas.vXY.size());
         EXPECT_EQ(Cv.vXY[0].R, 1.f);   // se setting of color above.
         EXPECT_EQ(Cv.vXY[0].G, 0.f);   // se setting of color above.
         EXPECT_EQ(Cv.vXY[0].B, 0.6f);  // se setting of color above.
      }
   }
}

//------------------------------------------------------------------------------
TEST(Canvas, WritePPMLong)
{
   char const *ptrFilename = "WritePPMLong.ppm";
   ww::canvas Canvas(10, 2);
   // NOTE: Write the same color to all pixels.
   for (size_t Idx = 0;           //<!
        Idx < Canvas.vXY.size();  //<!
        ++Idx)
   {
      Canvas.vXY[Idx] = ww::Color(1.f, 0.0f, 0.6f);
   }
   ww::WriteToPPM(Canvas, ptrFilename);

   // ---
   // NOTE: Now test that things are actually equal.
   // ---
   {
      std::shared_ptr<ww::canvas> ptrCanvas = ww::ReadFromPPM(ptrFilename);
      if (ptrCanvas)
      {
         ww::canvas &Cv = *ptrCanvas;
         EXPECT_EQ(ptrCanvas != 0, true);
         EXPECT_EQ(Cv.W, Canvas.W);
         EXPECT_EQ(Cv.H, Canvas.H);
         EXPECT_EQ(Cv.vXY.size(), Canvas.vXY.size());
         EXPECT_EQ(Cv.vXY[0].R, 1.f);   // se setting of color above.
         EXPECT_EQ(Cv.vXY[0].G, 0.f);   // se setting of color above.
         EXPECT_EQ(Cv.vXY[0].B, 0.6f);  // se setting of color above.
      }
   }
}

//------------------------------------------------------------------------------
TEST(Canvas, PPMHeader)
{
   ww::canvas Canvas(1500, 2000);

   size_t const CalcSizeHeader = (sizeof("P3\n") - 1) + std::to_string(Canvas.W).size() + (sizeof(" ") - 1) +
                                 std::to_string(Canvas.H).size() + (sizeof("\n255") - 1);
   std::string const Header = ww::PPMHeader(Canvas);

   EXPECT_EQ(CalcSizeHeader == Header.size(), true);
}

//------------------------------------------------------------------------------
TEST(Canvas, WritePPMLong2)
{
   char const *ptrFilename = "W.ppm";
   ww::canvas Canvas(100, 100);
   for (auto &Color : Canvas.vXY)
   {
      Color = ww::Color(1.f, 0.0f, 0.6f);
   }
   ww::WriteToPPM(Canvas, ptrFilename);

   // ---
   // NOTE: Now test that things are actually equal.
   // ---
   {
      std::shared_ptr<ww::canvas> ptrCanvas = ww::ReadFromPPM(ptrFilename);
      ww::canvas &Cv = *ptrCanvas;
      if (ptrCanvas)
      {
         EXPECT_EQ(ptrCanvas != 0, true);
         EXPECT_EQ(Cv.W, Canvas.W);
         EXPECT_EQ(Cv.H, Canvas.H);
         EXPECT_EQ(Cv.vXY.size(), Canvas.vXY.size());
         EXPECT_EQ(ww::Equal(Cv.vXY[0].R, 1.f), true);   // se setting of color above.
         EXPECT_EQ(ww::Equal(Cv.vXY[0].G, 0.f), true);   // se setting of color above.
         EXPECT_EQ(ww::Equal(Cv.vXY[0].B, 0.6f), true);  // se setting of color above.
      }
   }
}

//------------------------------------------------------------------------------
TEST(Canvas, ReadFromPPM)
{
   char const *ptrFilename = "ReadFromPPM.ppm";
   {
      ww::canvas Canvas(1, 2);
      for (auto &Color : Canvas.vXY)
      {
         Color = ww::Color(1.f, 0.0f, 0.6f);
      }
      ww::WriteToPPM(Canvas, ptrFilename);
   }
   {
      std::shared_ptr<ww::canvas> ptrCanvas = ww::ReadFromPPM(ptrFilename);
      EXPECT_EQ(ptrCanvas != 0, true);
   }
}
// ---
// NOTE: Test function for the tuples.
//     : Prints out expected results for some simple tests define in the
//     : raytracing book.
// ---
void RunTupleTest(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);

   int Count{};
   for (size_t Idx = 0;  //<!
        Idx < 20;        //<!
        ++Idx)
   {
      Count++;
      int const Result = RUN_ALL_TESTS();
      std::cout << "Testing result was " << (Result ? "Failure." : "Success.")  //<!
                << ". Test count " << Count << std::endl;
      if (Result != 0) break;
   }
}
#endif
