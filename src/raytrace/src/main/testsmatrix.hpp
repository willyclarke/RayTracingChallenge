/******************************************************************************
 * Filename : testsmatrix.hpp
 * Date     : 2018 Sep 15
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Google tests for chapter 3 of the Ray Tracing Challenge.
 ******************************************************************************/
#ifndef SRC_RAYTRACE_SRC_MAIN_TESTSMATRIX_HPP
#define SRC_RAYTRACE_SRC_MAIN_TESTSMATRIX_HPP

#include <datastructures.hpp>

#include <cmath>
#include <iostream>
#include <memory>  // for shared pointer.

#include "gtest/gtest.h"

namespace rtcch3
{
//------------------------------------------------------------------------------
TEST(Matrix, InitializationToZero)
{
  ww::matrix M{};
  EXPECT_EQ(M.Dimension, 4);

  float Sum{};
  for (size_t Row = 0;  ///<!
       Row < 4;         ///<!
       ++Row)
  {
    for (size_t Col = 0;  ///<!
         Col < 4;         ///<!
         ++Col)
    {
      Sum += ww::Get(M, Row, Col);
      EXPECT_EQ(Sum, 0.f);
    }
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, InitializationToValue)
{
  ww::matrix M{ww::tup{0.f, 1.f, 2.f, 3.f},    //!<
               ww::tup{4.f, 5.f, 6.f, 7.f},    //!<
               ww::tup{8.f, 9.f, 10.f, 11.f},  //!<
               ww::tup{12.f, 13.f, 14.f, 15.f}};

  // std::cout << M << std::endl;

  size_t const TupleSize = M.Dimension;
  EXPECT_EQ(TupleSize, 4);
  int Count{};
  for (size_t Row = 0;   ///<!
       Row < TupleSize;  ///<!
       ++Row)
  {
    for (size_t Col = 0;   ///<!
         Col < TupleSize;  ///<!
         ++Col, Count++)
    {
      EXPECT_EQ(ww::Get(M, Row, Col), Count);
    }
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, InitializationAndInspection)
{
  ww::matrix const M{ww::tup{1.f, 2.f, 3.f, 4.f},      //!<
                     ww::tup{5.5f, 6.5f, 7.5f, 8.5f},  //!<
                     ww::tup{9.f, 10.f, 11.f, 12.f},   //!<
                     ww::tup{13.5f, 14.5f, 15.5f, 16.5f}};

  // std::cout << M << std::endl;
  EXPECT_EQ(Get(M, 0, 0), 1.f);
  EXPECT_EQ(Get(M, 0, 3), 4.f);
  EXPECT_EQ(Get(M, 1, 0), 5.5f);
  EXPECT_EQ(Get(M, 1, 2), 7.5f);
  EXPECT_EQ(Get(M, 2, 2), 11.f);
  EXPECT_EQ(Get(M, 3, 0), 13.5f);
  EXPECT_EQ(Get(M, 3, 2), 15.5f);
}

//------------------------------------------------------------------------------
TEST(Matrix, TwoByTwo)
{
  ww::matrix M22 = Matrix22(ww::tup{1.f, 2.f}, ww::tup{3.f, 4.f});
  EXPECT_EQ(Get(M22, 0, 0), 1.f);
  EXPECT_EQ(Get(M22, 0, 1), 2.f);
  EXPECT_EQ(Get(M22, 1, 0), 3.f);
  EXPECT_EQ(Get(M22, 1, 1), 4.f);
  //std::cout << M22 << std::endl;
}

//------------------------------------------------------------------------------
TEST(Matrix, ThreeByThree)
{
  ww::matrix M33 = Matrix33(ww::tup{1.f, 2.f, 3.f}, ww::tup{4.f, 5.f, 6.f}, ww::tup{7.f, 8.f, 9.f});
  EXPECT_EQ(Get(M33, 0, 0), 1.f);
  EXPECT_EQ(Get(M33, 0, 1), 2.f);
  EXPECT_EQ(Get(M33, 0, 2), 3.f);
  EXPECT_EQ(Get(M33, 1, 0), 4.f);
  EXPECT_EQ(Get(M33, 1, 1), 5.f);
  EXPECT_EQ(Get(M33, 1, 2), 6.f);
  EXPECT_EQ(Get(M33, 2, 0), 7.f);
  EXPECT_EQ(Get(M33, 2, 1), 8.f);
  EXPECT_EQ(Get(M33, 2, 2), 9.f);
  //std::cout << M33 << std::endl;
  //M33.Dimension = 4;
  //std::cout << M33 << std::endl;
}
//------------------------------------------------------------------------------
void RunMatrixTest(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  int Count{};
  for (size_t Idx = 0;  //<!
       Idx < 1;         //<!
       ++Idx)
  {
    Count++;
    int const Result = RUN_ALL_TESTS();
    std::cout << "Testing result was " << (Result ? "Failure." : "Success.")  //<!
              << ". Test count " << Count << std::endl;
    if (Result != 0) break;
  }
}
}  // namespace rtcch3
#endif
