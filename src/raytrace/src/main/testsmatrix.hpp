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
  // std::cout << M22 << std::endl;
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
  // std::cout << M33 << std::endl;
  // M33.Dimension = 4;
  // std::cout << M33 << std::endl;
}
//------------------------------------------------------------------------------
TEST(Matrix, Multiply)
{
  ww::matrix const M{ww::tup{1.f, 2.f, 3.f, 4.f},  //!<
                     ww::tup{2.f, 4.f, 4.f, 2.f},  //!<
                     ww::tup{8.f, 6.f, 4.f, 1.f},  //!<
                     ww::tup{0.f, 0.f, 0.f, 1.f}};
  ww::tup const T{1.f, 2.f, 3.f, 1.f};
  ww::tup const Result = M * T;
  EXPECT_EQ(Result.C[0], 18.f);
  EXPECT_EQ(Result.C[1], 24.f);
  EXPECT_EQ(Result.C[2], 33.f);
  EXPECT_EQ(Result.C[3], 1.f);
}
//------------------------------------------------------------------------------
TEST(Matrix, Identity)
{
  {
    // ---
    // NOTE: Verify that matrix multiplied by Identity matrix gives the same matrix as result.
    // ---
    ww::matrix const M{ww::tup{0.f, 1.f, 2.f, 4.f},   //!<
                       ww::tup{1.f, 2.f, 4.f, 8.f},   //!<
                       ww::tup{2.f, 4.f, 8.f, 16.f},  //!<
                       ww::tup{4.f, 8.f, 16.f, 32.f}};
    ww::matrix const A = M * ww::I();

    EXPECT_EQ(ww::Equal(M, A), true);
  }

  // ---
  // NOTE: Test that multiplication of Identity matrix and a tuple gives the same tuple.
  // ---
  {
    ww::tup A{1.f, 2.f, 3.f, 4.f};
    ww::tup B = ww::I() * A;
    EXPECT_EQ(ww::Equal(A, B), true);
  }
}
//------------------------------------------------------------------------------
TEST(Matrix, Transpose)
{
  ww::matrix const M{
      ww::tup{0.f, 9.f, 3.f, 0.f},  //
      ww::tup{9.f, 8.f, 0.f, 8.f},  //
      ww::tup{1.f, 8.f, 5.f, 3.f},  //
      ww::tup{0.f, 0.f, 5.f, 8.f}   //
  };
  ww::matrix const MT = ww::Transpose(ww::Transpose(M));
  EXPECT_EQ(ww::Equal(M, MT), true);
  ww::matrix A = ww::I();
  ww::matrix AT = ww::Transpose(A);
  EXPECT_EQ(ww::Equal(AT, ww::I()), true);
}
//------------------------------------------------------------------------------
TEST(Matrix, Determinant22)
{
  ww::matrix M22 = ww::Matrix22(ww::tup{1.f, 5.f, 0.f, 0.f}, ww::tup{-3.f, 2.f, 0.f, 0.f});
  float Det = ww::Determinant22(M22);
  EXPECT_EQ(ww::Equal(Det, 17.f), true);
}
//------------------------------------------------------------------------------
TEST(Matrix, Submatrix33)
{
  {
    ww::matrix const M{
        ww::tup{-6.f, 1.f, 1.f, 6.f},  //
        ww::tup{-8.f, 5.f, 8.f, 6.f},  //
        ww::tup{-1.f, 0.f, 8.f, 2.f},  //
        ww::tup{-7.f, 1.f, -1.f, 1.f}  //
    };
    ww::matrix const TR = ww::Matrix33(ww::tup{-6.f, 1.f, 1.f, 0.f},  //
                                       ww::tup{-8.f, 5.f, 8.f, 0.f},  //
                                       ww::tup{-1.f, 0.f, 8.f, 0.f});

    ww::matrix const R = ww::SubMatrix(M, 3, 3);
    EXPECT_EQ(ww::Equal(R, TR), true);
  }
  {
    ww::matrix const M{
        ww::tup{-6.f, 1.f, 1.f, 6.f},  //
        ww::tup{-8.f, 5.f, 8.f, 6.f},  //
        ww::tup{-1.f, 0.f, 8.f, 2.f},  //
        ww::tup{-7.f, 1.f, -1.f, 1.f}  //
    };
    ww::matrix TR = ww::Matrix33(ww::tup{-6.f, 1.f, 6.f, 0.f},  //
                                 ww::tup{-8.f, 8.f, 6.f, 0.f},  //
                                 ww::tup{-7.f, -1.f, 1.f, 0.f});

    ww::matrix R = ww::SubMatrix(M, 2, 1);
    EXPECT_EQ(ww::Equal(R, TR), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Submatrix22)
{
  {
    ww::matrix const M = ww::Matrix33(ww::tup{1.f, 5.f, 0.f, 0.f},   //
                                      ww::tup{-3.f, 2.f, 7.f, 0.f},  //
                                      ww::tup{0.f, 6.f, -3.f, 0.f});
    ww::matrix const R = ww::SubMatrix(M, 0, 2);
    ww::matrix const TR = ww::Matrix22(ww::tup{-3.f, 2.f, 0.f, 0.f},  //
                                       ww::tup{0.f, 6.f, 0.f, 0.f});
    EXPECT_EQ(ww::Equal(R, TR), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Minorx33)
{
  ww::matrix const A = ww::Matrix33(ww::tup{3.f, 5.f, 0.f, 0.f},    //
                                    ww::tup{2.f, -1.f, -7.f, 0.f},  //
                                    ww::tup{6.f, -1.f, 5.f, 0.f});
  ww::matrix const B = ww::SubMatrix(A, 1, 0);
  ww::matrix const TR = ww::Matrix22(ww::tup{5.f, 0.f, 0.f, 0.f},  //
                                     ww::tup{-1.f, 5.f, 0.f, 0.f});
  EXPECT_EQ(ww::Equal(B, TR), true);

  float const CalculatedMinor = Minor(A, 1, 0);
  EXPECT_EQ(CalculatedMinor, 25.f);

  float const DeterminantB = Determinant22(B);
  EXPECT_EQ(DeterminantB, 25.f);
}

//------------------------------------------------------------------------------
TEST(Matrix, Cofactor)
{
  {
    ww::matrix const A = ww::Matrix33(ww::tup{3.f, 5.f, 0.f, 0.f},    //
                                      ww::tup{2.f, -1.f, -7.f, 0.f},  //
                                      ww::tup{6.f, -1.f, 5.f, 0.f});
    EXPECT_EQ(A.Dimension, 3);
    EXPECT_EQ(ww::Minor(A, 0, 0), -12.f);
    EXPECT_EQ(ww::Cofactor22(A, 0, 0), -12.f);
    EXPECT_EQ(ww::Minor(A, 1, 0), 25.f);
    EXPECT_EQ(ww::Cofactor22(A, 1, 0), -25.f);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Determinant33)
{
  {
    ww::matrix const A = ww::Matrix33(ww::tup{1.f, 2.f, 6.f, 0.f},    //
                                      ww::tup{-5.f, 8.f, -4.f, 0.f},  //
                                      ww::tup{2.f, 6.f, 4.f, 0.f});

    EXPECT_EQ(A.Dimension, 3);
    EXPECT_EQ(ww::Cofactor22(A, 0, 0), 56.f);
    EXPECT_EQ(ww::Cofactor22(A, 0, 1), 12.f);
    EXPECT_EQ(ww::Cofactor22(A, 0, 2), -46.f);
    float const Det = A.R0.C[0] * ww::Cofactor22(A, 0, 0) +  //
                      A.R0.C[1] * ww::Cofactor22(A, 0, 1) +  //
                      A.R0.C[2] * ww::Cofactor22(A, 0, 2);
    EXPECT_EQ(Det, -196);
    EXPECT_EQ(ww::Determinant33(A), -196.f);
    EXPECT_EQ(Det, ww::Determinant(A));
  }
  {
    ww::matrix const A = ww::Matrix44(ww::tup{-2.f, -8.f, 3.f, 5.f},  //
                                      ww::tup{-3.f, 1.f, 7.f, 3.f},   //
                                      ww::tup{1.f, 2.f, -9.f, 6.f},   //
                                      ww::tup{-6.f, 7.f, 7.f, -9.f});
    float const D = Determinant44(A);
    EXPECT_EQ(D, -4071.f);
    EXPECT_EQ(D, ww::Determinant(A));
    EXPECT_EQ(ww::Cofactor33(A, 0, 0), 690.f);
    EXPECT_EQ(ww::Cofactor33(A, 0, 1), 447.f);
    EXPECT_EQ(ww::Cofactor33(A, 0, 2), 210.f);
    EXPECT_EQ(ww::Cofactor33(A, 0, 3), 51.f);
  }
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
};  // namespace rtcch3
#endif
