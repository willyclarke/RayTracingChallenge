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
    EXPECT_EQ(ww::Cofactor33(A, 0, 0), -12.f);
    EXPECT_EQ(ww::Minor(A, 1, 0), 25.f);
    EXPECT_EQ(ww::Cofactor33(A, 1, 0), -25.f);
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
    EXPECT_EQ(ww::Cofactor33(A, 0, 0), 56.f);
    EXPECT_EQ(ww::Cofactor33(A, 0, 1), 12.f);
    EXPECT_EQ(ww::Cofactor33(A, 0, 2), -46.f);
    float const Det = A.R0.C[0] * ww::Cofactor33(A, 0, 0) +  //
                      A.R0.C[1] * ww::Cofactor33(A, 0, 1) +  //
                      A.R0.C[2] * ww::Cofactor33(A, 0, 2);
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
    EXPECT_EQ(ww::Cofactor44(A, 0, 0), 690.f);
    EXPECT_EQ(ww::Cofactor44(A, 0, 1), 447.f);
    EXPECT_EQ(ww::Cofactor44(A, 0, 2), 210.f);
    EXPECT_EQ(ww::Cofactor44(A, 0, 3), 51.f);
    EXPECT_EQ(ww::Cofactor(A, 0, 0), 690.f);
    EXPECT_EQ(ww::Cofactor(A, 0, 1), 447.f);
    EXPECT_EQ(ww::Cofactor(A, 0, 2), 210.f);
    EXPECT_EQ(ww::Cofactor(A, 0, 3), 51.f);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Inversion)
{
  {
    ww::matrix const A = ww::Matrix44(ww::tup{6.f, 4.f, 4.f, 4.f},    //
                                      ww::tup{5.f, 5.f, 7.f, 6.f},    //
                                      ww::tup{4.f, -9.f, 3.f, -7.f},  //
                                      ww::tup{9.f, 1.f, 7.f, -6.f});
    EXPECT_EQ(-2120.f, ww::Determinant(A));
    ww::is_invertible_return const Result = ww::IsInvertible(A);
    EXPECT_EQ(true, Result.IsInvertible);
  }  // namespace rtcch3
  {
    ww::matrix const A = ww::Matrix44(ww::tup{-4.f, 2.f, -2.f, -3.f},  //
                                      ww::tup{9.f, 6.f, 2.f, 6.f},     //
                                      ww::tup{0.f, -5.f, 1.f, -5.f},   //
                                      ww::tup{0.f, 0.f, 0.f, 0.f});
    EXPECT_EQ(0.f, ww::Determinant(A));
    ww::is_invertible_return const Result = ww::IsInvertible(A);
    EXPECT_EQ(0.f, Result.Determinant);
    EXPECT_EQ(false, Result.IsInvertible);
  }
  {
    ww::matrix const A = ww::Matrix44(ww::tup{-5.f, 2.f, 6.f, -8.f},  //
                                      ww::tup{1.f, -5.f, 1.f, 8.f},   //
                                      ww::tup{7.f, 7.f, -6.f, -7.f},  //
                                      ww::tup{1.f, -3.f, 7.f, 4.f});
    EXPECT_EQ(532.f, ww::Determinant(A));
    ww::is_invertible_return const Result = ww::IsInvertible(A);
    EXPECT_EQ(532.f, Result.Determinant);
    EXPECT_EQ(true, Result.IsInvertible);
    EXPECT_EQ(ww::Cofactor(A, 2, 3), -160.f);
    EXPECT_EQ(ww::Cofactor(A, 3, 2), 105.f);

    float const DetA = ww::Determinant(A);
    float const CfA23 = ww::Cofactor(A, 2, 3);
    float const CfA32 = ww::Cofactor(A, 3, 2);
    float const B32 = CfA23 / DetA;
    float const B23 = CfA32 / DetA;
    EXPECT_EQ(B32, -160.f / 532.f);
    EXPECT_EQ(B23, 105.f / 532.f);

    ww::matrix const IA = ww::Inverse(A);
    ww::matrix const CheckA = ww::Matrix44(ww::tup{0.21805, 0.45113, 0.24060, -0.04511},    //!<
                                           ww::tup{-0.80827, -1.45677, -0.44361, 0.52068},  //!<
                                           ww::tup{-0.07895, -0.22368, -0.05263, 0.19737},  //!<
                                           ww::tup{-0.52256, -0.81391, -0.30075, 0.30639});
    EXPECT_EQ(ww::Equal(IA, CheckA), true);
  }
  {
    ww::matrix const A = ww::Matrix44(ww::tup{8.f, -5.f, 9.f, 2.f},  //
                                      ww::tup{7.f, 5.f, 6.f, 1.f},   //
                                      ww::tup{-6.f, 0.f, 9.f, 6.f},  //
                                      ww::tup{-3.f, 0.f, -9.f, -4.f});
    ww::matrix const IA = ww::Inverse(A);
    ww::matrix const CheckA = ww::Matrix44(ww::tup{-0.15385, -0.15385, -0.28205, -0.53846},  //!<
                                           ww::tup{-0.07692, 0.12308, 0.02564, 0.03077},     //<!
                                           ww::tup{0.35897, 0.35897, 0.43590, 0.92308},      //!<
                                           ww::tup{-0.69231, -0.69231, -0.76923, -1.92308});
    EXPECT_EQ(ww::Equal(IA, CheckA), true);
  }
  {
    ww::matrix const A = ww::Matrix44(ww::tup{9.f, 3.f, 0.f, 9.f},      //
                                      ww::tup{-5.f, -2.f, -6.f, -3.f},  //
                                      ww::tup{-4.f, 9.f, 6.f, 4.f},     //
                                      ww::tup{-7.f, 6.f, 6.f, 2.f});
    ww::matrix const IA = ww::Inverse(A);
    ww::matrix const CheckA = ww::Matrix44({-0.04074, -0.07778, 0.14444, -0.22222},  //!<
                                           {-0.07778, 0.03333, 0.36667, -0.33333},   //!<
                                           {-0.02901, -0.14630, -0.10926, 0.12963},  //!<
                                           {0.17778, 0.06667, -0.26667, 0.33333});
    EXPECT_EQ(ww::Equal(IA, CheckA), true);
  }
  {
    ww::matrix const A = ww::Matrix44(ww::tup{3.f, -9.f, 7.f, 3.f},   //
                                      ww::tup{3.f, -8.f, 2.f, -9.f},  //
                                      ww::tup{-4.f, 4.f, 4.f, 1.f},   //
                                      ww::tup{-6.f, 5.f, -1.f, 1.f});

    ww::matrix const B = ww::Matrix44(ww::tup{8.f, 2.f, 2.f, 2.f},   //
                                      ww::tup{3.f, -1.f, 7.f, 0.f},  //
                                      ww::tup{7.f, 0.f, 5.f, 4.f},   //
                                      ww::tup{-6.f, -2.f, 0.f, 5.f});

    // NOTE: Test that multiplication with the inverse gives us back the original matrix.
    ww::matrix const C = A * B;  // --> C * B^-1 = A.
    ww::matrix const CheckA = C * ww::Inverse(B);
    EXPECT_EQ(ww::Equal(A, CheckA), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, PuttingItTogether)
{
  // NOTE: 1. Invert the identity matrix.
  {
    // NOTE: Invert the identity matrix
    ww::matrix const Im = ww::I();
    ww::matrix const ImInv = ww::Inverse(Im);
    // NOTE: Verify that the inverse of the identity matrix is the identity.
    EXPECT_EQ(ww::Equal(ImInv, Im), true);
    // std::cout << ImInv << std::endl;
  }

  // NOTE: 2. What do you get when you multiply a matrix with its inverse.
  {
    ww::matrix const B = ww::Matrix44(ww::tup{8.f, 2.f, 2.f, 2.f},   //
                                      ww::tup{3.f, -1.f, 7.f, 0.f},  //
                                      ww::tup{7.f, 0.f, 5.f, 4.f},   //
                                      ww::tup{-6.f, -2.f, 0.f, 5.f});
    ww::matrix const BInv = ww::Inverse(B);
    ww::matrix const C = B * BInv;
    EXPECT_EQ(ww::Equal(ww::I(), C), true);
    // std::cout << C << std::endl;
  }

  // NOTE: 3. Is there any difference ...
  {
    // NOTE: Check the difference between the inverse of the transpose and
    //      the transpose of the inverse.
    ww::matrix const A = ww::Matrix44(ww::tup{8.f, 2.f, 2.f, 2.f},   //
                                      ww::tup{3.f, -1.f, 7.f, 0.f},  //
                                      ww::tup{7.f, 0.f, 5.f, 4.f},   //
                                      ww::tup{-6.f, -2.f, 0.f, 5.f});
    ww::matrix const InvOfTranspose = ww::Inverse(ww::Transpose(A));
    ww::matrix const TransposeOfInverse = ww::Transpose(ww::Inverse(A));
    EXPECT_EQ(ww::Equal(InvOfTranspose, TransposeOfInverse), true);

    // std::cout << InvOfTranspose << std::endl;
    //{  -0.30000 ,   1.25185 ,   0.30741 ,   0.14074 }
    //{  -0.20000 ,   0.52593 ,   0.30370 ,  -0.02963 }
    //{   0.40000 ,  -1.23704 ,  -0.34815 ,  -0.01481 }
    //{  -0.20000 ,   0.48889 ,   0.15556 ,   0.15556 }

    // std::cout << TransposeOfInverse << std::endl;
    //{  -0.30000 ,   1.25185 ,   0.30741 ,   0.14074 }
    //{  -0.20000 ,   0.52593 ,   0.30370 ,  -0.02963 }
    //{   0.40000 ,  -1.23704 ,  -0.34815 ,  -0.01481 }
    //{  -0.20000 ,   0.48889 ,   0.15556 ,   0.15556 }
    // NOTE: So to conclude they are equal.
  }
  // NOTE: 4. Multiply matrix by tuple ...
  {
    ww::tup const T = ww::tup{0.1f, 0.2f, 0.3f, 0.4f};
    ww::matrix M = ww::I();
    // ww::Set(M, 0, 0, 0.1f); // set an element to a different number
    ww::Set(M, 1, 0, 10.f);  // set an element to a different number
    ww::Set(M, 2, 0, 10.f);  // set an element to a different number
    ww::tup const T2 = M * T;
    // std::cout << T2 << std::endl;
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Translation)
{
  // NOTE: A simple translation.
  // from X:-3 -> X:2 is a translation of 5 units.
  // from Y:4 -> Y:1 is a translation of -3 units.
  // from Z:5 -> Z:2 is a translation of 2 units.
  {
    ww::tup const P = ww::Point(-3.f, 4.f, 5.f);
    ww::tup const T = ww::Translation(5.f, -3.f, 2.f) * P;
    EXPECT_EQ(ww::Equal(T, ww::Point(2.f, 1.f, 7.f)), true);
  }
  // NOTE: move the translation in reverse by multiplication with the inverse
  {
    ww::matrix const Transform = ww::Translation(5.f, -3.f, 2.f);
    ww::matrix const InvTransf = ww::Inverse(Transform);
    ww::tup const P = ww::Point(-3.f, 4.f, 5.f);
    ww::tup const P2 = InvTransf * P;
    EXPECT_EQ(ww::Equal(P2, ww::Point(-8.f, 7.f, 3.f)), true);
  }
  // NOTE: Multiplying a vector with a translation should not change the vector, it is just
  //       an arrow with a different position.
  {
    ww::tup V = ww::Vector(-3.f, 4.f, 5.f);
    ww::tup const V2 = ww::Translation(5.f, -3.f, 2.f) * V;
    EXPECT_EQ(ww::Equal(V2, V), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Scaling)
{
  // NOTE: Scaling applied to a point.
  {
    ww::tup const P = ww::Point(-4.f, 6.f, 8.f);
    ww::tup const P2 = ww::Scale(2.f, 3.f, 4.f) * P;
    EXPECT_EQ(ww::Equal(P2, ww::Point(-8.f, 18.f, 32.f)), true);
  }
  // NOTE: Scaling applied to a vector.
  {
    ww::tup const V = ww::Vector(-4.f, 6.f, 8.f);
    ww::tup const V2 = ww::Scale(2.f, 3.f, 4.f) * V;
    EXPECT_EQ(ww::Equal(V2, ww::Vector(-8.f, 18.f, 32.f)), true);
  }
  // NOTE: Expect the tuple to "grow" in the opposite direction when scaling with the inverse.
  {
    ww::matrix const InvTransform = ww::Inverse(ww::Scale(2.f, 3.f, 4.f));
    ww::tup const V = InvTransform * ww::Vector(-4.f, 6.f, 8.f);
    EXPECT_EQ(ww::Equal(V, ww::Vector(-2.f, 2.f, 2.f)), true);
  }
  // NOTE: Test that reflection is the same a changing the sign of one axis.
  {
      ww::tup const P = ww::Scale(-1.f, 1.f, 1.f) * ww::Point(2.f, 3.f, 4.f);
      EXPECT_EQ(ww::Equal(P, ww::Point(-2.f, 3.f, 4.f)), true);
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
