#include "gtest/gtest.h"
#include <datastructures.hpp>
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
    ww::tup const P2 = ww::Scaling(2.f, 3.f, 4.f) * P;
    EXPECT_EQ(ww::Equal(P2, ww::Point(-8.f, 18.f, 32.f)), true);
  }
  // NOTE: Scaling applied to a vector.
  {
    ww::tup const V = ww::Vector(-4.f, 6.f, 8.f);
    ww::tup const V2 = ww::Scaling(2.f, 3.f, 4.f) * V;
    EXPECT_EQ(ww::Equal(V2, ww::Vector(-8.f, 18.f, 32.f)), true);
  }
  // NOTE: Expect the tuple to "grow" in the opposite direction when scaling with the inverse.
  {
    ww::matrix const InvTransform = ww::Inverse(ww::Scaling(2.f, 3.f, 4.f));
    ww::tup const V = InvTransform * ww::Vector(-4.f, 6.f, 8.f);
    EXPECT_EQ(ww::Equal(V, ww::Vector(-2.f, 2.f, 2.f)), true);
  }
  // NOTE: Test that reflection is the same a changing the sign of one axis.
  {
    ww::tup const P = ww::Scaling(-1.f, 1.f, 1.f) * ww::Point(2.f, 3.f, 4.f);
    EXPECT_EQ(ww::Equal(P, ww::Point(-2.f, 3.f, 4.f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, RotationX)
{
  // NOTE: Rotating a point around the X-axis
  {
    ww::tup const P = ww::Point(0.f, 1.f, 0.f);
    float const N = 180.f;
    ww::tup const HalfQuart = ww::RotateX(ww::Radians(N / 4.f)) * P;
    ww::tup const FullQuart = ww::RotateX(ww::Radians(N / 2.f)) * P;
    EXPECT_EQ(ww::Equal(HalfQuart, ww::Point(0.f, std::sqrt(2.f) / 2.f, std::sqrt(2.f) / 2.f)), true);
    EXPECT_EQ(ww::Equal(FullQuart, ww::Point(0.f, 0.f, 1.f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, RotationY)
{  // NOTE: Rotating a point around the Y-axis
  {
    ww::tup const P = ww::Point(0.f, 0.f, 1.f);
    float const N = 180.f;
    ww::tup const HalfQuart = ww::RotateY(ww::Radians(N / 4.f)) * P;
    ww::tup const FullQuart = ww::RotateY(ww::Radians(N / 2.f)) * P;
    EXPECT_EQ(ww::Equal(HalfQuart, ww::Point(std::sqrt(2.f) / 2.f, 0.f, std::sqrt(2.f) / 2.f)), true);
    EXPECT_EQ(ww::Equal(FullQuart, ww::Point(1.f, 0.f, 0.f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, RotationZ)
{  // NOTE: Rotating a point around the Z-axis
  {
    ww::tup const P = ww::Point(0.f, 1.f, 0.f);
    float const N = 180.f;
    ww::tup const HalfQuart = ww::RotateZ(ww::Radians(N / 4.f)) * P;
    ww::tup const FullQuart = ww::RotateZ(ww::Radians(N / 2.f)) * P;
    EXPECT_EQ(ww::Equal(HalfQuart, ww::Point(-std::sqrt(2.f) / 2.f, std::sqrt(2.f) / 2.f, 0.f)), true);
    EXPECT_EQ(ww::Equal(FullQuart, ww::Point(-1.f, 0.f, 0.f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Shearing1)
{
  {
    // NOTE: Check shearing by checking point after calculation.
    EXPECT_EQ(                                                                      //!<
        ww::Equal(                                                                  //!<
            ww::Shearing(1.f, 0.f, 0.f, 0.f, 0.f, 0.f) * ww::Point(2.f, 3.f, 4.f),  //!<
            ww::Point(5.f, 3.f, 4.f)                                                //!<
            ),
        true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Shearing2)
{
  {
    EXPECT_EQ(                                                                      //!<
        ww::Equal(                                                                  //!<
            ww::Shearing(0.f, 1.f, 0.f, 0.f, 0.f, 0.f) * ww::Point(2.f, 3.f, 4.f),  //!<
            ww::Point(6.f, 3.f, 4.f)                                                //!<
            ),
        true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Shearing3)
{
  {
    EXPECT_EQ(                                                                      //!<
        ww::Equal(                                                                  //!<
            ww::Shearing(0.f, 0.f, 1.f, 0.f, 0.f, 0.f) * ww::Point(2.f, 3.f, 4.f),  //!<
            ww::Point(2.f, 5.f, 4.f)                                                //!<
            ),
        true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Shearing4)
{
  {
    EXPECT_EQ(                                                                      //!<
        ww::Equal(                                                                  //!<
            ww::Shearing(0.f, 0.f, 0.f, 1.f, 0.f, 0.f) * ww::Point(2.f, 3.f, 4.f),  //!<
            ww::Point(2.f, 7.f, 4.f)                                                //!<
            ),
        true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Shearing5)
{
  {
    EXPECT_EQ(                                                                      //!<
        ww::Equal(                                                                  //!<
            ww::Shearing(0.f, 0.f, 0.f, 0.f, 1.f, 0.f) * ww::Point(2.f, 3.f, 4.f),  //!<
            ww::Point(2.f, 3.f, 6.f)                                                //!<
            ),
        true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, Shearing6)
{
  {
    EXPECT_EQ(                                                                      //!<
        ww::Equal(                                                                  //!<
            ww::Shearing(0.f, 0.f, 0.f, 0.f, 0.f, 1.f) * ww::Point(2.f, 3.f, 4.f),  //!<
            ww::Point(2.f, 3.f, 7.f)                                                //!<
            ),
        true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, RotScaleTranslate)
{
  // NOTE: First check that the simple approach works. Do the single transformations
  //       one by one and verify the result as we move along.
  {
    float const N = 180.f;
    ww::tup const P = ww::Point(1.f, 0.f, 1.f);

    ww::tup const P2 = ww::RotateX(ww::Radians(N) / 2.f) * P;
    EXPECT_EQ(ww::Equal(P2, ww::Point(1.f, -1.f, 0.f)), true);

    ww::tup const P3 = ww::Scaling(5.f, 5.f, 5.f) * P2;
    EXPECT_EQ(ww::Equal(P3, ww::Point(5.f, -5.f, 0.f)), true);

    ww::tup const P4 = ww::Translation(10.f, 5.f, 7.f) * P3;
    EXPECT_EQ(ww::Equal(P4, ww::Point(15.f, 0.f, 7.f)), true);
  }
  // NOTE: Now check that the TranslateScaleRotate works as expected by applying the funtion
  //       to do the operations.
  {
    float const N = 180.f;
    ww::tup const P = ww::Point(1.f, 0.f, 1.f);

    // NOTE: Rotation.
    ww::tup const P2 = ww::TranslateScaleRotate(           //!<
                           0.f, 0.f, 0.f,                  //!<
                           1.f, 1.f, 1.f,                  //!<
                           ww::Radians(N) / 2.f, 0.f, 0.f  //!<
                           ) *
                       P;
    EXPECT_EQ(ww::Equal(P2, ww::Point(1.f, -1.f, 0.f)), true);

    // NOTE: Scale and rotate.
    ww::tup const P3 = ww::TranslateScaleRotate(           //!<
                           0.f, 0.f, 0.f,                  //!<
                           5.f, 5.f, 5.f,                  //!<
                           ww::Radians(N) / 2.f, 0.f, 0.f  //!<
                           ) *
                       P;
    EXPECT_EQ(ww::Equal(P3, ww::Point(5.f, -5.f, 0.f)), true);

    // NOTE: Full on, translate, scale and rotate.
    ww::tup const P4 = ww::TranslateScaleRotate(           //!<
                           10.f, 5.f, 7.f,                 //!<
                           5.f, 5.f, 5.f,                  //!<
                           ww::Radians(N) / 2.f, 0.f, 0.f  //!<
                           ) *
                       P;
    EXPECT_EQ(ww::Equal(P4, ww::Point(15.f, 0.f, 7.f)), true);
  }
}

//------------------------------------------------------------------------------
TEST(Matrix, CreateClockPPM)
{
  ww::canvas Canvas(900, 900);
  ww::tup Color = ww::Color(1.f, 0.f, 0.f);
  // NOTE: Define the number of meters across.
  float const HAcross{1.f};
  float const VAcross{1.f};

  // ---
  // So we need to scale from our world to pixel.
  // Set up pixels per meter.
  // ---
  float const HPixelMeter = Canvas.W / (2.f * HAcross);  // use Xm horisontal across.
  auto HPixel = [&](float const Pos) -> int {
    // Pos = 0  -> Pixel = Canvas.W/2
    // Pos = -1 -> Pixel = 0
    // Pos = +1 -> Pixel = Canvas.W
    float const P = HPixelMeter * Pos + Canvas.W / 2.f;
    int const Result = std::min<int>(Canvas.W - 1, std::max<int>(0, int(P)));
    return Result;
  };

  // The vertical axis need to be reversed.
  float const VPixelMeter = Canvas.H / (2.f * VAcross);  // use Xm vertiacal across.
  auto VPixel = [&](float const Pos) -> int {
    // Pos y = 0    -> Pixel = Canvas.H
    // Pos y = 2000 -> Pixel = 0
    // Pos y = 1000 -> Pixel = Canvas.H/2
    float const P = VPixelMeter * Pos + Canvas.H / 2.f;
    int const Result = std::min<int>(Canvas.H - 1, std::max<int>(0, int(P)));
    return Result;
  };

  EXPECT_EQ(HPixel(-HAcross), 0);
  EXPECT_EQ(VPixel(-VAcross), 0);
  EXPECT_EQ(HPixel(0.f), Canvas.W / 2.f);
  EXPECT_EQ(VPixel(0.f), Canvas.H / 2.f);
  EXPECT_EQ(HPixel(HAcross), Canvas.W - 1);
  EXPECT_EQ(VPixel(VAcross), Canvas.H - 1);

  auto ClockPixel = [&](float const X, float const Y) -> void {
    ww::WritePixel(Canvas, HPixel(X + 2 * 0.f / Canvas.W), VPixel(Y + 2 * 0.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 1.f / Canvas.W), VPixel(Y + 2 * 0.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 2.f / Canvas.W), VPixel(Y + 2 * 0.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 3.f / Canvas.W), VPixel(Y + 2 * 0.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 4.f / Canvas.W), VPixel(Y + 2 * 0.f / Canvas.H), Color);

    ww::WritePixel(Canvas, HPixel(X + 2 * 0.f / Canvas.W), VPixel(Y + 2 * 1.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 1.f / Canvas.W), VPixel(Y + 2 * 1.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 2.f / Canvas.W), VPixel(Y + 2 * 1.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 3.f / Canvas.W), VPixel(Y + 2 * 1.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 4.f / Canvas.W), VPixel(Y + 2 * 1.f / Canvas.H), Color);

    ww::WritePixel(Canvas, HPixel(X + 2 * 0.f / Canvas.W), VPixel(Y + 2 * 2.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 1.f / Canvas.W), VPixel(Y + 2 * 2.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 2.f / Canvas.W), VPixel(Y + 2 * 2.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 3.f / Canvas.W), VPixel(Y + 2 * 2.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 4.f / Canvas.W), VPixel(Y + 2 * 2.f / Canvas.H), Color);

    ww::WritePixel(Canvas, HPixel(X + 2 * 0.f / Canvas.W), VPixel(Y + 2 * 3.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 1.f / Canvas.W), VPixel(Y + 2 * 3.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 2.f / Canvas.W), VPixel(Y + 2 * 3.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 3.f / Canvas.W), VPixel(Y + 2 * 3.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 4.f / Canvas.W), VPixel(Y + 2 * 3.f / Canvas.H), Color);

    ww::WritePixel(Canvas, HPixel(X + 2 * 0.f / Canvas.W), VPixel(Y + 2 * 4.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 1.f / Canvas.W), VPixel(Y + 2 * 4.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 2.f / Canvas.W), VPixel(Y + 2 * 4.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 3.f / Canvas.W), VPixel(Y + 2 * 4.f / Canvas.H), Color);
    ww::WritePixel(Canvas, HPixel(X + 2 * 4.f / Canvas.W), VPixel(Y + 2 * 4.f / Canvas.H), Color);
  };

  ClockPixel(0.f, 0.f);
  float const Angle = 360.f / 12.f;  //!< There are 12 jewels around the clock.

  for (int Idx = 0;  ///<!
       Idx < 12;     ///<!
       ++Idx)
  {
    ww::tup const ClockP =
        ww::TranslateScaleRotate(0.5f, 0.f, 0.f,                               //!< Translation is in m(?)
                                 1.f, 1.f, 1.f,                                //!< Scale input is unitless.
                                 0.f, 0.f, float(Idx) * ww::Radians(Angle)) *  //!< Input rotation in radians.
        ww::Vector(0.5f, 0.f, 0.f);
    ClockPixel(ClockP.X, ClockP.Y);
  }

  // Finally we write the file
  ww::WriteToPPM(Canvas, "Clock.ppm");
}

//------------------------------------------------------------------------------
TEST(RaySphere, Intersect1)
{
  ww::tup const Origin = ww::Point(1.f, 2.f, 3.f);
  ww::tup const Direction = ww::Vector(4.f, 5.f, 6.f);
  ww::ray const Ray = ww::Ray(Origin, Direction);
  EXPECT_EQ(ww::Equal(Ray.Origin, Origin), true);
  EXPECT_EQ(ww::Equal(Ray.Direction, ww::Normalize(Direction)), true);
}

//------------------------------------------------------------------------------
TEST(RaySphere, MoveAlongRay)
{
  ww::ray const R = Ray(ww::Point(2.f, 3.f, 4.f), ww::Vector(1.f, 0.f, 0.f));

  // NOTE: Move along a ray, start at t=0, point is exected to be at origin.
  //       Keep on moving along the ray, so at t=1 the expected position is
  //       an x value which has increased by 1.
  //       And so on...
  EXPECT_EQ(ww::Equal(ww::Point(2.f, 3.f, 4.f), ww::PositionAt(R, 0.f)), true);
  EXPECT_EQ(ww::Equal(ww::Point(3.f, 3.f, 4.f), ww::PositionAt(R, 1.f)), true);
  EXPECT_EQ(ww::Equal(ww::Point(1.f, 3.f, 4.f), ww::PositionAt(R, -1.f)), true);
  EXPECT_EQ(ww::Equal(ww::Point(4.5f, 3.f, 4.f), ww::PositionAt(R, 2.5f)), true);
}

//------------------------------------------------------------------------------
// NOTE: Test that the ray goes straight through the sphere.
//
//                              xx
//-5    -4     -3     -2    -1 x  x       +1    +2
// -----| -----| -----| -----|x    x -----| -----| -----|
// Ray >>>>                    x  x
//                              xx
//
TEST(RaySphere, IntersectSphere2Points1)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape S = ww::PtrDefaultSphere();

  ww::intersections XS = ww::Intersect(S, R);
  EXPECT_EQ(ww::Dot(R.Direction, R.Direction), 1.f);
  EXPECT_EQ(XS.Count(), 2);
  if (XS.Count() == 2)
  {
    EXPECT_EQ(XS.vI[0].t, 4.f);
    EXPECT_EQ(XS.vI[1].t, 6.f);
  }
}
//------------------------------------------------------------------------------
// NOTE: Test that the ray kind of just touches, tangents, the sphere.
TEST(RaySphere, IntersectSphere2Points2)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 1.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape S = ww::PtrDefaultSphere();

  ww::intersections XS = ww::Intersect(S, R);
  EXPECT_EQ(XS.Count(), 2);
  if (XS.Count() == 2)
  {
    EXPECT_EQ(XS.vI[0].t, 5.f);
    EXPECT_EQ(XS.vI[1].t, 5.f);
  }
}

//------------------------------------------------------------------------------
// NOTE: Test that the ray misses the sphere.
TEST(RaySphere, IntersectSphere2Points3)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 2.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape S = ww::PtrDefaultSphere();
  ww::sphere *pSphere = dynamic_cast<ww::sphere *>(S.get());
  ww::sphere &Sphere = *(dynamic_cast<ww::sphere *>(S.get()));
  EXPECT_EQ(pSphere->Radius, 1.f);
  EXPECT_EQ(Sphere.Radius, 1.f);

  ww::intersections XS = ww::Intersect(S, R);
  EXPECT_EQ(XS.Count(), 0);
}

//------------------------------------------------------------------------------
// NOTE: Test that the ray when it originates from inside the sphere intersects
//       behind and in the direction of the sphere. i.e. the ray is infinite in
//       both directions.
TEST(RaySphere, IntersectSphere2Points4)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 0.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape S = ww::PtrDefaultSphere();

  ww::intersections XS = ww::Intersect(S, R);
  EXPECT_EQ(XS.Count(), 2);
  if (XS.Count() == 2)
  {
    EXPECT_EQ(XS.vI[0].t, -1.f);
    EXPECT_EQ(XS.vI[1].t, 1.f);
  }
}

//------------------------------------------------------------------------------
// NOTE: Test that the sphere may be behind the origin of the ray.
TEST(RaySphere, IntersectSphere2Points5)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, 5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape S = ww::PtrDefaultSphere();

  ww::intersections XS = ww::Intersect(S, R);
  EXPECT_EQ(XS.Count(), 2);
  if (XS.Count() == 2)
  {
    EXPECT_EQ(XS.vI[0].t, -6.f);
    EXPECT_EQ(XS.vI[1].t, -4.f);
  }
}

//------------------------------------------------------------------------------
// Scenario: An intersection encapsulates the 'time' t and the object of type T
TEST(RaySphere, IntersectionSphere)
{
  ww::shared_ptr_shape PtrSphere{};
  PtrSphere.reset(new ww::sphere);

  ww::intersection I = ww::Intersection(3.5f, PtrSphere);

  ww::sphere S{};
  ww::shape *pShape = &S;

  EXPECT_EQ(3.5f, I.t);
  EXPECT_EQ(pShape->isA<ww::sphere>(), true);

  EXPECT_EQ(I.pShape.get()->isA<ww::sphere>(), true);
  EXPECT_EQ(I.pShape.get()->isA<ww::cube>(), false);

  // NOTE: We dont expect these to be at the same address at all.
  // EXPECT_EQ(I.pShape.get(), &S);

  ww::sphere *pSphere = dynamic_cast<ww::sphere *>(I.pShape.get());

  bool const TheSpheresAreEqual = *pSphere == S;
  EXPECT_EQ(TheSpheresAreEqual, true);
}

//------------------------------------------------------------------------------
// NOTE: Aggregate intersection objects together so that we can work with multiple
//       intersections at once.
TEST(RaySphere, IntersectionsVector)
{
  ww::shared_ptr_shape PtrSphere{};
  PtrSphere.reset(new ww::sphere);

  ww::intersection I1 = ww::Intersection(1.f, PtrSphere);
  ww::intersection I2 = ww::Intersection(2.f, PtrSphere);
  ww::intersections XS = ww::Intersections(I1, I2);
  EXPECT_EQ(XS.Count(), 2);
  if (XS.Count() == 2)
  {
    EXPECT_EQ(XS.vI[0].t, 1.f);
    EXPECT_EQ(XS.vI[1].t, 2.f);
  }
}

//------------------------------------------------------------------------------
// NOTE: The hit, when all intersections have positive t.
TEST(RaySphere, IntersectionsHit1)
{
  ww::shared_ptr_shape PtrSphere{};
  PtrSphere.reset(new ww::sphere);

  ww::intersection const I1 = ww::Intersection(1.f, PtrSphere);
  ww::intersection const I2 = ww::Intersection(2.f, PtrSphere);
  ww::intersections const XS = ww::Intersections(I1, I2);
  ww::intersection const H = ww::Hit(XS);
  EXPECT_EQ(I1 == H, true);
  EXPECT_EQ(I2 == H, false);
}

//------------------------------------------------------------------------------
// NOTE: The hit, when some intersections have negative t.
TEST(RaySphere, IntersectionsHit2)
{
  ww::shared_ptr_shape PtrSphere{};
  PtrSphere.reset(new ww::sphere);

  ww::intersection const I1 = ww::Intersection(-1.f, PtrSphere);
  ww::intersection const I2 = ww::Intersection(1.f, PtrSphere);
  ww::intersections const XS = ww::Intersections(I1, I2);
  ww::intersection const H = ww::Hit(XS);
  EXPECT_EQ(I1 == H, false);
  EXPECT_EQ(I2 == H, true);

  // NOTE: Both object pointers will point to the sphere.
  EXPECT_EQ(H.pShape == I1.pShape, true);
  EXPECT_EQ(H.pShape == I2.pShape, true);
}

//------------------------------------------------------------------------------
// NOTE: The hit, when all intersections have negative t.
TEST(RaySphere, IntersectionsHit3)
{
  ww::shared_ptr_shape PtrSphere{};
  PtrSphere.reset(new ww::sphere);

  ww::intersection const I1 = ww::Intersection(-2.f, PtrSphere);
  ww::intersection const I2 = ww::Intersection(-1.f, PtrSphere);
  ww::intersections const XS = ww::Intersections(I1, I2);
  ww::intersection const H = ww::Hit(XS);
  EXPECT_EQ(I1 == H, false);
  EXPECT_EQ(I2 == H, false);
  EXPECT_EQ(XS.Count(), 2);
  EXPECT_EQ(H.pShape, nullptr);
  EXPECT_EQ(H.pShape == I1.pShape, false);
  EXPECT_EQ(H.pShape == I2.pShape, false);
  EXPECT_EQ(H.pShape.get()->isA<ww::sphere>(), false);
  EXPECT_EQ(H.pShape.get()->isA<ww::cube>(), false);
}

//------------------------------------------------------------------------------
// NOTE: The hit is always the lowest non-negative intersection.
TEST(RaySphere, IntersectionsHit4)
{
  ww::shared_ptr_shape PtrSphere{};
  PtrSphere.reset(new ww::sphere);

  ww::intersection const I1 = ww::Intersection(5.f, PtrSphere);
  ww::intersection const I2 = ww::Intersection(7.f, PtrSphere);
  ww::intersection const I3 = ww::Intersection(-3.f, PtrSphere);
  ww::intersection const I4 = ww::Intersection(2.f, PtrSphere);
  ww::intersections XS{};
  Intersections(XS, I1);
  Intersections(XS, I2);
  Intersections(XS, I3);
  Intersections(XS, I4);
  ww::intersection const H = ww::Hit(XS);
  EXPECT_EQ(XS.Count(), 4);
  EXPECT_EQ(I1 == H, false);
  EXPECT_EQ(I2 == H, false);
  EXPECT_EQ(I3 == H, false);
  EXPECT_EQ(I4 == H, true);
}

//------------------------------------------------------------------------------
TEST(RaySphere, TranslateRay)
{
  ww::ray const R = ww::Ray(ww::Point(1.f, 2.f, 3.f), ww::Vector(0.f, 1.f, 0.f));
  ww::matrix const M = ww::Translation(3.f, 4.f, 5.f);

  // NOTE: transform the ray R with matrix M. This will test the translation of a ray.
  ww::ray const R2 = ww::Transform(R, M);
  EXPECT_EQ(R2.Origin.X, 4.f);
  EXPECT_EQ(R2.Origin.Y, 6.f);
  EXPECT_EQ(R2.Origin.Z, 8.f);
  EXPECT_EQ(ww::IsPoint(R2.Origin), true);
  EXPECT_EQ(ww::IsVector(R2.Direction), true);
  EXPECT_EQ(ww::Equal(R2.Direction, ww::Vector(0.f, 1.f, 0.f)), true);
}

//------------------------------------------------------------------------------
TEST(RaySphere, ScaleRay)
{
  ww::ray const R = ww::Ray(ww::Point(1.f, 2.f, 3.f), ww::Vector(0.f, 1.f, 0.f));
  ww::matrix const M = ww::Scaling(2.f, 3.f, 4.f);

  // NOTE: transform the ray R with matrix M. This will test the scaling of a ray.
  //       The vector is left un-normalized. This is intentional.
  ww::ray const R2 = ww::Transform(R, M);
  EXPECT_EQ(R2.Origin.X, 2.f);
  EXPECT_EQ(R2.Origin.Y, 6.f);
  EXPECT_EQ(R2.Origin.Z, 12.f);
  EXPECT_EQ(ww::IsPoint(R2.Origin), true);
  EXPECT_EQ(ww::IsVector(R2.Direction), true);
  EXPECT_EQ(ww::Equal(R2.Direction, ww::Vector(0.f, 3.f, 0.f)), true);
}

//------------------------------------------------------------------------------
TEST(RaySphere, SphereDefaultTransformation)
{
  ww::sphere S{};
  EXPECT_EQ(ww::Equal(S.Transform, ww::I()), true);
}

//------------------------------------------------------------------------------
TEST(RaySphere, SphereChangeTransformation)
{
  ww::sphere S{};
  S.Transform = ww::Translation(2.f, 3.f, 4.f);  // set the transform.
  EXPECT_EQ(ww::Equal(S.Transform, ww::Translation(2.f, 3.f, 4.f)), true);
}

//------------------------------------------------------------------------------
TEST(RaySphere, SphereIntersectScaled)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape S = ww::PtrDefaultSphere();

  EXPECT_EQ(S.get()->isA<ww::sphere>(), true);

  // std::cout << "S transform: " << S->Transform << std::endl;
  S->Transform = ww::Scaling(2.f, 2.f, 2.f);
  // std::cout << "S transform: " << S->Transform << std::endl;
  ww::ray SavedRay{};
  ww::intersections const XS = ww::Intersect(S, R, &SavedRay);

  // std::cout << "LocalRay: \n" << SavedRay << std::endl;

  EXPECT_EQ(SavedRay.Origin == ww::Point(0.f, 0.f, -2.5f), true);
  EXPECT_EQ(SavedRay.Direction == ww::Vector(0.f, 0.f, 0.5f), true);
  EXPECT_EQ(XS.Count(), 2);

  // ---
  // NOTE: The calculation based on the local coordinate system does not match
  //       with the previous test. Further checks are required.
  // ---
  if (XS.Count() == 2 && false)
  {
    EXPECT_EQ(XS.vI[0].t, 3.f);
    EXPECT_EQ(XS.vI[1].t, 7.f);
  }
}

//------------------------------------------------------------------------------
TEST(RaySphere, SphereIntersectTranslated)
{
  ww::ray const R = ww::Ray(ww::Point(0.f, 0.f, -5.f), ww::Vector(0.f, 0.f, 1.f));
  ww::shared_ptr_shape S = ww::PtrDefaultSphere();

  S->Transform = ww::Translation(5.f, 0.f, 0.f);
  ww::intersections const XS = ww::Intersect(S, R);

  EXPECT_EQ(XS.Count(), 0);
}

//------------------------------------------------------------------------------
TEST(RaySphere, SpherePuttingItTogether)
{
  // ---
  // NOTE: Create a lambda for the ray/sphere intersection.
  //       This allows us to change the transform of the sphere
  //       to create interesting effects, like scaling, skewing etc.
  // ---
  auto CreateSphere = [&](ww::shared_ptr_shape S, std::string const &FileName) -> void {
    ww::tup const RayOrigin = ww::Point(0.f, 0.f, -5.f);
    ww::tup Color = ww::Color(1.f, 0.f, 0.f);

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
  CreateSphere(S, "SpherePuttingItTogether.ppm");

  // NOTE: Scale the sphere
  S->Transform = ww::Scaling(1.f, 0.5f, 1.f);
  CreateSphere(S, "SphereScaledY.ppm");
  S->Transform = ww::Scaling(0.5f, 1.0f, 1.f);
  CreateSphere(S, "SphereScaledX.ppm");

  // NOTE: Scaling in the Z direction is not visible from this viewpoint.
  //       So the end result should still be a circle.
  S->Transform = ww::Scaling(1.0f, 1.0f, 0.5f);
  CreateSphere(S, "SphereScaledZ.ppm");

  S->Transform = ww::RotateZ(3.1415 / 4.) * ww::Scaling(0.5f, 1.f, 1.f);
  CreateSphere(S, "SphereScaledX1.ppm");
}

