/******************************************************************************
 * Filename : datastructures.hpp
 * Date     : 2018 Aug 25
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Data structures for My Raytracing Challenge.
 ******************************************************************************/
#ifndef COMMON_DATASTRUCTURES_HPP
#define COMMON_DATASTRUCTURES_HPP

#include <iomanip>  // for setw().
#include <iostream>
#include <limits>

namespace ww
{
struct tup
{
   float X{};
   float Y{};
   float Z{};
   float W{};  //!< is 1.0 when tuple is point and 0.0 when tuple is a vector.
};

union tup {
   //!< Ctor
   tup() : X{}, Y{}, Z{}, W{} {};
   tup(float IX, float IY, float IZ, float IW) : X{IX}, Y{IY}, Z{IZ}, W{IW} {};
   ~tup() {}
   tup(tup const &Other)
   {
      X = Other.X;
      Y = Other.Y;
      Z = Other.Z;
      W = Other.W;
   }
   struct
   {
      float X{};
      float Y{};
      float Z{};
      float W{};  //!< is 1.0 when tuple is point and 0.0 when tuple is a vector.
   };
   struct
   {
      float R;
      float G;
      float B;
      float I;  //!< Intensity is 1.0 at max and 0.0 at pitch black.
   };
};  // namespace ww

constexpr float EPSILON = 0.0000001;  // 1E27 * std::numeric_limits<float>::min();

// NOTE: Declarations.
tup Add(tup const &A, tup const &B);
tup Cross(tup const &A, tup const &B);

// ---
// NOTE: For an explanation of what the dot product actually is you can take a look
//     : at the description at
//     : http://betterexplained.com/articles/vector-calculus-understanding-the-dot-product .
// ---
float Dot(tup const &A, tup const &B);

// ---
// NOTE: For a discussion on how to do comparison with floating point number the following web
//     : site can be consulted:
//     https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
// ---
bool Equal(float const A, float const B);
bool Equal(tup const &A, tup const &B);
bool IsVector(tup const &Tup);
bool IsPoint(tup const &Tup);
float MagSquared(tup const &Tup);
float Mag(tup const &Tup);
tup Multiply(float const S, tup const &Tup);
tup Negate(tup const &Tup);
tup Normal(tup const &Tup);
tup Point(float A, float B, float C);
tup Sub(tup const &A, tup const &B);
tup Vector(float A, float B, float C);
};  // namespace ww

// ---
// NOTE: Declare the operator for inclusion elsewhere.
// ---
std::ostream &operator<<(std::ostream &stream, const ww::tup &T);
ww::tup operator+(ww::tup const &A, ww::tup const &B);
ww::tup operator-(ww::tup const &Tup);
ww::tup operator-(ww::tup const &A, ww::tup const &B);
ww::tup operator*(float const S, ww::tup const &Tup);
ww::tup operator*(ww::tup const &Tup, float const S);
ww::tup operator/(ww::tup const &Tup, float const S);
#endif
