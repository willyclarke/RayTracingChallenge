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

constexpr float EPSILON = 0.0000001;  // 1E27 * std::numeric_limits<float>::min();

// NOTE: Declarations.
tup Add(tup const &A, tup const &B);
tup Cross(tup const &A, tup const &B);
float Dot(tup const &A, tup const &B);
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
ww::tup operator-(ww::tup const &Tup);
ww::tup operator*(float const S, ww::tup const &Tup);
ww::tup operator*(ww::tup const &Tup, float const S);
ww::tup operator/(ww::tup const &Tup, float const S);
#endif
