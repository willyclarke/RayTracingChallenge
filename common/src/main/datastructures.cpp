/******************************************************************************
 * Filename : datastructures.cpp
 * Date     : 2018 Aug 26
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: This file contains the code that should be common to the
 *          : various applications in this project.
 ******************************************************************************/
#include "datastructures.hpp"

#include <cmath>
#include <iostream>

namespace ww
{
//------------------------------------------------------------------------------
tup Add(tup const &A, tup const &B)
{
   tup const Result = {A.X + B.X, A.Y + B.Y, A.Z + B.Z, A.W + B.W};
   return (Result);
}

//------------------------------------------------------------------------------
bool Equal(float const A, float const B)
{
   if (std::fabs(A - B) < EPSILON)
   {
      return true;
   }
   return false;
}

//------------------------------------------------------------------------------
bool Equal(tup const &A, tup const &B)
{
   bool const Result = Equal(A.X, B.X) &&  //<!
                       Equal(A.Y, B.Y) &&  //<!
                       Equal(A.Z, B.Z) &&  //<!
                       Equal(A.W, B.W);
   return (Result);
}

//------------------------------------------------------------------------------
bool IsPoint(tup const &Tup)
{
   bool Result{};
   Result = (Tup.W != 0);
   return (Result);
}

//------------------------------------------------------------------------------
bool IsVector(tup const &Tup)
{
   bool Result{};
   Result = !(Tup.W != 0);
   return (Result);
}

//------------------------------------------------------------------------------
tup Multiply(float const S, tup const &Tup)
{
   tup const Result{S * Tup.X, S * Tup.Y, S * Tup.Z, S * Tup.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Negate(tup const &Tup)
{
   tup const Result{-Tup.X, -Tup.Y, -Tup.Z, -Tup.W};
   return (Result);
}

//------------------------------------------------------------------------------
float MagSquared(tup const &Tup)
{
    float const Result = Tup.X * Tup.X + //<!
        Tup.Y * Tup.Y + //<!
        Tup.Z * Tup.Z + //<!
        Tup.W * Tup.W ; //<!
    return (Result);
}

//------------------------------------------------------------------------------
float Mag(tup const &Tup)
{
    return (std::sqrt(MagSquared(Tup)));
}

//------------------------------------------------------------------------------
tup Point(float A, float B, float C)
{
   tup Result{A, B, C, 1.f};
   return (Result);
}

//------------------------------------------------------------------------------
tup Sub(tup const &A, tup const &B)
{
   tup const Result = {A.X - B.X, A.Y - B.Y, A.Z - B.Z, A.W - B.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Vector(float A, float B, float C)
{
   tup Result{A, B, C, 0.f};
   return (Result);
}
};  // namespace ww

// ---
// NOTE: Stream operator
// ---
std::ostream &operator<<(std::ostream &stream, const ww::tup &T)
{
   stream << ((T.W != 0) ? "Point:" : "Vector:");
   stream << std::setw(4) << T.X         //<!
          << " " << std::setw(4) << T.Y  //<!
          << " " << std::setw(4) << T.Z  //<!
          << " " << std::setw(4) << T.W;
   return stream;
}

// ---
// NOTE: The Negate operator.
// ---
ww::tup operator-(ww::tup const &Tup) { return (ww::Negate(Tup)); }
ww::tup operator*(float const S, ww::tup const &Tup) { return (ww::Multiply(S, Tup)); }
ww::tup operator*(ww::tup const &Tup, float const S) { return (ww::Multiply(S, Tup)); }
// ---
// NOTE: Division operator does not check for divide by zero; Who cares?
// ---
ww::tup operator/(ww::tup const &Tup, float const S) { return (ww::Multiply(1.f / S, Tup)); }
