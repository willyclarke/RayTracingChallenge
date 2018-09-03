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
#include <iomanip>  // for std::setprecision
#include <iostream>

// ---
// NOTE: Stream operator
// ---
std::ostream &operator<<(std::ostream &stream, const ww::tup &T)
{
   // ---
   // NOTE: The width need to be big enough to hold a negative sign.
   // ---
   size_t const P{5};
   size_t const W{P + 5};
   stream << ((T.W != 0) ? "Point:" : "Vector:");
   stream << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.X  //<!
          << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.Y  //<!
          << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.Z  //<!
          << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.W;
   return stream;
}

namespace ww
{
//------------------------------------------------------------------------------
tup Add(tup const &A, tup const &B)
{
   tup const Result{A.X + B.X, A.Y + B.Y, A.Z + B.Z, A.W + B.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Cross(tup const &A, tup const &B)
{
   tup const Result = Vector(A.Y * B.Z - A.Z * B.Y, A.Z * B.X - A.X * B.Z, A.X * B.Y - A.Y * B.X);
   return (Result);
}

//------------------------------------------------------------------------------
float Dot(tup const &A, tup const &B)
{
   float const Result = A.X * B.X +  //!<
                        A.Y * B.Y +  //<!
                        A.Y * B.Y +  //<!
                        A.Y * B.Y;   //<!
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
float MagSquared(tup const &Tup)
{
   float const Result = Tup.X * Tup.X +  //<!
                        Tup.Y * Tup.Y +  //<!
                        Tup.Z * Tup.Z +  //<!
                        Tup.W * Tup.W;   //<!
   return (Result);
}

//------------------------------------------------------------------------------
float Mag(tup const &Tup) { return (std::sqrt(MagSquared(Tup))); }

//------------------------------------------------------------------------------
tup Negate(tup const &Tup)
{
   tup const Result{-Tup.X, -Tup.Y, -Tup.Z, -Tup.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Normal(tup const &Tup)
{
   tup const Result = Tup / Mag(Tup);
   return (Result);
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
// NOTE: The Negate operator.
// ---
ww::tup operator+(ww::tup const &A, ww::tup const &B) { return (ww::Add(A, B)); }
ww::tup operator-(ww::tup const &Tup) { return (ww::Negate(Tup)); }
ww::tup operator-(ww::tup const &A, ww::tup const &B) { return (ww::Sub(A, B)); }
ww::tup operator*(float const S, ww::tup const &Tup) { return (ww::Multiply(S, Tup)); }
ww::tup operator*(ww::tup const &Tup, float const S) { return (ww::Multiply(S, Tup)); }
// ---
// NOTE: Division operator does not check for divide by zero; Who cares?
// ---
ww::tup operator/(ww::tup const &Tup, float const S) { return (ww::Multiply(1.f / S, Tup)); }
